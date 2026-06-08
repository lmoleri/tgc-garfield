// Garfield++ TGC (Thin Gap Chamber) simulation
//
// Geometry : 10 anode wires (50 μm diameter, 1.8 mm pitch) at y = 0,
//            two grounded cathode planes at y = ±gap_cm.
//            y > 0 is the non-readout ("cathode_top") side;
//            y < 0 is the readout pad ("cathode") side.
// Gas      : Ar:CO2 70:30, 1 atm, 293.15 K.
// Source   : N = E/W primary electrons placed at a configurable distance from
//            the wire plane (signed: positive → readout pad side, y < 0).
//            N is computed from the photon energy and the gas W-value.
// Readout  : (1) all wires as a single "anode" channel,
//            (2) the bottom cathode as a "cathode" channel.
//            Configurable via readout.type:
//              "conductive" (default) – grounded plane at y = −gap_cm.
//              "resistive" – infinitely-thin resistive layer (at y = −gap_cm)
//              on an insulating substrate (Kapton ε_r=3.5 / FR4 ε_r=4.6),
//              with external conductive pads. No charge spreading is simulated;
//              the deposited charge remains at its landing point, but the local
//              surface potential decays to 0 V (grounded edges) with time
//              constant τ = ε₀ε_r ρ_s L²/(π²d). The Ramo weighting potential is
//              corrected for the dielectric layer (α = ε_r·gap/(d+ε_r·gap));
//              the delayed signal (both during drift and after collection) is
//              computed via ComponentUser::SetDelayedWeightingPotential.
//
// Each primary electron is transported by AvalancheMicroscopic.
// Induced signals are computed via Shockley-Ramo weighting fields.

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "nlohmann/json.hpp"

namespace fs = std::filesystem;

namespace {

using Garfield::AvalancheMicroscopic;
using Garfield::ComponentAnalyticField;
using Garfield::ComponentUser;
using Garfield::DriftLineRKF;
using Garfield::MediumMagboltz;
using Garfield::Sensor;
using json = nlohmann::json;


// ─── Configuration structs ────────────────────────────────────────────────────

struct GeometryConfig {
  double wirePitchCm    = 0.18;
  double wireDiamUm     = 50.0;
  double gapCm          = 0.14;
  int    nWires         = 10;
  double wireVoltageV   = 1900.0;
};

struct ReadoutConfig {
  std::string type                    = "conductive"; // "conductive" | "resistive"
  std::string insulatorMaterial       = "kapton";     // "kapton" | "fr4"
  double      insulatorThicknessUm    = 100.0;
  double      surfaceResistivityOhmSq = 500e3;
  bool        enableDelayedSignal     = true;  // if false, skip SetDelayedWeightingPotential
};

struct SourceConfig {
  double energyKeV                                          = 5.9;
  std::optional<std::vector<double>> fixedDistMm            =       // nullopt → uniform random over gap
      std::vector<double>{0.2, 0.5, 0.9, 1.2};
  std::optional<std::vector<double>> fixedXCmList;                  // nullopt → uniform random over wire span
};

struct GasConfig {
  std::string gas1           = "ar";      // first gas species (Magboltz name, lowercase)
  double      frac1          = 70.0;      // fraction of gas1 [%]; gas2 gets 100 − frac1
  std::string gas2           = "co2";     // second gas species
  std::string ionSpecies     = "co2";     // base name for IonMobility_{X}+_{X}.txt
  double temperatureK        = 293.15;
  double pressureTorr        = 760.0;
  bool   enablePenning       = true;
  int    nCollisions         = 10;
  double maxElectronEnergyEV = 2000.0;
  int    nFieldPoints        = 20;         // number of E-field grid points for Magboltz
  double eFieldMaxVcm        = 300000.0;  // upper E-field limit [V/cm] for the gas table
  double wValueEV            = 26.0;      // effective ionisation energy [eV per ion pair]
};

struct SimulationConfig {
  std::size_t nEvents          = 1000;
  std::size_t maxAvalancheSize = 500000;
  double      timeWindowNs     = 300.0;
  double      timeStepNs       = 0.5;
  bool        enableIonDrift   = true;
  bool        storeDriftLines  = false; // store intermediate e⁻ drift steps for 3D viz
};

// ─── 3D-visualisation display limits ─────────────────────────────────────────
constexpr std::size_t kMaxDispIonPaths = 100;  // ion drift paths saved per event
constexpr std::size_t kMaxDispCloudPts = 500;  // avalanche-cloud points saved per event

struct Config {
  GeometryConfig   geometry;
  ReadoutConfig    readout;
  SourceConfig     source;
  GasConfig        gas;
  SimulationConfig simulation;
};

// ─── Per-distance summary ─────────────────────────────────────────────────────

struct DistanceSummary {
  std::optional<double> distanceMm;           // nullopt = random per event
  std::optional<double> xPositionCm;          // nullopt = random per event
  std::size_t nEvents               = 0;
  std::size_t nInteracted           = 0;
  double      interactionFraction   = 0.;
  double      meanAnodeChargeFC     = 0.;
  double      rmsAnodeChargeFC      = 0.;
  double      semAnodeChargeFC      = 0.;
  double      meanCathodeChargeFC      = 0.;
  double      rmsCathodeChargeFC       = 0.;
  double      semCathodeChargeFC       = 0.;
  double      meanCathodeTopChargeFC   = 0.;
  double      rmsCathodeTopChargeFC    = 0.;
  double      semCathodeTopChargeFC    = 0.;
  double      meanChargeRatio       = 0.;
  double      rmsChargeRatio        = 0.;
  double      semChargeRatio        = 0.;
  double      meanPrimaryElectrons  = 0.;
  double      meanAvalancheSize     = 0.;
};

// ─── Utility ──────────────────────────────────────────────────────────────────

std::string FormatNumber(double v, int precision = 4) {
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(precision) << v;
  std::string s = ss.str();
  while (!s.empty() && s.back() == '0') s.pop_back();
  if (!s.empty() && s.back() == '.') s.pop_back();
  return s.empty() ? "0" : s;
}

std::string FileSafeNumber(double v) {
  std::string s = FormatNumber(v);
  std::replace(s.begin(), s.end(), '.', 'p');
  std::replace(s.begin(), s.end(), '-', 'm');
  return s;
}

std::string DeriveGasFileName(const GasConfig& g) {
  auto I = [](double v) {
    return std::to_string(static_cast<long long>(std::llround(v)));
  };
  const int    efKv  = static_cast<int>(std::llround(g.eFieldMaxVcm / 1000.0));
  const double frac2 = 100.0 - g.frac1;
  const std::string prefix = g.gas1 + I(g.frac1) + "_" + g.gas2 + "_" + I(frac2);
  return prefix
       + "_T"  + I(g.temperatureK)
       + "_P"  + I(g.pressureTorr)
       + "_Ee" + I(g.maxElectronEnergyEV)
       + "_Ef" + std::to_string(efKv) + "k"
       + "_n"  + std::to_string(g.nFieldPoints)
       + "_c"  + std::to_string(g.nCollisions)
       + (g.enablePenning ? "_pen" : "_nopen")
       + ".gas";
}

void EnsureDirectory(const fs::path& p) { fs::create_directories(p); }

/// Estimate the peak near-wire electric field using the Sauli 1977 capacitance
/// formula (eq. 2.3): E_peak = V / (r × (ln(pitch/(2π r)) + π gap/pitch)).
/// Returns V/cm.  Used to sanity-check gas.e_field_max_vcm.
double ComputePeakFieldVcm(const GeometryConfig& g) {
  const double rCm = g.wireDiamUm * 0.5e-4;   // wire radius [cm]
  const double cap = std::log(g.wirePitchCm / (2. * M_PI * rCm))
                     + M_PI * g.gapCm / g.wirePitchCm;
  return g.wireVoltageV / (rCm * cap);
}

template <typename T>
double Mean(const std::vector<T>& v) {
  if (v.empty()) return 0.;
  return std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
}

template <typename T>
double Rms(const std::vector<T>& v, double mean) {
  if (v.size() < 2) return 0.;
  double var = 0.;
  for (const auto& x : v) { double d = static_cast<double>(x) - mean; var += d * d; }
  return std::sqrt(var / static_cast<double>(v.size()));
}

double Sem(double rms, std::size_t n) {
  return n < 2 ? 0. : rms / std::sqrt(static_cast<double>(n));
}

// ─── JSON helpers ─────────────────────────────────────────────────────────────

[[noreturn]] void ThrowJsonTypeError(const std::initializer_list<std::string_view>& path,
                                     std::string_view expect) {
  std::string msg = "JSON error at '";
  bool first = true;
  for (auto p : path) { if (!first) msg += '.'; first = false; msg += std::string(p); }
  msg += "': expected " + std::string(expect) + ".";
  throw std::runtime_error(msg);
}

const json* FindMember(const json& obj, std::string_view key,
                       const std::initializer_list<std::string_view>& path) {
  if (!obj.is_object()) ThrowJsonTypeError(path, "an object");
  auto it = obj.find(std::string(key));
  return it == obj.end() ? nullptr : &(*it);
}

const json* FindSection(const json& obj, std::string_view key) {
  auto* p = FindMember(obj, key, {key});
  if (p && !p->is_object()) ThrowJsonTypeError({key}, "an object");
  return p;
}

double ReadDouble(const json& obj, std::string_view sec, std::string_view key, double fb) {
  auto* v = FindMember(obj, key, {sec, key});
  if (!v) return fb;
  if (!v->is_number()) ThrowJsonTypeError({sec, key}, "a number");
  return v->get<double>();
}

int ReadInt(const json& obj, std::string_view sec, std::string_view key, int fb) {
  auto* v = FindMember(obj, key, {sec, key});
  if (!v) return fb;
  if (!v->is_number_integer()) ThrowJsonTypeError({sec, key}, "an integer");
  return v->get<int>();
}

std::size_t ReadSizeT(const json& obj, std::string_view sec, std::string_view key, std::size_t fb) {
  auto* v = FindMember(obj, key, {sec, key});
  if (!v) return fb;
  if (!v->is_number_integer() && !v->is_number_unsigned())
    ThrowJsonTypeError({sec, key}, "a non-negative integer");
  auto val = v->get<long long>();
  if (val < 0) throw std::runtime_error("Expected non-negative integer at key '" + std::string(key) + "'");
  return static_cast<std::size_t>(val);
}

bool ReadBool(const json& obj, std::string_view sec, std::string_view key, bool fb) {
  auto* v = FindMember(obj, key, {sec, key});
  if (!v) return fb;
  if (!v->is_boolean()) ThrowJsonTypeError({sec, key}, "a boolean");
  return v->get<bool>();
}

std::string ReadString(const json& obj, std::string_view sec, std::string_view key,
                       const std::string& fb) {
  auto* v = FindMember(obj, key, {sec, key});
  if (!v) return fb;
  if (!v->is_string()) ThrowJsonTypeError({sec, key}, "a string");
  return v->get<std::string>();
}

std::vector<double> ReadDoubleArray(const json& obj, std::string_view sec, std::string_view key,
                                    const std::vector<double>& fb) {
  auto* v = FindMember(obj, key, {sec, key});
  if (!v) return fb;
  if (!v->is_array()) ThrowJsonTypeError({sec, key}, "an array of numbers");
  std::vector<double> result;
  result.reserve(v->size());
  for (std::size_t i = 0; i < v->size(); ++i) {
    const auto& item = v->at(i);
    if (!item.is_number()) ThrowJsonTypeError({sec, key, std::to_string(i)}, "a number");
    result.push_back(item.get<double>());
  }
  return result;
}

json ReadJsonFile(const fs::path& p) {
  std::ifstream s(p);
  if (!s) throw std::runtime_error("Cannot open JSON file: " + p.string());
  try { return json::parse(s); }
  catch (const json::parse_error& e) {
    throw std::runtime_error("JSON parse error in '" + p.string() + "': " + e.what());
  }
}

void WriteJsonFile(const fs::path& p, const json& payload) {
  std::ofstream s(p);
  if (!s) throw std::runtime_error("Cannot write JSON file: " + p.string());
  s << std::setw(2) << payload << '\n';
}

// ─── CLI ──────────────────────────────────────────────────────────────────────

struct CliOptions {
  fs::path configPath{"config/default_tgc.json"};
  fs::path outDir{"results"};
  std::optional<double> singleDistanceMm;
};

[[noreturn]] void PrintUsageAndExit(const char* prog, int code) {
  std::ostream& out = code == 0 ? std::cout : std::cerr;
  out << "Usage: " << prog << " [options]\n"
         "  --config <path>    JSON config file (default: config/default_tgc.json)\n"
         "  --out    <dir>     Output directory (default: results)\n"
         "  --distance <mm>    Run only this source distance (overrides config list)\n"
         "  --help             Show this message\n";
  std::exit(code);
}

CliOptions ParseCli(int argc, char* argv[]) {
  CliOptions opts;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--config") {
      if (i + 1 >= argc) PrintUsageAndExit(argv[0], 1);
      opts.configPath = argv[++i];
    } else if (arg == "--out") {
      if (i + 1 >= argc) PrintUsageAndExit(argv[0], 1);
      opts.outDir = argv[++i];
    } else if (arg == "--distance") {
      if (i + 1 >= argc) PrintUsageAndExit(argv[0], 1);
      opts.singleDistanceMm = std::stod(argv[++i]);
    } else if (arg == "--help" || arg == "-h") {
      PrintUsageAndExit(argv[0], 0);
    } else {
      throw std::runtime_error("Unknown argument: " + arg);
    }
  }
  return opts;
}

// ─── Config loading ───────────────────────────────────────────────────────────

Config LoadConfig(const fs::path& path) {
  if (!fs::exists(path))
    throw std::runtime_error("Configuration file not found: " + path.string());
  const json root = ReadJsonFile(path);
  if (!root.is_object()) ThrowJsonTypeError({"<root>"}, "a JSON object");

  Config cfg;

  if (const auto* g = FindSection(root, "geometry")) {
    cfg.geometry.wirePitchCm  = ReadDouble(*g, "geometry", "wire_pitch_cm",  cfg.geometry.wirePitchCm);
    cfg.geometry.wireDiamUm   = ReadDouble(*g, "geometry", "wire_diameter_um", cfg.geometry.wireDiamUm);
    cfg.geometry.gapCm        = ReadDouble(*g, "geometry", "gap_cm",         cfg.geometry.gapCm);
    cfg.geometry.nWires       = ReadInt   (*g, "geometry", "n_wires",         cfg.geometry.nWires);
    cfg.geometry.wireVoltageV = ReadDouble(*g, "geometry", "wire_voltage_V",  cfg.geometry.wireVoltageV);
  }

  if (const auto* r = FindSection(root, "readout")) {
    cfg.readout.type = ReadString(*r, "readout", "type", cfg.readout.type);
    cfg.readout.insulatorMaterial = ReadString(*r, "readout", "insulator_material",
                                               cfg.readout.insulatorMaterial);
    cfg.readout.insulatorThicknessUm    = ReadDouble(*r, "readout", "insulator_thickness_um",
                                                      cfg.readout.insulatorThicknessUm);
    cfg.readout.surfaceResistivityOhmSq = ReadDouble(*r, "readout", "surface_resistivity_ohm_sq",
                                                      cfg.readout.surfaceResistivityOhmSq);
    cfg.readout.enableDelayedSignal     = ReadBool  (*r, "readout", "enable_delayed_signal",
                                                      cfg.readout.enableDelayedSignal);
  }

  if (const auto* s = FindSection(root, "source")) {
    cfg.source.energyKeV    = ReadDouble(*s, "source", "energy_keV", cfg.source.energyKeV);
    {
      auto* dl = FindMember(*s, "source_distances_mm", {"source", "source_distances_mm"});
      if (dl && dl->is_array()) {
        cfg.source.fixedDistMm = std::vector<double>{};
        for (auto& v : *dl) cfg.source.fixedDistMm->push_back(v.get<double>());
      } else {
        cfg.source.fixedDistMm = std::nullopt;  // random per event
      }
    }
    auto* xpList = FindMember(*s, "x_positions_cm", {"source", "x_positions_cm"});
    if (xpList && xpList->is_array()) {
      cfg.source.fixedXCmList = std::vector<double>{};
      for (auto& v : *xpList)
        cfg.source.fixedXCmList->push_back(v.get<double>());
    } else {
      // backward compat: old scalar x_position_cm key
      auto* xp = FindMember(*s, "x_position_cm", {"source", "x_position_cm"});
      if (xp && xp->is_number())
        cfg.source.fixedXCmList = std::vector<double>{xp->get<double>()};
    }
  }

  if (const auto* g = FindSection(root, "gas")) {
    cfg.gas.gas1       = ReadString(*g, "gas", "gas1",               cfg.gas.gas1);
    cfg.gas.frac1      = ReadDouble(*g, "gas", "gas1_fraction_pct",  cfg.gas.frac1);
    cfg.gas.gas2       = ReadString(*g, "gas", "gas2",               cfg.gas.gas2);
    cfg.gas.ionSpecies = ReadString(*g, "gas", "ion_species",        cfg.gas.ionSpecies);
    cfg.gas.temperatureK  = ReadDouble(*g, "gas", "temperature_K",       cfg.gas.temperatureK);
    cfg.gas.pressureTorr  = ReadDouble(*g, "gas", "pressure_Torr",       cfg.gas.pressureTorr);
    cfg.gas.enablePenning = ReadBool  (*g, "gas", "enable_penning",       cfg.gas.enablePenning);
    cfg.gas.nCollisions         = ReadInt   (*g, "gas", "n_magboltz_collisions",  cfg.gas.nCollisions);
    cfg.gas.maxElectronEnergyEV = ReadDouble(*g, "gas", "max_electron_energy_eV", cfg.gas.maxElectronEnergyEV);
    cfg.gas.nFieldPoints        = ReadInt   (*g, "gas", "n_field_points",         cfg.gas.nFieldPoints);
    cfg.gas.eFieldMaxVcm        = ReadDouble(*g, "gas", "e_field_max_vcm",        cfg.gas.eFieldMaxVcm);
    cfg.gas.wValueEV            = ReadDouble(*g, "gas", "w_value_eV",             cfg.gas.wValueEV);
  }

  if (const auto* s = FindSection(root, "simulation")) {
    cfg.simulation.nEvents          = ReadSizeT (*s, "simulation", "n_events",          cfg.simulation.nEvents);
    cfg.simulation.maxAvalancheSize = ReadSizeT (*s, "simulation", "max_avalanche_size", cfg.simulation.maxAvalancheSize);
    cfg.simulation.timeWindowNs     = ReadDouble(*s, "simulation", "time_window_ns",     cfg.simulation.timeWindowNs);
    cfg.simulation.timeStepNs       = ReadDouble(*s, "simulation", "time_step_ns",       cfg.simulation.timeStepNs);
    cfg.simulation.enableIonDrift   = ReadBool  (*s, "simulation", "enable_ion_drift",   cfg.simulation.enableIonDrift);
    cfg.simulation.storeDriftLines  = ReadBool  (*s, "simulation", "store_drift_lines",  cfg.simulation.storeDriftLines);
  }

  if (cfg.readout.type != "conductive" && cfg.readout.type != "resistive")
    throw std::runtime_error("readout.type must be 'conductive' or 'resistive'");
  if (cfg.readout.type == "resistive") {
    if (cfg.readout.insulatorMaterial != "kapton" && cfg.readout.insulatorMaterial != "fr4")
      throw std::runtime_error("readout.insulator_material must be 'kapton' or 'fr4'");
    if (cfg.readout.insulatorThicknessUm <= 0.)
      throw std::runtime_error("readout.insulator_thickness_um must be positive");
    if (cfg.readout.surfaceResistivityOhmSq <= 0.)
      throw std::runtime_error("readout.surface_resistivity_ohm_sq must be positive");
  }

  if (cfg.geometry.wirePitchCm <= 0.)  throw std::runtime_error("geometry.wire_pitch_cm must be positive");
  if (cfg.geometry.wireDiamUm  <= 0.)  throw std::runtime_error("geometry.wire_diameter_um must be positive");
  if (cfg.geometry.gapCm       <= 0.)  throw std::runtime_error("geometry.gap_cm must be positive");
  if (cfg.geometry.nWires      <= 0)   throw std::runtime_error("geometry.n_wires must be positive");
  if (cfg.geometry.wireVoltageV<= 0.)  throw std::runtime_error("geometry.wire_voltage_V must be positive");
  if (cfg.source.fixedDistMm.has_value() && cfg.source.fixedDistMm->empty())
    throw std::runtime_error("source.source_distances_mm must not be empty when set");
  if (cfg.gas.frac1 <= 0. || cfg.gas.frac1 >= 100.)
    throw std::runtime_error("gas.gas1_fraction_pct must be in (0, 100)");
  if (cfg.gas.temperatureK     <= 0.)  throw std::runtime_error("gas.temperature_K must be positive");
  if (cfg.gas.pressureTorr     <= 0.)  throw std::runtime_error("gas.pressure_Torr must be positive");
  if (cfg.simulation.nEvents   == 0)   throw std::runtime_error("simulation.n_events must be at least 1");
  if (cfg.simulation.timeWindowNs <= 0.) throw std::runtime_error("simulation.time_window_ns must be positive");
  if (cfg.simulation.timeStepNs   <= 0.) throw std::runtime_error("simulation.time_step_ns must be positive");

  return cfg;
}

// ─── Gas setup ────────────────────────────────────────────────────────────────

static void ExportGasProps(MediumMagboltz& gas, const std::string& outPath,
                           const std::string& ionMobFile = "") {
  std::vector<double> efields, bfields, angles;
  gas.GetFieldGrid(efields, bfields, angles);

  std::ofstream f(outPath);
  if (!f) {
    std::cerr << "  Warning: could not write gas properties to " << outPath << "\n";
    return;
  }
  // Write ion mobility basename as a comment so the GUI can parse ion/gas names.
  if (!ionMobFile.empty()) {
    const auto sep = ionMobFile.find_last_of("/\\");
    const std::string base = (sep == std::string::npos)
                             ? ionMobFile : ionMobFile.substr(sep + 1);
    f << "# ion_mobility: " << base << "\n";
  }
  f << "e_field_Vcm,vd_cm_per_us,alpha_per_cm,eta_per_cm,"
       "dl_sqrtcm,dt_sqrtcm,v_ion_cm_per_us,mu_ion_cm2_per_Vus\n";

  for (double E : efields) {
    double vx = 0, vy = 0, vz = 0;
    gas.ElectronVelocity(E, 0, 0, 0, 0, 0, vx, vy, vz);
    double alpha = 0, eta = 0, dl = 0, dt = 0;
    gas.ElectronTownsend(E, 0, 0, 0, 0, 0, alpha);
    gas.ElectronAttachment(E, 0, 0, 0, 0, 0, eta);
    gas.ElectronDiffusion(E, 0, 0, 0, 0, 0, dl, dt);

    double v_ion = 0, mu_ion = 0;
    double vix = 0, viy = 0, viz = 0;
    if (gas.IonVelocity(E, 0, 0, 0, 0, 0, vix, viy, viz)) {
      v_ion  = std::abs(vix) * 1.e3;              // cm/ns → cm/μs
      mu_ion = (E > 0.) ? v_ion / E : 0.;         // cm/μs / (V/cm) = cm²/(V·μs)
    }

    f << std::scientific << std::setprecision(6)
      << E       << ","
      << vx * 1.e3 << ","   // cm/ns → cm/μs
      << alpha   << ","
      << eta     << ","
      << dl      << ","
      << dt      << ","
      << v_ion   << ","
      << mu_ion  << "\n";
  }
  std::cout << "  Gas properties exported to: " << outPath << "\n";
}

void SetupGas(MediumMagboltz& gas, const GasConfig& cfg) {
  gas.SetTemperature(cfg.temperatureK);
  gas.SetPressure(cfg.pressureTorr);

  const std::string gasFile = DeriveGasFileName(cfg);

  if (fs::exists(gasFile)) {
    std::cout << "  Loading gas table from: " << gasFile << "\n";
    gas.LoadGasFile(gasFile);
    // LoadGasFile restores the EFINAL ceiling stored in the file.  Override it
    // so the collision-rate table for AvalancheMicroscopic is pre-built to a
    // high enough energy to cover electrons near the wire (~160 kV/cm at 1900 V
    // can push electrons past 400 eV), preventing silent transport aborts.
    gas.SetMaxElectronEnergy(cfg.maxElectronEnergyEV);
  } else {
    std::cout << "  Gas file not found: " << gasFile << "\n"
              << "  Running Magboltz for " << cfg.nFieldPoints
              << " field points up to " << static_cast<int>(cfg.eFieldMaxVcm) << " V/cm"
              << " (first run: ~5 min for smoke grid, ~1-2 h for full grid)...\n";
    // Set energy ceiling before generation — electrons near the wire (>100 kV/cm
    // at 1900 V) reach energies well above Magboltz's default ~40 eV ceiling.
    gas.SetMaxElectronEnergy(cfg.maxElectronEnergyEV);
    // Logarithmically-spaced field grid from gentle drift region to near-wire avalanche.
    // At E>100 kV/cm the Magboltz SST/TOF tracks millions of avalanche electrons and
    // can take hours; reduce e_field_max_vcm / n_field_points for faster smoke runs.
    gas.SetFieldGrid(100., cfg.eFieldMaxVcm, cfg.nFieldPoints, /*logspacing=*/true);
    gas.GenerateGasTable(cfg.nCollisions, /*verbose=*/false);
    gas.WriteGasFile(gasFile);
    std::cout << "  Gas table saved to: " << gasFile << "\n";
  }

  if (cfg.enablePenning) {
    if (!gas.EnablePenningTransfer()) {
      std::cerr << "  Warning: Penning transfer could not be enabled.\n";
    } else {
      std::cout << "  Penning transfer enabled.\n";
    }
  }

  // Load ion mobility table for the configured ion species.
  // Garfield++ ships IonMobility_{X}+_{X}.txt for single-component gas X;
  // the species name is uppercased to match the filename convention.
  const char* garfieldInstall = std::getenv("GARFIELD_INSTALL");
  std::string loadedMob;   // basename recorded here for the props CSV comment
  if (garfieldInstall) {
    std::string ionUpper = cfg.ionSpecies;
    std::transform(ionUpper.begin(), ionUpper.end(), ionUpper.begin(), ::toupper);
    const std::string mobFile = std::string(garfieldInstall) +
                                "/share/Garfield/Data/IonMobility_"
                                + ionUpper + "+_" + ionUpper + ".txt";
    if (fs::exists(mobFile)) {
      gas.LoadIonMobility(mobFile);
      std::cout << "  " << ionUpper << "+ ion mobility loaded.\n";
      loadedMob = mobFile;
    } else {
      std::cerr << "  Warning: IonMobility_" << ionUpper << "+_" << ionUpper
                << ".txt not found at " << mobFile << "\n";
    }
  } else {
    std::cerr << "  Warning: GARFIELD_INSTALL not set; ion mobility not loaded.\n";
  }

  const std::string propsFile = gasFile.substr(0, gasFile.size() - 4) + "_props.csv";
  ExportGasProps(gas, propsFile, loadedMob);
}

// ─── Geometry and sensor setup ────────────────────────────────────────────────

void BuildGeometry(ComponentAnalyticField& cmp, MediumMagboltz& gas,
                   const GeometryConfig& geom) {
  const double wireDiamCm = geom.wireDiamUm * 1.e-4; // μm → cm

  cmp.SetMedium(&gas);

  // Cathode planes at ±gap_cm.
  // Bottom cathode (y = -gap) is the readout electrode labelled "cathode".
  // Top cathode (y = +gap) is ground only; not added as a readout electrode.
  cmp.AddPlaneY(-geom.gapCm, 0., "cathode");
  cmp.AddPlaneY(+geom.gapCm, 0., "cathode_top");

  // Anode wires at y = 0, centred at x = 0, all at +wireVoltageV.
  // All wires share the label "anode" so Sensor sums them as one channel.
  for (int i = 0; i < geom.nWires; ++i) {
    const double xw = (i - (geom.nWires - 1) / 2.) * geom.wirePitchCm;
    cmp.AddWire(xw, 0., wireDiamCm, geom.wireVoltageV, "anode");
  }
}

void SetupResistiveReadout(ComponentUser& cmp,
                           const ReadoutConfig& ro, const GeometryConfig& geom,
                           const SimulationConfig& sim) {
  const double epsR   = (ro.insulatorMaterial == "fr4") ? 4.6 : 3.5;
  const double dInsCm = ro.insulatorThicknessUm * 1e-4;        // μm → cm
  const double alpha  = epsR * geom.gapCm / (dInsCm + epsR * geom.gapCm);

  // τ = ε₀ ε_r ρ_s L² / (π² d_ins),  L = half wire-array width [m]
  const double eps0SI = 8.854e-12;                              // F/m
  const double dInsM  = ro.insulatorThicknessUm * 1e-6;        // μm → m
  const double L_m    = geom.nWires * geom.wirePitchCm * 0.5e-2; // cm → m
  const double tauNs  = eps0SI * epsR * ro.surfaceResistivityOhmSq
                        * (L_m * L_m) / (M_PI * M_PI * dInsM) * 1e9;

  const double gap = geom.gapCm;

  // Static weighting potential: W_s(y) = α (y + gap) / gap
  cmp.SetWeightingPotential(
      [alpha, gap](const double, const double y, const double) {
        return alpha * (y + gap) / gap;
      }, "cathode");

  // Static weighting field: E_w = −α/gap ŷ (uniform in y)
  cmp.SetWeightingField(
      [alpha, gap](const double, const double, const double,
                   double& wx, double& wy, double& wz) {
        wx = wz = 0.;
        wy = -alpha / gap;
      }, "cathode");

  if (ro.enableDelayedSignal) {
    // Delayed weighting potential: W_d(y,t) = W_s(y) × (exp(−t/τ) − 1)
    // Total W = W_s + W_d = W_s × exp(−t/τ) → decays to 0 as surface potential
    // at landing point is pulled to ground by the resistive sheet.
    cmp.SetDelayedWeightingPotential(
        [alpha, gap, tauNs](const double, const double y, const double,
                            const double t) {
          return alpha * (y + gap) / gap * (std::exp(-t / tauNs) - 1.0);
        }, "cathode");

    // Delayed weighting field: E_wd = −dW_d/dy = (α/gap)(1 − exp(−t/τ)) ŷ
    cmp.SetDelayedWeightingField(
        [alpha, gap, tauNs](const double, const double, const double,
                            const double t,
                            double& wx, double& wy, double& wz) {
          wx = wz = 0.;
          wy = -alpha / gap * (std::exp(-t / tauNs) - 1.0);
        }, "cathode");

    // Sample delayed signal over max(5τ, full time window), 200 points
    const double maxDelayNs = std::max(5.0 * tauNs, sim.timeWindowNs);
    constexpr std::size_t nDelayPts = 200;
    std::vector<double> dtimes(nDelayPts);
    for (std::size_t i = 0; i < nDelayPts; ++i)
      dtimes[i] = (i + 1) * maxDelayNs / nDelayPts;
    cmp.SetDelayedSignalTimes(dtimes);
  }

  std::cout << "  Resistive readout: α = " << alpha
            << ", τ = " << tauNs << " ns"
            << (ro.enableDelayedSignal ? "" : " (delayed signal disabled)")
            << "\n";
}

void SetupSensor(Sensor& sensor, ComponentAnalyticField& cmp,
                 ComponentUser* cmpReadout, const Config& cfg) {
  const auto& geom = cfg.geometry;
  const auto& sim  = cfg.simulation;

  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&cmp, "anode");       // all wires together
  // For resistive readout, use a ComponentUser with the dielectric-corrected
  // (and time-delayed) weighting potential for the cathode electrode.
  sensor.AddElectrode(cmpReadout ? static_cast<Garfield::Component*>(cmpReadout)
                                 : static_cast<Garfield::Component*>(&cmp),
                      "cathode");           // bottom cathode plane (readout)
  sensor.AddElectrode(&cmp, "cathode_top"); // top cathode plane (Ramo cross-check)
  if (cmpReadout && cfg.readout.enableDelayedSignal) sensor.EnableDelayedSignal();

  const std::size_t nBins =
      static_cast<std::size_t>(std::round(sim.timeWindowNs / sim.timeStepNs));
  sensor.SetTimeWindow(0., sim.timeStepNs, nBins);

  // Extend the active area one pitch beyond the outermost wire in x,
  // and 1 % of the gap beyond the cathode planes in y.
  const double xHalf    = ((geom.nWires - 1) / 2. + 1.) * geom.wirePitchCm;
  const double yMargin  = 0.01 * geom.gapCm;
  sensor.SetArea(-xHalf, -geom.gapCm - yMargin, -0.5,
                  xHalf,  geom.gapCm + yMargin,  0.5);
}

// ─── Per-distance simulation loop ─────────────────────────────────────────────

DistanceSummary RunDistancePoint(const Config& cfg,
                                 std::optional<double> distOptMm,
                                 Sensor& sensor,
                                 TDirectory* distDir,
                                 std::optional<double> fixedXCm = std::nullopt) {
  const auto& geom = cfg.geometry;
  const auto& sim  = cfg.simulation;

  // Half-span of the wire array for random x sampling
  const double xHalfWires = (geom.nWires - 1) / 2. * geom.wirePitchCm;
  // Half-gap for random distance sampling (mm)
  const double gapHalfMm  = geom.gapCm * 10.0 / 2.0;

  const std::size_t nBins =
      static_cast<std::size_t>(std::round(sim.timeWindowNs / sim.timeStepNs));

  // ── Histograms ──────────────────────────────────────────────────────────────
  TH1D hAnodeQ("h_anode_charge",
               "Induced charge on anode;Q_{anode} [fC];Events", 200, 0., 0.);
  TH1D hCathodeQ("h_cathode_charge",
                 "Induced charge on cathode;Q_{cathode} [fC];Events", 200, 0., 0.);
  TH1D hRatio("h_ratio_charge",
              "Charge ratio;Q_{cathode}/Q_{anode};Events", 100, 0., 2.);
  TH1D hNprimary("h_n_primary_electrons",
                 "Primary electrons per event;N_{e,primary};Events", 400, -0.5, 399.5);
  TH1D hAvalSize("h_avalanche_size",
                 "Total avalanche size;N_{e,total};Events", 200, 0., 0.);

  TH1D hCathodeTopQ("h_cathode_top_charge",
                   "Induced charge on cathode_top;Q_{cathode\_top} [fC];Events", 200, 0., 0.);

  TProfile pAnodeSignal("p_anode_signal",
                        "Mean anode signal;t [ns];#LTi_{anode}#GT [fC/ns]",
                        static_cast<int>(nBins), 0., sim.timeWindowNs);
  TProfile pCathodeSignal("p_cathode_signal",
                          "Mean cathode signal;t [ns];#LTi_{cathode}#GT [fC/ns]",
                          static_cast<int>(nBins), 0., sim.timeWindowNs);
  TProfile pCathodeTopSignal("p_cathode_top_signal",
                             "Mean cathode_top signal;t [ns];#LTi_{cathode\_top}#GT [fC/ns]",
                             static_cast<int>(nBins), 0., sim.timeWindowNs);

  for (TH1* h : std::initializer_list<TH1*>{
           &hAnodeQ, &hCathodeQ, &hCathodeTopQ, &hRatio,
           &hNprimary, &hAvalSize,
           &pAnodeSignal, &pCathodeSignal, &pCathodeTopSignal}) {
    h->SetDirectory(nullptr);
  }

  // ── Per-event signal tree ────────────────────────────────────────────────────
  TTree signalTree("t_signals", "Per-event signal waveforms");
  signalTree.SetDirectory(nullptr);
  std::vector<float> anodeSig(nBins, 0.f), cathodeSig(nBins, 0.f);
  int   evtId = 0;
  float evtQa = 0.f, evtQc = 0.f;
  signalTree.Branch("event",             &evtId, "event/I");
  signalTree.Branch("anode_charge_fC",   &evtQa, "anode_charge_fC/F");
  signalTree.Branch("cathode_charge_fC", &evtQc, "cathode_charge_fC/F");
  signalTree.Branch("anode",   &anodeSig);
  signalTree.Branch("cathode", &cathodeSig);

  // ── 3D track branches ────────────────────────────────────────────────────────
  // primary_x/y/z : points along the primary electron drift line (≥ 2 always;
  //                 many intermediate steps when store_drift_lines=true)
  // cloud_x/y/z   : start positions of secondary electron tracks (avalanche cloud)
  // ion_x/y/z     : ion drift paths for up to kMaxDispIonPaths ions, flattened
  // ion_npts       : number of drift-line points per stored ion path
  std::vector<float> primaryX, primaryY, primaryZ;
  std::vector<float> cloudX,   cloudY,   cloudZ;
  std::vector<float> ionX,     ionY,     ionZ;
  std::vector<int>   ionNpts;
  signalTree.Branch("primary_x", &primaryX);
  signalTree.Branch("primary_y", &primaryY);
  signalTree.Branch("primary_z", &primaryZ);
  signalTree.Branch("cloud_x",   &cloudX);
  signalTree.Branch("cloud_y",   &cloudY);
  signalTree.Branch("cloud_z",   &cloudZ);
  signalTree.Branch("ion_x",     &ionX);
  signalTree.Branch("ion_y",     &ionY);
  signalTree.Branch("ion_z",     &ionZ);
  signalTree.Branch("ion_npts",  &ionNpts);

  // ── Transport objects ────────────────────────────────────────────────────────
  // Average primary electrons: N = E_photon / W-value (e.g. 5900 eV / 26 eV ≈ 227)
  const int nPrimary = std::max(1,
      static_cast<int>(std::round(cfg.source.energyKeV * 1.e3 / cfg.gas.wValueEV)));

  AvalancheMicroscopic aval(&sensor);
  if (sim.maxAvalancheSize > 0) aval.EnableAvalancheSizeLimit(sim.maxAvalancheSize);
  // Enables storage of intermediate collision steps along each drift line.
  // Without this, GetNumberOfElectronDriftLinePoints returns 2 (start + end only).
  if (sim.storeDriftLines) aval.EnableDriftLines(true);

  // Ion drift: DriftLineRKF transports each positive ion from its creation
  // position and adds the Ramo-theorem induced current to the sensor.
  // Constructed only when ion transport is enabled to avoid needless overhead.
  std::optional<DriftLineRKF> ionDrift;
  if (sim.enableIonDrift) ionDrift.emplace(&sensor);

  // ── Accumulators for summary statistics ──────────────────────────────────────
  std::vector<double> anodeCharges, cathodeCharges, cathodeTopCharges, chargeRatios;
  std::vector<double> primaryCounts, avalancheSizes;
  anodeCharges.reserve(sim.nEvents);
  cathodeCharges.reserve(sim.nEvents);
  cathodeTopCharges.reserve(sim.nEvents);

  std::size_t nInteracted = 0;
  const std::size_t progressStep = std::max<std::size_t>(1, sim.nEvents / 10);
  const std::string distLabel = distOptMm.has_value()
      ? FormatNumber(*distOptMm) + " mm" : "random";

  // ── Event loop ───────────────────────────────────────────────────────────────
  for (std::size_t ev = 0; ev < sim.nEvents; ++ev) {
    sensor.ClearSignal();

    const double x0 = fixedXCm.has_value()
                          ? *fixedXCm
                          : gRandom->Uniform(-xHalfWires, xHalfWires);

    const double distMm   = distOptMm.has_value()
                                ? *distOptMm
                                : gRandom->Uniform(-gapHalfMm, gapHalfMm);
    const double sourceYCm = -distMm * 0.1;  // mm → cm, positive distance → readout side (y < 0)
    const double y0 = std::max(-geom.gapCm + 1.e-4,
                                std::min(geom.gapCm - 1.e-4, sourceYCm));

    // Transport one representative electron from (x0, y0, 0) and scale all results
    // by nPrimary.  This is exact for the mean Q_cathode/Q_anode ratio (the key
    // observable): all electrons start at the same position, and the Shockley-Ramo
    // weighting is linear in charge, so scaling by nPrimary gives the correct
    // expected total charge and ratio.  Running nPrimary independent avalanches
    // would only add ~1/sqrt(nPrimary) statistical noise on top of the much larger
    // single-avalanche (Polya) fluctuation — not worth the ~227× runtime cost.
    ++nInteracted;
    hNprimary.Fill(nPrimary);
    primaryCounts.push_back(static_cast<double>(nPrimary));

    aval.AvalancheElectron(x0, y0, 0., 0., 0.1); // 0.1 eV ≈ thermal
    int ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);
    const int totalAvalElectrons = ne * nPrimary;

    hAvalSize.Fill(static_cast<double>(totalAvalElectrons));
    avalancheSizes.push_back(static_cast<double>(totalAvalElectrons));

    // ── 3D track data ────────────────────────────────────────────────────────
    // nEp is hoisted here so both the cloud loop and the ion-drift loop reuse it.
    const std::size_t nEp = aval.GetNumberOfElectronEndpoints();

    // Primary electron drift line (track 0).
    // GetNumberOfElectronDriftLinePoints(0) returns ≥ 2 (start + end) always;
    // with storeDriftLines=true it returns all intermediate collision steps.
    primaryX.clear(); primaryY.clear(); primaryZ.clear();
    {
      const std::size_t nPts = aval.GetNumberOfElectronDriftLinePoints(0);
      primaryX.reserve(nPts); primaryY.reserve(nPts); primaryZ.reserve(nPts);
      for (std::size_t ip = 0; ip < nPts; ++ip) {
        double px, py, pz, pt;
        aval.GetElectronDriftLinePoint(px, py, pz, pt, ip, /*track=*/0);
        primaryX.push_back(static_cast<float>(px));
        primaryY.push_back(static_cast<float>(py));
        primaryZ.push_back(static_cast<float>(pz));
      }
    }

    // Avalanche cloud: start positions of secondary electron tracks (tracks 1…nEp-1).
    // These cluster near the wire and visualise the avalanche extent.
    cloudX.clear(); cloudY.clear(); cloudZ.clear();
    {
      const std::size_t nSec   = nEp > 0 ? nEp - 1 : 0;
      const std::size_t stride = (nSec > kMaxDispCloudPts && kMaxDispCloudPts > 0)
                                  ? nSec / kMaxDispCloudPts : 1;
      cloudX.reserve(std::min(nSec, kMaxDispCloudPts));
      for (std::size_t i = 1; i < nEp; i += stride) {
        double x0c, y0c, z0c, t0c, e0c, x1c, y1c, z1c, t1c, e1c; int stc;
        aval.GetElectronEndpoint(i, x0c, y0c, z0c, t0c, e0c,
                                    x1c, y1c, z1c, t1c, e1c, stc);
        cloudX.push_back(static_cast<float>(x0c));
        cloudY.push_back(static_cast<float>(y0c));
        cloudZ.push_back(static_cast<float>(z0c));
      }
    }

    // Drift every positive ion from the position where it was created.
    // AvalancheMicroscopic records each electron track; the start point of
    // track 0 is the primary photoionisation ion, and the start points of
    // tracks 1…n are the positions of the avalanche ions from ionising
    // collisions near the wire.  DriftLineRKF::DriftIon() transports each ion
    // to the cathode and adds its Ramo-theorem induced current to the sensor;
    // contributions outside the sensor time window are clipped automatically.
    // For the first kMaxDispIonPaths ions the full drift-line path is also
    // extracted immediately after each DriftIon call (before the next call
    // overwrites DriftLineRKF's internal state).
    ionX.clear(); ionY.clear(); ionZ.clear(); ionNpts.clear();
    if (sim.enableIonDrift) {
      for (std::size_t i = 0; i < nEp; ++i) {
        double xi0, yi0, zi0, ti0, ei0;
        double xi1, yi1, zi1, ti1, ei1;
        int st;
        aval.GetElectronEndpoint(i, xi0, yi0, zi0, ti0, ei0,
                                    xi1, yi1, zi1, ti1, ei1, st);
        ionDrift->DriftIon(xi0, yi0, zi0, ti0);

        if (i < kMaxDispIonPaths) {
          const std::size_t nPts = ionDrift->GetNumberOfDriftLinePoints();
          ionNpts.push_back(static_cast<int>(nPts));
          for (std::size_t ip = 0; ip < nPts; ++ip) {
            double ix, iy, iz, it;
            ionDrift->GetDriftLinePoint(ip, ix, iy, iz, it);
            ionX.push_back(static_cast<float>(ix));
            ionY.push_back(static_cast<float>(iy));
            ionZ.push_back(static_cast<float>(iz));
          }
        }
      }
    }

    // Integrate the binned induced-current signal to obtain charge.
    // GetSignal returns the induced current in fC/ns; multiplying by the bin
    // width (ns) and summing gives charge in fC.  Scale by nPrimary.
    //
    // Sign convention (Shockley-Ramo):
    //   anode   : electrons moving toward wire + ions moving away from wire
    //             → both give net integral NEGATIVE in Garfield++
    //             → negate to obtain the conventionally positive collected charge
    //   cathode : dominated by ions drifting toward the readout pad
    //             → net integral is POSITIVE
    //
    // NOTE: Sensor::GetInducedCharge() uses a separate per-electrode "charge"
    // accumulator that AvalancheMicroscopic does not populate (it only calls
    // AddSignal into the time-binned arrays), so GetInducedCharge always returns
    // zero here.  The manual integral below is the correct approach.
    double rawAnode = 0., rawCathode = 0., rawCathodeTop = 0.;
    for (std::size_t k = 0; k < nBins; ++k) {
      const double sigA = sensor.GetSignal("anode",       k);
      const double sigC = sensor.GetSignal("cathode",     k);
      const double sigT = sensor.GetSignal("cathode_top", k);
      rawAnode      += sigA;
      rawCathode    += sigC;
      rawCathodeTop += sigT;
      const double t = (static_cast<double>(k) + 0.5) * sim.timeStepNs;
      pAnodeSignal.Fill(t,      sigA * nPrimary);
      pCathodeSignal.Fill(t,    sigC * nPrimary);
      pCathodeTopSignal.Fill(t, sigT * nPrimary);
      anodeSig[k]   = static_cast<float>(sigA * nPrimary);
      cathodeSig[k] = static_cast<float>(sigC * nPrimary);
    }

    const double qAnode      = -rawAnode      * sim.timeStepNs * nPrimary; // [fC]
    const double qCathode    =  rawCathode    * sim.timeStepNs * nPrimary; // [fC]
    const double qCathodeTop =  rawCathodeTop * sim.timeStepNs * nPrimary; // [fC]

    evtId = static_cast<int>(ev);
    evtQa = static_cast<float>(qAnode);
    evtQc = static_cast<float>(qCathode);
    signalTree.Fill();

    hAnodeQ.Fill(qAnode);
    hCathodeQ.Fill(qCathode);
    hCathodeTopQ.Fill(qCathodeTop);
    anodeCharges.push_back(qAnode);
    cathodeCharges.push_back(qCathode);
    cathodeTopCharges.push_back(qCathodeTop);

    if (qAnode > 0.) {
      const double ratio = qCathode / qAnode;
      hRatio.Fill(ratio);
      chargeRatios.push_back(ratio);
    }

    if ((ev + 1) % progressStep == 0 || ev + 1 == sim.nEvents) {
      std::cout << "  dist=" << distLabel << ": "
                << (ev + 1) << "/" << sim.nEvents
                << " events processed\n";
    }
  }

  // ── Write histograms ─────────────────────────────────────────────────────────
  if (distDir) {
    distDir->cd();
    hAnodeQ.Write("h_anode_charge");
    hCathodeQ.Write("h_cathode_charge");
    hCathodeTopQ.Write("h_cathode_top_charge");
    hRatio.Write("h_ratio_charge");
    hNprimary.Write("h_n_primary_electrons");
    hAvalSize.Write("h_avalanche_size");
    pAnodeSignal.Write("p_anode_signal");
    pCathodeSignal.Write("p_cathode_signal");
    pCathodeTopSignal.Write("p_cathode_top_signal");
    signalTree.Write("t_signals");
  }

  // ── Build summary ─────────────────────────────────────────────────────────────
  DistanceSummary s;
  s.distanceMm          = distOptMm;  // nullopt when random per-event
  s.nEvents             = sim.nEvents;
  s.nInteracted         = nInteracted;
  s.interactionFraction = sim.nEvents > 0
                              ? static_cast<double>(nInteracted) / static_cast<double>(sim.nEvents)
                              : 0.;

  s.meanAnodeChargeFC   = Mean(anodeCharges);
  s.rmsAnodeChargeFC    = Rms(anodeCharges, s.meanAnodeChargeFC);
  s.semAnodeChargeFC    = Sem(s.rmsAnodeChargeFC, anodeCharges.size());

  s.meanCathodeChargeFC = Mean(cathodeCharges);
  s.rmsCathodeChargeFC  = Rms(cathodeCharges, s.meanCathodeChargeFC);
  s.semCathodeChargeFC  = Sem(s.rmsCathodeChargeFC, cathodeCharges.size());

  s.meanCathodeTopChargeFC = Mean(cathodeTopCharges);
  s.rmsCathodeTopChargeFC  = Rms(cathodeTopCharges, s.meanCathodeTopChargeFC);
  s.semCathodeTopChargeFC  = Sem(s.rmsCathodeTopChargeFC, cathodeTopCharges.size());

  s.meanChargeRatio     = Mean(chargeRatios);
  s.rmsChargeRatio      = Rms(chargeRatios, s.meanChargeRatio);
  s.semChargeRatio      = Sem(s.rmsChargeRatio, chargeRatios.size());

  s.meanPrimaryElectrons = Mean(primaryCounts);
  s.meanAvalancheSize    = Mean(avalancheSizes);

  return s;
}

// ─── Summary graphs ───────────────────────────────────────────────────────────

void WriteSummaryGraphs(const std::vector<DistanceSummary>& sums,
                        TDirectory* summaryDir, const fs::path& pngPath) {
  if (sums.empty()) return;
  const std::size_t n = sums.size();
  std::vector<double> x(n), xe(n, 0.);
  std::vector<double> qa(n), qaE(n), qc(n), qcE(n), qt(n), qtE(n), rat(n), ratE(n);

  for (std::size_t i = 0; i < n; ++i) {
    x[i]    = sums[i].distanceMm.value_or(static_cast<double>(i));
    qa[i]   = sums[i].meanAnodeChargeFC;      qaE[i]  = sums[i].semAnodeChargeFC;
    qc[i]   = sums[i].meanCathodeChargeFC;    qcE[i]  = sums[i].semCathodeChargeFC;
    qt[i]   = sums[i].meanCathodeTopChargeFC; qtE[i]  = sums[i].semCathodeTopChargeFC;
    rat[i]  = sums[i].meanChargeRatio;        ratE[i] = sums[i].semChargeRatio;
  }

  auto MakeGraph = [&](const char* name, const char* title,
                       const std::vector<double>& y, const std::vector<double>& ye,
                       int marker) {
    TGraphErrors g(static_cast<int>(n), x.data(), y.data(), xe.data(), ye.data());
    g.SetName(name);
    g.SetTitle(title);
    g.SetMarkerStyle(marker);
    g.SetLineWidth(2);
    if (summaryDir) { summaryDir->cd(); g.Write(); }
    return g;
  };

  auto gAnode      = MakeGraph("g_anode_charge",
    "Mean anode charge;Source distance from wire plane [mm];Q_{anode} [fC]", qa, qaE, 20);
  auto gCathode    = MakeGraph("g_cathode_charge",
    "Mean cathode charge;Source distance from wire plane [mm];Q_{cathode} [fC]", qc, qcE, 21);
  auto gCathodeTop = MakeGraph("g_cathode_top_charge",
    "Mean cathode_top charge;Source distance from wire plane [mm];Q_{cathode\_top} [fC]",
    qt, qtE, 22);
  auto gRatio      = MakeGraph("g_charge_ratio",
    "Charge ratio;Source distance from wire plane [mm];Q_{cathode}/Q_{anode}", rat, ratE, 23);

  TCanvas canvas("c_tgc_summary", "TGC summary", 2400, 500);
  canvas.Divide(4, 1);
  canvas.cd(1); gAnode.Draw("APL");
  canvas.cd(2); gCathode.Draw("APL");
  canvas.cd(3); gCathodeTop.Draw("APL");
  canvas.cd(4); gRatio.Draw("APL");
  EnsureDirectory(pngPath.parent_path());
  canvas.SaveAs(pngPath.string().c_str());
}

// ─── CSV summary ─────────────────────────────────────────────────────────────

void WriteSummaryCsv(const fs::path& path, const std::vector<DistanceSummary>& sums) {
  std::ofstream f(path);
  if (!f) throw std::runtime_error("Cannot write CSV: " + path.string());

  f << "source_distance_mm,x_position_cm,n_events,n_interacted,interaction_fraction,"
       "mean_anode_charge_fC,rms_anode_charge_fC,sem_anode_charge_fC,"
       "mean_cathode_charge_fC,rms_cathode_charge_fC,sem_cathode_charge_fC,"
       "mean_cathode_top_charge_fC,rms_cathode_top_charge_fC,sem_cathode_top_charge_fC,"
       "mean_charge_ratio,rms_charge_ratio,sem_charge_ratio,"
       "mean_primary_electrons,mean_avalanche_size\n";

  f << std::fixed << std::setprecision(6);
  for (const auto& s : sums) {
    if (s.distanceMm) f << *s.distanceMm; else f << "random";
    f << ',';
    if (s.xPositionCm.has_value()) f << *s.xPositionCm;
    f << ','
      << s.nEvents                   << ','
      << s.nInteracted               << ','
      << s.interactionFraction       << ','
      << s.meanAnodeChargeFC         << ','
      << s.rmsAnodeChargeFC          << ','
      << s.semAnodeChargeFC          << ','
      << s.meanCathodeChargeFC       << ','
      << s.rmsCathodeChargeFC        << ','
      << s.semCathodeChargeFC        << ','
      << s.meanCathodeTopChargeFC    << ','
      << s.rmsCathodeTopChargeFC     << ','
      << s.semCathodeTopChargeFC     << ','
      << s.meanChargeRatio           << ','
      << s.rmsChargeRatio            << ','
      << s.semChargeRatio            << ','
      << s.meanPrimaryElectrons      << ','
      << s.meanAvalancheSize         << '\n';
  }
}

// ─── Config echo ──────────────────────────────────────────────────────────────

json ConfigToJson(const Config& cfg) {
  json jSrc = {
    {"energy_keV", cfg.source.energyKeV}
  };
  jSrc["source_distances_mm"] = cfg.source.fixedDistMm.has_value()
                                     ? json(*cfg.source.fixedDistMm)
                                     : json(nullptr);
  jSrc["x_positions_cm"] = cfg.source.fixedXCmList.has_value()
                               ? json(*cfg.source.fixedXCmList)
                               : json(nullptr);
  return {
    {"geometry", {
      {"wire_pitch_cm",    cfg.geometry.wirePitchCm},
      {"wire_diameter_um", cfg.geometry.wireDiamUm},
      {"gap_cm",           cfg.geometry.gapCm},
      {"n_wires",          cfg.geometry.nWires},
      {"wire_voltage_V",   cfg.geometry.wireVoltageV}
    }},
    {"readout", {
      {"type",                       cfg.readout.type},
      {"insulator_material",         cfg.readout.insulatorMaterial},
      {"insulator_thickness_um",     cfg.readout.insulatorThicknessUm},
      {"surface_resistivity_ohm_sq", cfg.readout.surfaceResistivityOhmSq},
      {"enable_delayed_signal",      cfg.readout.enableDelayedSignal}
    }},
    {"source", jSrc},
    {"gas", {
      {"gas1",                   cfg.gas.gas1},
      {"gas1_fraction_pct",      cfg.gas.frac1},
      {"gas2",                   cfg.gas.gas2},
      {"ion_species",            cfg.gas.ionSpecies},
      {"temperature_K",          cfg.gas.temperatureK},
      {"pressure_Torr",          cfg.gas.pressureTorr},
      {"enable_penning",         cfg.gas.enablePenning},
      {"n_magboltz_collisions",  cfg.gas.nCollisions},
      {"max_electron_energy_eV", cfg.gas.maxElectronEnergyEV},
      {"n_field_points",         cfg.gas.nFieldPoints},
      {"e_field_max_vcm",        cfg.gas.eFieldMaxVcm},
      {"w_value_eV",             cfg.gas.wValueEV}
    }},
    {"simulation", {
      {"n_events",           cfg.simulation.nEvents},
      {"max_avalanche_size", cfg.simulation.maxAvalancheSize},
      {"time_window_ns",     cfg.simulation.timeWindowNs},
      {"time_step_ns",       cfg.simulation.timeStepNs},
      {"enable_ion_drift",   cfg.simulation.enableIonDrift},
      {"store_drift_lines",  cfg.simulation.storeDriftLines}
    }}
  };
}

std::string BuildRunFolderName(const Config& cfg) {
  std::ostringstream ss;
  ss << "V" << static_cast<int>(cfg.geometry.wireVoltageV)
     << "V__n" << cfg.simulation.nEvents;
  return ss.str();
}

} // namespace

// ─── main ─────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
  try {
    // Flush cout after every write so the GUI log panel streams in real time
    // even when stdout is a pipe (pipes switch cout to full buffering by default).
    std::cout << std::unitbuf;

    gROOT->SetBatch(true);
    gStyle->SetOptStat(1110);
    TH1::AddDirectory(false);
    TH1::StatOverflows(true);
    gRandom->SetSeed(0);

    const auto opts = ParseCli(argc, argv);
    Config cfg = LoadConfig(opts.configPath);

    if (opts.singleDistanceMm)
      cfg.source.fixedDistMm = std::vector<double>{*opts.singleDistanceMm};

    const fs::path runDir = opts.outDir / BuildRunFolderName(cfg);
    EnsureDirectory(runDir);

    std::cout << "TGC Garfield++ simulation\n"
              << "  config  : " << opts.configPath << "\n"
              << "  output  : " << runDir << "\n"
              << "  geometry: " << cfg.geometry.nWires << " wires, "
              << cfg.geometry.wirePitchCm * 10. << " mm pitch, "
              << cfg.geometry.wireDiamUm << " μm diameter, "
              << cfg.geometry.gapCm * 10. << " mm gap"
              << ", E_peak ≈ " << static_cast<int>(ComputePeakFieldVcm(cfg.geometry) / 1000.) << " kV/cm\n"
              << "  voltage : " << cfg.geometry.wireVoltageV << " V (wires), 0 V (cathodes)\n"
              << "  readout : " << cfg.readout.type
              << (cfg.readout.type == "resistive"
                    ? (" (" + cfg.readout.insulatorMaterial
                       + ", " + std::to_string(static_cast<int>(cfg.readout.insulatorThicknessUm)) + " μm"
                       + ", " + std::to_string(static_cast<int>(cfg.readout.surfaceResistivityOhmSq / 1000)) + " kΩ/sq"
                       + (cfg.readout.enableDelayedSignal ? ", delayed signal on" : ", delayed signal off") + ")")
                    : std::string())
              << "\n"
              << "  gas     : " << cfg.gas.gas1 << ":" << cfg.gas.gas2
              << " " << static_cast<int>(cfg.gas.frac1) << ":"
              << static_cast<int>(100. - cfg.gas.frac1) << ", "
              << cfg.gas.temperatureK << " K, "
              << cfg.gas.pressureTorr << " Torr\n"
              << "  source  : " << cfg.source.energyKeV << " keV, "
              << (cfg.source.fixedDistMm.has_value()
                    ? std::to_string(cfg.source.fixedDistMm->size()) + " distance point(s)"
                    : "random distance")
              << (cfg.source.fixedXCmList.has_value()
                    ? (", " + std::to_string(cfg.source.fixedXCmList->size()) + " x-position(s)")
                    : std::string(", x random"))
              << "\n"
              << "  events  : " << cfg.simulation.nEvents << " per point\n";

    // Sanity-check gas.e_field_max_vcm against the estimated peak near-wire field.
    {
      const double ePeakVcm = ComputePeakFieldVcm(cfg.geometry);
      if (cfg.gas.eFieldMaxVcm < ePeakVcm) {
        std::cerr << "  WARNING: gas.e_field_max_vcm ("
                  << static_cast<int>(cfg.gas.eFieldMaxVcm / 1000.) << " kV/cm) is below"
                  << " the estimated peak near-wire field ("
                  << static_cast<int>(ePeakVcm / 1000.) << " kV/cm).\n"
                  << "  Magboltz will extrapolate — increase e_field_max_vcm.\n";
      } else if (cfg.gas.eFieldMaxVcm < 1.5 * ePeakVcm) {
        std::cout << "  Note: e_field_max_vcm is only "
                  << std::fixed << std::setprecision(1)
                  << cfg.gas.eFieldMaxVcm / ePeakVcm
                  << "× E_peak — recommend ≥ 1.5×.\n"
                  << std::defaultfloat;
      }
    }

    // Gas is shared across all distance points
    std::cout << "\nSetting up gas...\n";
    MediumMagboltz gas(cfg.gas.gas1, cfg.gas.frac1,
                       cfg.gas.gas2, 100. - cfg.gas.frac1);
    SetupGas(gas, cfg.gas);

    // Geometry and sensor are shared across all distance points
    ComponentUser cmpReadout;
    if (cfg.readout.type == "resistive") {
      std::cout << "\nSetting up resistive readout...\n";
      SetupResistiveReadout(cmpReadout, cfg.readout, cfg.geometry, cfg.simulation);
    }

    ComponentAnalyticField cmp;
    BuildGeometry(cmp, gas, cfg.geometry);

    Sensor sensor;
    SetupSensor(sensor, cmp,
                cfg.readout.type == "resistive" ? &cmpReadout : nullptr,
                cfg);

    // ROOT output
    TFile rootFile((runDir / "tgc_sim.root").string().c_str(), "RECREATE");
    if (rootFile.IsZombie())
      throw std::runtime_error("Failed to create ROOT file in " + runDir.string());

    TDirectory* summaryDir = rootFile.mkdir("summary");

    std::vector<DistanceSummary> allSummaries;

    // Build flat list of (optional) distances to iterate:
    //   nullopt  → pick uniform random y per event  (directory name: "dist_rnd[_xYYmm]")
    //   value    → fixed distance for all events in that directory
    std::vector<std::optional<double>> dList;
    if (cfg.source.fixedDistMm.has_value() && !cfg.source.fixedDistMm->empty()) {
      for (double d : *cfg.source.fixedDistMm)
        dList.push_back(d);
    } else {
      dList.push_back(std::nullopt);  // random per event
    }

    // Build flat list of (optional) x-positions to iterate:
    //   nullopt  → pick uniform random x per event  (no _x suffix in directory name)
    //   value    → fixed x for all events in that directory
    std::vector<std::optional<double>> xList;
    if (cfg.source.fixedXCmList.has_value() && !cfg.source.fixedXCmList->empty()) {
      for (double x : *cfg.source.fixedXCmList)
        xList.push_back(x);
    } else {
      xList.push_back(std::nullopt);  // random per event
    }

    for (const auto& dOpt : dList) {
      for (const auto& xOpt : xList) {
        std::string label = dOpt.has_value() ? FormatNumber(*dOpt) + " mm" : "random";
        if (xOpt.has_value())
          label += "  x=" + FormatNumber(*xOpt * 10.0) + " mm";
        std::cout << "\n--- Source distance: " << label << " ---\n";

        std::string tag = dOpt.has_value()
            ? "dist_" + FileSafeNumber(*dOpt) + "mm"
            : "dist_rnd";
        if (xOpt.has_value())
          tag += "_x" + FileSafeNumber(*xOpt * 10.0) + "mm";  // cm → mm for name
        TDirectory* distDir = rootFile.mkdir(tag.c_str());
        if (!distDir) throw std::runtime_error("Failed to create ROOT dir: " + tag);

        DistanceSummary summary = RunDistancePoint(cfg, dOpt, sensor, distDir, xOpt);
        summary.xPositionCm = xOpt;
        allSummaries.push_back(summary);

        std::cout << "  ⟨Q_anode⟩       = " << FormatNumber(summary.meanAnodeChargeFC)      << " fC"
                  << "  ±" << FormatNumber(summary.semAnodeChargeFC)      << " (SEM)\n"
                  << "  ⟨Q_cathode⟩     = " << FormatNumber(summary.meanCathodeChargeFC)   << " fC"
                  << "  ±" << FormatNumber(summary.semCathodeChargeFC)    << " (SEM)\n"
                  << "  ⟨Q_cathode_top⟩ = " << FormatNumber(summary.meanCathodeTopChargeFC) << " fC"
                  << "  ±" << FormatNumber(summary.semCathodeTopChargeFC) << " (SEM)\n"
                  << "  Ramo check: |Q_anode| = "
                  << FormatNumber(std::abs(summary.meanAnodeChargeFC))
                  << "  |Q_cathode|+|Q_top| = "
                  << FormatNumber(std::abs(summary.meanCathodeChargeFC)
                                  + std::abs(summary.meanCathodeTopChargeFC)) << " fC\n"
                  << "  ⟨ratio⟩         = " << FormatNumber(summary.meanChargeRatio) << "\n"
                  << "  interaction fraction: "
                  << FormatNumber(summary.interactionFraction * 100., 2) << "%\n"
                  << "  ⟨avalanche size⟩: "
                  << FormatNumber(summary.meanAvalancheSize, 0) << " electrons\n";
      }
    }

    WriteSummaryGraphs(allSummaries, summaryDir, runDir / "summary" / "tgc_summary.png");
    rootFile.Write();
    rootFile.Close();

    WriteSummaryCsv(runDir / "summary.csv", allSummaries);
    WriteJsonFile(runDir / "run_config.json", ConfigToJson(cfg));

    std::cout << "\nDone. Results written to " << runDir << "\n";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
