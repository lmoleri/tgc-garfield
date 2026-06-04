// Garfield++ TGC (Thin Gap Chamber) simulation
//
// Geometry : 10 anode wires (50 μm diameter, 1.8 mm pitch) at y = 0,
//            two grounded cathode planes at y = ±gap_cm.
//            y > 0 is the non-readout ("cathode_top") side;
//            y < 0 is the readout pad ("cathode") side.
// Gas      : Ar:CO2 70:30, 1 atm, 293.15 K.
// Source   : N = E/W primary electrons placed at a configurable distance from
//            the wire plane (signed: positive → cathode_top side).
//            N is computed from the photon energy and the gas W-value.
// Readout  : (1) all wires as a single "anode" channel,
//            (2) the bottom cathode plane as a "cathode" channel.
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
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "nlohmann/json.hpp"

namespace fs = std::filesystem;

namespace {

using Garfield::AvalancheMicroscopic;
using Garfield::ComponentAnalyticField;
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

struct SourceConfig {
  double energyKeV                     = 5.9;
  std::vector<double> distancesMm      = {0.2, 0.5, 0.9, 1.2};
  std::optional<double> fixedXCm;      // nullopt → uniform random over wire span
};

struct GasConfig {
  double temperatureK        = 293.15;
  double pressureTorr        = 760.0;
  std::string gasFile        = "ar_70_co2_30_e400.gas";
  bool   enablePenning       = true;
  int    nCollisions         = 10;
  double maxElectronEnergyEV = 400.0;
  int    nFieldPoints        = 20;         // number of E-field grid points for Magboltz
  double eFieldMaxVcm        = 300000.0;  // upper E-field limit [V/cm] for the gas table
  double wValueEV            = 26.0;      // effective ionisation energy [eV per ion pair]
};

struct SimulationConfig {
  std::size_t nEvents          = 1000;
  std::size_t maxAvalancheSize = 500000;
  double      timeWindowNs     = 300.0;
  double      timeStepNs       = 0.5;
};

struct Config {
  GeometryConfig   geometry;
  SourceConfig     source;
  GasConfig        gas;
  SimulationConfig simulation;
};

// ─── Per-distance summary ─────────────────────────────────────────────────────

struct DistanceSummary {
  double      distanceMm            = 0.;
  std::size_t nEvents               = 0;
  std::size_t nInteracted           = 0;
  double      interactionFraction   = 0.;
  double      meanAnodeChargeFC     = 0.;
  double      rmsAnodeChargeFC      = 0.;
  double      semAnodeChargeFC      = 0.;
  double      meanCathodeChargeFC   = 0.;
  double      rmsCathodeChargeFC    = 0.;
  double      semCathodeChargeFC    = 0.;
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

void EnsureDirectory(const fs::path& p) { fs::create_directories(p); }

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

std::string ReadString(const json& obj, std::string_view sec, std::string_view key, const std::string& fb) {
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

  if (const auto* s = FindSection(root, "source")) {
    cfg.source.energyKeV    = ReadDouble(*s, "source", "energy_keV", cfg.source.energyKeV);
    cfg.source.distancesMm  = ReadDoubleArray(*s, "source", "source_distances_mm", cfg.source.distancesMm);
    auto* xp = FindMember(*s, "x_position_cm", {"source", "x_position_cm"});
    if (xp && !xp->is_null()) {
      if (!xp->is_number()) ThrowJsonTypeError({"source", "x_position_cm"}, "a number or null");
      cfg.source.fixedXCm = xp->get<double>();
    }
  }

  if (const auto* g = FindSection(root, "gas")) {
    cfg.gas.temperatureK  = ReadDouble(*g, "gas", "temperature_K",       cfg.gas.temperatureK);
    cfg.gas.pressureTorr  = ReadDouble(*g, "gas", "pressure_Torr",       cfg.gas.pressureTorr);
    cfg.gas.gasFile       = ReadString(*g, "gas", "gas_file",             cfg.gas.gasFile);
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
  }

  if (cfg.geometry.wirePitchCm <= 0.)  throw std::runtime_error("geometry.wire_pitch_cm must be positive");
  if (cfg.geometry.wireDiamUm  <= 0.)  throw std::runtime_error("geometry.wire_diameter_um must be positive");
  if (cfg.geometry.gapCm       <= 0.)  throw std::runtime_error("geometry.gap_cm must be positive");
  if (cfg.geometry.nWires      <= 0)   throw std::runtime_error("geometry.n_wires must be positive");
  if (cfg.geometry.wireVoltageV<= 0.)  throw std::runtime_error("geometry.wire_voltage_V must be positive");
  if (cfg.source.distancesMm.empty())  throw std::runtime_error("source.source_distances_mm must not be empty");
  if (cfg.gas.temperatureK     <= 0.)  throw std::runtime_error("gas.temperature_K must be positive");
  if (cfg.gas.pressureTorr     <= 0.)  throw std::runtime_error("gas.pressure_Torr must be positive");
  if (cfg.simulation.nEvents   == 0)   throw std::runtime_error("simulation.n_events must be at least 1");
  if (cfg.simulation.timeWindowNs <= 0.) throw std::runtime_error("simulation.time_window_ns must be positive");
  if (cfg.simulation.timeStepNs   <= 0.) throw std::runtime_error("simulation.time_step_ns must be positive");

  return cfg;
}

// ─── Gas setup ────────────────────────────────────────────────────────────────

void SetupGas(MediumMagboltz& gas, const GasConfig& cfg) {
  gas.SetTemperature(cfg.temperatureK);
  gas.SetPressure(cfg.pressureTorr);

  if (fs::exists(cfg.gasFile)) {
    std::cout << "  Loading gas table from: " << cfg.gasFile << "\n";
    gas.LoadGasFile(cfg.gasFile);
    // LoadGasFile restores the EFINAL ceiling stored in the file.  Override it
    // so the collision-rate table for AvalancheMicroscopic is pre-built to a
    // high enough energy to cover electrons near the wire (~160 kV/cm at 1900 V
    // can push electrons past 400 eV), preventing silent transport aborts.
    gas.SetMaxElectronEnergy(cfg.maxElectronEnergyEV);
  } else {
    std::cout << "  Gas file not found: " << cfg.gasFile << "\n"
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
    gas.WriteGasFile(cfg.gasFile);
    std::cout << "  Gas table saved to: " << cfg.gasFile << "\n";
  }

  if (cfg.enablePenning) {
    if (!gas.EnablePenningTransfer()) {
      std::cerr << "  Warning: Penning transfer could not be enabled.\n";
    } else {
      std::cout << "  Penning transfer enabled.\n";
    }
  }

  // CO2+ is the dominant drifting ion in Ar:CO2 mixtures.  After the initial
  // photoelectric absorption, Ar+ rapidly charge-transfers to CO2+ because
  // the CO2 ionisation potential (13.78 eV) is lower than Ar (15.76 eV).
  // Garfield++ accepts a single positive-ion mobility table per gas object.
  const char* garfieldInstall = std::getenv("GARFIELD_INSTALL");
  if (garfieldInstall) {
    const std::string mobFile = std::string(garfieldInstall) +
                                "/share/Garfield/Data/IonMobility_CO2+_CO2.txt";
    if (fs::exists(mobFile)) {
      gas.LoadIonMobility(mobFile);
      std::cout << "  CO2+ ion mobility loaded.\n";
    } else {
      std::cerr << "  Warning: IonMobility_CO2+_CO2.txt not found at " << mobFile << "\n";
    }
  } else {
    std::cerr << "  Warning: GARFIELD_INSTALL not set; ion mobility not loaded.\n";
  }
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

void SetupSensor(Sensor& sensor, ComponentAnalyticField& cmp, const Config& cfg) {
  const auto& geom = cfg.geometry;
  const auto& sim  = cfg.simulation;

  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&cmp, "anode");    // all wires together
  sensor.AddElectrode(&cmp, "cathode");  // bottom cathode plane

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
                                 double sourceDistanceMm,
                                 Sensor& sensor,
                                 TDirectory* distDir) {
  const auto& geom = cfg.geometry;
  const auto& sim  = cfg.simulation;

  const double sourceYCm = sourceDistanceMm * 0.1; // mm → cm
  // Clamp to strictly inside the gas gap (away from cathode surfaces)
  const double y0 = std::max(-geom.gapCm + 1.e-4,
                              std::min(geom.gapCm - 1.e-4, sourceYCm));
  // Half-span of the wire array for random x sampling
  const double xHalfWires = (geom.nWires - 1) / 2. * geom.wirePitchCm;

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

  TProfile pAnodeSignal("p_anode_signal",
                        "Mean anode signal;t [ns];#LTi_{anode}#GT [fC/ns]",
                        static_cast<int>(nBins), 0., sim.timeWindowNs);
  TProfile pCathodeSignal("p_cathode_signal",
                          "Mean cathode signal;t [ns];#LTi_{cathode}#GT [fC/ns]",
                          static_cast<int>(nBins), 0., sim.timeWindowNs);

  for (TH1* h : std::initializer_list<TH1*>{
           &hAnodeQ, &hCathodeQ, &hRatio,
           &hNprimary, &hAvalSize,
           &pAnodeSignal, &pCathodeSignal}) {
    h->SetDirectory(nullptr);
  }

  // ── Transport objects ────────────────────────────────────────────────────────
  // Average primary electrons: N = E_photon / W-value (e.g. 5900 eV / 26 eV ≈ 227)
  const int nPrimary = std::max(1,
      static_cast<int>(std::round(cfg.source.energyKeV * 1.e3 / cfg.gas.wValueEV)));

  AvalancheMicroscopic aval(&sensor);
  if (sim.maxAvalancheSize > 0) aval.EnableAvalancheSizeLimit(sim.maxAvalancheSize);

  // ── Accumulators for summary statistics ──────────────────────────────────────
  std::vector<double> anodeCharges, cathodeCharges, chargeRatios;
  std::vector<double> primaryCounts, avalancheSizes;
  anodeCharges.reserve(sim.nEvents);
  cathodeCharges.reserve(sim.nEvents);

  std::size_t nInteracted = 0;
  const std::size_t progressStep = std::max<std::size_t>(1, sim.nEvents / 10);

  // ── Event loop ───────────────────────────────────────────────────────────────
  for (std::size_t ev = 0; ev < sim.nEvents; ++ev) {
    sensor.ClearSignal();

    const double x0 = cfg.source.fixedXCm.has_value()
                          ? *cfg.source.fixedXCm
                          : gRandom->Uniform(-xHalfWires, xHalfWires);

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

    // Integrate the binned induced-current signal to obtain charge.
    // GetSignal returns the induced current in fC/ns; multiplying by the bin
    // width (ns) and summing gives charge in fC.  Scale by nPrimary.
    //
    // Sign convention (Shockley-Ramo, confirmed by inspection of ROOT output):
    //   anode   : dominated by electron collection → net integral is NEGATIVE
    //             → negate to obtain the conventionally positive collected charge
    //   cathode : dominated by ion drift toward pad → net integral is POSITIVE
    //
    // NOTE: Sensor::GetInducedCharge() uses a separate per-electrode "charge"
    // accumulator that AvalancheMicroscopic does not populate (it only calls
    // AddSignal into the time-binned arrays), so GetInducedCharge always returns
    // zero here.  The manual integral below is the correct approach.
    double rawAnode = 0., rawCathode = 0.;
    for (std::size_t k = 0; k < nBins; ++k) {
      const double sigA = sensor.GetSignal("anode",   k);
      const double sigC = sensor.GetSignal("cathode", k);
      rawAnode   += sigA;
      rawCathode += sigC;
      const double t = (static_cast<double>(k) + 0.5) * sim.timeStepNs;
      pAnodeSignal.Fill(t,   sigA * nPrimary);
      pCathodeSignal.Fill(t, sigC * nPrimary);
    }

    const double qAnode   = -rawAnode   * sim.timeStepNs * nPrimary; // [fC]
    const double qCathode =  rawCathode * sim.timeStepNs * nPrimary; // [fC]

    hAnodeQ.Fill(qAnode);
    hCathodeQ.Fill(qCathode);
    anodeCharges.push_back(qAnode);
    cathodeCharges.push_back(qCathode);

    if (qAnode > 0.) {
      const double ratio = qCathode / qAnode;
      hRatio.Fill(ratio);
      chargeRatios.push_back(ratio);
    }

    if ((ev + 1) % progressStep == 0 || ev + 1 == sim.nEvents) {
      std::cout << "  dist=" << FormatNumber(sourceDistanceMm) << " mm: "
                << (ev + 1) << "/" << sim.nEvents
                << " events processed\n";
    }
  }

  // ── Write histograms ─────────────────────────────────────────────────────────
  if (distDir) {
    distDir->cd();
    hAnodeQ.Write("h_anode_charge");
    hCathodeQ.Write("h_cathode_charge");
    hRatio.Write("h_ratio_charge");
    hNprimary.Write("h_n_primary_electrons");
    hAvalSize.Write("h_avalanche_size");
    pAnodeSignal.Write("p_anode_signal");
    pCathodeSignal.Write("p_cathode_signal");
  }

  // ── Build summary ─────────────────────────────────────────────────────────────
  DistanceSummary s;
  s.distanceMm          = sourceDistanceMm;
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
  std::vector<double> qa(n), qaE(n), qc(n), qcE(n), rat(n), ratE(n);

  for (std::size_t i = 0; i < n; ++i) {
    x[i]    = sums[i].distanceMm;
    qa[i]   = sums[i].meanAnodeChargeFC;   qaE[i]  = sums[i].semAnodeChargeFC;
    qc[i]   = sums[i].meanCathodeChargeFC; qcE[i]  = sums[i].semCathodeChargeFC;
    rat[i]  = sums[i].meanChargeRatio;     ratE[i] = sums[i].semChargeRatio;
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

  auto gAnode   = MakeGraph("g_anode_charge",
    "Mean anode charge;Source distance from wire plane [mm];Q_{anode} [fC]", qa, qaE, 20);
  auto gCathode = MakeGraph("g_cathode_charge",
    "Mean cathode charge;Source distance from wire plane [mm];Q_{cathode} [fC]", qc, qcE, 21);
  auto gRatio   = MakeGraph("g_charge_ratio",
    "Charge ratio;Source distance from wire plane [mm];Q_{cathode}/Q_{anode}", rat, ratE, 22);

  TCanvas canvas("c_tgc_summary", "TGC summary", 1800, 500);
  canvas.Divide(3, 1);
  canvas.cd(1); gAnode.Draw("APL");
  canvas.cd(2); gCathode.Draw("APL");
  canvas.cd(3); gRatio.Draw("APL");
  EnsureDirectory(pngPath.parent_path());
  canvas.SaveAs(pngPath.string().c_str());
}

// ─── CSV summary ─────────────────────────────────────────────────────────────

void WriteSummaryCsv(const fs::path& path, const std::vector<DistanceSummary>& sums) {
  std::ofstream f(path);
  if (!f) throw std::runtime_error("Cannot write CSV: " + path.string());

  f << "source_distance_mm,n_events,n_interacted,interaction_fraction,"
       "mean_anode_charge_fC,rms_anode_charge_fC,sem_anode_charge_fC,"
       "mean_cathode_charge_fC,rms_cathode_charge_fC,sem_cathode_charge_fC,"
       "mean_charge_ratio,rms_charge_ratio,sem_charge_ratio,"
       "mean_primary_electrons,mean_avalanche_size\n";

  f << std::fixed << std::setprecision(6);
  for (const auto& s : sums) {
    f << s.distanceMm         << ','
      << s.nEvents            << ','
      << s.nInteracted        << ','
      << s.interactionFraction<< ','
      << s.meanAnodeChargeFC  << ','
      << s.rmsAnodeChargeFC   << ','
      << s.semAnodeChargeFC   << ','
      << s.meanCathodeChargeFC<< ','
      << s.rmsCathodeChargeFC << ','
      << s.semCathodeChargeFC << ','
      << s.meanChargeRatio    << ','
      << s.rmsChargeRatio     << ','
      << s.semChargeRatio     << ','
      << s.meanPrimaryElectrons << ','
      << s.meanAvalancheSize  << '\n';
  }
}

// ─── Config echo ──────────────────────────────────────────────────────────────

json ConfigToJson(const Config& cfg) {
  json jSrc = {
    {"energy_keV",           cfg.source.energyKeV},
    {"source_distances_mm",  cfg.source.distancesMm}
  };
  jSrc["x_position_cm"] = cfg.source.fixedXCm.has_value()
                               ? json(*cfg.source.fixedXCm)
                               : json(nullptr);
  return {
    {"geometry", {
      {"wire_pitch_cm",    cfg.geometry.wirePitchCm},
      {"wire_diameter_um", cfg.geometry.wireDiamUm},
      {"gap_cm",           cfg.geometry.gapCm},
      {"n_wires",          cfg.geometry.nWires},
      {"wire_voltage_V",   cfg.geometry.wireVoltageV}
    }},
    {"source", jSrc},
    {"gas", {
      {"temperature_K",         cfg.gas.temperatureK},
      {"pressure_Torr",         cfg.gas.pressureTorr},
      {"gas_file",              cfg.gas.gasFile},
      {"enable_penning",        cfg.gas.enablePenning},
      {"n_magboltz_collisions",  cfg.gas.nCollisions},
      {"max_electron_energy_eV", cfg.gas.maxElectronEnergyEV},
      {"n_field_points",         cfg.gas.nFieldPoints},
      {"e_field_max_vcm",        cfg.gas.eFieldMaxVcm},
      {"w_value_eV",             cfg.gas.wValueEV}
    }},
    {"simulation", {
      {"n_events",          cfg.simulation.nEvents},
      {"max_avalanche_size",cfg.simulation.maxAvalancheSize},
      {"time_window_ns",    cfg.simulation.timeWindowNs},
      {"time_step_ns",      cfg.simulation.timeStepNs}
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
      cfg.source.distancesMm = {*opts.singleDistanceMm};

    const fs::path runDir = opts.outDir / BuildRunFolderName(cfg);
    EnsureDirectory(runDir);

    std::cout << "TGC Garfield++ simulation\n"
              << "  config  : " << opts.configPath << "\n"
              << "  output  : " << runDir << "\n"
              << "  geometry: " << cfg.geometry.nWires << " wires, "
              << cfg.geometry.wirePitchCm * 10. << " mm pitch, "
              << cfg.geometry.wireDiamUm << " μm diameter, "
              << cfg.geometry.gapCm * 10. << " mm gap\n"
              << "  voltage : " << cfg.geometry.wireVoltageV << " V (wires), 0 V (cathodes)\n"
              << "  gas     : Ar:CO2 70:30, " << cfg.gas.temperatureK << " K, "
              << cfg.gas.pressureTorr << " Torr\n"
              << "  source  : " << cfg.source.energyKeV << " keV, "
              << cfg.source.distancesMm.size() << " distance point(s)\n"
              << "  events  : " << cfg.simulation.nEvents << " per point\n";

    // Gas is shared across all distance points
    std::cout << "\nSetting up gas...\n";
    MediumMagboltz gas("ar", 70., "co2", 30.);
    SetupGas(gas, cfg.gas);

    // Geometry and sensor are shared across all distance points
    ComponentAnalyticField cmp;
    BuildGeometry(cmp, gas, cfg.geometry);

    Sensor sensor;
    SetupSensor(sensor, cmp, cfg);

    // ROOT output
    TFile rootFile((runDir / "tgc_sim.root").string().c_str(), "RECREATE");
    if (rootFile.IsZombie())
      throw std::runtime_error("Failed to create ROOT file in " + runDir.string());

    TDirectory* summaryDir = rootFile.mkdir("summary");

    std::vector<DistanceSummary> allSummaries;

    for (const double distMm : cfg.source.distancesMm) {
      std::cout << "\n--- Source distance: " << FormatNumber(distMm) << " mm ---\n";

      const std::string tag = "dist_" + FileSafeNumber(distMm) + "mm";
      TDirectory* distDir = rootFile.mkdir(tag.c_str());
      if (!distDir) throw std::runtime_error("Failed to create ROOT dir: " + tag);

      DistanceSummary summary = RunDistancePoint(cfg, distMm, sensor, distDir);
      allSummaries.push_back(summary);

      std::cout << "  ⟨Q_anode⟩   = " << FormatNumber(summary.meanAnodeChargeFC)   << " fC"
                << "  ±" << FormatNumber(summary.semAnodeChargeFC)   << " (SEM)\n"
                << "  ⟨Q_cathode⟩ = " << FormatNumber(summary.meanCathodeChargeFC) << " fC"
                << "  ±" << FormatNumber(summary.semCathodeChargeFC) << " (SEM)\n"
                << "  ⟨ratio⟩     = " << FormatNumber(summary.meanChargeRatio)     << "\n"
                << "  interaction fraction: "
                << FormatNumber(summary.interactionFraction * 100., 2) << "%\n"
                << "  ⟨avalanche size⟩: "
                << FormatNumber(summary.meanAvalancheSize, 0) << " electrons\n";
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
