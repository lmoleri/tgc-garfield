# TGC TODO

## Open

- [x] **Expose all gas config parameters in GUI**
  Three parameters in `GasConfig` are currently hidden in `ConfigPanel._hidden_gas`
  (`gui/app.py`) and are not editable in the UI — only via JSON file:
  - `max_electron_energy_eV` (default 2000.0)
  - `n_field_points` (default 20)
  - `e_field_max_vcm` (default 300000.0)
  Add spinboxes for these in the Gas section of `ConfigPanel`, and connect them to
  `_update_gas_file_label` (all three affect the derived `.gas` filename).

- [ ] **Investigate: no visible difference between source distances +0.7 mm and −0.7 mm**
  By the sign convention (`source_distances_mm` positive → readout-cathode side y < 0;
  negative → cathode_top side y > 0), ±0.7 mm place the primary electron on opposite
  sides of the wire plane. Observed result: Q_cathode and charge ratio are
  indistinguishable — this is not physical.
  Expected physics: avalanche ions are created along the electron's full drift path toward
  the wire. An electron arriving from the readout side (y < 0) creates ions on the
  readout-cathode side → those ions drift to the readout cathode → larger Q_cathode and
  different waveform shape. An electron arriving from the cathode_top side (y > 0) creates
  ions on the cathode_top side → those ions drift away from the readout cathode → smaller
  Q_cathode (or negative Ramo contribution from cathode_top-bound ions).
  Suggested investigation:
  - Print (xi0, yi0) ion-creation positions from `GetElectronEndpoints` to check whether
    the expected y-sign asymmetry in ion positions is present.
  - Compare ±1.2 mm (electrons starting close to each cathode) where the effect is largest.
