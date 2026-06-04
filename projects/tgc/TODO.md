# TGC TODO

## Open

- [ ] **Expose all gas config parameters in GUI**
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
  sides of the wire plane. Observed result: Q_anode, Q_cathode, and charge ratio are
  indistinguishable.
  Leading hypothesis: the two positions are equidistant from the wire (0.7 mm each), so
  drift time and avalanche are identical. The only asymmetry is which cathode receives the
  primary photoionisation ion (~1 elementary charge vs ~10^4 avalanche electrons), which
  is below statistical noise at any practical event count. Need to verify this explanation
  (e.g. run with `enable_ion_drift = false` to isolate the electron signal, or compare
  asymmetric distances such as +0.2 mm vs −0.2 mm where drift-time differences are
  largest relative to the gap).
