# Compact frozen data products

The [Zenodo data version](https://doi.org/10.5281/zenodo.24604508) contains the 17 requested cubes and derived maps, including the separately named Helioviewer reconstruction.

This directory contains derived products, not raw STEREO/SECCHI FITS observations.

[January COR1 data cubes, X-front maps, control intervals, EUVI time series and reproduction commands](cor1-release/README.md)

## Event and phase products

- `all12_events.csv` — frozen 12-event audit table.
- `event_phase_jitter_table.csv` — event-level phase / phase-jitter diagnostics.
- `event_phase_jitter_results.json` — structured output of the phase-jitter analysis.
- `height_matrix_no_BH.csv` — height-resolved common-shift diagnostics using the frozen raw-p convention.
- `resonator_candidate_metrics.csv` — compact low-coronal candidate metrics.

## Transport, geometry, and transition products

- `level0_nonlinear_transport_results.json` — Level-0-derived nonlinear transport diagnostic summary.
- `dynamic_geometry_gate_results.json` — dynamic geometry / cusp-gate diagnostic summary.
- `habbal_standing_mhd_results.json` — stationary-transition and MHD-branch necessary-condition screen.
- `ridge_onset_results.json` — event/ridge onset summary products.

## Expansion and map products

- `spatial_expanding_results.json` — expansion-aware Event 9 / Event 12 results.
- `spatial_expanding_pair_null.csv` — shifted-pair control table.
- `spatial_expanding_posteriors.npz` — compact expansion posterior arrays.
- `two_event_maps.npz` — compact Event 9 / Event 12 polar-map arrays used in the frozen analysis.

The JSON/CSV/NPZ files contain no credentials, private correspondence, or raw mission FITS data.
