# Compact frozen data products

This directory contains compact derived data products associated with the released analysis. Raw STEREO/SECCHI FITS observations are not redistributed.

## Event and phase products

- `all12_events.csv` — frozen table for the 12-event sample.
- `event_phase_jitter_table.csv` — event-level phase and phase-jitter diagnostics.
- `event_phase_jitter_results.json` — structured output of the phase-jitter analysis.
- `height_matrix_no_BH.csv` — height-resolved common-shift diagnostics using the frozen raw-p convention.
- `resonator_candidate_metrics.csv` — compact low-coronal candidate metrics.

## Transport, geometry, and transition products

- `level0_nonlinear_transport_results.json` — Level-0-derived nonlinear-transport diagnostic summary.
- `dynamic_geometry_gate_results.json` — dynamic geometry and cusp-gate diagnostic summary.
- `habbal_standing_mhd_results.json` — stationary-transition and MHD-branch necessary-condition screen.
- `ridge_onset_results.json` — event and ridge-onset summary products.

## Expansion and map products

- `spatial_expanding_results.json` — expansion-aware Event 9 and Event 12 results.
- `spatial_expanding_pair_null.csv` — shifted-pair control table.
- `spatial_expanding_posteriors.npz` — compact expansion-posterior arrays.
- `two_event_maps.npz` — compact Event 9 and Event 12 polar-map arrays used in the frozen analysis.

## Reproducibility note

The CSV, JSON, and NPZ files preserve compact numerical products used to inspect and reproduce the reported event-level diagnostics.

Some products are derived from bias/exposure-normalized Level-0 SECCHI observations. These products are retained as diagnostic results. Full image-level confirmation requires the corresponding public STEREO/SECCHI observations and the appropriate SECCHI_PREP Level-1 calibration and background-processing workflow.

See `../REPRODUCIBILITY.md` for the full reproducibility boundary and data-processing requirements.