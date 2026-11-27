# Reproducibility boundary

This archive intentionally distinguishes two levels of reproducibility.

## 1. Compact frozen products

The files in `data/` preserve the numerical products used to audit the reported analysis, including:

- the frozen 12-event table (`all12_events.csv`);
- event phase-jitter outputs (`event_phase_jitter_table.csv`, `event_phase_jitter_results.json`);
- height-localization/common-shift results (`height_matrix_no_BH.csv`);
- low-coronal candidate metrics (`resonator_candidate_metrics.csv`);
- nonlinear transport and dynamic-geometry summaries;
- spatial-expansion posteriors and shifted-pair controls;
- compact Event 9 / Event 12 map arrays.

These products can be inspected without downloading raw mission data.

## 2. Full image-level rerun

The [COR1 release guide](data/cor1-release/README.md) lists the released inputs and commands for rerunning the cell, node-order, standing-transition, EUVI, COR2 and WISPR diagnostics. A raw-image run uses public mission data.

For a full mission-level rerun, use the public STEREO/SECCHI observations and the appropriate calibration environment. In particular, the independent Level-0-derived polarization diagnostic is not a substitute for a complete `SECCHI_PREP` Level-1 calibration/background chain.

## Statistical convention

The archived products report the raw empirical probabilities used in this analysis. Files explicitly marked `no_BH` retain that frozen convention; no Benjamini–Hochberg correction has been retroactively applied to the archived products.
