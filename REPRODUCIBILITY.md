# Reproducibility boundary

This repository distinguishes two complementary levels of reproducibility.

## 1. Compact frozen products

The files in `data/` preserve the numerical products associated with the reported analysis, including:

- the frozen 12-event table (`all12_events.csv`);
- event phase-jitter outputs (`event_phase_jitter_table.csv`, `event_phase_jitter_results.json`);
- height-localization and common-shift results (`height_matrix_no_BH.csv`);
- low-coronal candidate metrics (`resonator_candidate_metrics.csv`);
- nonlinear-transport and dynamic-geometry summaries;
- spatial-expansion posteriors and shifted-pair controls;
- compact Event 9 and Event 12 map arrays.

These products can be inspected and compared without downloading the raw mission observations.

## 2. Full image-level reproduction

The programs under `analysis/frozen_programs/` are the frozen scientific scripts associated with the released analysis. Some scripts require intermediate image arrays or public STEREO/SECCHI observations that are not redistributed in this repository.

This release is therefore not intended as a one-command raw-FITS-to-publication pipeline. It preserves the scientific analysis code, compact numerical products, and the information needed to trace the reported diagnostics.

### Level-0-derived diagnostics

Some COR1 morphology, geometry, polarization, and transport diagnostics were obtained from bias/exposure-normalized Level-0 SECCHI observations and derived intermediate products.

These Level-0-derived products are preserved because they represent the analysis state from which the reported diagnostic results were obtained.

They should be interpreted as diagnostic image products rather than as absolute mission-calibrated pB measurements.

### Level-1 confirmation

A full publication-grade image-level confirmation requires the corresponding public STEREO/SECCHI observations and the appropriate SECCHI calibration environment, including a complete `SECCHI_PREP` Level-1 calibration and background-processing workflow where applicable.

Calibration, background treatment, image registration, and related preprocessing choices can affect measured amplitudes and morphology. For this reason, the Level-0-derived diagnostic products and a future Level-1 confirmation should be treated as distinct stages of the reproducibility chain.

## Statistical convention

The archived products report the raw empirical probabilities used in this analysis.

Files explicitly marked `no_BH` retain that frozen convention. No Benjamini–Hochberg correction has been retroactively applied to the released numerical products.