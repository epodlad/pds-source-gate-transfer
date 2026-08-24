# PDS source–gate–transfer analysis code

Research software and compact frozen products associated with the manuscript:

**Where Do Periodic Density Structures Acquire Coherence? A Low-Coronal Resonator Candidate and a 12-Event STEREO/SECCHI Transfer Test**

Author: **Olena Podladchikova**

## Scope

This repository preserves the analysis code used for the 12-event STEREO-A/SECCHI PDS study. The scientific programs under `analysis/frozen_programs/` are preserved unchanged from the recovered source-and-code package. The public release cleans the repository structure, removes manuscript/editorial working files, documents the reproducibility boundary, and completes the Python dependency list without changing the frozen scientific analysis.

The code covers:

- event-by-event phase and phase-jitter diagnostics;
- nonlinear travel-time tests;
- Level-0 COR1 diagnostic reconstruction;
- expanding-ridge and streamer-width sensitivity;
- X-/diamond-like ridge-node tests for Events 9 and 12;
- dynamic cusp/current-sheet geometry tests;
- stationary-transition / MHD-branch necessary-condition screens;
- the frozen 12-event coronal reassessment.

## What is included

`analysis/frozen_programs/` contains the frozen scientific scripts.

`analysis/specifications/` contains the corresponding frozen test specifications/protocols.

`data/` contains compact CSV, JSON, and NPZ products used to document and reproduce the reported numerical diagnostics without redistributing raw mission FITS files.

## What is not included

Raw STEREO/SECCHI observations are not redistributed. Several frozen scripts require intermediate image cubes or archive retrieval products that are not present in this compact release. These scripts are preserved for method provenance and for a full rerun in the documented SECCHI analysis environment.

The compact data products are sufficient to inspect the reported event-level quantities, phase-jitter results, height-localization diagnostics, expansion tests, and map/posterior products. A mission-level confirmation of the Level-0-derived diagnostic requires a complete `SECCHI_PREP` Level-1 calibration/background workflow.

## Python environment

Recommended: Python 3.11 or newer.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\\Scripts\\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Run the archive validation checks with:

```bash
python tests/validate_archive.py
```

## Reproducibility

See `REPRODUCIBILITY.md` for the distinction between compact-product reproduction and a full image-level rerun. See `CODE_PROVENANCE.md` for the provenance of the frozen scripts.

## Citation

Citation metadata are provided in `CITATION.cff`. The Zenodo DOI can be added after the public v1.0.0 release is archived.

## License

This software is released under the **BSD 3-Clause License**. See `LICENSE`.
