# PDS source–gate–transfer analysis code

Research software and compact data products associated with the manuscript:

**Where Do Periodic Density Structures Acquire Coherence? A Low-Coronal Resonator Candidate and a 12-Event STEREO/SECCHI Transfer Test**

Author: **Olena Podladchikova**

## Scope

This repository provides the analysis code and compact frozen data products for the 12-event STEREO-A/SECCHI PDS study. The scientific programs under `analysis/frozen_programs/` implement the event-phase and phase-jitter, nonlinear-transport, expansion-aware geometry, dynamic-geometry, and stationary-transition diagnostics used in the study. The accompanying files document the derived products, software dependencies, and reproducibility requirements.

The code covers:

- event-by-event phase and phase-jitter diagnostics;
- nonlinear travel-time tests;
- Level-0 COR1 diagnostic reconstruction;
- expanding-ridge and streamer-width sensitivity;
- X-/diamond-like ridge-node tests for Events 9 and 12;
- dynamic cusp/current-sheet geometry tests;
- stationary-transition and MHD-branch necessary-condition screens;
- the frozen 12-event coronal reassessment.

## Repository contents

`analysis/frozen_programs/` contains the scientific analysis scripts.

`analysis/specifications/` contains the corresponding analysis specifications and test protocols.

`data/` contains compact CSV, JSON, and NPZ products supporting the reported numerical diagnostics without redistributing raw mission FITS observations.

`tests/` contains archive validation checks for the released software and compact products.

## Data and full image-level reproduction

Raw STEREO/SECCHI observations are not redistributed and remain available from the mission archive.

Some analysis scripts require intermediate image cubes or archive-derived products that are not included in this compact repository. A full image-level rerun therefore requires the corresponding public STEREO/SECCHI observations and the appropriate SECCHI calibration environment.

The included compact products support inspection and reproduction of the reported event-level quantities, phase-jitter results, height-localization diagnostics, expansion tests, and map/posterior products.

A mission-level confirmation of the Level-0-derived diagnostic requires a complete `SECCHI_PREP` Level-1 calibration and background-processing workflow.

## Python environment

Recommended: Python 3.11 or newer.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Run the archive validation checks with:

```bash
python tests/validate_archive.py
```

## Reproducibility

See `REPRODUCIBILITY.md` for details on reproduction from the compact data products and the requirements for a full image-level rerun.

See `CODE_PROVENANCE.md` for information about the frozen scientific scripts included in this release.

## Citation

If you use this software or its derived products in scientific work, please cite the software release and the associated manuscript.

Citation metadata are provided in `CITATION.cff`.

## License

This software is released under the **BSD 3-Clause License**. See `LICENSE`.