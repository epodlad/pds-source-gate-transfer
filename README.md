# PDS source–gate–transfer analysis code

Research software and compact data products associated with the manuscript:

**Where Do Periodic Density Structures Acquire Coherence? A Low-Coronal Resonator Candidate and a 12-Event STEREO/SECCHI Transfer Test**

Author: **Olena Podladchikova**

## Scope

Periodic density structures (PDS) are recurring density variations observed in the solar wind. This code tests where their regular timing may arise and how it is transferred through the corona, using a sample of 12 STEREO-A/SECCHI events.

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

The v1.0.0 software archive contains compact products. The [January 2008 reconstruction and test results](https://github.com/epodlad/pds-source-gate-transfer/tree/recovery/euvi-20261118/recovery) include recovered COR1 maps, an executable 12-event phase test, and the [EUVI frame inventory](https://github.com/epodlad/pds-source-gate-transfer/blob/recovery/euvi-20261118/recovery/euvi_20080111_14/euvi_candidate_inventory.csv). The exact original EUVI frame selection and statistical extraction remain open. The reconstruction covers 11 January 2008 at 12:00 UTC to 14 January at 00:00 UTC.

| STEREO-A observation | Data provenance | Mission dataset DOI |
|---|---|---|
| EUVI 171 and 195 Å | The reconstruction reads individual Level 0.5 brightness FITS from the [NASA GSFC STEREO Science Center](https://stereo-ssc.nascom.nasa.gov/pub/ins_data/secchi/L0/a/img/euvi/). The inventory lists the individual file URLs. | [10.48322/260z-fs61](https://doi.org/10.48322/260z-fs61) |
| COR1 polarization triplets | The original nonlinear-transport code reads the [NRL SECCHI post-flight sequence archive](https://secchi.nrl.navy.mil/postflight/cor1/L0/a/seq/). The processed-image sensitivity comparison used Helioviewer JP2 images. | [10.48322/9fer-5m90](https://doi.org/10.48322/9fer-5m90) |
| COR2 | Event anchors are preserved in the [event table](data/event_phase_jitter_table.csv). The original frame list and retrieval mirror remain to be identified. | [10.48322/qfy4-7q03](https://doi.org/10.48322/qfy4-7q03) |

The original compact software release is archived on [Zenodo, 10.5281/zenodo.22087646](https://doi.org/10.5281/zenodo.22087646). The [AASTeX data-availability text](manuscript/data_availability.tex) includes the versioned software DOI, mission dataset citations and the archive locations.

Some analysis scripts require intermediate image cubes or archive-derived products that are not included in this compact repository. A full image-level rerun therefore requires the corresponding public STEREO/SECCHI observations and the appropriate SECCHI calibration environment.

The included compact products support inspection and reproduction of the reported event-level quantities, phase-jitter results, height-localization diagnostics, expansion tests, and map/posterior products.

A mission-level confirmation of the Level-0-derived diagnostic requires a complete `SECCHI_PREP` Level-1 calibration and background-processing workflow.

## Python environment

Recommended: Python 3.11 or newer. Download or clone this repository, open a terminal in its folder (the one containing `requirements.txt`), and create a Python environment.

### macOS or Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### Windows — Command Prompt

```bat
py -3 -m venv .venv
.venv\Scripts\activate.bat
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### Check the compact archive

In the activated environment on either system, run:

```bash
python tests/validate_archive.py
```

This checks the archived programs and compact data products. It does not rerun the full analysis from raw mission images; those requirements are described below.

## Reproducibility

See `REPRODUCIBILITY.md` for details on reproduction from the compact data products and the requirements for a full image-level rerun.

See `CODE_PROVENANCE.md` for information about the frozen scientific scripts included in this release.

## Citation

If you use this software or its derived products in scientific work, please cite the archived software release and the associated manuscript.

**Zenodo DOI (v1.0.0):** https://doi.org/10.5281/zenodo.22087646

Citation metadata are provided in `CITATION.cff`.

## License

This software is released under the **BSD 3-Clause License**. See `LICENSE`.