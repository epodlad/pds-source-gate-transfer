# Reproducing the January 2008 observational analysis

This folder accompanies *Where Do Periodic Density Structures Acquire Coherence? A Low-Coronal Resonator Candidate and a 12-Event STEREO/SECCHI Transfer Test*. It adds recovered intermediate arrays to the compact software release. The archive version of the manuscript used for figure checks is the August 24 preprint.

## Verified calculations

| Calculation | Starting data | Result |
|---|---|---|
| COR1 event-phase and jitter tests | Frozen time–height maps, 12 COR2 event times and the original reference model | The original program reproduced the complete archived result JSON, including each event, the null sample count and the global tests. |
| Manuscript figures and tables | Numerical inputs in the repository’s data directory | Seven figures matched the archived renders pixel for pixel at 108 dpi and in extracted text. Both LaTeX tables matched byte for byte. |
| Expansion test for events 9 and 12 | The 131 archived shifted control pairs | Two controls met the measured-width threshold. The corrected count gives (2+1)/(131+1) = 0.022727. |
| EUVI low-coronal brightness | Archived EUVI result tables and their ten multiple-testing families | All 120 Benjamini–Hochberg values were recalculated exactly. The nominal band-power probabilities at 171 Å and 1.10–1.20 solar radii were 0.034931 and 0.027865. Their adjusted values were 0.253050 and 0.334378. The spectral maximum was at 128.21 min. The separate circular-shift phase statistic gave p = 0.911111. |

Run from the repository root with the scientific Python dependencies installed. Choose a new output directory.

```bash
python recovery/article1_20080111_14/reproduce.py article1_rerun
python recovery/euvi_20080111_14/check_results.py
```

The two NPZ files here preserve the paired Level-0-derived COR1 brightness and polarization maps and the separate processed-image sensitivity maps. The phase program uses the first file. The archived phase calculation and the EUVI band-power probabilities have different input requirements.

## Original inputs still sought

The historical ridge-onset program imports `estimate_path14` from `test_cor1_nonradial_source_20080111_14.py` and `trace_euvi_path` from `test_euvi195_lower_resonator_20080111_14.py`. Its input names include `euvi171_polar_cube_10min.npz`, `euvi195_polar_cube_10min.npz`, `cor1_pb_fitpol_cube_15min.npz` and `cor2_blind_tracks.npz`. The archived summary specifies 357 EUVI 171 Å frames and 358 EUVI 195 Å frames. The exact selection, sampling geometry, original EUVI time series and original null-generation implementation are still needed for an exact raw-to-probability comparison. The COR1 primary probabilities 0.0084 and 0.0248 also remain unverified from the original image inputs.

## Observation and software provenance

| Product | Source used or identified | Dataset |
|---|---|---|
| EUVI-A FITS used for this reconstruction | [NASA GSFC STEREO Science Center, Ahead EUVI brightness frames in the L0 tree](https://stereo-ssc.nascom.nasa.gov/pub/ins_data/secchi/L0/a/img/euvi/). The mission identifies these files as Level 0.5. Each candidate has a direct archive URL in the adjacent EUVI inventory. | [STEREO-A/SECCHI/EUVI FITS, 10.48322/260z-fs61](https://doi.org/10.48322/260z-fs61) |
| COR1-A exposure triplets | The frozen Level-0 nonlinear-transport program requests exposures from the [NRL SECCHI post-flight sequence archive](https://secchi.nrl.navy.mil/postflight/cor1/L0/a/seq/). The separate processed-image test uses Helioviewer JP2. | [STEREO-A/SECCHI/COR1 FITS, 10.48322/9fer-5m90](https://doi.org/10.48322/9fer-5m90) |
| COR2-A | The recovered event tables preserve the anchors. The exact mirror and full identities of the original COR2 FITS have yet to be established. | [STEREO-A/SECCHI/COR2 FITS, 10.48322/qfy4-7q03](https://doi.org/10.48322/qfy4-7q03) |
| Archived software v1.0.0 | The versioned original compact release. Recovery additions on this branch postdate the deposition. | [Zenodo, 10.5281/zenodo.22087646](https://doi.org/10.5281/zenodo.22087646) |

A local brightness excess and alignment of downstream event phases are distinct measurements. These computations leave the physical origin of the density structures open.
