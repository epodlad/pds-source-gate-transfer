# January 2008 EUVI east-limb reconstruction

The two original diagnostic PNGs show the traced path for each channel against an axis labelled “Code PA (deg)”. We digitized the centre of the white line at 21 radii between 1.10 and 1.20 solar radii and retained its projected line thickness in the path CSVs. The two original PNGs are included with the code. Run prepare_geometry.py to regenerate the static-map inputs. Running trace_spectrum.py with the previously released sector NPZ files regenerates the path tables and the diagnostic for the earlier angle convention. The path samples and the digitized static-map ranks in this folder are sufficient to repeat the full-frame calculation.

We compared the archived static brightness pattern with every verified observation after mapping the code angle α to helioprojective coordinates. For angle θ we use x = r sin θ, y = r cos θ (x positive toward solar west). Among eight cardinal rotations/reflections, θ = 90° − α gave the strongest morphology agreement in both channels. This maps α ≃ 180° to the east limb. At each of 19 radii (1.05–1.41 solar radii) the test computes the angular Spearman correlation of the time-median observed image with the colorbar-derived archival brightness ranks. The published summary gives the median of the 19 correlations. This is a descriptive registration check. No period or time-series statistic was used to choose the orientation.

## Inputs and execution

Download the six full-frame FITS ZIP assets and verified_frames.csv from the [January 2008 observation release](https://github.com/epodlad/pds-source-gate-transfer/releases/tag/euvi-20080111-14-reconstruction-v1). Keep the ZIPs intact in one directory. Python with NumPy and SciPy is sufficient for the following command.

```bash
OPENBLAS_NUM_THREADS=1 python full_frame_check.py verified_frames.csv /path/to/FITS_ZIPs output --save-cubes
```

The program verifies all 716 exposure hashes and their observation times before reporting success. Its eastern cubes use 71 radii (1.00–1.70 solar radii) and 281 code angles (100°–240°), with a separate UTC time and filename for each of 358 samples per channel. Values are bilinearly sampled exposure-normalized detector DN/s. Off-detector values are NaN. No background subtraction, bad-pixel correction or SECCHI_PREP Level-1 photometric calibration is applied. The original 171 Å report used 357 frames. The six source archives are unchanged. The previously released 140°–200° sectors used a different mapping to solar coordinates.

The time series for each channel is extracted directly from individual FITS using its digitized path, averaging 21 equally weighted angular samples across ±2.5° and then averaging the 21 radii. The path need not fall on the coarser cube grid. For a separate spectral sensitivity test we interpolate two missing 10-minute slots, remove a linear trend, and apply a Hann window. The statistic is the fraction of positive-frequency power in 80–130 minutes. Each of 20,000 stationary Gaussian AR(1) simulations undergoes the same missing-slot handling, detrending and windowing. Its lag-one coefficient is estimated from adjacent observed slots after detrending. This is a new conditional null test. The original EUVI extraction, red-noise and filter-matched algorithms, exact frame selection and their primary probabilities remain unresolved. The geometry agreement does not by itself establish a periodic source or phase linkage with COR2.

The [original software and compact tables](https://doi.org/10.5281/zenodo.22087646) and this reconstruction have distinct provenance. Original FITS are from the NASA GSFC STEREO Science Center, [STEREO-A/SECCHI/EUVI FITS](https://doi.org/10.48322/260z-fs61). Credit the STEREO/SECCHI team.
