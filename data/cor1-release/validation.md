# Validation of the COR1 streamer analysis

The [clean GitHub Actions run](https://github.com/epodlad/pds-source-gate-transfer/actions/runs/46061626668) downloaded the public release, checked the published SHA-256 sums and reran eight frozen programs. All eight numerical outputs matched their archived references within the tolerances in replay.py. The [release guide](README.md) gives the commands and input files.

BFF stands for Bandpass Frame Filtering. The panel-b code applies Gaussian smooths with full widths at half maximum of 75 and 618 min to each pixel time series, subtracts the long smoothing from the short one, and normalizes the residual. These widths are filter parameters. [Alzate and Di Matteo](https://doi.org/10.3847/1538-4357/ae07dd)

## X-front contrast and ridge-width sensitivity

Both tests use the same 131 shifted pairs, with the later anchor 360 min after the earlier anchor. A control must meet the test statistic at both events in at least two reductions. The contrast test measures the coherent X contrast. The ridge test measures the posterior probability for pattern speeds between 25 and 300 km/s. The p values below are raw empirical values with the plus-one correction.

| Test | Image reductions | Control pairs matching or exceeding the observed pair | Raw p |
| --- | --- | ---: | ---: |
| X-front contrast | base-60, base-120, NRGF-60 | 1/131 | 2/132 = 0.0152 |
| Figure 8b, fixed ridge width | BFF, base-60, NRGF-60 | 0/131 | 1/132 = 0.0076 |
| Figure 8b, width proportional to radius | BFF, base-60, NRGF-60 | 0/131 | 1/132 = 0.0076 |
| Figure 8b, measured width | BFF, base-60, NRGF-60 | 2/131 | 3/132 = 0.0227 |

In minutes since 2008-01-11 00:00 UT, the control anchors are (2100, 2460) for X-front contrast and (1500, 1860) and (3210, 3570) for measured width. The [archived Figure 8b table](../spatial_expanding_pair_null.csv) contains a row for every pair.

## Stationary transition and node order

The frozen stationary-transition test gives raw cross-validation p values of 0.3600 for acceleration and 0.3333 for deceleration. One of 12 ridges is identified. Zero of 12 show both an identified ridge and a phase kink. Polarized brightness satisfies the compression quality criterion at 1/4 nodes in event 9 and 1/3 in event 12. Across 1.6–3.0 solar radii the width coverage is 0.161 in total brightness and 0.321 in polarized brightness.

The data do not meet the frozen criteria for a stationary, event-linked transition. A moving or reforming transition remains possible. The activation-order test gives the following model-aligned times for the primary BFF reduction on 13 January 2008.

| Event | Increasing radius (solar radii) | Times at those radii (UT) | Reductions passing all quality checks |
| --- | --- | --- | ---: |
| 9 | 1.775, 2.125, 2.500, 2.825 | 09:09:01, 10:34:36, 12:24:15, 14:14:51 | 0/3 |
| 12 | 1.950, 2.350, 2.700 | 19:30:58, 20:32:20, 20:27:55 | 0/3 |

Event 9 is ordered outward in 3/3 reductions and event 12 in 2/3. Neither event passes the full activation-order criteria. The formal score probability reaches 1/132, but the failed quality checks preclude a detection. This test supplies no p value for a repeated activation interval. See node_times.md in the release for all times, the timing grid and the quality checks.

## Exploratory streamer-width measurement

The following values use the median over 60 h (240 frames). There are 21 sampled radii between 1.80 and 2.30 solar radii. The table gives the fraction of valid widths and the smallest measured values' radii, in solar radii. The last column assumes a circular section and an isothermal temperature of 1 MK. No minimum has valid measurements at both adjacent radii, so these data do not locate a throat.

| Half-window | Brightness | Valid radii | Radius at minimum width | Radius at minimum effective area |
| ---: | --- | ---: | ---: | ---: |
| 12° | tB | 1/21 (0.048) | 2.300 | 2.300 |
| 12° | pB | 5/21 (0.238) | 2.175 | 2.175 |
| 16° | tB | 2/21 (0.095) | 2.275 | 2.275 |
| 16° | pB | 4/21 (0.190) | 1.975 | 1.975 |
| 20° | tB | 3/21 (0.143) | 2.250 | 2.250 |
| 20° | pB | 2/21 (0.095) | 1.850 | 2.275 |

The same analysis using all 352 frames gives 2–3 valid radii in tB and 3–4 in pB. The complete profiles and both time windows are in width_profiles.csv and width_summary.csv.

## Reconstruction of the archived Helioviewer input

The Helioviewer JPEG2000 images provide processed 8-bit COR1-A total brightness without absolute photometric calibration. The reconstructed archive input has 352 time steps, 77 radii and 71 position angles. The original sector cube's checksum is unknown. The eight arrays below were compared against the original **derived transport map**. The difference is reconstructed minus archived on common finite samples. Pearson r is rounded to six decimals.

| Array | Finite pairs | Pearson r | Median difference | Max absolute difference |
| --- | ---: | ---: | ---: | ---: |
| minutes | 351 | 1.000000 | 0 | 0 |
| radii | 58 | 1.000000 | 0 | 0 |
| zmap | 20358 | 1.000000 | 0 | 0 |
| raw_bff | 20358 | 1.000000 | 0 | 0 |
| lags | 62 | 1.000000 | 0 | 0 |
| event_stack | 3596 | 1.000000 | 0 | 0 |
| traced_path_pa | 77 | 1.000000 | 0 | 0 |
| event_times | 12 | 1.000000 | 0 | 0 |

All array dimensions and finite-data masks agree. The archived transport map SHA-256 is 8d06e2a4fdc64765f0407763ea95bfbe320cf68bb43223efee9c0526e3edaadb. The regenerated sector SHA-256 is aa716c71f373c04fbcea13aff0e3d41f33ceb1d5f0d2cdcbbb4b1cac471b20e0. See the [machine-readable comparison](helioviewer_reconstruction_verification.json).

The upstream frozen peak finder was also run on the separate Level-0 total- and polarized-brightness maps. It recovered the same primary peak radii for events 9 and 12. The node-order and expansion programs subsequently use those fixed nominal radii. The panel-b probabilities rely on both brightness products. The Helioviewer total-brightness sector alone cannot supply that test.

Suggested sentence for the manuscript: “Reconstructing the archived COR1 input from Helioviewer JPEG2000 images reproduced all eight arrays of the original transport map exactly.”
