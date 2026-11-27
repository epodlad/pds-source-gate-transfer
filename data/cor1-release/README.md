# STEREO-A streamer maps and control intervals

The [Zenodo data version](https://doi.org/10.5281/zenodo.24604508) contains the 17 requested cubes and derived maps, including the separately named Helioviewer reconstruction. [Validation results](validation.md) · [GitHub inventory](github_manifest.csv) · [SHA-256 for every GitHub asset](github_checksums.sha256) · [SHA-256 for every Zenodo file](zenodo_checksums.sha256). The GitHub release retains the original attachments and adds cor1_helioviewer_reconstruction_audit.zip, with the comparison code, JSON and supplement checksums. The two Helioviewer sector filenames on GitHub have identical bytes.

[Download the release assets](https://github.com/epodlad/pds-source-gate-transfer/releases/tag/data-cubes-20080111-14-v1). The [file inventory](manifest.csv) lists the scientific contents, producer, reader, size and SHA-256. The inventory inside the revision ZIP lists its members. The source observations are from the STEREO/SECCHI team. [COR1 mission data](https://doi.org/10.48322/9fer-5m90) · [EUVI mission data](https://doi.org/10.48322/260z-fs61)

From a fresh checkout of the tagged repository, install the requirements and run:

```sh
gh release download data-cubes-20080111-14-v1 --dir assets
(cd assets && sha256sum --check checksums.sha256)
python -m pip install -r requirements.txt
python analysis/reproduction/replay.py assets ../cor1-quick
python analysis/reproduction/replay.py assets ../cor1-all --all
python analysis/reproduction/pair_counts.py data/spatial_expanding_pair_null.csv assets/pds_20080113_two_event_representation_pair_null.csv
python analysis/reproduction/width_aeff.py assets ../width-replay
unzip -q assets/revision_code_and_tables.zip -d revision
```

The quick run rebuilds the stationary-transition and node-order JSONs. The full run also generates both X-front maps, the alternative representations, width-law maps, control tables and figures. New output directories preserve previous runs. The frozen programs are used unchanged. Comparisons allow absolute numerical differences of 1e-10 and relative differences of 1e-8. Each array's ordering and each categorical result are preserved in the JSON comparison.

| Observable | Data, command or output |
| --- | --- |
| EUVI 171 and 195 Å brightness along the eastern path, 1.10–1.20 solar radii | [East-limb cubes and average series](https://github.com/epodlad/pds-source-gate-transfer/releases/tag/euvi-20080111-14-east-limb-v1) and [exposure inventory](https://github.com/epodlad/pds-source-gate-transfer/releases/tag/euvi-20080111-14-reconstruction-v1). The asset stroke_171_timeseries.csv contains the nine narrow and nine wide extractions. stroke_171_spectra.csv gives their band peaks. fits_path_timeseries.csv and fits_path_spectra.csv contain both wavelengths. Use the revision EUVI code with the source inputs and FITS archives named in fetch_fits.py. |
| Historical narrow-strip period and phase | In euvi171_lower_resonator_results.json, select metrics entry brightness_1.10_1.20R. The archived period is 128.2143 min, delay 107.25375 min, R=0.066414, circular-shift p=0.911111 (21 peaks). The exact original numerical phase-driving series was unavailable in the located archive. The separately labelled digitized_archived_filtered_curves.csv preserves points recovered from the original plotted curves. The image reconstruction and nine shifted strips have their own measurements. |
| Projected speed in COR1 L0/seq | In revision/cor1_l0 run the commands below. The three peaks in both ideal-polarizer proxies occur at 18:20, 18:50 and 19:20 UT on 13 January, at 2.15, 2.20 and 2.25 solar radii. The projected endpoint speed is 19.325 km/s conditional on one front. The source manifest covers 1,056 exposures in 352 triplets. expected/front_summary.json stores the original result. |
| Independent COR2 detection | In revision/AAS80409_revision/cor2_recovery run the commands below. Inspect results/matching_counts.csv for the separate comparison. Only the 21 original event times and sampled scores are available. The historical continuous COR2 detection score was not recovered. |
| Candidate X fronts and alternative image reductions | The full run writes event09_2d_xfront, two_candidate_2d_xfront and two_event_representation_robustness below results/pds_20080111_14_nonlinear_level0. Individual frozen-template raw p: 0.1389 and 0.1458. The original six-hour pair test gives 0/131. The alternative X-contrast test gives 1/131, p=0.0152. |
| Expansion sensitivity, Figure fig8_expansion_sensitivity panel b | pair_counts.py recomputes the panel from the 131 rows in data/spatial_expanding_pair_null.csv. For fixed, spherical (width proportional to r) and measured width the primary counts are 0, 0, 2. Raw p=(k+1)/132. The three reductions are BFF, base60 and NRGF60. The X-contrast test instead uses base60, base120 and NRGF60 and compares median coherent X contrast. Both tests use exactly the same six-hour control pairs. |
| Stationary transition and activation times | See standing_transition.md and node_times.md in the assets. The primary model-aligned node peaks and the native exposures answer different questions. A score-null p conditional on failed quality checks cannot establish propagation. |
| Isothermal effective-area sensitivity | width_aeff.py is separate exploratory code. It reports FWHM across 21 sampled radii from 1.80 to 2.30 solar radii using half-windows of 12°, 16° and 20°. Projected linear width is r times the angular FWHM in radians. Assuming a circular section, A=πW²/4. At 1 MK use c_T²=k_B T/(0.6 m_p). The geometric factor and c_T are model inputs. The minimum tabulated A_eff=A exp(GM/(r c_T²)) is among valid bins only. Neither the width nor A_eff has a sampled minimum bracketed by valid immediate neighbors in this test. |
| WISPR ray intersections | In revision/AAS80409_revision/nozzle_test/fan_checks run python recheck.py. The program recalculates offsets and spreads from archived intersection coordinates and vertex geometry. It does not reselect rays in the original exposures. |

COR1 L0/seq uses the public mission archive. The FITS can be downloaded to a local computer with enough storage.

```sh
cd revision/cor1_l0
python fetch_cor1.py inputs --manifest raw_manifest.csv
python reduce_native.py inputs raw_manifest.csv path_cor1.csv checks
python compare_l0.py
python extra_checks.py
python overview.py . figures
python plot_level0.py
```

For the independent COR2 score, start in revision/AAS80409_revision/cor2_recovery. The scripts reconstruct the independently defined score from public images.

```sh
mkdir -p rerun
python archive.py fits-cache rerun/frames.csv
python archive.py fits-cache rerun/padding.csv --padding
python extract.py fits-cache rerun/frames.csv rerun/padding.csv rerun/sector.npz
python compare.py rerun/sector.npz ../../COR1_height_recovery/source/ridge_onset_results.json rerun/frames.csv rerun/results
```

### Reconstruction of the archived Helioviewer input

The separately named cor1a_hv_polar_sector_15min_reconstruction.npz contains processed, photometrically uncalibrated 8-bit COR1-A total brightness from Helioviewer JPEG2000, source ID 28. The frozen program makes the sector cube from the archive's original times and sampling. The archived transport maps retain their original bytes. The original sector-cube checksum is unavailable.

To run from the reconstructed input, start from a repository checkout with the release downloaded to assets.

```sh
mkdir -p hv-rerun/results/pds_20080111_14_nonlinear
cp analysis/frozen_programs/analyze_pds_nonlinear_transport_20080111_14.py hv-rerun/
cp assets/cor1a_hv_polar_sector_15min_reconstruction.npz hv-rerun/results/pds_20080111_14_nonlinear/cor1a_hv_polar_sector_15min.npz
python hv-rerun/analyze_pds_nonlinear_transport_20080111_14.py
python assets/compare_transport_maps.py assets/pds_20080111_14_nonlinear_transport_maps.npz hv-rerun/results/pds_20080111_14_nonlinear/pds_20080111_14_nonlinear_transport_maps.npz
```

Each of the eight derived map arrays reproduces all finite values and missing pixels. The accompanying helioviewer_reconstruction_verification.json gives Pearson correlations, median signed differences and maximum absolute differences. The event-9 and event-12 node radii and Figure 8b require a separate Level-0 analysis with both total and polarized brightness. The same JSON lists the Level-0 rerun results.
