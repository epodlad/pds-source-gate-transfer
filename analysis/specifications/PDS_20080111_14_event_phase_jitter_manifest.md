# Reproducibility package: event-dependent phase jitter

This package contains the event-by-event phase-jitter analysis for the 12
independently selected COR2 anchors from 11--14 January 2008.

## Main outputs

- `PDS_20080111_14_event_phase_jitter_no_BH.md` -- result, full 12-event
  table, physical interpretation and calibration limitation.
- `pds_20080111_14_event_phase_jitter_results.json` -- event metrics and
  family-level common-shift tests.
- `pds_20080111_14_event_phase_jitter_table.csv` -- flat event table.
- `pds_20080111_14_event_phase_curves.png` -- separate tB/pB phase curves for
  every event.
- `pds_20080111_14_event_phase_summary.png` -- phase heatmap, global offsets
  and individual raw p values.
- `pds_20080111_14_event_phase_jitter_curves.npz` -- numeric phase curves and
  peak amplitudes.

## Reproduction inputs

- `analyze_pds_event_phase_jitter_20080111_14.py` -- primary analysis.
- `analyze_pds_nonlinear_transport_20080111_14.py` -- shared frozen event and
  kinematic-model definitions.
- `pds_20080111_14_level0_nonlinear_transport_maps.npz` -- compact Level-0
  BFF time--height maps.
- `pds_20080111_14_level0_nonlinear_transport_results.json` -- frozen common
  reference model and prior held-out test.

All reported p values are raw empirical p values.  No BH adjustment is used.
The direct Level-0 polarization reduction is not an absolute SECCHI_PREP pB
calibration, so no density compression ratio is inferred.
