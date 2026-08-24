# Spatially expanding-ridge pilot: COR1-A events #9 and #12

## Status and purpose

This specification is frozen before calculation.  It is a targeted sensitivity
test, not a new discovery scan.  It asks whether allowing an already identified
X/ridge pattern to broaden spatially changes the event #9 -> #12 conclusion.

The original fixed-template result and the previous expansion-aware timing
result remain unchanged and are reported beside this pilot.

## Frozen events, products, and reference

- Event #9: 2008-01-13 15:37:30 UT; nodes 1.775, 2.125, 2.500, 2.825 R_sun;
  decelerating 100 -> 30 km/s model with change radius 2.20 R_sun.
- Event #12: 2008-01-13 21:37:30 UT; nodes 1.950, 2.350, 2.700 R_sun;
  decelerating 100 -> 25 km/s model with change radius 2.90 R_sun.
- Representations: temporal BFF, base-difference 60 min, and NRGF
  base-difference 60 min.
- The pair separation remains exactly 360 min.
- The same 131 shifted early/late pair positions are used as the empirical
  reference.  All p values are raw; no BH correction is applied.

## Spatial matched filters

The existing opposite-slope template is retained.  Only its radial thickness
is changed.  The inner one-sigma thickness is fixed at 0.025 R_sun, one radial
sample of the working map.  At radius r:

1. `fixed`: sigma_r(r) = 0.025 R_sun.
2. `spherical`: sigma_r(r) = 0.025 (r/r_inner) R_sun, corresponding to
   A(r)/A_inner = (r/r_inner)^2.
3. `measured`: sigma_r(r) = 0.025 sqrt[A(r)/A_inner] R_sun, where
   A is the pB projected-width proxy [r sigma_PA]^2 sampled on the frozen
   travel-time path.  Missing measurements fall back to the spherical law.

The measured node area ratios are linearly interpolated in radius and are not
forced to be monotonic.  Thickness is clipped to 0.025--0.075 R_sun.

Each arm response is a Gaussian-weighted mean across +/-3 sigma_r normal to the
arm (radial-normal approximation).  Background sidebands use the same kernel at
both signs and are centred at max(0.10 R_sun, 3 sigma_r) from the arm.  The arm
PA range, slope grid, vertex range, nodes, and tB/pB coherent-minimum rule remain
the same as in the frozen X-front test.

For each event, representation, and width law, the vertex/slope combination is
selected only at the real zero-offset event map from the already frozen grids.
That geometry is then held fixed for all timing offsets and all 131 controls.
This makes the pilot exploratory but prevents re-fitting each control.

## Timing and transport

- Offset grid: -60 to +60 min in 15-min steps.
- Per-node tB and pB response curves are robust-standardized independently,
  converted to softmax probabilities, and combined by geometric mean.
- Timing broadening is sigma_t = clip[15 sqrt(A/A_inner), 15, 30] min for the
  spherical and measured models; sigma_t=15 min for the fixed model.
- Exact enumeration of all node-offset combinations is used.
- Report P(outward order), P(all segment speeds in 25--300 km/s), and
  P(common residual phase within 30 min), with COR2 anchor jitter +/-15 min.

The 25--300 km/s range is used only as a probability integral, never as a hard
event rejection.

## Decision and interpretation

The pilot changes the physical conclusion only if the measured-width result is
stable in at least two of three representations and is not reproduced by the
same shifted-pair reference.  The comparison of fixed, spherical, and measured
filters is treated as a sensitivity bracket, not as three independent tests.

This analysis can support an expanding or reforming morphology/transport
interpretation.  It cannot establish a density jump, Rankine--Hugoniot closure,
magnetosonic Mach number, or an MHD shock branch.
