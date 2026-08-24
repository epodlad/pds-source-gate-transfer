# Physics-aware expanding-front sensitivity: events #9 and #12

This specification is fixed after the strict node-activation test and is
therefore an explicitly labelled physics-motivated sensitivity analysis, not
an independent preregistered confirmation.

## Scientific question

Does the already frozen #9 -> #12 X/diamond-like pair remain compatible with
outward pattern transport when the timing inference accounts for radial
dilution, projected streamer expansion, front broadening, and the finite
COR1/COR2 cadence?

## Fixed targets and image geometry

- Event #9: 2008-01-13 15:37:30 UT; nodes 1.775, 2.125, 2.500, 2.825 R_sun;
  frozen 100 -> 30 km/s path with transition at 2.20 R_sun.
- Event #12: 2008-01-13 21:37:30 UT; nodes 1.950, 2.350, 2.700 R_sun;
  frozen 100 -> 25 km/s path with transition at 2.90 R_sun.
- COR2 anchor at 3.0 R_sun and pair separation 360 min.
- Previously selected PA path, node radii, line vertices and line slopes.
- Representations: primary BFF, base60, and nrgf60.
- No event, node, PA, slope, vertex, propagation family, speed, break radius,
  or pair separation is reselected.

## Expansion and timing uncertainty

1. The projected streamer area is estimated from the calibrated pB transverse
   second-moment width along the frozen event path:
   `A(r) proportional to [r * sigma_PA(r)]^2`.
2. If a node lacks a valid measured area, use the spherical fallback
   `A(r)/A(r0) = [r/r0]^2`.
3. The timing-broadening kernel is
   `sigma_t(r) = clip[15 min * sqrt(A(r)/A(r0)), 15 min, 30 min]`.
   It represents cadence plus projected expansion/front broadening; it is not
   interpreted as a measured thermodynamic pulse width.
4. The COR2 anchor time is marginalized uniformly over -15, 0, +15 min.
5. Activation offsets are sampled from -60 to +60 min in 15-min steps.

## Product timing posterior

At every frozen node and for each representation:

- evaluate the fixed opposite-slope X response versus activation offset;
- robust-standardize the nine-point tB and pB curves separately;
- convert each standardized curve to a unit-temperature softmax probability;
- convolve each probability with the node's expansion kernel;
- combine tB and pB symmetrically using their geometric mean and renormalize.

Every node posterior is normalized independently. Thus no equality of
amplitude with radius is imposed, and spherical/super-radial brightness
dilution is treated as a nuisance rather than evidence against transport.

## Exact sequence probability

Enumerate every combination of node offsets (9^4 for event #9 and 9^3 for
event #12) and the three COR2 anchor jitters. For every combination calculate:

- absolute activation times;
- all adjacent-node speeds;
- the final-node-to-3-R_sun speed;
- whether times are strictly outward ordered;
- whether every speed is within the unchanged physical band 25--300 km/s;
- whether the residual node phases fit within a 30-min range.

Report posterior probabilities for outward ordering, full transport-band
compatibility, and common residual phase. No hard decision is made from one
peak bin.

## Shifted-pair reference

- Repeat the complete calculation at the already fixed 131 control pairs,
  preserving the 360-min separation.
- A representation-level control pair exceeds the real pair when both its
  event-like transport probabilities equal or exceed the corresponding real
  event probabilities.
- The primary control statistic requires exceedance in at least two of the
  three representations.
- Empirical p values use `(k+1)/(N+1)` and are reported raw, without BH.

## Interpretation rule

An event is descriptively expanding-transport compatible in a representation
when its posterior transport probability is at least 0.50. The pair is called
processing-robust expanding-transport compatible only when both events meet
that threshold in at least two representations and the primary >=2-of-3
shifted-pair raw p is below 0.05.

If transport probability is at least 0.50 but common residual-phase
probability is below 0.50, label the pattern moving/reforming rather than a
single rigid front. If ordering probability is at least 0.50 but transport
probability is below 0.50, label it ordered but speed-indeterminate.

This analysis cannot establish density compression, Rankine--Hugoniot
consistency, Mach number, or slow/fast/compound MHD branch. Formal
SECCHI_PREP Level-1 confirmation remains required for a shock claim.
