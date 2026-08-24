# Frozen representation-robustness test for events #9 and #12

This protocol is written before rendering or inspecting the new base-difference
and NRGF maps.

## Fixed events

- Event #9: 2008-01-13 15:37:30 UT; pB deceleration 100 -> 30 km/s with break
  at 2.20 R_sun; nodes 1.775, 2.125, 2.500, 2.825 R_sun.
- Event #12: 2008-01-13 21:37:30 UT; tB deceleration 100 -> 25 km/s with break
  at 2.90 R_sun; nodes 1.950, 2.350, 2.700 R_sun.
- Ordered-pair separation: exactly 360 min.

No additional global time shift, height-dependent phase shift, PA recentering,
node selection, or change of the 360-min separation is allowed.

## Three frozen representations

The input remains the existing bias/exposure-normalized COR1-A Level-0
tangential fitpol sector cube.  The following alternatives are constructed
independently for tB and pB:

1. `base60`: image minus the image 60 min earlier.
2. `base120`: image minus the image 120 min earlier.
3. `nrgf60`: NRGF at each time/radius, defined by subtraction of the PA mean
   and division by the PA standard deviation, minus the NRGF image 60 min
   earlier.

Every difference cube is robustly standardized in time at each radius/PA pixel.
The event map samples the exact COR2 anchor plus the frozen propagation delay
at every radius, using only linear interpolation on the 15-min cadence.

## Frozen 2-D statistic

- PA offsets: -12 to +12 degrees around the previously traced static pB ray.
- X arms: |PA offset|=2--8 degrees.
- Vertex tolerance: +/-0.075 R_sun.
- Absolute slope grid: 0.025, 0.0375, 0.0500, 0.0625, 0.0750
  R_sun/degree.
- Parallel sidebands: +/-0.10 R_sun.
- At a node, both opposite diagonal contrasts must occur in both tB and pB;
  the coherent node value is the minimum of the four required responses.
- Event statistic in each representation: median coherent response across all
  frozen nodes.  Positive-node counts are reported independently.

## Ordered-pair null

The #9-like template is evaluated at an early control anchor and the #12-like
template exactly 360 min later.  The pair is shifted through every valid
30-min control position.

- Representation-specific exceedance: both the early and late control scores
  equal or exceed the corresponding real-event scores.
- Primary robustness exceedance: the shifted pair exceeds the real pair in at
  least two of the three frozen representations.
- Raw empirical p: `(exceedance count + 1)/(valid shifted pairs + 1)`.

## Interpretation rule

Representation robustness requires primary pair p<0.05 and, for each real
event, at least 75% positive coherent nodes in at least two of the three
representations.  This is an independent-processing check, not an independent
data set and not SECCHI_PREP Level-1 validation.  All p values are raw; no BH
adjustment is used.
