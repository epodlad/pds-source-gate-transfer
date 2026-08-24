# Frozen geometry-dependent Habbal and dynamic cusp-gate screen

Date frozen: 2026-08-21, before running the new geometry/event-lock tests.
All reported probabilities are raw empirical probabilities.  No BH correction
is applied.

## Fixed data and events

1. The same 12 January 2008 COR2 anchors, the same traced COR1 path, and the
   same calibrated Level-0-derived total-B and signed-pB cubes used in the
   frozen event-phase analysis are retained.
2. The independent 2008-07-26 23:30 UT event and both calibrated July
   backgrounds (080723 and 080802) are retained for sensitivity checks.
3. PDS are treated as advected density structures.  Recurrence is not assumed
   to require an oscillatory clock.
4. Event times are projected back from 3.0 R_sun with a fixed 200 km/s speed.
   Fixed 100 and 300 km/s projections are sensitivity tests.  Pattern speed is
   not treated as plasma bulk speed.

## Projected geometry observables

At every time and radius, a transverse profile within +/-12 degrees of the
frozen path is formed.  The baseline and noise come from the +/-8--12 degree
sidebands.  Positive excess over the baseline supplies:

- transverse centroid displacement (axis proxy);
- second-moment angular width;
- integrated excess brightness (density-column proxy).

A point is retained only when its peak excess exceeds twice the sideband MAD.
The physical-width proxy is `r * sigma_theta`, and
`A_proxy=(r*sigma_theta)^2`.  This is a projected emissivity geometry, not a
magnetic flux-tube area.

The primary Habbal geometry metric is the maximum positive
`d ln(A_proxy)/dr` over 1.9--2.9 R_sun at the back-projected event time.  The
associated radius and the maximum absolute axis shear are secondary metrics.

## Test 1: geometry dependence of PDS visibility

Across the fixed 12 January events, Spearman correlations are computed between
the event geometry metrics and the already frozen common-path ridge score.
Two-sided label permutations give raw p values.  A max-|rho| permutation over
the declared geometry metrics is also reported as an omnibus diagnostic.

Strong geometry dependence requires:

1. positive correlation of rapid area divergence with the common-path ridge
   score at raw p<0.05 for the 200 km/s projection;
2. the same correlation sign at both 100 and 300 km/s, with p<0.10 in at
   least one sensitivity speed;
3. no reversal between total-B and pB where signed-pB geometry is measurable.

Failure means geometry control is not detected; it does not show that magnetic
geometry is physically irrelevant.

## Test 2: moving or reforming cusp/current-sheet gate

For log area, centroid displacement, and log integrated excess, an event
change is the difference between the post interval (+15 to +60 min) and the
pre interval (-60 to -15 min), evaluated along each fixed back-projected path.
The primary statistic is the median absolute standardized change across the
12 events as a function of radius.

The null shifts the complete 12-event train together by all allowed 15-min
offsets, excluding shifts within +/-180 min of the true train.  This preserves
event spacing and slow evolution.  Raw p values are reported at 2.10 and 2.90
R_sun and for the maximum over 1.9--2.9 R_sun; the scan p includes the radial
look-elsewhere effect.

Strong dynamic-gate support requires at least two of the three observables to
have scan p<0.05, with peak radii agreeing within 0.15 R_sun, and the result to
retain the same radius family in pB or a fixed-speed sensitivity run.

## Test 3: periodic geometry modulation

At the fixed radii 2.10, 2.50, 2.80, and 2.90 R_sun, each geometry time series
is linearly detrended and its maximum Fourier power in the predeclared
80--130 min band is compared with 4000 fitted AR(1) red-noise surrogates.
The maximum inside the band is the statistic, so its raw p already includes
the within-band frequency scan.  No across-radius or across-observable BH
correction is applied.

Periodic gate support requires a raw p<0.05 peak with periods agreeing within
15 min in integrated excess and at least one geometry observable, reproduced
in total-B/pB or in both July backgrounds.  Isolated single-channel peaks are
reported but not interpreted as a clock.

## MHD-branch interpretation

- A density-column pulse alone does not distinguish slow, fast, or compound
  transitions.
- Slow-shock identification additionally requires the slow-Mach crossing,
  sub-Alfvenic normal flow, tangential-field decrease, compression/heating,
  and Rankine--Hugoniot consistency.
- Fast-shock identification requires a fast-Mach crossing and tangential-field
  increase, in addition to compression/heating and jump consistency.
- Coupled density, width, and axis changes can support a reforming/compound
  gate phenomenology, but imaging without B, T, and bulk velocity cannot
  identify the MHD branch.
- Compact slow modes remain possible sub-streamer modulators; this analysis
  does not reinstate a full-streamer standing slow-mode claim or literal
  Laval/shock-diamond interpretation.

## Post-extraction quality clarification

The January moment maps retain only 46--51% of all time--radius points after
the frozen 2-MAD profile-quality gate.  Long-gap interpolation is therefore
forbidden.  An event geometry value may use linear interpolation only when
the two valid brackets are separated by at most 30 min; otherwise the nearest
valid value is used only when it lies within +/-30 min.  Dynamic pre/post
tests allow interpolation only across gaps of at most 30 min.  The stricter
periodic test still requires at least 70% direct coverage.  This clarification
was added after inspecting coverage and before accepting any positive result.

## User-requested single-event exploratory reading

After the frozen population tests were completed, the user explicitly asked
whether a physically interesting example could be retained without requiring
the full-series threshold.  The strongest already identified individual ridge
may therefore be described as a *single-event candidate*.  This post-primary
reading must retain the individual raw probabilities, product failures, and
compression caveat and cannot change the population decisions above.
