# Frozen test: Habbal standing transition and MHD shock branch

Date frozen: 2026-08-21.  All probabilities are raw empirical p values; no
Benjamini--Hochberg correction is used.

## Fixed events and representations

1. Primary ensemble: the 12 previously fixed COR2-A anchors from
   2008-01-11 14:37:30 through 2008-01-13 21:37:30 UT at PA=174.5 deg.
2. Fixed January morphology candidates: event #9 (2008-01-13 15:37:30 UT)
   and event #12 (2008-01-13 21:37:30 UT).  No additional January event may
   be promoted after this protocol is frozen.
3. Independent event: 2008-07-26 23:30 UT on the previously traced curved
   COR1 path ending at PA=107 deg in COR2.
4. Primary calibrated brightness observables are signed pB and total-B.
   Level-0-derived maps are diagnostic only.

## Habbal standing-transition hypothesis

A stationary transition is identified only if all observable necessary
conditions pass:

1. **Predictive nonlinear path:** a one-break acceleration or deceleration
   model has held-out raw p<0.05 and a raw gain p<0.05 over constant speed.
2. **Repeatable height:** the selected break radii in the two six-event folds
   agree within 0.20 R_sun.
3. **Event support:** at least four of the 12 events have individually
   identified ridges (raw p<0.05), and at least three of those have a
   height-dependent phase kink with raw p<0.05 within +/-0.15 R_sun of one
   common transition radius.
4. **Compression support:** both fixed January candidates have at least 75%
   valid pB nodes and a median apparent column compression 1<C_col<4 at the
   same transition family.
5. **Nozzle proxy:** an independently measured time-median streamer-width
   expansion proxy has a reproducible localized gradient/extremum within
   +/-0.20 R_sun of the transition radius in both total-B and pB.  Width is a
   projected emissivity proxy, not a direct magnetic flux-tube area.

Failure of any condition rejects a detected *stationary* Habbal transition
for this data set.  It does not reject a moving/reforming shock, a
time-dependent throat, reconnection release, or a visibility threshold.

## Slow/fast/compound MHD branch screen

For temperature T, density n_e, magnetic field B, shock-normal angle
theta_Bn, upstream normal flow u_n, and front speed V_sh,n, compute

    c_f,s^2 = 0.5 * [(v_A^2 + c_sound^2)
              +/- sqrt((v_A^2 + c_sound^2)^2
                       - 4 v_A^2 c_sound^2 cos^2(theta_Bn))]

and U_n=abs(u_n-V_sh,n).  The screen uses predeclared broad coronal ranges,
not fitted best values:

- T = 0.8--1.8 MK;
- electron density at 2--3 R_sun = 1e5--1e7 cm^-3;
- B = 0.01--1.0 G;
- theta_Bn = 0--85 deg;
- u_n = 20--400 km/s;
- V_sh,n = 0--400 km/s.

The ranges are a sensitivity envelope.  They cannot replace event-specific
spectroscopy or a magnetic-field reconstruction.

Branch-compatible parameter samples obey:

- **slow:** U_n>c_slow, U_n<v_A,n and U_n<c_fast;
- **fast:** U_n>c_fast;
- **sub-slow:** U_n<c_slow;
- **intermediate/compound-compatible:** U_n crosses the normal Alfvén
  characteristic or the imaging signatures cannot be assigned uniquely to a
  pure slow/fast branch.

Remote identification requires more than a characteristic crossing:

- density/temperature increase for either compressive shock;
- tangential magnetic-field decrease for slow/switch-off shocks;
- tangential magnetic-field increase for fast shocks;
- Rankine--Hugoniot consistency in the front frame.

Because pB measures line-of-sight electron column density and provides no
direct B or T, the present data can yield branch compatibility or exclusion,
not a confirmed slow/fast shock, unless the imaging-only standing-transition
and compression criteria also pass and external B/T/u constraints are added.

## Decisions

- A low raw p for a combined ridge-plus-morphology statistic cannot rescue a
  nonsignificant X statistic or invalid compression.
- Pattern-alignment speed is never substituted silently for plasma bulk speed.
- The January and July events are reported separately; their probabilities
  are not pooled.

## Post-render quality clarification

Added after the first diagnostic render, before interpreting the result: a
nominal streamer-width gradient peak is not accepted as a nozzle proxy when
the FWHM is measurable over less than 75% of the predeclared 1.6--3.0 R_sun
interval, or when the maximum lies within 0.075 R_sun of either search
boundary.  Both the nominal location and the failed quality flag are retained
in the output so this clarification cannot manufacture a positive result.

## Supplementary fixed-edge diagnostic

After the primary necessary-condition screen, a supplementary direct
stationarity diagnostic is added and labelled as such.  Each time-dependent
brightness profile along the frozen path is detrended by a cubic radial
profile, smoothed over 0.05 R_sun, differentiated, and robustly standardized.
The null circularly shifts each whole gradient profile in radius.  The shared
COR1 inner-edge response is excluded by restricting the solar search to
1.9--2.9 R_sun.  Raw p values are reported both for the maximum in that
interval and at the already frozen 2.10, 2.80, or 2.90 R_sun candidates.  This
diagnostic is necessary but not sufficient: a fixed emissivity or
instrumental edge is not by itself a shock and cannot override failed
compression or Rankine--Hugoniot conditions.
