# Frozen coronal-physics reassessment of all 12 COR2-anchored events

## Purpose

This is a post-extraction reassessment.  It does not reprocess images, move
ridges, select new nodes, or recompute p values.  It changes only the physical
interpretation of already frozen measurements after allowing for an expanding,
turbulent coronal streamer rather than a rigid constant-width template.

The original strict tables and raw p values remain the audit baseline.

## Frozen inputs

- Twelve strong COR2 anchors from 2008-01-11--14.
- Strict matched tB/pB radial-peak counts and cell-statistic raw p values.
- Best and common-ridge raw p values.
- Frozen best kinematic family and speeds.
- tB/pB phase-offset difference and phase-height diagnostics.
- Previously measured projected streamer-width/area diagnostics where valid.
- The separate spatially expanding-ridge pilot for events #9 and #12, including
  the fixed 360-min pair reference with 131 shifted pairs.
- All p values are raw; no BH correction is introduced in this reassessment.

## Coronal corrections

1. Constant cross-section is not required.  Where a valid pB width proxy exists,
   projected A(r) is used qualitatively; otherwise spherical A proportional to
   r^2 is the conservative fallback.
2. Constant amplitude is not required.  Density/brightness dilution is
   compatible with a broad continuity prior n v A approximately constant and
   line-of-sight Thomson weighting.
3. Exact phase preservation is not required.  A tB/pB offset difference of
   0--30 min is called compatible with the adopted timing uncertainty; 45--60
   min is marginal; larger differences remain unresolved/incoherent.
4. A fitted speed of 20, 25, or 300 km/s is not treated as a hard boundary
   failure.  The frozen value is retained and interpreted as lying near the
   edge of a broad uncertainty interval.
5. None of these corrections creates a new ridge or a new statistical
   detection.

## Event categories

The categories are descriptive evidence tiers, not additional significance
tests.

### A. Expansion-stable ordered candidate

Assigned only to event #9 or #12 when the separately frozen measured-width
pilot gives P(outward order) >= 0.50 in at least two of three representations
and the primary fixed-pair reference has raw p < 0.05.

### B. Multi-ridge coronal-compatible morphology

At least three strict matched tB/pB peaks and at least one pre-existing raw
screening value <= 0.10 among strict cell-statistic p, best-ridge p, or
common-ridge p.  The 0.10 level is an exploratory screening label, not a claim
of statistical detection.

### C. Weak multi-ridge morphology

At least three strict matched tB/pB peaks, but none of the three pre-existing
screening p values is <= 0.10.

### D. Isolated or incomplete chain

Fewer than three strict matched tB/pB peaks.  A single ridge or phase hint does
not establish a cross-height PDS chain.

Category A takes precedence over B--D.  The correction can move a strict
failure to a compatibility/inconclusive label, but cannot promote any event to
a confirmed shock or confirmed clock.

## Separate physical axes

Every event is reported on five axes:

- radial morphology;
- kinematic/transport compatibility;
- phase or source-clock compatibility;
- geometry/nozzle compatibility;
- shock-specific evidence.

The final shock-specific value is `not established` for all events because no
event has event-specific B, T, bulk normal velocity, density jump, and valid
Rankine--Hugoniot closure.  The global 120-min low-coronal clock also remains
not established event by event.
