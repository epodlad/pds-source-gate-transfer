# EUVI reconstruction checkpoint (2026-11-18)

This is a new acquisition/reconstruction of STEREO-A EUVI observations from
2008-01-11 12:00 UTC through 2008-01-14 00:00 UTC. It does not replace the
frozen published products and has not reproduced the original EUVI p-values.

The original report records 357 frames at 171 A and 358 at 195 A, approximately
10 minutes apart. The original frame manifest, cube builder and statistical
script have not been recovered. This reconstruction selects archive `n4euA`
filenames at minutes 00/10/20/30/40/50 plus 30 seconds for 195 A, and minutes
01/11/21/31/41/51 at 00 seconds for 171 A. Every channel is checked from the
FITS header. Both candidate lists have 358 frames. No frame will be removed
merely to force agreement with the old count.

The code preserves original full-resolution FITS bytes and their SHA256 hashes.
It also constructs NEW polar sector cubes on an explicitly specified grid:
101 radii from 1.0 to 1.5 solar radii and 241 angles from 140 to 200 degrees.
The original sector grid is unknown. Angle convention is the legacy code
convention x=+r*sin(PA), y=+r*cos(PA), not an asserted standard solar PA.
Header WCS maps solar coordinates to detector pixels; interpolation is linear.
Off-detector samples are NaN. Brightness is divided only by EXPTIME. Bias,
flat-field, degradation correction, SECCHI_PREP, background removal, cosmic-ray
removal, time interpolation and gamma/noise filtering are NOT applied.
These cubes are diagnostic DN/s, not calibrated radiance or Level 1.

Run from a directory containing the inventory CSV:

```
python rebuild.py --inventory euvi_candidate_inventory.csv --output euvi_rebuilt --workers 8
```

Dependencies: Python, NumPy, SciPy, Astropy. Existing validated sector results
are resumed. Keep the full output directory while downloading. The resulting
`progress.json` lists every selected observation and any errors. All 358 frames
per channel are candidates until the original selection is recovered. A complete
file count alone is not scientific validation of geometry or source association.

Official archive:
https://stereo-ssc.nascom.nasa.gov/pub/ins_data/secchi/L0/a/img/euvi/
The mission describes this L0 directory as Level 0.5 FITS data.

This checkpoint contains recovery code and archive inventory. Full raw FITS and
rebuilt cubes must be distributed with a final manifest and explicit coverage;
they must not be described as already published here unless actually uploaded.
