"""Compare a clean rebuild with the independent run from the released FITS."""
import csv
import json
import sys
from pathlib import Path
import numpy as np

here=Path(__file__).resolve().parent
out=Path(sys.argv[1])
actual=json.loads((out/'full_frame_results.json').read_text())
expected=json.loads((here/'full_frame_results.json').read_text())
for wave in ('171','195'):
    ref=list(csv.DictReader((here/f'east_brightness_{wave}.csv').open()))
    now=list(csv.DictReader((out/f'east_brightness_{wave}.csv').open()))
    assert len(ref)==len(now)==actual[wave]['frames']==358
    assert all((r['utc'],r['filename'])==(n['utc'],n['filename']) for r,n in zip(ref,now))
    np.testing.assert_allclose([float(r['mean_dn_per_s']) for r in ref],
                               [float(n['mean_dn_per_s']) for n in now],rtol=1e-9,atol=1e-9)
    assert max(actual[wave]['geometry'],key=actual[wave]['geometry'].get)=='-1*angle+90'
    np.testing.assert_allclose(actual[wave]['geometry']['-1*angle+90'],expected[wave]['geometry']['-1*angle+90'],atol=.001)
    with np.load(out/f'euvi{wave}_east_limb_diagnostic.npz') as c:
        assert c['intensity_dn_per_s'].shape==(358,71,281)
        assert c['filename'].tolist()==[r['filename'] for r in ref]
        assert np.isfinite(c['intensity_dn_per_s'][:,10:21,150:180]).all()
print('Verified both series, east-limb registration, exposure names and cube dimensions.')
