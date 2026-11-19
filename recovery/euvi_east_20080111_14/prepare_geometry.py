"""Extract archival static-image ranks on the fixed comparison grid from the original plots."""
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from PIL import Image
root=Path(__file__).resolve().parent
for w,x0,x1,barx in [(171,122.5,1914.5,2040),(195,138.5,1880.5,2006)]:
    a=np.array(Image.open(root/f'euvi{w}_lower_resonator_test.png').convert('RGB'))
    assert a.shape[:2]==(2700,2160)
    radii=np.arange(1.05,1.411,.02)
    angles=np.arange(106.,237.,1.)
    tick_y=[530.5,450.5,371.5,292.5,212.5,133.5,54.5]
    yy=np.rint(np.polyval(np.polyfit(np.arange(11,18)/10,tick_y,1),radii)).astype(int)
    xx=np.rint(x0+(angles-100)*(x1-x0)/140).astype(int)
    palette=cKDTree(a[55:600,barx,:].astype(float))
    error,row=palette.query(a[yy[:,None],xx[None,:]].astype(float))
    np.savez_compressed(root/f'archived_image_{w}.npz',radii=radii,angles=angles,valid=(error<4),ranks=-row)
