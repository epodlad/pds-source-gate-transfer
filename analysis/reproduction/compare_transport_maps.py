"""Compare a regenerated COR1 transport map to its archive on the common finite pixels."""
import argparse
import json
from pathlib import Path

import numpy as np


def compare_arrays(reference, reconstruction):
    reference = np.asarray(reference, dtype=np.float64)
    reconstruction = np.asarray(reconstruction, dtype=np.float64)
    if reference.shape != reconstruction.shape:
        raise ValueError(f"Shape mismatch: {reference.shape} vs {reconstruction.shape}")

    valid_ref = np.isfinite(reference)
    valid_new = np.isfinite(reconstruction)
    common = valid_ref & valid_new
    a, b = reference[common], reconstruction[common]
    diff = b - a
    corr = None
    if a.size > 1 and np.ptp(a) > 0 and np.ptp(b) > 0:
        corr = float(np.corrcoef(a, b)[0, 1])
    return {
        "shape": list(reference.shape),
        "finite_pairs": int(a.size),
        "finite_mask_equal": bool(np.array_equal(valid_ref, valid_new)),
        "arrays_equal_including_nan": bool(np.array_equal(reference, reconstruction, equal_nan=True)),
        "pearson_r": corr,
        "median_signed_difference": float(np.median(diff)) if a.size else None,
        "max_abs_difference": float(np.max(np.abs(diff))) if a.size else None,
    }


def compare_maps(reference_path, reconstruction_path):
    with np.load(reference_path, allow_pickle=False) as archived, np.load(reconstruction_path, allow_pickle=False) as regenerated:
        if set(archived.files) != set(regenerated.files):
            raise ValueError("Input files have different array names")
        return {key: compare_arrays(archived[key], regenerated[key]) for key in archived.files}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archived", type=Path)
    parser.add_argument("reconstructed", type=Path)
    args = parser.parse_args()
    print(json.dumps(compare_maps(args.archived, args.reconstructed), indent=2, allow_nan=False))
