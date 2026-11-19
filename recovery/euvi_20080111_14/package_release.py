"""Verify the exact FITS inventory and assemble per-day release assets."""
import argparse
import csv
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import shutil
import zipfile

import numpy as np


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def package(root, manifest, destination):
    with manifest.open(newline="") as f:
        expected = list(csv.DictReader(f))
    results = json.loads((root / "progress.json").read_text())
    by_name = {r["filename"]: r for r in expected}
    if len(expected) != 716 or len(by_name) != 716 or len(results) != 716:
        raise ValueError("Expected 716 distinct exposures.")
    if {r["filename"] for r in results} != by_name.keys():
        raise ValueError("Frame sets differ.")
    grouped = defaultdict(list)
    for row in results:
        name = row["filename"]
        ref = by_name[name]
        if Path(name).name != name or row["status"] != "ok":
            raise ValueError(f"Invalid frame: {name}")
        if (str(row["channel"]), row["date_obs"], row["url"], row["sha256"]) != (
            ref["channel"], ref["date_obs"], ref["url"], ref["sha256"]
        ):
            raise ValueError(f"Frame metadata differ: {name}")
        source = root / "raw" / name
        if source.stat().st_size != int(ref["size_bytes"]) or digest(source) != ref["sha256"]:
            raise ValueError(f"FITS hash or size differs: {name}")
        grouped[(name[:8], row["channel"])].append(source)
    destination.mkdir(parents=True, exist_ok=False)
    for (day, channel), frames in sorted(grouped.items()):
        archive = destination / f"stereo_a_euvi_{channel}_{day}_fits.zip"
        with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED, compresslevel=1) as output:
            for frame in sorted(frames):
                output.write(frame, frame.name)
        if archive.stat().st_size >= 2**31:
            raise ValueError(f"Asset exceeds 2 GiB: {archive}")
    for channel in (171, 195):
        cube = root / f"euvi{channel}_reconstructed_sector.npz"
        ordered = sorted((r for r in expected if int(r["channel"]) == channel), key=lambda r: r["date_obs"])
        with np.load(cube, allow_pickle=False) as z:
            if len(ordered) != 358 or z["intensity_dn_per_s"].shape != (358, 101, 241):
                raise ValueError(f"Unexpected cube dimensions: {channel}")
            if list(z["filenames"]) != [r["filename"] for r in ordered]:
                raise ValueError(f"Cube frame order differs: {channel}")
            if list(z["time_utc"]) != [r["date_obs"] for r in ordered]:
                raise ValueError(f"Cube time axis differs: {channel}")
        shutil.copy2(cube, destination)
    shutil.copy2(root / "validation.json", destination)
    shutil.copy2(manifest, destination / "verified_frames.csv")
    files = sorted(destination.iterdir())
    (destination / "SHA256SUMS.txt").write_text("".join(f"{digest(f)}  {f.name}\n" for f in files))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    package(args.root, args.manifest, args.destination)
