# Code provenance

The scientific scripts in `analysis/frozen_programs/` represent the frozen analysis programs used for the released study.

The scientific algorithms, numerical constants, thresholds, event selections, model grids, and statistical logic are preserved unchanged in this release.

The repository also includes compact derived data products, analysis specifications, dependency information, and validation tools needed to document and inspect the released analysis.

File integrity can be checked using the SHA-256 hashes provided in `checksums.sha256`.

The Python dependencies required by the frozen programs are listed in `requirements.txt`. `Pillow` is included because one analysis script imports the `PIL` image-processing library.

This release is intended to preserve the analysis state associated with the reported scientific results and to support their traceability and reproducibility.