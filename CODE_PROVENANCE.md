# Code provenance

The scientific scripts in `analysis/frozen_programs/` were copied unchanged from the recovered package `PDS_Streamer_ApJ_Source_and_Code.zip`.

For this release-candidate cleanup:

- scientific algorithms were not rewritten;
- numerical constants, thresholds, frozen event selections, and model grids were not altered;
- manuscript drafts, editor/referee working notes, submission checklists, compiled manuscript PDFs, old figure PDFs, and Overleaf-specific instructions were excluded;
- public-release documentation and validation files were added;
- `Pillow` was added to `requirements.txt` because one frozen analysis script imports `PIL`.

The purpose of this cleanup is to make the archive suitable for GitHub/Zenodo without changing the frozen scientific analysis.
