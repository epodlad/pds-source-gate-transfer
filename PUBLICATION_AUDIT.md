# Publication audit

Audit performed on the recovered `PDS_Streamer_ApJ_Source_and_Code.zip` before preparing the public software release.

## Author metadata

The public software metadata identify **Olena Podladchikova** as the author. No email address, home address, account identifier, ORCID guess, or other unnecessary personal metadata are included in the software release.

## Material removed from the public-code candidate

The recovered working package also contained manuscript-production and internal workflow material. These files were not copied into this release candidate because they are not needed to understand or archive the scientific code:

- manuscript source and compiled manuscript PDF;
- Overleaf-specific instructions;
- submission/editorial checklist;
- editor/referee response working note;
- internal physical-story notes;
- old generated manuscript figures and LaTeX tables;
- generated prose reports that duplicate the compact machine-readable products.

This also avoids exposing stale manuscript titles or pre-submission wording in the software archive.

## Scientific code

The 13 Python programs under `analysis/frozen_programs/` are preserved unchanged from the recovered source-and-code package. Their scientific constants, frozen event selections, thresholds, model grids, and statistical logic were not edited during cleanup.

The scripts parse successfully as Python. The archive validation test confirms that the compact event table contains 12 events and that the key JSON/CSV/NPZ products are readable.

## Dependency audit

The recovered package listed NumPy, Matplotlib, SciPy, and Astropy. One frozen script also imports `PIL`, so `Pillow>=10.0` has been added to the public `requirements.txt`. No scientific code change was required.

## Privacy / secret scan

No API keys, passwords, access tokens, local user paths, Gmail addresses, private correspondence, or ChatGPT/OpenAI references were found in the copied scientific scripts, specifications, or compact data products.

## Reproducibility boundary

The frozen scripts are not all directly runnable from the compact archive alone because some expect intermediate image cubes or archive-retrieved observations that are intentionally not redistributed. This is stated explicitly in `README.md` and `REPRODUCIBILITY.md` rather than implying a one-command raw-data reproduction.

## Public-release status

The scientific package is cleared for public release under the BSD 3-Clause License. After GitHub/Zenodo publication, add the final repository URL and Zenodo DOI to the citation metadata/README.
