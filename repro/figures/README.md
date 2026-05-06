# Per-figure reproducibility notes

This directory is the index for manuscript figure regeneration notes. The lightweight source repository records the expected organization and smoke-testable code paths; large figure inputs, controlled-access metadata, and generated full-resolution outputs should be archived in the final public release bundle when redistribution is permitted.

Recommended structure for each public-release figure directory:

- `README.md`: exact command(s), input artifact IDs, software environment, and expected output files
- `inputs.tsv`: input artifact IDs matching `repro/manifest.tsv`
- `expected_outputs.tsv`: output filenames, formats, and checksums when stable

Current release-candidate limitation: the Git repository contains smoke tests and release metadata, but not all large benchmark matrices or per-figure source data. Those assets must be synchronized with the public Zenodo/GitHub release before final journal upload if the journal requires full public reproducibility at submission.
