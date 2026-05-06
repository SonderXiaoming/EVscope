# EMapper MAQC/qPCR benchmark notes

This directory documents the benchmark comparison used in the manuscript text and supplementary figure legends. Large input matrices and generated benchmark outputs are not committed to this lightweight source repository; they should be archived with the final public release bundle when redistribution is permitted.

Current release-candidate scope:

- Public RNA-seq accessions used for lightweight benchmark examples: `SRR31350808`-`SRR31350811`.
- The qPCR comparison is reported as a descriptive concordance/sanity check, not as a ground-truth accuracy benchmark.
- Exact public matrix filenames, checksums, and access locations must be recorded in `repro/manifest.tsv` when the final Zenodo/GitHub release is created.

Minimum expected public-release artifacts when finalized:

- input matrix description (`inputs.tsv` or equivalent)
- benchmark command script or notebook
- processed concordance summary table
- figure-generation commands for Supplementary Figs. S4-S6
- checksums for redistributable inputs and outputs
