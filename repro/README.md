# EVscope reproducibility index

This directory records release-level reproducibility metadata for EVscope. Large reference files and sequencing data are not stored in Git; use the DOI/accessions in `manifest.tsv`.

Minimum release checks:

```bash
bash tests/smoke/run_smoke.sh
sha256sum -c repro/checksums.sha256
```

Benchmark scripts and per-figure commands should be added under `repro/benchmarks/` and `repro/figures/` as they are finalized.
