# Data preparation report: high100

**Date:** 2026-09-24
**Mode:** Full execution
**Overall status:** Data-ready with warnings

Phase 0 found no high100 FASTA or manifest. The held-out 40 assemblies and the low75 assemblies were already on disk and were not downloaded again.

Phase 1 used NCBI Datasets CLI 18.36.0 on 200 accessions pinned by `examples/high100/select_accessions.py` from RefSeq `assembly_summary.txt` (bacteria, archaea, viral, fungi, protozoa; fetched 2026-09-24) and the local NCBI taxdump. The zip is `data/raw/ncbi_high100.zip`, 495682943 bytes, a valid zip with 200 genomic FASTA members and `assembly_data_report.jsonl`. Uncompressed sequence in the pin sums to 1615673298 bp. A first write of that zip on the project disk exited 1 and left no file. The successful download was written under `/tmp` and copied to `data/raw`.

Phase 2 audit: `docs/dataset-audit-high100.md`, status Ready with warnings. Biological checks passed: genus prefix, two species per family, `pair_mode` matches the genus taxids, distinct sequence digests. Each domain has 13 same-genus families and 12 cross-genus families.

Phase 3: no new P0. MEGAHIT, Kraken2, Kaiju, and Samovar remain outside `environment.yml`, as already recorded for the held-out example. The home directory quota blocked Snakemake's default cache; the simulation was rerun with the cache on local disk. That is an environment limit, not a change to the example.

`examples/half100half` reuses the same FASTA files. It does not download a second package.

Family-rank scores are in `examples/high100/RESULTS.md` and `examples/half100half/RESULTS.md`. No resolver was changed.
