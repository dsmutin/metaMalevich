# Data preparation report: low75

**Date:** 2026-09-24
**Mode:** Full execution
**Overall status:** Data-ready with warnings

Phase 0 found no low75 FASTA or manifest. The held-out 40 assemblies were already valid and were not downloaded again.

Phase 1 used NCBI Datasets CLI 18.36.0 on 150 accessions pinned from RefSeq `assembly_summary.txt` (Last-Modified 2026-09-23) and the local NCBI taxdump. The zip is `data/raw/ncbi_low75.zip`, 110 MB, reported as a valid data package. Uncompressed sequence in the pin sums to 397184062 bp. That is under 1 GB, not a 50 GB download.

Phase 2 audit: `docs/dataset-audit-low75.md`, status Ready with warnings. Biological checks passed: genus prefix, two species per genus, distinct sequence digests.

Phase 3: no new P0. The tool profile is unchanged. MEGAHIT, Kraken2, Kaiju, and Samovar remain outside `environment.yml`, as already recorded for the held-out example.

The species-rank comparison is disjoint by design. See `examples/low75/RESULTS.md` and `examples/low75half/RESULTS.md`.
