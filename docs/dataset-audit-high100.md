# Dataset audit: high100

**Date:** 2026-09-24
**Overall status:** Ready with warnings

## Inventory

| File | Path | Notes |
|------|------|-------|
| Pinned accessions | examples/high100/accessions.tsv | 200 rows, 100 families |
| Half table | examples/half100half/accessions.tsv | 100 rows, 50 families |
| FASTA | data/raw/fasta/ | 200 files extracted from data/raw/ncbi_high100.zip |
| Assembly report | data/raw/high100_assembly_data_report.jsonl | from that zip |
| Manifest | data/manifests/high100_download_manifest.tsv | 200 rows, validation pass |
| Checksums | data/checksums/high100_checksums.txt | zip plus each FASTA |

The held-out 40-genome manifest and the low75 manifest were not rewritten. half100half does not download a second zip.

## Checks

Each domain (bacteria, archaea, viral, eukaryota) has 25 families: 13 `same_genus` and 12 `different_genera`. Twenty-five is odd, so the split is 13 and 12. Each family has two species taxids. Sim and db sequence digests differ. The assembly-report organism name starts with the genus in the table. `pair_mode` matches whether the two genus taxids are equal.

Eukaryotes are fungi and protozoa only. The selector keeps the smallest complete or chromosome family pairs. Those 25 pairs sum to 1208254752 bp. The smallest pair is 4749303 bp and the largest is 82423493 bp.

## Warnings

The read budget is 100000 paired fragments across 100 genomes, the same budget as the twenty-genome example. Most families will be missing from the assembly. Scores use family rank, so a call of the database species counts when it is in the simulated genome's family. The library default used by samovar10 stays species.

half100half keeps 50 families because 100 is even. It is every other family in table order.
