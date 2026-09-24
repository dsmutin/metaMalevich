# Dataset audit: low75

**Date:** 2026-09-24
**Overall status:** Ready with warnings

## Inventory

| File | Path | Notes |
|------|------|-------|
| Pinned accessions | examples/low75/accessions.tsv | 150 rows, 75 genera |
| Half table | examples/low75half/accessions.tsv | 76 rows, 38 genera |
| FASTA | data/raw/fasta/ | 150 files extracted from data/raw/ncbi_low75.zip (110 MB, datasets 18.36.0, valid data package) |
| Manifest | data/manifests/low75_download_manifest.tsv | 150 rows, validation pass |
| Checksums | data/checksums/low75_checksums.txt | zip plus each FASTA |

The held-out 40-genome manifest was not rewritten.

## Checks

Each genus has two different species taxids. Sim and db sequence digests differ. The assembly-report organism name starts with the genus in the table. Domains are 25 bacteria, 25 archaea, and 25 viruses.

## Warnings

Species-rank profiles cannot match the simulation, because the Kraken2 and Kaiju databases contain the other species of the genus. That is the design, not a missing file. The read budget is 100000 paired fragments across 75 genomes, the same budget as the twenty-genome example.

low75half keeps 38 genera because 75 is odd. It is every other genus in table order, not a rounded half of 37.5.
