# Dataset Audit Report

**Date:** 2026-09-24
**Auditor:** data-preparation orchestrator
**Overall status:** Ready with warnings

## 1. Data inventory

| File | Path | Rows × Cols | Notes |
|------|------|-------------|-------|
| Pinned accessions | examples/heldout_genera/accessions.tsv | 40 × 11 | 20 sim / 20 db |
| FASTA | data/raw/fasta/ | 40 files | one NCBI assembly each |
| Manifest | data/manifests/download_manifest.tsv | 40 × 11 | SHA-256 recorded |

## 2. Metadata completeness

Every accession has a genus, a role, and an organism name from the NCBI datasets summary.
`eco_coli` sim `GCF_000952955.1` had no tax_id in the summary used to pin the table; the download report supplies it.
Two RefSeq references (`Streptococcus agalactiae`, `Streptococcus thermophilus`) had no strain name in the summary.

## 3. Sample identifiers

Unique assembly accessions: yes. Sim keys after dropping the version do not collide.

## 4. Class balance

Four genera, five pairs each. Four of the five Escherichia pairs are not E. coli.
Shigella has two flexneri pairs (species 623 lineage and flexneri 2a taxid 42897). Species-rank scores add those read counts.

## 5. Duplicates and missing values

No duplicate accessions. Pair sequence digests differ, so a GenBank twin of a RefSeq assembly was not kept.

## 6. Outliers

Assembly length is the NCBI `total_sequence_length` recorded in the accession table (about 1.7–7.1 Mb).
Read depth is not a property of these assemblies; the simulated sample is 100000 paired fragments.

## 7. Sequencing layout

Not applicable to the downloaded assemblies. InSilicoSeq HiSeq paired reads are created later by Samovar.

## 8. Repository consistency

All 40 rows are NCBI Assembly accessions retrieved with the datasets CLI. Organism genus matches the request.
`Pseudomonas stutzeri` was not used: the NCBI summary for that taxid returned `Stutzerimonas stutzeri`.

## 9. Decision

Ready with warnings for the held-out simulation. The warnings are missing strain names on two references and the shared Shigella flexneri species rank.
