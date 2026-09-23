# Data Preparation Report

**Date:** 2026-09-24
**Mode:** Full execution
**Overall status:** Data-ready with warnings

## Executive summary

The requested dataset is 40 NCBI assemblies: four genera, five strains each, two assemblies per strain. One assembly of each pair was simulated into a 100000-fragment paired sample; the other was used for the Kraken2 and Kaiju databases. Download, checksums, and biological checks passed. The dataset audit is Ready with warnings. The project audit has no P0 item. Genomes and derived reads stay in gitignored `data/` and `examples/heldout_genera/work/`.

## Requested datasets

| Identifier | Source requested | Required |
|------------|------------------|----------|
| 40 RefSeq assemblies in `examples/heldout_genera/accessions.tsv` | NCBI Datasets genome assemblies | yes |

The accessions were chosen from NCBI Datasets summaries and written into that table before download. They are not invented stand-ins.

---

## Phase 0 — Project state detection

### Existing artifacts found

| Artifact | Path | Status |
|----------|------|--------|
| Download manifest | data/manifests/download_manifest.tsv | present — 40 rows |
| Acquisition report | data/manifests/acquisition_report.md | present |
| Dataset audit | docs/dataset-audit.md | present — Ready with warnings |
| Project audit | docs/project-audit.md | written in this phase |
| Download log | data/logs/download.log | present |
| Checksums | data/checksums/checksums.txt | present — 41 lines (zip plus 40 FASTA files) |

### Per-dataset state

| ID | Present | Technical validation | Biological validation | Decision |
|----|---------|---------------------|----------------------|----------|
| 40 accessions in `accessions.tsv` | yes, `data/raw/fasta/*.fna` (40 files) | pass, recorded in the manifest | pass, pair sequence digests differ | skip re-download |

### Skipped actions (Phase 0)

| Action | Reason |
|--------|--------|
| Second `datasets download` | Zip and 40 FASTA files were already present and the manifest records validation pass |

### Detected issues in existing data

| Issue | Severity | Details |
|-------|----------|---------|
| Two Streptococcus references have no strain name | warn | Stated in `docs/dataset-audit.md` |
| Two Shigella flexneri pairs share a species rank | warn | Species-rank scores add those read counts. Documented in the dataset audit and in `RESULTS.md`. |

---

## Phase 1 — Data acquisition (@get-data)

**Invoked:** Yes, earlier in this preparation. Not repeated.

### Summary

| Metric | Value |
|--------|-------|
| Datasets requested for download | 40 |
| Downloaded | 40 |
| Skipped (already valid) | 0 on the first run; the resume skips a valid zip |
| Failed | 0 |

### Repositories and methods

| Accession | Repository | Method | Status |
|-----------|------------|--------|--------|
| 40 accessions listed in `examples/heldout_genera/accessions.tsv` | NCBI Assembly via Datasets CLI 18.36.0 | `datasets download genome accession` | success |

`download_url` in the manifest is empty. The datasets summary used here has no per-file FTP URL. The zip SHA-256 is in `data/checksums/checksums.txt`.

### Failures (download vs validation)

| Accession | Failure type | Error | Blocks pipeline |
|-----------|--------------|-------|-----------------|
| — | — | — | — |

**Subordinate report:** [data/manifests/acquisition_report.md](../data/manifests/acquisition_report.md)

---

## Phase 2 — Dataset audit (@dataset-auditor)

**Invoked:** Yes

### Summary

**Audit status:** Ready with warnings

| Check area | Result |
|------------|--------|
| Metadata completeness | warn — `eco_coli` sim taxid filled from the download report; two Streptococcus references have no strain name |
| Paired-end consistency | pass for the simulated sample: 100000 R1 and 100000 R2 records |
| Sequencing depth | by design — 100000 paired fragments across 20 genomes, not fold coverage |
| Batch confounding | not applicable — one simulated sample |

**Subordinate report:** [docs/dataset-audit.md](dataset-audit.md)

### Critical issues

| Issue | Impact | Recommendation |
|-------|--------|----------------|
| — | — | — |

---

## Phase 3 — Scientific project audit (@project-auditor)

**Invoked:** Yes

### Summary

**Publication / analysis readiness:** N/A (tool profile). No P0.

| Dimension | Rating |
|-----------|--------|
| Reproducibility | Fair |
| Data organization | Good |
| Documentation | Good |

**Subordinate report:** [docs/project-audit.md](project-audit.md)

### Blocking issues (P0)

| Issue | Evidence |
|-------|----------|
| — | — |

P1, not blocking data use: MEGAHIT, minimap2, Kraken2, Kaiju, and Samovar are not pins in `environment.yml`. Versions and the Samovar remote are in `examples/heldout_genera/README.md`.

---

## Consolidated findings

### Detected issues (all phases)

| # | Phase | Severity | Issue | Status |
|---|-------|----------|-------|--------|
| 1 | 2 | warn | Two reference assemblies have no strain name | open, recorded |
| 2 | 2 | warn | Two flexneri pairs merge at species 623 | open, recorded; scores add the read counts |
| 3 | 3 | P1 | Example tools are outside `environment.yml` | open |

### Unresolved problems

- Re-running the example on a new machine still needs Samovar, Kraken2, Kaiju, MEGAHIT, and minimap2 in addition to `environment.yml`.
- The simulated sample is shallow relative to the 20 genomes. That is the design in `run_example.py` (`TOTAL_READS = 100000`), not a missing file.

### Recommended next steps

1. Keep using the pinned accessions and the gitignored `data/` tree.
2. Treat `examples/heldout_genera/RESULTS.md` as the comparison record. Do not retune resolvers from the 24-edge graph unless a later sample shows a different error.

---

## Success criteria checklist

| Criterion | Met | Evidence |
|-----------|-----|----------|
| All datasets acquired or available | yes | 40 FASTA files and 40 manifest rows |
| Technical validation passed | yes | manifest validation status and `data/checksums/checksums.txt` |
| Biological validation passed | yes | acquisition report and dataset audit: pair digests differ |
| Dataset audit successful | yes, with warnings | `docs/dataset-audit.md` status Ready with warnings |
| No critical project audit blockers | yes | `docs/project-audit.md` has no P0 |
| Reports and manifests generated | yes | paths in the table above |
| Ready for downstream analyses | yes, with the warnings above | comparison already run; see `examples/heldout_genera/RESULTS.md` |

**Data-ready:** Yes, with warnings. The warnings do not remove any of the 40 assemblies.

---

## Referenced reports (preserved, not replaced)

| Report | Path |
|--------|------|
| Acquisition | data/manifests/acquisition_report.md |
| Download manifest | data/manifests/download_manifest.tsv |
| Dataset audit | docs/dataset-audit.md |
| Project audit | docs/project-audit.md |
| Download log | data/logs/download.log |
| Example scores | examples/heldout_genera/RESULTS.md |
