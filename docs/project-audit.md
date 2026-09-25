# Project Audit Report

**Date:** 2026-09-24
**Repository:** metamalevich (`vaegbin_improved`)
**Branch / commit:** master. Held-out example `f1a9514`. Leakage module `82969a3` (no `VERSION` change). FASTG scoring `f38fcb6` sets `VERSION` to 0.2.1.
**Auditor mode:** Read-only inspection; this report is the written deliverable
**Harness:** tool
**Properties in scope:** scientific-integrity, validation-first, reproducibility, missing-data-policy, english-docs, no-push, follow-contributing

`method_decision`, SLURM, and artifact-registry are disabled. They are not treated as missing. `@verify-methods` was not run.

## Executive summary

- The package is a conda tool with English docs, pinned library dependencies in `environment.yml`, and a mandatory pytest suite.
- The held-out genera example pins 40 NCBI assemblies and keeps genomes, reads, and indexes out of git (`data/`, `examples/heldout_genera/work/`).
- Comparison scores exist in the example `RESULTS.md`. The assembly graph used for those scores is the k141 FASTG (6439 nodes, 24 edges), exported to CFA and CDBG.
- No resolver change was made from this run. Neighbour methods do not improve L1 and node F1 together on that graph.
- MEGAHIT, minimap2, Kraken2, Kaiju, and Samovar are required to re-run the example and are not pins in `environment.yml`.
- No P0 blocker for using the prepared assemblies or the recorded scores.

## Overall scores

| Dimension | Rating | Summary |
|-----------|--------|---------|
| Project completeness | Good | Library, tests, samovar10 bench, and the held-out example are in the tree. Raw genomes are local and gitignored. |
| Documentation | Good | README, example README, and `examples/heldout_genera/RESULTS.md` state how to run and what was measured. |
| Reproducibility | Fair | Library pins and seed 42 are recorded. Example binaries are versioned in the example README, not in `environment.yml`. |
| Methodological rigor | Good | Scores cite `work/metrics.tsv`. A hypothesis is not called better unless the recorded rule (both L1 and macro F1) holds. |
| Publication readiness | N/A | Tool profile. |

## 1. Project completeness

| Component | Status | Evidence |
|-----------|--------|----------|
| Raw assemblies | Present locally, not in git | `data/raw/fasta/` has 40 `.fna` files; `.gitignore` ignores `/data/` |
| Manifest and checksums | Present | `data/manifests/download_manifest.tsv` (40 rows), `data/checksums/checksums.txt` (41 lines) |
| Dataset audit | Present | `docs/dataset-audit.md`, status Ready with warnings |
| Example runner | Present | `examples/heldout_genera/run_example.py` |
| Measured comparison | Present | `examples/heldout_genera/RESULTS.md`; tables in gitignored `work/` |
| Library environment | Present | `environment.yml` |
| Manuscript drafts | Out of scope | paper profile is disabled |

## 2. Missing documentation

| Document | Gap | Impact | Priority |
|----------|-----|--------|----------|
| method-decision.md | Absent | Disabled by the harness. Not a gap. | — |
| docs/artifact-registry.md | Absent | Disabled. Not a gap. | — |

## 3. Reproducibility issues

| Issue | Evidence | Impact | Priority |
|-------|----------|--------|----------|
| Example tools are not in `environment.yml` | `environment.yml` lists Python, pytest, numpy, pyyaml, pandas, altair, and `gxx_linux-64=16.2.0`. MEGAHIT 1.2.9, minimap2 2.31, Kraken2 2.0.7-beta, Kaiju 1.10.1, and Samovar are named in `examples/heldout_genera/README.md` only. | A fresh conda env from `environment.yml` does not assemble or classify this example. | P1 |
| Samovar is an external binary | README records `git@github.com:ctlab/samovar.git`. The run log calls a Samovar binary outside this repository. | Another machine needs its own Samovar checkout. | P1 |
| Seed and read count are fixed in the example | `run_example.py` uses `TOTAL_READS = 100000` and `SEED = 42`. | The simulation design is recoverable from the script. | — |

## 4. Inconsistent methods

**Source of truth:** not applicable. `method_decision` is disabled, so methods were not reconstructed with `@verify-methods`.

| Conflict | From verify-methods | Audit note | Priority |
|----------|---------------------|------------|----------|
| — | not run | `RESULTS.md` does not claim a graph hypothesis beat hard colouring on both L1 and node F1. | — |

## 5. Outdated software

| Tool | Pinned version | Concern | Evidence basis | Priority |
|------|----------------|---------|----------------|----------|
| — | — | No tool was judged outdated. Versions in the example README were taken from the binaries used for this run. | example README | — |

## 6. Duplicated functionality

| Function | Locations | Recommendation |
|----------|-----------|----------------|
| Species rollup | `examples/heldout_genera/community.py` `species_taxon` and the library taxonomy helpers used by the bench | Leave them separate. The example rollup is the scoring rule for this simulation and is tested in `tests/test_heldout_community.py`. | P3 |

## 7. Incomplete analyses

| Analysis | State | Evidence | Priority |
|----------|-------|----------|----------|
| Held-out compare | Done for this sample | `RESULTS.md` and `figures/heldout_abundance.html` | — |
| Resolver retune | Not done | `RESULTS.md`: 24 edges, and coverage weighting does not improve L1 and node F1 together | — |
| Decaying leakage | Committed, not scored on this sample | `82969a3` adds `src/metamalevich/leakage.py`. It is not a row in `work/metrics.tsv`. `src/metamalevich/flip.py` is untracked and was not used for the comparison. | P2 |

## 8. Code quality notes

| Area | Finding | Location | Priority |
|------|---------|----------|----------|
| ToCUMG failure is logged and scoring continues | `except Exception` writes the error to `work/commands.log` and returns | `examples/heldout_genera/run_example.py` `_export_assembly_tocumg` | P2 |
| Empty FASTG is regenerated | Size 0 is treated as missing, and only `k<int>.contigs.fa` is passed to `contig2fastg` | `stage_assemble`, `community.highest_megahit_contigs` | — |

## 9. Statistical methodology

Not reconstructed. `method_decision` is off.

| Check | Pass/Fail | Cite |
|-------|-----------|------|
| Effect sizes in the example | Pass | `RESULTS.md` reports L1, R², presence F1, and node F1 from `work/metrics.tsv` |
| Multiple-testing correction | Not applicable | No hypothesis test was run |

## 10. Figures and outputs

| Check | Status | Notes |
|-------|--------|-------|
| Held-out abundance chart | Present | `examples/heldout_genera/figures/heldout_abundance.html` and `.json` (Altair) |
| samovar10 error charts | Present | `figures/samovar10_species_error.html`, `figures/samovar10_enterococcus.html` |
| Publication figure review | N/A | publication-figures was not requested |

## 11. Tracking artifacts

| File | Status | Notes |
|------|--------|-------|
| todo.md | Feature checklist | Held-out genera item is checked only after the comparison and the mandatory tests. |
| VERSION | 0.2.1 at `f38fcb6` | Patch for the FASTG selection fix. `82969a3` did not bump `VERSION`. |

## 12. Publication readiness

**Rating:** N/A

**Rationale:** The harness profile is `tool`.

**Blocking items for this audit:** none at P0.

## Findings

| Priority | Issue | Evidence | Recommended action |
|----------|-------|----------|--------------------|
| P1 | Example binaries are outside the conda pin file | `environment.yml` versus `examples/heldout_genera/README.md` | Pin MEGAHIT and minimap2 when the example is meant to install from the env file, or keep the README install line as the source of those versions. |
| P1 | Samovar is not an env dependency | example README remote `git@github.com:ctlab/samovar.git` | Document the checkout as a prerequisite, which the README already does. |
| P2 | Leakage is committed and not in the held-out metrics | `82969a3`, `src/metamalevich/leakage.py`; `work/metrics.tsv` has no leakage row. `src/metamalevich/flip.py` is untracked. | Do not describe either as a result of this comparison. `82969a3` also left `VERSION` at 0.2.0. |
| P2 | ToCUMG errors do not fail the compare stage | `_export_assembly_tocumg` | This run did export (`export.json`, `topology_unchanged` true). A later failure would only appear in `commands.log`. |
