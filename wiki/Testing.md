# Testing

## Test data

| Path | Role |
|------|------|
| `examples/toy/data/` | toy CLI output |
| `tests/test_colouring.py` | taxonomy, multi-label colours, toy L1 |
| `samovar10/` | closed community, 10 genomes, whole-contig Kraken2, 399 contigs |
| `samovar10_ont1b/` | same community layout with the ONT Kraken2 table |
| `intermediate/{dataset}/` | Kraken counts, kNN edges, CFA, CDBG (not committed) |
| `benchmark/{hypothesis}/{dataset}/` | profile, node assignments, metrics |

`strong100` has a Flye GFA, but segment ids (`edge_*`) are not the ground-truth contig ids. It is not a scored dataset. Illumina directories have read-level Kraken2 only.

## Integrative testing

1. Mandatory tests: `pytest -m mandatory`
2. CLI write path: `tests/test_integration.py`
3. Toy: `python examples/toy/run.py`
4. Benchmark: `metamalevich bench --data-root . --benchmark benchmark`
5. Full suite: GitHub Action `full-tests.yml` on release or `workflow_dispatch`

Required CI on every push or pull request: `required-tests.yml` (mandatory only), conda from `environment.yml`.
