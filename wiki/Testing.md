# Testing

## Test data

| Path | Role |
|------|------|
| `examples/toy/data/` | Toy CLI JSON. Three nodes; neighbours correct a bad hard call |
| `tests/test_colouring.py` | Taxonomy, explicit colour operations, toy L1 |
| `tests/test_reprofile.py` | Unclassified bases counted once; partition rules |
| `tests/test_kmer.py` | 4-mer graph builder limits |
| `tests/test_bridge.py` | ToCUMG export keeps topology |
| `tests/test_evaluate.py` | L1 / F1 and the "both must improve" rule |
| `samovar10/` | Closed community, 10 genomes, whole-contig Kraken2, 399 contigs. Not in git |
| `intermediate/{dataset}/` | Kraken counts, kNN edges, CFA, CDBG. Not in git |
| `benchmark/{hypothesis}/samovar10/` | Profile, node assignments, metrics |

`strong100` has a Flye GFA. Segment ids (`edge_*`) are not the ground-truth contig ids. It is not a scored dataset. Illumina directories have read-level Kraken2 only. `samovar10_ont1b` matches `samovar10` byte for byte on the FASTA, the whole-contig Kraken2 output, and the ground truth.

## Integrative testing

1. Mandatory tests: `pytest -m mandatory`
2. CLI write path: `tests/test_integration.py`
3. Toy: `python examples/toy/run.py`
4. Benchmark, when `samovar10/` is on disk: `python -m metamalevich bench --data-root . --benchmark benchmark`
5. Full suite: GitHub Action `full-tests.yml` on release or `workflow_dispatch`

Required CI on every push or pull request to `master` or `main`: `required-tests.yml` (mandatory marker only), conda from `environment.yml`. The default branch of the code repository is `master`.

## What the mandatory suite does not do

It does not download Kraken indexes, does not assemble reads, and does not require the `samovar10` bundle. The community scores on [Benchmark](Benchmark) are a separate run.
