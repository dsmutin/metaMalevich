# Contracts

Change a row only together with tests. The specification's invariants are I1–I8.

## Implementation table

| Contract | Implementation |
|----------|----------------|
| Single version source (`VERSION`) | `metamalevich.__version__` reads `VERSION` |
| CLI `--version` and JSON run | `metamalevich.cli.main` |
| Pipeline result keys `status`, `ok`, `input_path` | `metamalevich.baseline.run_pipeline` |
| Conda-only install (`environment.yml`) | env file + `PYTHONPATH=src` |
| Mandatory pytest on every commit | `pytest -m mandatory` |
| Optional pytest on release / manual | `pytest -m optional` |
| Toy example runs the tool | `examples/toy/run.py` |
| English docs on public APIs | module docstrings |
| No `git push` unless asked | rule `no-push` (the wiki remote is separate) |
| I1 topology unchanged by TCA | `bridge.export_tocumg` checks the graph after `colour_cfa` |
| I2 multi-label evidence kept until resolution | `evidence.EvidenceGraph.apply`; empty `layer_id` cannot `replace` |
| I3 node and edge colours are separate | `node_colours` / `edge_colours`; edge colour is min of endpoints |
| I4 resolvers return distributions | `resolve.resolve` |
| I5 no forced linearization | Contigs stay nodes. No path is required to emit a profile |
| I6 provenance | `colour_layer.yaml` and the `metrics.yaml` manifest |
| I7 taxonomy version | Kraken report path stored as the taxonomy version |
| I8 unique, shared, ambiguous, unclassified | `reprofile.profile_nodes` |
| Named aggregations | `aggregate.aggregate`: `count`, `weighted_count`, `majority`, `lca`, `probability_sum` |
| Better than the initial colouring | L1 down and macro F1 up (`bench.compare_to_initial`) |
| Toy graph beats the hard colouring on L1 | `toy.run_toy` |

## Still open relative to the specification

| Contract | Status |
|----------|--------|
| Read reassignment (`assigned_reads`) | Column exists and is 0 |
| Path consistency | Replaced by neighbour agreement. Not a path score |
| Uncertainty intervals | Not computed |
| Evaluation at phylum through genus | Scored rank is species |
| Baselines B1–B8 beyond contig Kraken2 | Not run |
| GCN and path-aware resolvers | Not implemented |

Posterior weights are not stored in the MetaMetro uint8 colour mask. The mask is the set of taxa with positive weight. Weights stay in `intermediate/{dataset}/node_evidence.tsv` and in each hypothesis profile.
