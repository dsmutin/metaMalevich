# Contracts

Change a row only together with tests.

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
| No `git push` unless the human asks | rule `no-push` |
| Taxonomy ids, parents, ranks, source, version | `taxonomy.parse_kraken_report` |
| TCA keeps multi-label evidence until an explicit colour operation | `evidence.EvidenceGraph.apply` |
| Aggregations are named (`count`, `weighted_count`, `majority`, `lca`, `probability_sum`) | `aggregate.aggregate` |
| Node and edge colours are separate tables | `EvidenceGraph.node_colours` / `edge_colours` |
| Resolvers return distributions | `resolve.resolve` |
| Profiler reports unique, shared, ambiguous, and unclassified bases | `reprofile.profile_nodes` |
| CFA to ToCUMG uses MetaMetro and does not change topology | `bridge.export_tocumg` |
| Toy graph reprofile beats the hard colouring on L1 | `toy.run_toy` |

## Pipeline

```
TCA (Kraken2 k-mer counts on contigs, plus minimizer-Jaccard edges)
  → ToCUMG (MetaMetro CFA, then cfa_to_cdbg)
  → Resolve
  → Reprofile
```

Posterior weights are not stored in the MetaMetro `uint8` colour mask. The mask is the set of taxa with positive weight. Weights stay in `intermediate/{dataset}/node_evidence.tsv` and in each hypothesis profile.
