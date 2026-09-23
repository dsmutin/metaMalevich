# Pipeline

The specification separates four contracts. A new resolver does not change colouring, and the profiler does not change the resolver.

```text
Kraken2 whole-contig output
        │
        ▼
cpp/kraken_count  →  species rollup
        │
assembly FASTA    →  cpp/kmer_knn (4-mer cosine, top-k 8, min_sim 0.15)
        │
        ▼
node colours and edge colours (separate tables)
        │
        ▼
MetaMetro colour_cfa → cfa_to_cdbg
        │
        ▼
resolver (full distribution per node)
        │
        ▼
length-weighted profile
```

Complex TCA→ToCUMG conversions from the specification are not implemented. The built-in MetaMetro path above is the one that runs.

## Implementation order

From section 40 of the specification. A row is filled only when that stage is in this tree.

| Step | Spec | In this tree |
|------|------|----------------|
| 1 | Taxonomy model | `taxonomy.parse_kraken_report`: id, parent, rank, name, source, version |
| 2 | TCA contract | `evidence.EvidenceGraph.apply` with explicit `replace` / `merge` / `intersect` / `subtract` |
| 3 | Node and edge colours | Separate tables. Multi-label until a resolver runs |
| 4 | ToCUMG | MetaMetro CFA, then CDBG. Posterior weights stay outside the uint8 mask |
| 5 | TCA → ToCUMG | `bridge.export_tocumg`. Topology is checked unchanged |
| 6 | Deterministic resolvers | Hard call, `probability_sum`, LCA, gated neighbour |
| 7 | Bayesian resolver | `bayesian_edge` |
| 8 | Read-to-graph reassignment | Not implemented. `assigned_reads` is 0 |
| 9 | Graph-aware profiler | `reprofile.profile_nodes`, weighted by node length |
| 10 | GCN | Not started. The specification says not to start here |
| 11 | Path-aware resolver | Not implemented. The score that is recorded is neighbour agreement |
| 12 | Benchmark | `samovar10` hypotheses under `benchmark/{hypothesis}/` |

## What a bench writes

Under `intermediate/{dataset}/` (gitignored, reused when the input signature matches):

| File | Contents |
|------|----------|
| `kraken_counts.tsv` | Per-contig `taxid` counts from the k-mer column |
| `kraken_calls.tsv` | Classifier call (not the k-mer argmax) |
| `knn_edges.tsv` | `edge_id`, `source`, `target`, `orientation`, `weight` |
| `node_evidence.tsv` | Species weight and the raw species count |
| `edge_evidence.tsv` | Renormalized minimum of the two endpoint distributions |
| `colour_layer.yaml` | Layer id, source, method, taxonomy version |
| `tocumg/` | CFA and CDBG |

Under `benchmark/{hypothesis}/{dataset}/`:

| File | Contents |
|------|----------|
| `profile.tsv` | Columns in section 29 of the specification |
| `node_assignments.tsv` | Argmax taxon, probability, truth taxon, length |
| `metrics.yaml` | Scores plus the reproducibility manifest (section 36) |
| `summary.tsv` | One-line scores |

`benchmark/summary/` holds `metrics.tsv`, `versus_initial.yaml`, and Altair `l1`, `f1`, and `abundance` HTML/JSON.

## Not in this version

Specified, and not built:

| Spec | Why it is absent |
|------|------------------|
| Read reassignment and Bracken-style read profiles | No read-to-contig map in the scored bundle |
| Path-aware Bayesian resolution, local subgraphs, hybrid ML | Stopped at the deterministic milestone |
| GCN, GraphSAGE, GAT | Section 40: do not start with the GCN |
| Tier 1 besides contig Kraken2 (Bracken, MetaPhlAn) | Those outputs are not the scored comparison |
| Tier 2 (GraphBin2, MetaCoAG, Bifrost, GGCAT) | Not run |
| Posterior intervals, bootstrap, graph perturbation | Point estimates plus the node distribution are what is stored |
| Ranks other than species in the scored F1 | Evaluation rolls labels to species |
| `strong100` | Flye `edge_*` ids are not the `contig_*` ground truth |
| Illumina Kraken bundles | Read-level output, no whole-contig file |

`samovar10_ont1b` is not a second dataset: its FASTA, whole-contig Kraken2 output, and ground truth match `samovar10`.
