# Audit of Kraken2 node colouring

## metaMalevich remote

`https://github.com/dsmutin/metaMalevich` was created on 2026-09-23 and is an empty Git repository (`size: 0`). The GitHub contents API returns 409 for the default branch, and the wiki clone has no pages yet. There is no Kraken2 colouring source in that remote.

## Local artefacts in this work tree

The colouring products live next to the assemblies:

| Path | What it is |
|---|---|
| `samovar10/kraken2/whole_contig.output` | Kraken2 calls plus the per-k-mer taxid counts for each contig |
| `samovar10/kraken2/whole_contig.report` | Taxonomy tree used as the hierarchy (names, ranks, parents by indentation) |
| `samovar10/taxonomy/strategy_a_whole_contig_probs.csv` | One row per contig, probability fixed at 1.0 on the Kraken call |
| `samovar10/taxonomy/strategy_b_reads_to_contig_probs.csv` | Multi-label weights whose ids are `edge_*`. The file has no source or target column |
| `samovar10/graph/graph_meta.txt` | `n_nodes=399`, `n_edges=6384`, `top_k=8`, `min_sim=0.15` |
| `samovar10_ont1b/` | Same layout as `samovar10`, including whole-contig Kraken2 |
| `strong100/assembly_graph.gfa` | Flye repeat graph. Segment names are `edge_*` |
| `strong100/ground_truth/ground_truth.csv` | Labels on `contig_*`, which are not the Flye segment ids |
| `strong100/taxonomy/strategy_a_whole_contig_probs.csv` | Hard Kraken labels on `contig_*` |
| `strong100/taxonomy/strategy_c_blobtools_probs.csv` | Hard labels on `contig_*` |

`strategy_a` does not keep a distribution. Several `samovar10` calls are strains under the genome's species (for example taxid 1206105 under 1351, Enterococcus faecalis). The Kraken report also shows Enterococcus contigs split across 1351, 2774749, and 2774762.

## What this tool uses

Node colours are the k-mer counts in `whole_contig.output`, summed in `cpp/kraken_count` and rolled to species with the report hierarchy. That is the `probability_sum` aggregation. The hard Kraken call, rolled to the same species rank, is the `initial_colouring` hypothesis.

The kNN edge list that `graph_meta.txt` describes was not stored, and `mean_weight` is not a Jaccard (it is 2.69 on `samovar10` and 24.6 on `samovar10_ont1b`). A 21-mer minimizer Jaccard at `min_sim=0.15` kept only 82 edges, because disjoint 100 kb contigs from one genome share few exact 21-mers. The rebuilt graph is therefore a canonical 4-mer cosine graph with the recorded `top_k=8` and `min_sim=0.15`. New edge ids are not joined to `strategy_b`, because that table has no endpoints. Edge colours are the renormalized minimum of the two endpoint species distributions, so an edge can carry fewer taxa than either node.

`strong100` is not in the scored benchmark. Its assembly-graph nodes and its ground-truth contig ids are different sequences (lengths do not match), and this tree has no path map between them.

`samovar10_ont1b` assembly FASTA, whole-contig Kraken2 output, and ground truth are byte-identical to `samovar10`. The second benchmark folder repeats the same run.

## Measured profiles

Species-rank scores on `samovar10` (399 contigs, 10 genomes). Lower L1 is better.

| Hypothesis | L1 | Macro F1 | Accuracy |
|---|---:|---:|---:|
| initial_colouring | 0.1586 | 0.913 | 0.915 |
| probability_sum | 0.1441 | 0.918 | 0.925 |
| gated_neighbour | 0.1378 | 0.878 | 0.917 |
| bayesian_edge | 0.1594 | 0.886 | 0.912 |
| lca | 1.608 | 0.594 | 0.193 |

`probability_sum` improves both abundance and species F1 relative to the hard Kraken call. `gated_neighbour` has the lowest L1 and raises path consistency from 0.884 to 0.975, and it lowers macro F1. Relative abundance of Enterococcus sp. DIV2432 (truth 0.081) is 0.031 for the hard call, 0.038 for `probability_sum`, and 0.020 for `gated_neighbour`. Composition neighbours of those contigs are other Enterococcus contigs already labelled faecalis, so neighbour mixing amplifies that swap. The soft species-level k-mer profile is the reprofile that moves DIV2432 toward the truth.

Illumina bundles under `samovar10_illumina_10M` and `samovar10_illumina_100m` have read-level Kraken2 output and no whole-contig output, so they are not coloured here.

## MetaMetro

CFA colouring and CFA-to-CDBG compaction are MetaMetro's `colour_cfa` and `cfa_to_cdbg`. Colour sets in that schema are integer masks. Posterior weights stay in the hypothesis TSVs. The observed multi-label mask is written to `intermediate/{dataset}/tocumg/`.
