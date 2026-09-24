# low75 comparison

All 17 assembly hypotheses were scored on the k141 FASTG (2685 nodes, 24 edges), plus Kraken2 and Kaiju on the reads. The database species is never the simulated species, so species-rank profiles are disjoint.

| method | weight | l1 | presence_f1 | node_f1 | node_accuracy |
| --- | --- | ---: | ---: | ---: | ---: |
| kraken2 | read_count | 2 | 0 |  |  |
| kaiju | read_count | 2 | 0 |  |  |
| every assembly hypothesis | node_length | 2 | 0 | 0 | 0 |

`every assembly hypothesis` is `initial_colouring`, `probability_sum`, `lca`, `gated_neighbour`, `bayesian_edge`, `edge_union`, `leakage`, `leakage_flipped`, `label_drop`, `label_drop_flipped`, `logistic_drop`, `logistic_drop_flipped`, `confident_lock`, `unanimous_rescue`, `cut_then_leak`, `genus_plurality`, and `margin_mixture`. Bray-Curtis is 1. Pearson R² is 0.002 to 0.004. Coverage-weighted rows are in `work/metrics.tsv` and are also L1 2. No resolver was adopted.

Rolling the same argmax calls up to genus, from the assignment tables and the Kraken2 taxdump, separates `lca` from the rest. This is not the library species metric.

| method | genus node accuracy | genus base accuracy | genus L1 |
| --- | ---: | ---: | ---: |
| initial_colouring | 0.2171 | 0.2336 | 1.5261 |
| probability_sum and the other non-LCA hypotheses | 0.2179 | 0.2340 | 1.5240 |
| lca | 0.1058 | 0.1159 | 1.7835 |

Kraken2 left 76157 of 100000 reads unclassified. Of the remaining reads, 23401 were placed in one of the 75 pinned genera. None was placed on a simulated species taxid. The largest species-rank Kraken fraction in `work/comparison.tsv` is unclassified 0.76599. The next calls are database-side species such as Acinetobacter calcoaceticus (0.04933), not the simulated species from that genus.

MEGAHIT wrote 2684 contigs, 1350101 bp, N50 465 bp. The k141 FASTG has 2685 nodes and 24 edges. The coloured CFA and CDBG have those counts, 35 colours, and `topology_unchanged` true.

The read budget is the same 100000 paired fragments used for twenty genomes in `examples/heldout_genera`, spread across 75 genomes.
