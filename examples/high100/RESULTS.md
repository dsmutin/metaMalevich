# high100 comparison

Scores roll to family. A database species in the same family as the simulated genome counts as that family. One hundred families were simulated (25 bacteria, 25 archaea, 25 viruses, 25 eukaryotes). Thirteen families in each domain use two species of one genus. Twelve use two genera.

All 17 assembly hypotheses were scored on the k141 FASTG (989 nodes, 4 edges), plus Kraken2 and Kaiju on the reads.

| method | weight | l1 | presence_f1 | presence_precision | presence_recall | node_f1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| kraken2 | read_count | 1.6323 | 0.901099 | 1 | 0.82 |  |
| kaiju | read_count | 1.2605 | 0.95288 | 1 | 0.91 |  |
| initial_colouring | node_length | 1.71456 | 0.387097 | 1 | 0.24 | 0.220142 |
| probability_sum and the other non-LCA hypotheses | node_length | 1.71376 | 0.387097 | 1 | 0.24 | 0.220142 |
| lca | node_length | 1.87903 | 0.165138 | 1 | 0.09 | 0.106753 |

`probability_sum and the other non-LCA hypotheses` is `probability_sum`, `gated_neighbour`, `bayesian_edge`, `edge_union`, `leakage`, `leakage_flipped`, `label_drop`, `label_drop_flipped`, `logistic_drop`, `logistic_drop_flipped`, `confident_lock`, `unanimous_rescue`, `cut_then_leak`, `genus_plurality`, and `margin_mixture`. Their presence F1 and node F1 match. `label_drop` and `logistic_drop` differ from `probability_sum` only past the printed L1 (Pearson 0.0128154 versus 0.0127056). Coverage-weighted rows are in `work/metrics.tsv`. They do not change which hypothesis is ahead, and their L1 values are higher (1.91501 for the non-LCA group, 1.9764 for `lca`).

`probability_sum` improves L1 by 0.0008 over hard colouring and does not improve presence F1 or node F1. No resolver was adopted.

Kraken2 classified 19017 of 100000 read pairs (19.02%). After the family roll, the unclassified fraction in `work/comparison.tsv` is 0.81317. Every non-zero Kraken2 family is one of the 100 simulated families (82 families, precision 1). Kaiju's unclassified fraction is 0.61366, with 91 families and precision 1.

MEGAHIT wrote 988 contigs, 520172 bp, N50 428 bp. The k141 FASTG has 989 nodes and 4 edges. The coloured CFA and CDBG have 989 nodes, 4 edges, 24 colours, and `topology_unchanged` true. Twenty-four families have assembled bases. The other 76 families are absent from a length-weighted contig profile, which is why assembly recall is 0.24 while the read classifiers reach 0.82 and 0.91.

The read budget is 100000 paired fragments, seed 42, the same budget as the twenty-genome example, spread across 100 genomes. Eukaryotic pairs in the pin range from about 4.7 Mb to about 82 Mb.

## Read fill from the assembly graph

The table above weights contig length. Filling an unclassified read from the contig it maps to is a different score. The k141 FASTG has 4 edges, so neighbour agreement barely moves Kaiju family L1, from 1.2605 to 1.24614. The k21 graph has 149792 nodes and 25674 edges. Its neighbour fill reaches 1.23956. The 4-mer graph on the final contigs reaches 1.2089. At 1000000 reads the same Kaiju baselines are 1.26107, the k141 neighbour fill is 1.13009, and the 4-mer graph is 0.835244 (`examples/high100_enriched/RESULTS.md`).

| method | family L1 | unclassified read fraction |
| --- | ---: | ---: |
| kraken2 | 1.6323 | 0.81317 |
| kaiju | 1.2605 | 0.61366 |
| kraken2_k141_neighbour | 1.61056 | 0.79923 |
| kaiju_k141_neighbour | 1.24614 | 0.60173 |
| kraken2_k21_neighbour | 1.6111 | 0.8022 |
| kaiju_k21_neighbour | 1.23956 | 0.60137 |
| kraken2_4mer_graph | 1.56398 | 0.64707 |
| kaiju_4mer_graph | 1.2089 | 0.47124 |
