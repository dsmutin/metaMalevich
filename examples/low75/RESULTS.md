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

## Genus rank, read baselines, and the 4-mer graph

On the k141 FASTG, 0.753 of bases have no Kraken genus and 0.013 are the wrong genus. All 24 overlap edges join contigs of the same true genus, and none joins an unclassified contig to a labelled one, so neighbour smoothing cannot move those bases.

`composition_genus_graph` labels each unclassified contig with the nearest Kraken-labelled contig under canonical 4-mer cosine. That raises the correct-genus base share from 0.234 to 0.688 and the wrong-genus share to 0.312, because every contig is forced onto a genus. The read profile does not do that. An unclassified Kraken2 or Kaiju read that maps to a contig inherits the contig genus; a read that already has a genus keeps it.

| method | genus L1 | unclassified read fraction |
| --- | ---: | ---: |
| kraken2 | 1.54014 | 0.76599 |
| kaiju | 1.01624 | 0.49412 |
| kraken2_4mer_graph | 1.20352 | 0.52264 |
| kaiju_4mer_graph | 0.78304 | 0.31074 |

Charts and the half-strain rows are in `ds/`. Twenty-nine of 75 genera have no assembled bases, so a length-weighted contig profile cannot represent them. The table above is read-weighted.

## Initial MEGAHIT graph (k=21)

`contig2fastg 21` on `k21.contigs.fa` is the first MEGAHIT graph. It has 122002 nodes and 26914 edges. The final k141 graph has 2685 nodes and 24 edges. Of the k21 edges, 0.775 join two unlabelled contigs. Among edges with one labelled end, 750 match the unlabelled contig's true genus and 1010 do not. Copying a neighbour genus when 80% of neighbour length agrees moves correct bases only from 0.181 to 0.183.

| method | genus L1 | unclassified read fraction |
| --- | ---: | ---: |
| kraken2_k21_neighbour | 1.46336 | 0.72692 |
| kaiju_k21_neighbour | 0.9688 | 0.46864 |

Requiring 4-mer cosine at least 0.5 on the same edges leaves Kaiju L1 at 0.96896. That is a smaller gain than the 4-mer graph on the final contigs. Tables are in `ds/initial_graph_*.tsv` and `ds/edge_4mer.tsv`.

A logistic model trained on `low75half` and applied here uses 4-mer cosine, length, MEGAHIT coverage, and Kraken k-mer support to pick one labelled neighbour. Training accuracy on the other community is 0.703. Held-out Kaiju genus L1 is 0.96896 and Kraken2 is 1.46358. Correct bases stay near 0.183. The features do not yet separate a true neighbour from a false one well enough to beat the 4-mer graph. Rows are in `ds/colour_features.tsv`.
