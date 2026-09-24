# low75half comparison

Thirty-eight of the seventy-five low75 genera. All 17 assembly hypotheses were scored on the k141 FASTG (4119 nodes, 296 edges). The database species is still not the simulated species.

| method | weight | l1 | presence_f1 | node_f1 | node_accuracy |
| --- | --- | ---: | ---: | ---: | ---: |
| kraken2 | read_count | 2 | 0 |  |  |
| kaiju | read_count | 2 | 0 |  |  |
| every assembly hypothesis | node_length | 2 | 0 | 0 | 0 |

The assembly list is the same as `examples/low75/RESULTS.md`. Bray-Curtis is 1. Kraken2 R² is 0.00463717. Coverage-weighted rows in `work/metrics.tsv` are also L1 2.

Genus rollup of the same argmax calls:

| method | genus node accuracy | genus base accuracy | genus L1 |
| --- | ---: | ---: | ---: |
| initial_colouring | 0.2731 | 0.2599 | 1.4800 |
| probability_sum and the other non-LCA hypotheses | 0.2736 | 0.2601 | 1.4781 |
| lca | 0.1440 | 0.1273 | 1.7704 |

MEGAHIT wrote 4090 contigs, 1881332 bp, N50 429 bp. The k141 FASTG has 4119 nodes and 296 edges. The coloured CFA and CDBG have those counts, 27 colours, and `topology_unchanged` true.

The read budget is again 100000 paired fragments, redrawn for these 38 genomes with seed 42.

## Genus rank and the 4-mer graph

Kraken leaves 0.724 of assembled bases without a genus. Wrong-genus bases are 0.016. The 296 overlap edges do not connect those unclassified contigs to a labelled contig. Nearest-neighbour transfer on canonical 4-mers raises the correct-genus base share from 0.260 to 0.827.

Read-weighted genus L1, with unclassified reads filled from the mapped contig:

| method | genus L1 |
| --- | ---: |
| kraken2 | 1.56126 |
| kaiju | 1.1342 |
| kraken2_4mer_graph | 1.09142 |
| kaiju_4mer_graph | 0.75986 |

The same table is in `examples/low75/ds/genus_metrics.tsv`.

## Initial MEGAHIT graph (k=21)

`contig2fastg 21` wrote 108039 nodes and 64134 edges. The final graph has 4119 nodes and 296 edges. Of the k21 edges, 0.867 join two unlabelled contigs. One-label edges match the true genus 1374 times and miss it 692 times. Neighbour agreement raises correct bases from 0.223 to 0.228.

| method | genus L1 | unclassified read fraction |
| --- | ---: | ---: |
| kraken2_k21_neighbour | 1.37908 | 0.68474 |
| kaiju_k21_neighbour | 1.01604 | 0.4983 |

A 4-mer cosine filter of 0.5 on those edges leaves Kaiju L1 at 1.0163. The final-contig 4-mer graph is still the larger gain.

The same logistic model, trained on `low75` and tested here, has training accuracy 0.609. Held-out Kaiju genus L1 is 1.01812 and Kraken2 is 1.38288. Correct bases are 0.225.
