# half100half comparison

Every other high100 family (50 of 100), scored at family rank. The read draw is a new seed-42 allocation of 100000 paired fragments. FASTA files are the high100 download.

All 17 assembly hypotheses were scored on the k141 FASTG (2532 nodes, 0 edges), plus Kraken2 and Kaiju on the reads.

| method | weight | l1 | presence_f1 | presence_precision | presence_recall | node_f1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| kraken2 | read_count | 1.67738 | 0.876404 | 1 | 0.78 |  |
| kaiju | read_count | 1.34148 | 0.93617 | 1 | 0.88 |  |
| initial_colouring | node_length | 1.71942 | 0.717949 | 1 | 0.56 | 0.229508 |
| gated_neighbour, bayesian_edge | node_length | 1.71795 | 0.717949 | 1 | 0.56 | 0.232201 |
| probability_sum | node_length | 1.71901 | 0.717949 | 1 | 0.56 | 0.232201 |
| genus_plurality | node_length | 1.71929 | 0.717949 | 1 | 0.56 | 0.232201 |
| lca | node_length | 1.92501 | 0.305085 | 1 | 0.18 | 0.0632523 |

The other non-LCA hypotheses match `probability_sum` on the printed L1, presence F1, and node F1, except `label_drop` and `logistic_drop` (L1 1.71902). Coverage-weighted rows are in `work/metrics.tsv`. The lowest coverage-weighted L1 is 1.87727, still above the length-weighted rows, with the same presence F1.

`gated_neighbour` and `bayesian_edge` lower L1 by 0.00147 and raise node F1 from 0.229508 to 0.232201. Presence F1 does not change. The FASTG has no edges, so that gap is the k-mer argmax against the Kraken call. No resolver was changed.

Kraken2 classified 16371 of 100000 read pairs (16.37%). After the family roll, the unclassified fraction in `work/comparison.tsv` is 0.83809. All 39 non-zero Kraken2 families are simulated families. Kaiju's unclassified fraction is 0.66362, with 44 families and precision 1.

MEGAHIT wrote 2532 contigs, 1050150 bp, N50 381 bp. The k141 FASTG has 2532 nodes and 0 edges. The coloured CFA and CDBG have those counts, 28 colours, and `topology_unchanged` true. Assembly recall is 0.56 (28 of 50 families). The read classifiers reach 0.78 and 0.88.
