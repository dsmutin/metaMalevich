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
