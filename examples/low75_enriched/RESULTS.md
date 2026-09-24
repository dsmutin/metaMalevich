# low75_enriched

Same 75 genera as `examples/low75`, with 1000000 paired fragments instead of 100000. Scores are genus L1 on reads. The assembly graphs are MEGAHIT `contig2fastg` output: k141 is `assembly.fastg`, k21 is `initial_k21.fastg`.

| method | genus L1 | unclassified read fraction |
| --- | ---: | ---: |
| kraken2 | 1.54117 | 0.766586 |
| kaiju | 1.0121 | 0.492109 |
| kraken2_k141_neighbour | 1.05219 | 0.517713 |
| kaiju_k141_neighbour | 0.785012 | 0.374918 |
| kraken2_k21_neighbour | 1.13088 | 0.555312 |
| kaiju_k21_neighbour | 0.831512 | 0.395649 |
| kraken2_4mer_graph | 0.696196 | 0.171638 |
| kaiju_4mer_graph | 0.544312 | 0.111869 |

The k141 graph has 51535 nodes and 1292 edges. Neighbour agreement on that graph lowers Kaiju genus L1 from 1.0121 to 0.785012. At 100000 reads the same rule could not move labels: that FASTG had 24 edges. The k21 graph has 527112 nodes and 668818 edges. Its neighbour fill reaches Kaiju L1 0.831512, which is worse than the final-graph fill. The 4-mer graph on the final contigs is the lowest L1.

Of the 1292 final-graph edges, 1286 join contigs of the same true genus and 6 join different genera. Edges with Kraken labels on both ends agree 498 times and disagree 34 times. Edges with one label match the unlabelled contig's true genus 112 times and miss it 26 times. The shallow k21 graph had the opposite rescue balance, 750 matches against 1010 misses. Rows are in `ds/k141_edges.tsv`.

On the 100000-read community, Kaiju plus the 4-mer graph was 0.78304. Here the same rule is 0.544312.
