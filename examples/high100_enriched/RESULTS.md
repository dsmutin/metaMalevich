# high100_enriched

Same 100 families as `examples/high100`, with 1000000 paired fragments instead of 100000. Scores are family L1 on reads. The assembly graphs are MEGAHIT `contig2fastg` output.

| method | family L1 | unclassified read fraction |
| --- | ---: | ---: |
| kraken2 | 1.63308 | 0.814037 |
| kaiju | 1.26107 | 0.613819 |
| kraken2_k141_neighbour | 1.40679 | 0.687671 |
| kaiju_k141_neighbour | 1.13009 | 0.533185 |
| kraken2_k21_neighbour | 1.51839 | 0.753645 |
| kaiju_k21_neighbour | 1.20249 | 0.579648 |
| kraken2_4mer_graph | 1.03465 | 0.349426 |
| kaiju_4mer_graph | 0.835244 | 0.257421 |

The k141 graph has 58137 nodes and 620 edges. Neighbour agreement lowers Kaiju family L1 from 1.26107 to 1.13009 and Kraken2 from 1.63308 to 1.40679. The k21 graph has 876751 nodes and 823034 edges. Its neighbour fill is weaker than the final graph: Kaiju L1 1.20249. The 4-mer graph on the final contigs is the lowest L1.

Of those 620 edges, 614 join contigs of the same true family and 6 join different families. Edges with Kraken labels on both ends agree 206 times and disagree 16 times. Edges with one label match the unlabelled contig's true family 96 times and miss it 18 times. Rows are in `ds/k141_edges.tsv`.

A cosine gate of 0.5 or 0.7 does not change the 4-mer L1. At 0.95 Kaiju family L1 rises from 0.835244 to 1.03135. Letting the overlap graph override the 4-mer label moves that L1 from 0.835244 to 0.835936.
