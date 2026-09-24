# half100half_enriched

Every other high100 family, at 1000000 paired fragments instead of 100000. Scores are family L1 on reads. The graphs are MEGAHIT `contig2fastg` output.

| method | family L1 | unclassified read fraction |
| --- | ---: | ---: |
| kraken2 | 1.67935 | 0.838908 |
| kaiju | 1.3402 | 0.663717 |
| kraken2_k141_neighbour | 1.05056 | 0.501546 |
| kaiju_k141_neighbour | 0.92092 | 0.426942 |
| kraken2_k21_neighbour | 1.41503 | 0.697784 |
| kaiju_k21_neighbour | 1.20726 | 0.584535 |
| kraken2_4mer_graph | 0.713726 | 0.156717 |
| kaiju_4mer_graph | 0.635004 | 0.118426 |

The k141 graph has 31394 nodes and 1238 edges. Neighbour agreement lowers Kaiju family L1 from 1.3402 to 0.92092 and Kraken2 from 1.67935 to 1.05056. The k21 graph has 569445 nodes and 659448 edges. Its neighbour fill is weaker than the final graph: Kaiju L1 1.20726. The 4-mer graph on the final contigs is the lowest L1.
