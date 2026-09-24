# low75half_enriched

Every other low75 genus, at 1000000 paired fragments instead of 100000. Scores are genus L1 on reads. The graphs are MEGAHIT `contig2fastg` output.

| method | genus L1 | unclassified read fraction |
| --- | ---: | ---: |
| kraken2 | 1.56177 | 0.776901 |
| kaiju | 1.13852 | 0.560632 |
| kraken2_k141_neighbour | 0.88787 | 0.423337 |
| kaiju_k141_neighbour | 0.717654 | 0.331246 |
| kraken2_k21_neighbour | 1.14558 | 0.563693 |
| kaiju_k21_neighbour | 0.911108 | 0.439982 |
| kraken2_4mer_graph | 0.48872 | 0.07285 |
| kaiju_4mer_graph | 0.42235 | 0.044543 |

The k141 graph has 28455 nodes and 6578 edges. At 100000 reads the same community had 4119 nodes and 296 edges, and those edges did not connect an unclassified contig to a labelled one. Neighbour agreement now lowers Kaiju genus L1 from 1.13852 to 0.717654. The k21 graph has 356823 nodes and 558598 edges. Its neighbour fill is weaker: Kaiju L1 0.911108. The 4-mer graph on the final contigs is the lowest L1. At 100000 reads that rule was 0.75986.
