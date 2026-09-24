# Half-strain comparison

Numbers are from this run. `work/` is not in git. Tables are copied from `work/metrics.tsv` and from `work/classify/kraken2.output`.

## Design

Ten strains, every other pair in `examples/heldout_genera/accessions.tsv`. One assembly of each pair was simulated. The other ten assemblies are the only Kraken2 and Kaiju sequences. Read counts are `random.Random(42).lognormvariate(0, 1.5)` over these ten genomes, summing to 100000 paired fragments:

| accession | reads |
| --- | ---: |
| GCF_003363755.1 | 42234 |
| GCF_000952955.1 | 20458 |
| GCF_001558215.1 | 18351 |
| GCF_000816985.1 | 7270 |
| GCF_001549955.1 | 4501 |
| GCF_001518855.1 | 2532 |
| GCF_000743015.1 | 1914 |
| GCF_002900365.1 | 1478 |
| GCF_002741615.1 | 721 |
| GCF_000730425.1 | 541 |

R1 and R2 each contain 100000 records, and the R1 count for each genome key equals the table. Samovar's snakefile still expects `iss/initial/1_full_R1.fastq`. InSilicoSeq wrote `sample_full_R1.fastq`. The runner keeps those files.

## Assembly graph

MEGAHIT `final.contigs.fa`: 10071 contigs, 4770815 bp, minimum 223 bp, maximum 6167 bp, average 473 bp, N50 460 bp.

The coloured graph is `contig2fastg 141` on `k141.contigs.fa`: 10072 forward nodes, 56 edges. MetaMetro export wrote a CFA and a CDBG with those counts, 10 colours, and `topology_unchanged` true.

## Scores

Truth is the species read fraction. Taxon 0 in a classifier profile includes calls above species.

| method | weight | l1 | r2 | presence_f1 | node_f1 | node_accuracy |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| kraken2 | read_count | 0.3851 | 0.717358 | 1 |  |  |
| kaiju | read_count | 0.50994 | 0.541619 | 1 |  |  |
| initial_colouring | node_length | 0.383564 | 0.853606 | 1 | 0.49223 | 0.7917 |
| probability_sum | node_length | 0.406141 | 0.846613 | 1 | 0.487334 | 0.792097 |
| lca | node_length | 1.09451 | 0.0885802 | 1 | 0.348321 | 0.477065 |
| gated_neighbour | node_length | 0.408551 | 0.844416 | 1 | 0.487507 | 0.791898 |
| bayesian_edge | node_length | 0.407528 | 0.845372 | 1 | 0.488553 | 0.793288 |
| edge_union | node_length | 0.407331 | 0.845515 | 1 | 0.487507 | 0.791898 |
| leakage | node_length | 0.406972 | 0.845817 | 1 | 0.487509 | 0.792196 |
| leakage_flipped | node_length | 0.407209 | 0.84546 | 1 | 0.487507 | 0.791898 |
| label_drop | node_length | 0.404144 | 0.847518 | 1 | 0.487334 | 0.792097 |
| label_drop_flipped | node_length | 0.407331 | 0.845515 | 1 | 0.487507 | 0.791898 |
| logistic_drop | node_length | 0.379475 | 0.858672 | 1 | 0.487895 | 0.792593 |
| logistic_drop_flipped | node_length | 0.407398 | 0.845503 | 1 | 0.487507 | 0.791898 |

Coverage weighting raises L1. `initial_colouring_coverage` is 0.455632 against 0.383564 for length weighting. It was not adopted.

`logistic_drop` has the lowest length-weighted L1 of the graph rows (0.379475 against 0.383564 for hard colouring). Its node F1 is 0.487895, below hard colouring. The other graph rows stay about 0.02 worse on L1 than hard colouring. With 56 edges on 10072 nodes, neighbour updates do not improve both scores. No resolver parameter was changed.

## Where the estimates fail

Kraken2 report line `U` is 10798 of 100000 reads. After the species roll, taxon 0 is 0.18181: those reads plus 7383 calls above species (6488 family, 895 genus). Every simulated species has a non-zero Kraken2 and Kaiju fraction, and every assembly hypothesis does too (presence F1 is 1).

The largest read-fraction misses are Pseudomonas protegens (truth 0.42234, Kraken2 0.30966, Kaiju 0.33582) and Escherichia coli (truth 0.20458, Kraken2 0.14185, Kaiju 0.07835). On the length-weighted graph, coli is 0.14487 and protegens is 0.39878. Streptococcus mutans is over-represented in the assembly profile (truth 0.18351, hard colouring 0.25703). Pseudomonas aeruginosa is under-represented (truth 0.07270, hard colouring 0.01936).

Chart: `figures/heldout_abundance.html`.
