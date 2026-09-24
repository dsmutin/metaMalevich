# low75 comparison

The database species is never the simulated species. Species-rank scores against the simulated read fractions are therefore disjoint. `work/metrics.tsv` records L1 2, Bray-Curtis 1, and presence F1 0 for Kraken2, Kaiju, and every graph hypothesis, including `initial_colouring`. Node F1 and node accuracy on the assembly are 0. R² is about 0.002 to 0.004. No resolver changes that, and none was adopted.

Kraken2 left 76157 of 100000 reads unclassified. Of the remaining reads, 23401 were placed in one of the 75 pinned genera. None was placed on a simulated species taxid. The largest species-rank Kraken fraction in `work/comparison.tsv` is unclassified 0.76599. The next calls are database-side species such as Acinetobacter calcoaceticus (0.04933), not the simulated species from that genus.

MEGAHIT wrote 2684 contigs, 1350101 bp, N50 465 bp. The k141 FASTG has 2685 nodes and 24 edges. The coloured CFA and CDBG have those counts, 35 colours, and `topology_unchanged` true.

The read budget is the same 100000 paired fragments used for twenty genomes in `examples/heldout_genera`, spread across 75 genomes.
