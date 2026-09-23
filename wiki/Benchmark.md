# Benchmark

Scored dataset: `samovar10`. 399 contigs, 10 genomes, whole-contig Kraken2, ground truth in `ground_truth/ground_truth.csv`. On this bundle `genome_abundance` matches each genome's share of contig bases, so L1 compares two base-fraction profiles. Coverage is not in the score.

`samovar10_ont1b` repeats those inputs and is not run as a second dataset. `strong100` is not scored: Flye segment ids and the ground-truth contig ids are different sequences. Illumina directories have read-level Kraken2 only.

## Rule for a better reprofile

`compare_to_initial` sets `beats_initial` only when L1 is lower and macro F1 is higher than `initial_colouring`. A lower L1 alone is not enough.

Species-rank scores from `benchmark/summary/metrics.tsv`. Lower L1 is better. Macro F1 counts truth species that the profile never predicts as F1 0.

| Hypothesis | L1 | Macro F1 | Accuracy | Neighbour agreement | Unclassified fraction | Beats initial |
|------------|---:|---------:|---------:|--------------------:|----------------------:|:--------------|
| `initial_colouring` | 0.1586 | 0.913 | 0.915 | 0.884 | 0.0154 | — |
| `probability_sum` | 0.1441 | 0.918 | 0.925 | 0.893 | 0.0102 | yes |
| `gated_neighbour` | 0.1378 | 0.878 | 0.917 | 0.975 | 0.0017 | no |
| `bayesian_edge` | 0.1594 | 0.886 | 0.912 | 0.934 | 0.0053 | no |
| `lca` | 1.608 | 0.237 | 0.193 | 0.148 | 0.804 | no |

`probability_sum` is the reprofile that improves both abundance and species F1. `gated_neighbour` has the lowest L1 and the highest neighbour agreement. Neighbour agreement is the quantity neighbour smoothing directly increases, so it is not an independent check. Macro F1 falls.

Enterococcus sp. DIV2432 (taxid 2774762) has truth relative abundance 0.081. The hard call gives 0.031, `probability_sum` gives 0.038, and `gated_neighbour` gives 0.020. Composition neighbours of those contigs are other Enterococcus contigs already labelled *E. faecalis* (1351), so mixing amplifies that swap instead of correcting it.

`lca` and `bayesian_edge` stay in the bench as negative controls.

## Plots

Altair files in `benchmark/summary/`:

| File | Content |
|------|---------|
| `l1.html` | L1 by hypothesis |
| `f1.html` | Macro F1 by hypothesis |
| `abundance.html` | Estimated vs truth relative abundance |

## Manifest

Each `metrics.yaml` records the section 36 fields that this run actually has: input graph, 4-mer cosine method, k = 4, assembler note, TCA method, taxonomy report path, classifier `kraken2`, resolver name and parameters, profiling method, and software version. There is no random seed: the resolvers are deterministic. The manifest does not invent a NCBI taxonomy release or a read-mapper version.
