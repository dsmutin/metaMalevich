# Held-out genera comparison

Numbers below are from this run. `work/` is not in git. The tables are copied from `work/metrics.tsv` and from counts on `work/classify/kraken2.output`, `work/reprofile/contigs.paf`, and `work/reprofile/kraken_calls.tsv`.

## What was built

MEGAHIT 1.2.9 (`--presets meta-sensitive`) wrote `final.contigs.fa`: 6438 contigs, 2632196 bp, minimum 295 bp, maximum 5100 bp, average 408 bp, N50 390 bp (`work/megahit/final.log`).

The coloured graph is not that file. `megahit_toolkit contig2fastg 141` on `k141.contigs.fa` produced the FASTG. Forward records only: 6439 nodes, 2632341 bp, 24 overlap edges. MetaMetro `export_tocumg` wrote a CFA and a CDBG with those 6439 nodes and 24 edges, 18 colours, and `topology_unchanged` true (`work/reprofile/tocumg/export.json`). The CFA metadata field `graph_type` is the library's fixed label. This graph is `megahit_fastg`.

## Scores against the species read fractions

Truth is the lognormal read fraction rolled to species. The two Shigella flexneri pairs add to taxid 623. Kraken2 and Kaiju profiles roll the same way: a call at genus or coarser becomes taxon 0, so taxon 0 is larger than the classifier's own unclassified line.

| method | weight | l1 | r2 | presence_f1 | node_f1 | node_accuracy |
| --- | --- | --- | --- | --- | --- | --- |
| kraken2 | read_count | 0.7058 | 0.0786308 | 1 |  |  |
| kaiju | read_count | 0.7126 | 0.131816 | 1 |  |  |
| initial_colouring | node_length | 0.629774 | 0.233402 | 0.972973 | 0.408121 | 0.473676 |
| probability_sum | node_length | 0.643788 | 0.251634 | 0.972973 | 0.399307 | 0.474297 |
| lca | node_length | 1.55331 | 0.00604007 | 0.882353 | 0.220458 | 0.21463 |
| gated_neighbour | node_length | 0.642243 | 0.253658 | 0.972973 | 0.399332 | 0.474453 |
| bayesian_edge | node_length | 0.642794 | 0.254571 | 0.972973 | 0.398747 | 0.472278 |

Coverage-weighted rows (`length_times_multi`, MEGAHIT `_cov_` on the node id) are in `work/metrics.tsv`. The lowest L1 in that file is `bayesian_edge_coverage` at 0.619239. Its node F1 is 0.398747, below `initial_colouring`. `initial_colouring_coverage` L1 is 0.63391, worse than length weighting. Coverage weight is not adopted as the library profile.

On samovar10, `probability_sum` is the only hypothesis that improves both abundance L1 and species macro F1 relative to hard colouring. On this graph it does not: L1 rises from 0.629774 to 0.643788 and node F1 falls from 0.408121 to 0.399307. `gated_neighbour` and `bayesian_edge` move L1 by about 0.01 on 24 edges among 6439 nodes. No resolver parameter was changed.

## Where the estimates fail

Kraken2's report line `U` is 20766 of 100000 reads (0.20766). After the species roll, taxon 0 is 0.31262: those 20766 reads plus 10496 calls above species (6649 genus, 3839 family, 8 class). Every simulated species has a non-zero Kraken2 and Kaiju read fraction. Presence F1 is 1 for both.

The largest read-fraction errors are unclassified mass and Pseudomonas fluorescens (taxid 294): truth 0.28199, Kraken2 0.09015, Kaiju 0.1349.

Every FASTG node has a minimap2 `asm20` hit to a simulated genome. 1243818 of 2632341 bp have a Kraken species call equal to that genome's species (0.4725 of bases). The largest base blocks are fluorescens contigs:

| mapped species | Kraken species call | bases |
| --- | --- | --- |
| Pseudomonas fluorescens | unclassified | 413162 |
| Pseudomonas fluorescens | Pseudomonas fluorescens | 400769 |
| Pseudomonas fluorescens | Pseudomonas protegens | 188276 |
| Pseudomonas fluorescens | Pseudomonas putida | 69631 |

Fluorescens is 1147236 bp of the graph (0.4358 of mapped bases) against a read fraction of 0.28199. Length weighting therefore cannot recover the read fraction even on the contigs that are labelled fluorescens: `initial_colouring` estimates 0.15537805.

No node mapped to Streptococcus thermophilus. Its read fraction is 0.00212 (212 reads). Kraken2 estimates 0.002 and Kaiju 0.00203. Every assembly hypothesis estimates 0, which is the presence miss (recall 0.947368). Streptococcus mutans also has no mapped node, but `initial_colouring` still estimates 0.00013714029 from calls on contigs mapped to other genomes.

Charts: `figures/heldout_abundance.html` and `figures/samovar10_species_error.html`.
