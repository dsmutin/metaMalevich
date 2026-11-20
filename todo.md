# metamalevich features

Current expectation of the implementation. Check a box only after mandatory tests pass and a score exists. A checked row is in this tree. An empty box is specified or tried, and is not a result.

`benchmark/*/samovar10/metrics.yaml` records `software_version: 0.4.1`. The file `VERSION` is newer. Those committed numbers were not regenerated at the current version. Example scores live in each `examples/*/RESULTS.md`. Re-running an example needs the gitignored genomes and the tools named in that example's README.

## Baseline: annotators alone

- [x] Kraken2 and Kaiju on the reads, scored on the example communities against the simulation abundance table
- [x] `initial_colouring`: one-hot of the Kraken contig call, rolled to species
- [x] `probability_sum`: species k-mer fractions, no graph
- [x] `lca`: taxa at or above 5% of the k-mer mass, then their LCA
- [x] Colourings come from MetaMetro ToCUMG namespaces (`--colouring`); decaying leakage is the MetaMetro update
- [x] `solve` builds a MetaMetro bench with `kraken2`, `composition_kmeans`, and `decaying`, and runs `--resolver` (default `genus_plurality`; `gated_neighbour`, `probability_sum`, `hard`) on the `kraken2` taxa, checked against an NCBI TAXDUMP (#254). It writes per-unitig assignments and a length-weighted profile. A MetaMetro Kraken2 pin counts as Kraken2 input (`phage_10`), and its pinned TAXDUMP is used when none is given; without a given or pinned input a bench is skipped with a reason. Full CI runs it for `bubble_strain_2` and `synthetic_reads_2`, which are skipped because their toy sequences have no Kraken2 pin, and benchbuilds `phage_10` (`phage_species_5_x10`). The job fetches the pinned phage genomes with `metametro benchfetch phage_10`. The job installs MEGAHIT 1.2.9 and SAMOVAR, so `phage_10` is assembled and solved from the MetaMetro Kraken2 and TAXDUMP pins, and the job checks that its row is built and solved. Without those programs `phage_10` records its contract and the errors MetaMetro raised, and the run is not ok (#275, #233)
- [ ] `solve` metrics against gold (#280 open). `metrics` now compares the resolved argmax with gold taxon ids exactly, with no rank rollup, so a genus or below-species call never matches a species gold label (`phage_10`: 0 on every metric, since its calls are genus or coarser, or below species). Without gold, or with gold names that are not taxon ids, `metrics` is null with a reason. Rank-aware scoring is #1173

## Colour the graph, then add to that colouring

Node colours are the species k-mer distribution. Edge colours are the renormalized minimum of the two endpoints. MetaMetro stores the mask (`colour_cfa`, then `cfa_to_cdbg`). Posterior weights stay in the hypothesis tables.

- [ ] Bare annotators plus coverage on the node. Coverage is unused in the profile. `examples/half_strains/RESULTS.md` records a coverage-weighted L1 and does not adopt it. There is no coverage colour hypothesis in `bench.HYPOTHESES`.
- [x] Edge colour, and resolvers that read it: `edge_union`, `bayesian_edge`, `gated_neighbour`
- [x] Node k-mer composition. `samovar10` rebuilds a canonical 4-mer cosine kNN (`top_k=8`, `min_sim=0.15`). `composition_genus` copies a genus label across that space. Example communities keep the MEGAHIT FASTG as the assembler edge list; a 4-mer graph computed from those contigs is a second graph.
- [x] Decaying signal from neighbouring nodes: `leakage` (decay 0.5, 4 iterations)
- [x] Line-graph flip (nodes and edges swap): `flip.py`, then `leakage_flipped`, `label_drop_flipped`, `logistic_drop_flipped`
- [x] Logistic filter: `logistic_drop` fits pseudo-labels from neighbour agreement on the same graph. Ground truth is not a feature.
- [x] Further resolvers already scored: `label_drop`, `confident_lock`, `unanimous_rescue`, `cut_then_leak`, `genus_plurality`, `margin_mixture`

## Reprofile

The profile is the length-weighted sum of the resolved distributions. `assigned_reads` is 0.

- [x] Bayesian update on the edge colours: `bayesian_edge`
- [x] Machine learning in the library: logistic label drop only (`label_drop.logistic_label_drop`)
- [ ] GCN, GraphSAGE, and GAT. The pipeline page records them as not started.
- [ ] Path-aware resolver, local subgraphs, and read-to-graph reassignment

## Checks that are not met yet

- [ ] Abundance scores read the taxon id assigned to each read before simulation (`samovar10/reads/simulation_manifest.tsv` or the example abundance table). The `samovar10` bench scores `ground_truth/ground_truth.csv` and does not open the read files. Example node F1 uses minimap2.
- [ ] The simulated taxon id stays out of a GCN graph and its colours. No GCN is built. `examples/low75/ds/colour_features.py` fits a logistic model on minimap-derived genus truth.
