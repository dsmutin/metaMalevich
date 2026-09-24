# metamalevich features

Check a box only after mandatory tests pass.

- [x] Taxonomy, TCA colour layers, and explicit aggregations
- [x] Deterministic and edge-aware resolvers, length-weighted reprofile
- [x] MetaMetro CFA to ToCUMG export
- [x] Benchmark hypotheses on samovar10
- [x] Skip samovar10_ont1b; its inputs duplicate samovar10
- [x] Keep raw k-mer counts on colour rows written from the evidence graph
- [x] Roll chart truth abundances to species
- [x] Choose the LCA by tree depth rather than rank code
- [x] Refuse an all-pairs 4-mer graph above a few hundred contigs
- [x] State that profiles weight contig bases and do not use coverage
- [x] Regenerate the samovar10 benchmark after the metric fixes
- [x] Count unpredicted truth species as F1 0
- [x] Record 4-mer graph provenance without machine paths
- [x] Rebuild cached graphs and Kraken counts when inputs change
- [x] Run outside this source tree: pinned deps, compiler, no machine path
- [x] Reject mismatched contig ids and malformed Kraken rows
- [x] Treat a hypothesis as better only when L1 and macro F1 both improve
- [x] Record resolver parameters from the functions that use them
- [x] Report composition-neighbour agreement instead of path consistency
- [x] Test reprofiling, the k-mer graph, and ToCUMG topology
- [x] Held-out genera example: simulate, classify, assemble, and compare
- [x] Half-strain example: same pipeline on every other downloaded strain
- [x] Edge-graph flip, decaying leakage, and label drop on samovar10
- [x] Score the same hypotheses on held-out genera and the half-strain subset
- [x] Confident lock, unanimous rescue, edge cut, genus plurality, and margin mixture
- [x] Genus-rank 4-mer graph on low75 and low75half
- [x] 10×-depth low75_enriched and high100_enriched scored at genus and family against Kraken2 and Kaiju
- [ ] 10×-depth low75half_enriched and half100half_enriched, same scores after assembly
- [x] Learned colouring on the dense shallow graph: coverage and k-mer support do not separate true neighbours; at 10× depth the final overlap edges are same-taxon and neighbour fill plus 4-mers beat Kaiju
- [x] low75 and low75half: 25 genera each of bacteria, archaea, and viruses
- [x] high100 and half100half: 25 families each of bacteria, archaea, viruses, and small eukaryotes, scored at family rank
