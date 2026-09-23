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
- [ ] Half-strain example: same pipeline on every other downloaded strain
