# Colouring

TCA stores evidence. It does not pick a winning taxon.

The scored colouring reads Kraken2 whole-contig output: the classified taxid in column 3, and the `taxid:count` k-mer pairs. `cpp/kraken_count` sums those pairs. Pairs that are not numeric are skipped.

## Species rollup

K-mers at species or below (finer than `S`, including strains) roll to rank `S`. Genus-or-coarser k-mers and taxid 0 are not spread onto child species. After that drop, the remaining counts are normalized. That distribution is `probability_sum`.

A node with no species-level k-mers keeps the Kraken call when that call itself rolls to species. Otherwise the node is taxon 0.

The hard initial colouring is a one-hot of the same rolled call. It matches the shipped `strategy_a` table (probability fixed at 1), not the raw k-mer argmax. On `samovar10` the k-mer column is mostly unclassified and genus, so rolling those counts onto species would mark almost every base unclassified.

## Colour record

Each colour has `taxon_id`, `evidence_type`, `evidence_count`, `weight`, `score`, `confidence`, and `source`. A layer records `layer_id`, `source`, `method`, `database`, `database_version`, `parameters`, `taxonomy_version`, and `timestamp`.

Writing a layer names the operation: `replace`, `merge`, `intersect`, or `subtract`. `replace` on a non-empty table with an empty `layer_id` is refused. Aggregations that are named in code: `count`, `weighted_count`, `majority`, `lca`, `probability_sum`.

## Nodes and edges are separate

Node colours come from the species k-mer distribution above.

Edge colours are the renormalized minimum of the two endpoint distributions, so an edge can carry fewer taxa than either node. They are not a copy of a node colour. The original `strategy_b` table cannot be joined: its ids are `edge_*` and it has no source or target.

## Taxonomy

`taxonomy.parse_kraken_report` reads the Kraken report indentation as the parent stack. Each taxon keeps id, parent, rank, name, source, and version. The version stored for a bench is the report path inside the dataset bundle, not a NCBI release tag (the bundle does not ship one).

LCA of a set of taxa is the finest shared ancestor. Rank order treats `S1` as finer than `S`.
