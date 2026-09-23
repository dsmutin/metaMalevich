# Assembly graph

`samovar10/graph/graph_meta.txt` records `n_nodes=399`, `n_edges=6384`, `top_k=8`, `min_sim=0.15`. The edge list was not stored, and `mean_weight` is not a Jaccard (2.69 on `samovar10`). New edge ids are not the `strategy_b` ids.

## What is built

`cpp/kmer_knn` builds a canonical 4-mer frequency vector per contig and links each node to its top 8 neighbours with cosine similarity at least 0.15. Edges are symmetrized: both directions are kept, and the stored weight is the max of the two. The recorded graph has 399 nodes and 4450 directed edges (`benchmark/*/samovar10/metrics.yaml`).

A 21-mer minimizer Jaccard at the same `min_sim` kept 82 edges. Disjoint 100 kb contigs from one genome share few exact 21-mers, so that graph does not connect the contigs the composition graph connects.

Similarity is still all-pairs inside that top-k. `--max-nodes` defaults to 2000 and the builder exits above it. The error text says the scan is meant for a few hundred contigs, because the comparison is quadratic.

## Cache

`intermediate/{dataset}/knn_edges.tsv` is reused while a signature of the FASTA and the graph parameters still matches. Changing the FASTA or the parameters rebuilds the edges and the Kraken counts. Delete the intermediate directory to force a rebuild.

## What this graph is not

It is not the Flye repeat graph, not a compacted de Bruijn graph, and not the missing kNN whose weights were stored only as a mean. Contig sequences are still the nodes. MetaMetro then compacts that coloured graph with `cfa_to_cdbg` without changing which contigs exist when there is no overlap to merge. GC and entropy are written as CFA floats. Coverage is not an input.
