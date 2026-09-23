# Resolution

Every resolver returns a distribution per node. Hard labels are produced only by the method that says so. Resolvers do not edit the evidence maps they read.

Parameters below are `resolve.RESOLVER_PARAMETERS`, which the bench writes into each `metrics.yaml`.

| Hypothesis | Graph | Parameters |
|------------|-------|------------|
| `initial_colouring` | no | One-hot of the Kraken call rolled to species |
| `probability_sum` | no | Species k-mer fractions. Genus and unclassified mass is dropped, not pushed onto children |
| `lca` | no | Taxa at or above 5% of the k-mer mass, then one node at their LCA, rolled to species. Coarser calls become taxon 0 |
| `gated_neighbour` | yes | 3 iterations. Confident 0.8, consensus 0.75, uncertain mix 0.7, conflict mix 0.5 |
| `bayesian_edge` | yes | 3 iterations. Neighbour weight 1, edge weight 1, floor `1e-6`. A taxon absent from the node is added only when the neighbour vote or the edge argmax is at least 0.75 |

## Gated neighbour

A node whose maximum probability is at least 0.8 keeps its distribution, unless the neighbour vote reaches 0.75 on a different taxon. Then it takes a 0.5 mix toward that vote. A node below 0.8 takes a 0.7 mix toward the neighbour vote.

## Edge Bayesian

The node score for a taxon already in its own support is the product of its own weight, the neighbour vote, and the incident edge colour, with a floor so one zero does not erase the others. New taxa are not imported through that floor.

## What was not run

The specification also lists majority, weighted majority, path-aware search, local subgraphs, a GCN, and a hybrid with that GCN. Majority is implemented as an aggregation, not as a scored hypothesis. The others are not in the bench. See [Pipeline](Pipeline).

On the toy graph (three nodes; the middle node is 0.6/0.4 toward the wrong species and both neighbours are the true species), neighbour resolution lowers L1 and the CLI requires `beats_initial`. That toy is not `samovar10`. The community result is on [Benchmark](Benchmark).
