# Reprofiling

The profile is the length-weighted posterior. It does not replace a soft distribution with a hard label. Coverage is not used. Relative abundance is the share of contig bases, so unclassified bases stay in the denominator.

`assigned_reads` is 0. The scored bundle has no read-to-graph map, so read reassignment from the specification is not performed.

## Node partition

| Class | Rule |
|-------|------|
| ambiguous | maximum probability below 0.5 |
| shared | maximum at least 0.5 and the second taxon at least 0.1 |
| unique | maximum at least 0.5 and the second taxon below 0.1 |
| unclassified | empty distribution, or taxon 0 |

Shared nodes contribute to each taxon in proportion to the posterior. Unclassified mass is counted once: the taxon-0 portion and the residual `(1 - assigned) * length` are the same bases.

## Profile columns

`taxon_id`, `rank`, `name`, `estimated_abundance`, `relative_abundance`, `assigned_bases`, `assigned_reads`, `unique_bases`, `shared_bases`, `ambiguous_bases`, `confidence`.

`estimated_abundance` is assigned bases. The summary also reports `unclassified_bases` and `unclassified_fraction`.

## Uncertainty that is stored

The node distribution and the unique/shared/ambiguous/unclassified base counts. Posterior intervals, bootstrap intervals, and graph-perturbation intervals from the specification are not computed.
