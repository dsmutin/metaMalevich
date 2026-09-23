# metamalevich wiki

Improve metagenomic taxonomic re-profiling with coloured assembly graphs.

**Status:** in development. Version: see repository file `VERSION`.

## Starting interconnections

```
Kraken2 whole-contig output
        │
        ▼
cpp/kraken_count  →  species rollup (taxonomy.py, aggregate.py)
        │
assembly contigs  →  cpp/kmer_knn  →  edge colours (resolve.edge_support)
        │
        ▼
MetaMetro colour_cfa → cfa_to_cdbg (ToCUMG)
        │
        ▼
resolvers: hard call, probability_sum, LCA, gated neighbour, edge Bayesian
        │
        ▼
length-weighted profile (reprofile.py) → benchmark/{hypothesis}/{dataset}/
```

## Pages

- [Contracts](Contracts)
- [Testing](Testing)
