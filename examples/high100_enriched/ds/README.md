# Family scores for high100_enriched

`score_family.py` writes `family_metrics.tsv`. Neighbour rows copy a family along MEGAHIT overlap edges when 80% of neighbour length agrees. The 4-mer rows use the nearest Kraken-labelled final contig. `k141_edges.py` counts whether those overlap edges join the same true family.
