# Genus errors on low75 and low75half

`genus_graph.py` scores Kraken2 and Kaiju read profiles at genus, then fills unclassified reads from contigs labelled by `composition_genus_graph`. Set `TAXDUMP` to `nodes.dmp` (or its directory) and `SAMOVAR_SRC` to the Samovar checkout. Parents are walked with `NCBITaxonomyParser.get_ancestor_by_rank`. If either path is missing, the script stops.

```bash
python examples/low75/ds/genus_graph.py
```

Charts: `genus_l1` and `genus_error_share`. The half-strain directory receives the same tables.

`initial_graph.py` scores neighbour transfer on the k=21 MEGAHIT FASTG. `edge_4mer.py` keeps a k21 neighbour only when canonical 4-mer cosine is at least 0.5. `colour_features.py` fits a logistic model on one community and applies it to the other. Tables are copied under `examples/low75half/ds/`.
