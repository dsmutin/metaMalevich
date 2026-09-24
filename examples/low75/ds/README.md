# Genus errors on low75 and low75half

`genus_graph.py` scores Kraken2 and Kaiju read profiles at genus, then fills unclassified reads from contigs labelled by `composition_genus_graph`. The taxdump is `/mnt/tank/scratch/partition-metagenomics/databases/taxdump/nodes.dmp`, walked with Samovar `NCBITaxonomyParser.get_ancestor_by_rank`.

```bash
python examples/low75/ds/genus_graph.py
```

Charts: `genus_l1` and `genus_error_share`. The half-strain directory receives the same tables.
