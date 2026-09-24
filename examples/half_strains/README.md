# Half of the held-out strains

Same simulation, databases, assembly, and reprofiling as `examples/heldout_genera`, on ten of the twenty downloaded strains. Genomes are not downloaded again and are not committed. Reads, indexes, and the assembly stay in `work/`.

## Which strains

`accessions.tsv` is every other pair in the order pairs first appear in `examples/heldout_genera/accessions.tsv`. That is ten strains and twenty assemblies. Each kept pair still puts one assembly in the metagenome and the other in the database.

| pair_id | genus |
|---------|-------|
| eco_albertii | Escherichia |
| eco_marmotae | Escherichia |
| eco_coli | Escherichia |
| shi_sonnei | Shigella |
| shi_dysenteriae | Shigella |
| pse_aeruginosa | Pseudomonas |
| pse_fluorescens | Pseudomonas |
| pse_protegens | Pseudomonas |
| str_pyogenes | Streptococcus |
| str_mutans | Streptococcus |

Escherichia and Pseudomonas contribute three pairs each. Shigella and Streptococcus contribute two. The other ten strains stay in the full example only.

## Run

Use the same programs and versions as the full example. From the repository root, after `examples/heldout_genera/run_example.py --stage download`:

```bash
python examples/half_strains/run_example.py
```

Read allocation is again `random.Random(42).lognormvariate(0, 1.5)` and 100000 paired fragments. The draw is for these ten genomes, not a slice of the twenty-genome abundance table. Kraken2, Kaiju (`fetch_missing=False` on the pinned database FASTA), MEGAHIT `meta-sensitive`, the k141 FASTG, and the current resolvers are unchanged. Measured scores are in `RESULTS.md`.
