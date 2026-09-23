# Held-out genera

One simulated metagenome from twenty assemblies. The other twenty assemblies, one different genome from each pair, are the only sequences in the Kraken2 and Kaiju databases. Nothing under `data/` or `work/` is committed.

## What is pinned

`accessions.tsv` lists 40 NCBI Assembly accessions: five pairs in each of Escherichia, Shigella, Pseudomonas, and Streptococcus. Four Escherichia pairs are not E. coli. A pair is two assemblies of the same species with different sequence lengths. Where NCBI had a RefSeq reference, that assembly is the database genome and the other strain is simulated. `Pseudomonas stutzeri` is not included: the NCBI summary for that taxid returned `Stutzerimonas stutzeri`.

The read allocation is `random.Random.lognormvariate(0, 1.5)` with seed 42, then rounded so the 20 genomes sum to 100000 and none is zero. `100000` is the Samovar `--total_reads` value: paired fragments, not fold coverage. Across these genome sizes that sample is well below 1× for most genomes.

## Run

Install this repository (`environment.yml`) and Samovar (`git@github.com:ctlab/samovar.git`). The run also needs these programs on `PATH`:

| Program | Version used while preparing this example |
|---------|-----------------------------------------------|
| NCBI datasets | 18.36.0 |
| MEGAHIT | 1.2.9 (bioconda build `haf24da9_8`) |
| minimap2 | 2.31 (bioconda build `h118bc1c_0`), used only for contig-to-genome truth |
| Kraken2 | 2.0.7-beta |
| Kaiju | 1.10.1 |
| InSilicoSeq `iss` | 2.0.1 |
| snakemake | 9.12.0 |

MEGAHIT is not part of the base `environment.yml`. The build above was installed with:

```bash
conda install -n metamalevich -c conda-forge -c bioconda --freeze-installed -y megahit=1.2.9
```

From the repository root, with those programs on `PATH`:

```bash
python examples/heldout_genera/run_example.py
```

Stages resume. `--stage download` only fetches and validates assemblies. Samovar `generate` writes the InSilicoSeq pipeline; the script then runs `work/iss/.generate/generate.sh`. Samovar `build --type kraken2 --no-example-omit` indexes the database half. Its snakemake preprocess is given `mutation_rate: 0`; the Kraken2 index is built from the unmodified FASTA directory. `samovar build --type kaiju` downloads the latest RefSeq proteome for each taxid, which is not the pinned assembly, so Kaiju is built with Samovar's `add_database_kaiju(..., fetch_missing=False)` on those FASTA files (6-frame translation only).

Measured scores and the error breakdown are in `RESULTS.md`. MEGAHIT uses the documented `--presets meta-sensitive` and `--keep-tmp-files`. The assembly graph is `megahit_toolkit contig2fastg` on the highest `k<int>.contigs.fa` (not `k<int>.final.contigs.fa`, which the toolkit writes as an empty FASTG). Reverse-complement FASTG records are dropped. The example does not replace a missing assembly graph with a 4-mer kNN. When MetaMetro imports, the same graph is written as a coloured CFA and CDBG under `work/reprofile/tocumg/`. The CFA metadata field `graph_type` is the library's fixed label; the example records the graph as `megahit_fastg`.

## Scores

`work/metrics.tsv` compares Kraken2, Kaiju, and the current metaMalevich hypotheses to the simulated species read fractions. `r2` is the square of the Pearson correlation of those fractions. `presence_f1` is the set F1 of species with a non-zero fraction. Node F1 is filled when `minimap2` is on `PATH` (`-x asm20` when that preset is in the help text).

`plot_errors.py` writes Altair charts to `figures/`. The samovar10 charts use the committed benchmark table and do not need this download.
