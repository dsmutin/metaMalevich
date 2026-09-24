# metaMalevich <img src="data/metaMalevich.png" align="right" width="150" alt="metaMalevich logo">
### Taxonomic re-profiling on a coloured assembly graph

[![version](https://img.shields.io/badge/dynamic/yaml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fdsmutin%2FmetaMalevich%2Fmaster%2FVERSION&query=%24&label=version&color=blue)](VERSION)
[![required tests](https://img.shields.io/github/actions/workflow/status/dsmutin/metaMalevich/required-tests.yml?branch=master&label=required%20tests)](https://github.com/dsmutin/metaMalevich/actions/workflows/required-tests.yml)
[![full tests](https://img.shields.io/github/actions/workflow/status/dsmutin/metaMalevich/full-tests.yml?branch=master&label=full%20tests)](https://github.com/dsmutin/metaMalevich/actions/workflows/full-tests.yml)
[![license](https://img.shields.io/github/license/dsmutin/metaMalevich)](LICENSE)
[![conda](https://img.shields.io/badge/conda-environment.yml-44A833?logo=anaconda&logoColor=white)](environment.yml)

metaMalevich re-profiles a metagenome after assembly. A k-mer classifier labels contigs, but those labels stay multi-label colours on the assembly graph until a named resolver turns them into one distribution per contig. The profile is the length-weighted sum of those distributions. The aim is a community profile closer to the simulated truth than the hard classifier call, at the taxonomic rank the community was built to test.

The Python package and the conda environment are named `metamalevich`. Concepts and the scored tables are in the [wiki](https://github.com/dsmutin/metaMalevich/wiki).

## How it works

1. Read whole-contig Kraken2 output. `cpp/kraken_count` sums `taxid:count` pairs. Counts at species or below roll to species. Genus-or-coarser k-mers and taxid 0 are dropped, not pushed onto child species. The hard colouring is a one-hot of the rolled Kraken call.
2. Build the graph the resolver is allowed to use. The `samovar10` bench rebuilds a canonical 4-mer cosine k-nearest-neighbour graph (`top_k=8`, `min_sim=0.15`) with `cpp/kmer_knn`. Example communities use the MEGAHIT assembly graph (`megahit_toolkit contig2fastg` on the kept `k<int>.contigs.fa`). That FASTG is the assembler edge list. A second graph, when one is built, is computed from those contig sequences and is labelled as such.
3. Store node colours and edge colours separately. An edge colour is the renormalized minimum of the two endpoint distributions. [MetaMetro](https://github.com/dsmutin/MetaMetro) writes the coloured graph as a CFA and then a CDBG. Topology is checked after colouring.
4. Run a named resolver. Each one returns a full distribution per node. Graph-free methods are the hard call, `probability_sum`, and LCA. Graph methods include gated neighbour smoothing, an edge Bayesian update, decaying leakage, label drop, and a small set of lock, rescue, and mixture rules. Parameters are recorded with the hypothesis.
5. Write a length-weighted profile. Relative abundance is the share of contig bases. Coverage is not an input. Nodes are split into unique, shared, ambiguous, and unclassified. On `samovar10`, a hypothesis beats the hard colouring only when L1 falls and macro F1 rises. Example communities report the same hypotheses at the rank in the table below.

```text
Kraken2 whole-contig output
        │
        ▼
species k-mer colours on contig nodes
        │
assembly graph (MEGAHIT FASTG, or 4-mer cosine kNN on samovar10)
        │
        ▼
edge colours + MetaMetro CFA / CDBG
        │
        ▼
named resolver → length-weighted profile
```

## Install

Conda is the supported install. Python dependencies and `gxx_linux-64` are pinned in `environment.yml`. The compiler builds `cpp/kraken_count` and `cpp/kmer_knn` on first use.

```bash
git clone https://github.com/dsmutin/metaMalevich
cd metaMalevich
conda env create -f environment.yml
conda activate metamalevich
```

`environment.yml` sets `PYTHONPATH=src`. Run commands from the repository root. ToCUMG export needs a MetaMetro checkout:

```bash
git clone --depth 1 https://github.com/dsmutin/MetaMetro external/MetaMetro
```

The bench also accepts `METAMETRO_SRC`. Details: [How to install](https://github.com/dsmutin/metaMalevich/wiki/How-to-install).

## Usage

```bash
python -m metamalevich --version
python -m metamalevich
python examples/toy/run.py
python -m metamalevich bench --data-root . --benchmark benchmark
```

With no command, the CLI runs a three-node toy and prints JSON. `bench` reads Kraken2 output already stored for `samovar10`, builds the 4-mer graph, colours nodes and edges, and writes each hypothesis under `benchmark/{hypothesis}/{dataset}/`. Reusable counts and edges go to `./intermediate`. `samovar10_ont1b` repeats the same FASTA, Kraken2 output, and ground truth, so it is not a second dataset. On `samovar10`, `genome_abundance` matches the base share, so the benchmark scores how bases are classified.

## Examples

Genomes, reads, and indexes are not committed. Each directory's README is the protocol. Scores are in that directory's `RESULTS.md`, at the rank the design can still match.

| Example | Score rank | Design |
|---------|------------|--------|
| `examples/toy` | — | Three nodes. Neighbour mixing must beat the hard colouring on L1. |
| `samovar10`, `examples/heldout_genera`, `examples/half_strains` | species | The simulated species is in the reference. |
| `examples/low75`, `examples/low75half`, `examples/low75_enriched`, `examples/low75half_enriched` | genus | The database species is a different species of the same genus. Enriched runs use the same genera at 10× read depth. |
| `examples/high100`, `examples/half100half`, `examples/high100_enriched`, `examples/half100half_enriched` | family | Communities are paired inside a family, including pairs from different genera. Enriched runs use the same families at 10× read depth. |

`examples/heldout_genera/` builds twenty simulated assemblies and uses twenty held-out assemblies as the Kraken2 and Kaiju databases. `examples/half_strains/` repeats that pipeline on every other strain. `examples/low75/` uses 25 genera each of bacteria, archaea, and viruses. `examples/low75half/` keeps every other of those genera. `examples/high100/` uses 25 families each of bacteria, archaea, viruses, and small eukaryotes, two genomes per family. `examples/half100half/` keeps every other of those families.

## Tests

```bash
pytest -m mandatory    # every commit
pytest                 # mandatory and optional
```

## License

MIT. See [LICENSE](LICENSE).

Cite Kraken 2 for the classifier (Wood, Lu, and Langmead, Genome Biology, 2019, doi:10.1186/s13059-019-1891-0) and [MetaMetro](https://github.com/dsmutin/MetaMetro) for the coloured-graph containers.
