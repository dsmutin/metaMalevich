# metamalevich

[![version](https://img.shields.io/badge/dynamic/yaml?url=https%3A%2F%2Fraw.githubusercontent.com%2Fdsmutin%2FmetaMalevich%2Fmaster%2FVERSION&query=%24&label=version&color=blue)](VERSION)
[![required tests](https://img.shields.io/github/actions/workflow/status/dsmutin/metaMalevich/required-tests.yml?branch=master&label=required%20tests)](https://github.com/dsmutin/metaMalevich/actions/workflows/required-tests.yml)
[![full tests](https://img.shields.io/github/actions/workflow/status/dsmutin/metaMalevich/full-tests.yml?branch=master&label=full%20tests)](https://github.com/dsmutin/metaMalevich/actions/workflows/full-tests.yml)
[![warning](https://img.shields.io/badge/warning-in%20development-yellow)](https://shields.io/badges/static-badge)

Improve metagenomic taxonomic re-profiling with coloured assembly graphs

**Warning: in development.** Interfaces may change. See `VERSION` (single source of truth).

## Install

Conda is the only supported install:

```bash
conda env create -f environment.yml
conda activate metamalevich
```

`environment.yml` sets `PYTHONPATH=src`. Do not publish a pip-first install path.

## Usage

```bash
python -m metamalevich --version
python -m metamalevich
python examples/toy/run.py
python examples/heldout_genera/run_example.py
python examples/half_strains/run_example.py
python -m metamalevich bench --data-root . --benchmark benchmark
```

`examples/heldout_genera/` is a from-scratch community: twenty simulated assemblies and twenty held-out assemblies used as the Kraken2 and Kaiju databases. `examples/half_strains/` repeats that pipeline on every other strain (ten pairs) and reuses the same downloaded FASTA files. Genomes, reads, and indexes are not committed. See each directory's README.

`bench` reads Kraken2 whole-contig output already stored for `samovar10`, builds a canonical 4-mer cosine k-nearest-neighbour graph (`top_k=8`, `min_sim=0.15`) with `cpp/kmer_knn`, colours nodes and edges, and writes each hypothesis under `benchmark/{hypothesis}/{dataset}/`. `samovar10_ont1b` repeats the same FASTA, Kraken2 output, and ground truth, so it is not a second dataset. Relative abundance is the share of contig bases. Coverage is not used. On `samovar10`, `genome_abundance` matches that base share, so the benchmark scores how bases are classified. Reusable counts and edges go to `./intermediate`. MetaMetro performs the CFA to ToCUMG step. Set `METAMETRO_SRC`, or clone MetaMetro to `external/MetaMetro`.

## Tests

```bash
pytest -m mandatory    # every commit
pytest                 # mandatory + optional (release / manual CI)
```

## License

MIT. See [CONTRIBUTING.md](CONTRIBUTING.md).
