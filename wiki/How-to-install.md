# How to install

Python 3.12. Conda is the only supported install. There is no pip-first path.

## Core install

```bash
git clone https://github.com/dsmutin/metaMalevich
cd metaMalevich
conda env create -f environment.yml
conda activate metamalevich
python -m metamalevich --version
python examples/toy/run.py
```

`environment.yml` sets `PYTHONPATH=src`, so commands are run from the repository root. The env pins Python, pytest, numpy, pyyaml, pandas, Altair, and `gxx_linux-64`. The compiler builds `cpp/kraken_count` and `cpp/kmer_knn` on first use. Binaries are gitignored.

Check

```bash
python -c "import metamalevich; print(metamalevich.__version__)"
pytest -m mandatory
```

### What the env does

| Piece | Role |
|-------|------|
| `PYTHONPATH=src` | Import `metamalevich` without installing a console script |
| `python -m metamalevich` | CLI (`--version`, toy JSON, `bench` / `pipeline`) |
| `gxx_linux-64` | Compile the Kraken counter and the 4-mer graph |
| Altair | `benchmark/summary/*.html` |

## MetaMetro

ToCUMG export calls MetaMetro `colour_cfa` then `cfa_to_cdbg`. The bench looks for a directory that contains `metametro/__init__.py`:

| Order | Path |
|-------|------|
| 1 | `$METAMETRO_SRC` |
| 2 | `external/MetaMetro/src` under the repository root |
| 3 | `external/MetaMetro/src` under the current directory |

```bash
git clone --depth 1 https://github.com/dsmutin/MetaMetro external/MetaMetro
```

`external/` is gitignored. The bench fails if none of those trees import.

## Bench command

Needs a dataset bundle next to the repo (not part of the git tree):

```bash
python -m metamalevich bench --data-root . --benchmark benchmark
python -m metamalevich bench --data-root . --dataset samovar10 --intermediate intermediate
```

| Flag | Default | Meaning |
|------|---------|---------|
| `--data-root` | `.` | Directory that contains `samovar10/` |
| `--dataset` | every name in `DATASETS` | Repeat to select one bundle |
| `--benchmark` | `benchmark` | Hypothesis profiles and the summary |
| `--intermediate` | `intermediate` | Reused Kraken counts, edges, CFA, CDBG |
| `-o` | stdout | JSON summary path |

`pipeline` is the same entry as `bench`. With no command, the CLI runs the toy and prints JSON.

## Not installed by this recipe

Kraken2, Bracken, Flye, and the assemblers named in the specification are not dependencies. The bench reads Kraken2 output that is already on disk. Coverage tools are not used: relative abundance is the share of contig bases.
