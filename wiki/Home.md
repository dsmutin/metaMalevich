# **metamalevich**

*Taxonomic re-profiling on a coloured assembly graph*

metamalevich projects Kraken2 k-mer evidence onto contig nodes and composition edges, keeps multi-label colours until a named resolver runs, and writes a length-weighted species profile. MetaMetro stores the coloured graph. The scientific method is this repository.

Version is the file `VERSION` in the code repository (currently read by the package; do not copy the number into other files).

```bash
conda env create -f environment.yml
conda activate metamalevich
python -m metamalevich --version
python examples/toy/run.py
```

Source tree: [github.com/dsmutin/metaMalevich](https://github.com/dsmutin/metaMalevich). Cite Kraken 2 for the classifier (Wood, Lu, and Langmead, Genome Biology, 2019, doi:10.1186/s13059-019-1891-0) and [MetaMetro](https://github.com/dsmutin/MetaMetro) for the graph containers.

---

## Where to look

### Setup

| Page | What is there |
|------|----------------|
| [How to install](How-to-install) | Conda env, `PYTHONPATH=src`, C++ compiler, MetaMetro checkout, CLI. |
| [Testing](Testing) | Mandatory pytest, toy, `samovar10` bench, what is not a scored dataset. |

### Pipeline

| Page | What is there |
|------|----------------|
| [Pipeline](Pipeline) | Spec stages 1–12: which are in this tree, which are not. |
| [Colouring](Colouring) | TCA: Kraken2 counts, species rollup, colour layers, node vs edge colours. |
| [Assembly graph](Assembly-graph) | Canonical 4-mer cosine kNN. Why 21-mer Jaccard was not used. |
| [Resolution](Resolution) | Hard call, `probability_sum`, LCA, gated neighbour, edge Bayesian. |
| [Reprofiling](Reprofiling) | Length-weighted profile and the unique / shared / ambiguous / unclassified split. |

### Results

| Page | What is there |
|------|----------------|
| [Benchmark](Benchmark) | `samovar10` scores. A hypothesis beats the hard colouring only when L1 and macro F1 both improve. |
| [Contracts](Contracts) | Invariants and the code that implements each contract. |

---

## Typical path

1. [How to install](How-to-install)
2. `python examples/toy/run.py` — three nodes; the neighbour mix must beat the hard colouring on L1
3. `python -m metamalevich bench --data-root . --benchmark benchmark` when `samovar10/` is present
4. [Benchmark](Benchmark) for the recorded scores
5. [Pipeline](Pipeline) for the stages from the specification that are still absent (read reassignment, path resolver, GCN, tier-2/3 baselines)
