# Contributing to metamalevich

All development follows this guide. Project rules require it (`follow-contributing`).

English is required for every public function, class, module, CLI flag, data file, and user-facing document. Other languages or missing documentation are not allowed.

## Architecture

```
CLI (metamalevich.cli)
  → baseline pipeline (metamalevich.baseline.run_pipeline)
      → JSON result {status, ok, input_path}

tests/          mandatory vs optional pytest
examples/toy/   end-to-end run of the current (baseline) tool
cite/           BibTeX for integrated third-party tools
agents/         portable rules and skills (any IDE)
```

The wiki is not in this tree. It is [dsmutin/metaMalevich.wiki](https://github.com/dsmutin/metaMalevich/wiki), cloned beside this repository as `metaMalevich.wiki`.

Replace baseline bodies with real implementations. Keep the documented return keys until you change the contract on the wiki page `Contracts` and the tests together.

## Testing architecture

| Kind | Marker | Command | When |
|------|--------|---------|------|
| Required | `mandatory` | `pytest -m mandatory` | every commit; GitHub Action `required-tests` |
| Optional | `optional` | `pytest` (all) | release, workflow_dispatch, or an edit of Action `full-tests` |
| Light examples | — | `examples/*/run.py` and `vignettes/*` | full CI. `examples/toy/run.py` is the one in this tree |
| MetaMetro solve | — | `python -m metamalevich solve --bench bubble_strain_2 --bench synthetic_reads_2 --bench phage_10` | full CI. In-process benches are solved. `phage_10` is assembled only when its programs and genome FASTA files are already present |
| Communities | — | `examples/*/run_example.py` | not in CI; they download, classify, and assemble |

Do not mark a contract test `optional`. Optional tests are slow, extra, or nice-to-have.

After **any new feature**, run the **mandatory** suite (and toy if the CLI changed) before you stop.

## Benchmark truth

Abundance and read-level scores use the taxon id assigned to each read before the metagenome was generated. That record is the simulation design (an abundance table, or `samovar10/reads/simulation_manifest.tsv`). A mock taxon, a Kraken2 or Kaiju call, or a taxon inferred after assembly (minimap2, blobtools) is not the truth. If that pre-generation label is missing, stop.

The same id may be attached to a contig for a node score. The attachment does not invent a taxon.

## No simulated-taxon leakage

The simulated taxon id is a score label. It stays out of the graph a resolver or a GCN trains or predicts on: node features, edge features, node colours, and edge colours. A table written after inference may store the id beside the prediction. That table is not an input.

On a MetaMetro CGT, colours stay on the colour mask and labels stay on `node_labels`. Colours are not copied into the feature matrix. `examples/low75/ds/colour_features.py` fits a logistic model on minimap-derived genus truth. Do not copy that supervision into the graph colouring or into a GCN.

## Feature checklist (`todo.md`)

Track work in `todo.md` (checkboxes). One line per feature or fix. Check it off only when mandatory tests pass. This is a **feature list**, not a `/do` analysis graph.

## Versioning

Edit **only** `VERSION`. Everything else reads it.

Starting value: `0.0.1`.

| Change | Bump |
|--------|------|
| New feature | **minor** (`0.0.1` → `0.1.0`) |
| Fix or update of an existing feature | **patch** (`0.1.0` → `0.1.1`) |
| Release | **major** (`0.1.1` → `1.0.0`) |

## Install (conda only)

```bash
git clone --recurse-submodules https://github.com/dsmutin/metaMalevich
cd metaMalevich
conda env create -f environment.yml
conda activate metamalevich
```

`environment.yml` puts `external/MetaMetro/src` on `PYTHONPATH`. That submodule is [dsmutin/MetaMetro](https://github.com/dsmutin/MetaMetro). Benches and colouring methods are implemented there. Example scripts that walk NCBI read `TAXDUMP` and `SAMOVAR_SRC`. If either is unset, they stop.

## GitHub

Never `git push` unless the human explicitly asks. CI runs on GitHub after they push.

## Benchmarks

Do not add a new benchmark or pinned community in this repository. Add it in MetaMetro (`metametro benchbuild`) and open a pull request there. `examples/*/run_example.py` only forwards to that generator.

Do not hard-code a machine path. Use the MetaMetro bench directory, a CLI argument, or `METAMETRO_SRC`. If it is missing, stop.

Do not mock a benchmark graph, a taxonomy label, or a metric. A tiny graph that exists only inside a test file stays in that test. Example assembly graphs are MEGAHIT FASTG files, not substitutes.

Do not copy an evaluation target into graph features, colours, or any file a model reads as input. Ground truth stays in `ground_truth/` and is used only by the scorer. A table written after inference may store the id beside the prediction. That table is not an input.

## Colourings

Do not add a new colouring method in this repository. Add it in MetaMetro and open a pull request there. Resolvers (gated neighbour, decaying leakage as inference) stay here; they consume colours, they do not invent a palette.

Do not mock a colouring. Kraken2 and Kaiju colours come from those tools on a MetaMetro graph, or from a MetaMetro bench that already applied them.

Check that a selected colouring exists on the current ToCUMG (`manifest.yaml` / `metametro benchbuild --list-colourings`). If it is missing, stop.

Use MetaMetro to pick layers:

```bash
metamalevich --metametro-bench path/to/bench --colouring kraken2 --colouring decaying
```

Strong100 and the samovar10 bundles are prebuilt. They are not generators to copy into MetaMetro.

## Citations

Add a `.bib` entry in `cite/` only for tools this package actually integrates. Do not invent papers.
