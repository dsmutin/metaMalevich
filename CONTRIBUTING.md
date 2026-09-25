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
conda env create -f environment.yml
conda activate metamalevich
```

## GitHub

Never `git push` unless the human explicitly asks. CI runs on GitHub after they push.

## Benchmarks

Do not add a new benchmark or pinned community in this repository. Add it in MetaMetro (`metametro benchbuild`) and open a pull request there. `examples/*/run_example.py` only forwards to that generator.

Do not hard-code a machine path. Use the MetaMetro bench directory, a CLI argument, or `METAMETRO_SRC`. If it is missing, stop.

Do not mock a benchmark graph, a taxonomy label, or a metric. A tiny graph that exists only inside a test file stays in that test. Example assembly graphs are MEGAHIT FASTG files, not substitutes.

Do not copy an evaluation target into graph features, colours, or any file a model reads as input. Ground truth stays in `ground_truth/` and is used only by the scorer. A table written after inference may store the id beside the prediction. That table is not an input.

Strong100 and the samovar10 bundles are prebuilt. They are not generators to copy into MetaMetro.

## Citations

Add a `.bib` entry in `cite/` only for tools this package actually integrates. Do not invent papers.
