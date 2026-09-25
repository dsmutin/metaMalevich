"""Command-line entry for metamalevich."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

from metamalevich import __version__
from metamalevich.baseline import run_pipeline


def main(argv: list[str] | None = None) -> int:
    """Parse CLI arguments and run the toy pipeline or a benchmark."""
    parser = argparse.ArgumentParser(
        prog="metamalevich",
        description="Improve metagenomic taxonomic re-profiling with coloured assembly graphs",
    )
    parser.add_argument("--version", action="store_true", help="print version and exit")
    parser.add_argument("-o", "--output", default="-", help="output path or - for stdout")
    parser.add_argument(
        "command",
        nargs="?",
        choices=["bench", "pipeline", "solve"],
        help="bench scores samovar10; solve builds MetaMetro benches and resolves their colours",
    )
    parser.add_argument("--data-root", default=".", help="directory that contains samovar10")
    parser.add_argument("--dataset", action="append", dest="datasets", help="dataset name; repeat to run several")
    parser.add_argument("--benchmark", default="benchmark", help="directory for hypothesis outputs")
    parser.add_argument(
        "--intermediate",
        default="intermediate",
        help="directory for reusable Kraken counts, k-mer edges, and ToCUMG files",
    )
    parser.add_argument(
        "--metametro-bench",
        type=Path,
        default=None,
        help="MetaMetro benchbuild directory whose ToCUMG this run should read",
    )
    parser.add_argument(
        "--colouring",
        action="append",
        dest="colourings",
        help="ToCUMG colouring or namespace to keep; repeat. Default: all layers on the bench",
    )
    parser.add_argument(
        "--bench",
        action="append",
        dest="benches",
        help="MetaMetro bench name for solve; repeat. Legacy names such as phage_10 are accepted",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Assemble an external MetaMetro bench. Phage benches also assemble when samovar, MEGAHIT, and the genome FASTA files are already present",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=None,
        help="Directory for solve benchbuild output. Default: a temporary directory",
    )
    args = parser.parse_args(argv)
    if args.version:
        print(__version__)
        return 0
    if args.command == "solve":
        return _solve(args)
    if args.metametro_bench is not None:
        from metamalevich.tocumg import load_selected

        graph = load_selected(args.metametro_bench, args.colourings)
        namespaces = sorted({row["namespace"] for row in graph.colors or []})
        summary = {
            "status": "ok",
            "ok": True,
            "input_path": str(args.metametro_bench),
            "n_unitigs": len(graph.unitigs),
            "n_links": len(graph.links),
            "namespaces": namespaces,
            "colourings": list(args.colourings or []),
        }
        _emit(json.dumps(summary, indent=2), args.output)
        return 0
    if args.command in {"bench", "pipeline"}:
        from metamalevich.bench import DATASETS, run_benchmark

        names = args.datasets or list(DATASETS)
        summary = run_benchmark(
            Path(args.data_root),
            names,
            Path(args.benchmark),
            intermediate=Path(args.intermediate),
        )
        text = json.dumps({"status": "ok", "ok": True, "input_path": args.data_root, **summary}, indent=2, default=str)
        _emit(text, args.output)
        return 0
    result = run_pipeline(None)
    _emit(json.dumps(result, indent=2), args.output)
    return 0


def _solve(args: argparse.Namespace) -> int:
    """Build the requested MetaMetro benches and print the solve summary."""
    from metametro.errors import ContractError

    from metamalevich.solve import solve_benches

    if not args.benches:
        print("solve requires at least one --bench name", file=sys.stderr)
        return 2
    try:
        if args.workdir is None:
            with tempfile.TemporaryDirectory(prefix="metamalevich-solve-") as tmp:
                summary = solve_benches(
                    args.benches, Path(tmp), args.colourings, execute=args.execute
                )
        else:
            summary = solve_benches(
                args.benches, args.workdir, args.colourings, execute=args.execute
            )
    except (ValueError, ContractError) as exc:
        print(exc, file=sys.stderr)
        return 1
    _emit(json.dumps(summary, indent=2), args.output)
    return 0 if summary["ok"] else 1


def _emit(text: str, output: str) -> None:
    if output in {"", "-"}:
        print(text)
    else:
        with open(output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
