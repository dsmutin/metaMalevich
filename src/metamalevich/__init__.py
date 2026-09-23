"""metamalevich: Improve metagenomic taxonomic re-profiling with coloured assembly graphs"""

from __future__ import annotations

from pathlib import Path


def _version() -> str:
    """Return the package version from VERSION, then from installed metadata."""
    candidates = (
        Path(__file__).resolve().parents[2] / "VERSION",
        Path(__file__).resolve().parent / "VERSION",
    )
    for path in candidates:
        if path.is_file():
            return path.read_text(encoding="utf-8").strip()
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version("metamalevich")
    except PackageNotFoundError as exc:
        raise RuntimeError("VERSION file and package metadata are both missing") from exc


__version__ = _version()
