"""Figure provenance sidecar writer.

Every figure produced by a :class:`~ema.viz.base.VizStrategy` can have a
companion ``<stem>.meta.json`` file written next to it.  The sidecar records
the parameters that produced the figure so a researcher can reproduce or
audit it months later.

Usage example::

    from ema.viz._meta import write_figure_meta
    paths = save_matplotlib(fig, basepath)
    write_figure_meta(basepath, {
        "viz_strategy": "umap_matplotlib",
        "color_key": "leiden",
        "n_cells": 4200,
    })
"""
from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Version resolution — reads ema.__version__ once and caches it.
# ---------------------------------------------------------------------------

def _get_peakatail_version() -> str:
    try:
        import importlib.metadata
        return importlib.metadata.version("peakatail")
    except Exception:
        pass
    try:
        from ema import __version__  # type: ignore[attr-defined]
        return str(__version__)
    except Exception:
        pass
    return "unknown"


_VERSION_CACHE: str | None = None


def _version() -> str:
    global _VERSION_CACHE
    if _VERSION_CACHE is None:
        _VERSION_CACHE = _get_peakatail_version()
    return _VERSION_CACHE


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def write_figure_meta(basepath: Path, meta: dict[str, Any]) -> Path:
    """Write a provenance sidecar JSON next to *basepath*.

    The sidecar path is ``basepath.with_suffix(".meta.json")``.  The function
    always returns the sidecar path; if the write fails a warning is logged and
    the caller receives the path anyway (no exception is raised).

    Args:
        basepath: The base path of the figure file (no extension, or any
            extension — the suffix is replaced with ``.meta.json``).
        meta: Caller-supplied plot-specific provenance keys.  These are merged
            into the fixed schema keys (``figure_name``, ``viz_strategy``,
            ``generated_at``, ``peakatail_version``).  Caller keys take
            precedence for the plot-specific section; the four fixed keys are
            never overwritten by the caller.

    Returns:
        Path to the sidecar file (``<basepath>.meta.json``).
    """
    sidecar = basepath.with_suffix(".meta.json")
    payload: dict[str, Any] = {
        "figure_name": basepath.stem,
        "viz_strategy": meta.get("viz_strategy", "unknown"),
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
        "peakatail_version": _version(),
    }
    # Merge caller keys — fixed keys above win for the four reserved names.
    for k, v in meta.items():
        if k not in payload:
            payload[k] = v

    try:
        sidecar.parent.mkdir(parents=True, exist_ok=True)
        sidecar.write_text(json.dumps(payload, indent=2, default=str))
    except Exception as exc:  # pragma: no cover
        log.warning("write_figure_meta: could not write %s: %s", sidecar, exc)

    return sidecar


def write_figures_index(figs_dir: Path, command: str | None = None) -> Path | None:
    """Walk ``figs_dir``, collect every ``*.meta.json`` and write a manifest.

    Produces two files in ``figs_dir``:

      * ``figures_INDEX.json`` — machine-readable list ``[{file, meta}, ...]``
      * ``figures_INDEX.md``   — human-readable summary grouped by plot type

    The manifest gives a researcher a single file to read to understand what
    every figure in the directory is for, without opening 66 sidecars.

    Args:
        figs_dir: Directory containing rendered figures + their ``.meta.json``
            sidecars.  Returns ``None`` if the directory doesn't exist or is
            empty.
        command: Optional name of the producing CLI subcommand (e.g.
            ``"ema switch diff"``) — included in the markdown header.

    Returns:
        Path to the markdown file, or ``None`` on failure / empty dir.
    """
    figs_dir = Path(figs_dir)
    if not figs_dir.exists():
        return None

    entries: list[dict[str, Any]] = []
    for meta_path in sorted(figs_dir.glob("*.meta.json")):
        try:
            meta = json.loads(meta_path.read_text())
        except Exception as exc:
            log.warning("write_figures_index: bad sidecar %s: %s", meta_path, exc)
            continue
        # Collect the sibling artefacts (png/svg/html) that share this stem.
        stem = meta_path.with_suffix("").stem  # strip .meta then .json
        if stem.endswith(".meta"):
            stem = stem[:-5]
        siblings = sorted(
            p.name for p in figs_dir.iterdir()
            if p.is_file()
            and p.name.startswith(stem + ".")
            and not p.name.endswith(".meta.json")
        )
        entries.append({"stem": stem, "files": siblings, "meta": meta})

    if not entries:
        return None

    json_path = figs_dir / "figures_INDEX.json"
    md_path = figs_dir / "figures_INDEX.md"
    try:
        json_path.write_text(json.dumps(entries, indent=2, default=str))
    except Exception as exc:  # pragma: no cover
        log.warning("write_figures_index: could not write %s: %s", json_path, exc)

    # Group by plot_type for the markdown summary.
    by_type: dict[str, list[dict[str, Any]]] = {}
    for e in entries:
        plot_type = (
            e["meta"].get("viz_strategy", "unknown")
            .replace("_matplotlib", "")
            .replace("_plotly", "")
            .replace("_scanpy", "")
        )
        by_type.setdefault(plot_type, []).append(e)

    lines: list[str] = []
    lines.append(f"# Figures index")
    if command:
        lines.append(f"\nProduced by `{command}`.")
    lines.append(f"\nDirectory: `{figs_dir}`\n")
    lines.append(f"Total figures: {len(entries)} across {len(by_type)} plot type(s).\n")

    for plot_type in sorted(by_type.keys()):
        rows = by_type[plot_type]
        lines.append(f"\n## `{plot_type}` ({len(rows)})\n")
        # Pull a representative description from the first entry if present.
        desc = rows[0]["meta"].get("description")
        if desc:
            lines.append(f"_{desc}_\n")
        # Concise table per entry.
        for e in rows:
            m = e["meta"]
            tag_parts: list[str] = []
            if "cluster1" in m and "cluster2" in m:
                tag_parts.append(f"cluster {m['cluster1']} vs {m['cluster2']}")
            if "dataset_id" in m:
                tag_parts.append(f"dataset={m['dataset_id']}")
            if "n_tested" in m:
                tag_parts.append(f"n_tested={m['n_tested']}")
            if "n_significant" in m:
                tag_parts.append(f"sig={m['n_significant']}")
            if "n_observations" in m:
                tag_parts.append(f"n_obs={m['n_observations']}")
            if "n_clusters" in m:
                tag_parts.append(f"clusters={m['n_clusters']}")
            # gene_track-specific tags (gene_id + n_pas + n_clusters_rendered)
            if "gene_id" in m:
                tag_parts.append(f"gene={m['gene_id']}")
            if "n_pas" in m:
                tag_parts.append(f"n_pas={m['n_pas']}")
            if "n_clusters_rendered" in m:
                tag_parts.append(f"clusters_rendered={m['n_clusters_rendered']}")
            if "n_isoforms" in m and m["n_isoforms"]:
                tag_parts.append(f"isoforms={m['n_isoforms']}")
            tags = " · ".join(tag_parts) if tag_parts else "—"
            lines.append(f"- **{e['stem']}** ({len(e['files'])} files) — {tags}")
        lines.append("")

    try:
        md_path.write_text("\n".join(lines))
    except Exception as exc:  # pragma: no cover
        log.warning("write_figures_index: could not write %s: %s", md_path, exc)
        return None

    return md_path
