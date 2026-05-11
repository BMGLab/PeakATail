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
