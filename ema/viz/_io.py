"""IO helpers: dual-format save for matplotlib + plotly figures."""
from __future__ import annotations

from pathlib import Path
from typing import Any


def save_matplotlib(fig: Any, basepath: Path,
                    formats: tuple[str, ...] = ("png", "svg")) -> list[Path]:
    """Save a matplotlib Figure to PNG + SVG. Returns paths written."""
    basepath.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        p = basepath.with_suffix(f".{fmt}")
        fig.savefig(p, dpi=300, bbox_inches="tight", format=fmt)
        paths.append(p)
    try:
        import matplotlib.pyplot as plt
        plt.close(fig)
    except Exception:
        pass
    return paths


def save_plotly(fig: Any, basepath: Path,
                formats: tuple[str, ...] = ("html", "svg")) -> list[Path]:
    """Save a plotly Figure to HTML (interactive) and SVG (via kaleido)."""
    basepath.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        p = basepath.with_suffix(f".{fmt}")
        if fmt == "html":
            fig.write_html(p, include_plotlyjs="cdn")
        else:
            fig.write_image(p, format=fmt)
        paths.append(p)
    return paths
