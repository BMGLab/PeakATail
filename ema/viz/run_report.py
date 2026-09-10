"""Composite HTML run report for PeakATail pipeline outputs.

Generates a single self-contained HTML file that links to (and embeds previews
of) all figures produced during a pipeline run.  The report is assembled from
``out_dir`` by scanning for known file patterns — it gracefully skips missing
sections rather than crashing.

Public API::

    generate_run_report(out_dir: Path, output_html: Path) -> Path

The function never raises: all errors are caught and logged.  This preserves
pipeline integrity even when the report assembly fails.
"""
from __future__ import annotations

import datetime
import json
import logging
import os
from pathlib import Path
from typing import NamedTuple

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Internal data types
# ---------------------------------------------------------------------------

class _FigFile(NamedTuple):
    """A single figure file discovered during scan."""
    path: Path
    rel_path: str   # relative to out_dir (for display)
    section: str    # "tier1", "tier2", "tier3", "other"
    dataset: str    # dataset_id, or "" for cross-dataset/top-level files
    kind: str       # "umap", "cluster_sizes", "pas_overlap", …


# ---------------------------------------------------------------------------
# Pattern catalogue
# ---------------------------------------------------------------------------

# Mapping from stem-keyword to (section, kind)
_STEM_PATTERNS: list[tuple[str, str, str]] = [
    # Tier 1 (per-dataset)
    ("umap_",            "tier1", "umap"),
    ("clusters_",        "tier1", "cluster_sizes"),
    ("peak_qc",          "tier1", "peak_qc"),
    # Tier 2 (switch / diff)
    ("volcano",          "tier2", "volcano"),
    ("diff_agreement",   "tier2", "diff_agreement"),
    ("pdui",             "tier2", "pdui_distribution"),
    ("length_shift",     "tier2", "length_shifts"),
    ("sankey",           "tier2", "cluster_match_sankey"),
    ("match_conf",       "tier2", "match_confidence"),
    # Tier 3 (cross-dataset)
    ("pas_overlap",      "tier3", "pas_overlap"),
    ("atlas_snap",       "tier3", "atlas_snap_diag"),
    # Tier 4 (perf)
    ("tile_timing",      "tier4", "tile_timing"),
    ("resource_timeline","tier4", "resource_timeline"),
]

_SECTION_TITLES = {
    "tier1": "Tier 1 — Per-dataset QC & clustering",
    "tier2": "Tier 2 — Differential APA (switch outputs)",
    "tier3": "Tier 3 — Multi-sample cross-dataset",
    "tier4": "Tier 4 — Performance & resource usage",
    "other": "Other figures",
}

_PRIORITY_EXT = (".html", ".png", ".svg")


# ---------------------------------------------------------------------------
# Scanner
# ---------------------------------------------------------------------------

def _classify_file(f: Path, out_dir: Path) -> _FigFile | None:
    """Return a _FigFile for a plot file, or None if it should be ignored."""
    if f.suffix.lower() not in _PRIORITY_EXT:
        return None
    stem = f.stem.lower()
    section = "other"
    kind = stem
    for keyword, sec, knd in _STEM_PATTERNS:
        if keyword in stem:
            section = sec
            kind = knd
            break

    # Infer dataset from path.  New layout: 0X_<stage>/<ds>/figures/<file>.
    # Legacy layout: per_dataset/<ds>/figures/<file> (kept for old runs).
    import re as _re
    _STAGE_PAT = _re.compile(r"^\d{2}_")
    dataset = ""
    parts = f.parts
    # Try new layout: look for a numbered stage dir immediately above the ds dir.
    for i, part in enumerate(parts):
        if _STAGE_PAT.match(part) and i + 2 < len(parts) and parts[i + 2] == "figures":
            dataset = parts[i + 1]
            break
    # Fallback: legacy per_dataset/<ds>/... layout (for runs produced before
    # the 0X_<stage>/<ds>/ migration).
    if not dataset:
        try:
            pd_idx = parts.index("per_dataset")
            if pd_idx + 1 < len(parts):
                dataset = parts[pd_idx + 1]
        except ValueError:
            pass

    try:
        rel = str(f.relative_to(out_dir))
    except ValueError:
        rel = str(f)

    return _FigFile(path=f, rel_path=rel, section=section, dataset=dataset, kind=kind)


def _scan_figures(out_dir: Path) -> list[_FigFile]:
    """Walk out_dir and return all discovered figure files, sorted."""
    figs: list[_FigFile] = []
    for f in sorted(out_dir.rglob("*")):
        if not f.is_file():
            continue
        fig = _classify_file(f, out_dir)
        if fig is not None:
            figs.append(fig)
    return figs


# ---------------------------------------------------------------------------
# HTML generation helpers
# ---------------------------------------------------------------------------

def _read_run_config(out_dir: Path) -> dict:
    """Load run_config.json from out_dir, return {} on any error."""
    cfg_path = out_dir / "run_config.json"
    if cfg_path.exists():
        try:
            return json.loads(cfg_path.read_text())
        except Exception:
            pass
    return {}


def _rel_to_report(fig_path: Path, report_html: Path) -> str:
    """Relative path from the report HTML file to a figure file."""
    try:
        return os.path.relpath(fig_path, report_html.parent)
    except ValueError:
        return str(fig_path)


def _figure_card(fig: _FigFile, report_html: Path) -> str:
    """Return an HTML snippet for a single figure."""
    rel = _rel_to_report(fig.path, report_html)
    label = f"{fig.dataset}/{fig.kind}" if fig.dataset else fig.kind
    ext = fig.path.suffix.lower()

    if ext == ".html":
        return (
            f'<div class="fig-card">'
            f'<a href="{rel}" target="_blank">'
            f'<div class="fig-label">{label}</div>'
            f'<div class="fig-type html-link">Interactive HTML</div>'
            f'</a>'
            f'</div>\n'
        )
    elif ext == ".png":
        return (
            f'<div class="fig-card">'
            f'<a href="{rel}" target="_blank">'
            f'<img src="{rel}" alt="{label}" loading="lazy"/>'
            f'<div class="fig-label">{label}</div>'
            f'</a>'
            f'</div>\n'
        )
    else:
        return (
            f'<div class="fig-card">'
            f'<a href="{rel}" target="_blank">'
            f'<div class="fig-label">{label}</div>'
            f'<div class="fig-type">{ext.lstrip(".")}</div>'
            f'</a>'
            f'</div>\n'
        )


def _section_html(
    section_id: str,
    figs: list[_FigFile],
    report_html: Path,
) -> str:
    """Return HTML for one section block."""
    if not figs:
        return ""
    title = _SECTION_TITLES.get(section_id, section_id.title())
    cards = "".join(_figure_card(f, report_html) for f in figs)
    return (
        f'<section id="{section_id}">\n'
        f'  <h2>{title}</h2>\n'
        f'  <div class="fig-grid">\n{cards}  </div>\n'
        f'</section>\n'
    )


def _appendix_html(all_figs: list[_FigFile], report_html: Path) -> str:
    """Appendix: alphabetical listing of every figure file with a link."""
    rows = "".join(
        f'<li><a href="{_rel_to_report(f.path, report_html)}" target="_blank">'
        f'{f.rel_path}</a></li>\n'
        for f in sorted(all_figs, key=lambda x: x.rel_path)
    )
    return (
        f'<section id="appendix-files">\n'
        f'  <h2>Files appendix</h2>\n'
        f'  <ul class="file-list">\n{rows}  </ul>\n'
        f'</section>\n'
    )


# ---------------------------------------------------------------------------
# CSS + HTML template
# ---------------------------------------------------------------------------

_CSS = """
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
       margin: 0; padding: 0; background: #f8f9fa; color: #212529; }
header { background: #2d4e8c; color: white; padding: 20px 40px; }
header h1 { margin: 0; font-size: 1.6em; }
header .meta { font-size: 0.85em; margin-top: 6px; opacity: 0.85; }
main { max-width: 1400px; margin: 20px auto; padding: 0 20px; }
section { background: white; border-radius: 8px; padding: 20px 24px;
          margin-bottom: 20px; box-shadow: 0 1px 4px rgba(0,0,0,.08); }
section h2 { margin-top: 0; font-size: 1.15em; color: #2d4e8c;
             border-bottom: 2px solid #e9ecef; padding-bottom: 8px; }
.fig-grid { display: flex; flex-wrap: wrap; gap: 12px; margin-top: 12px; }
.fig-card { border: 1px solid #dee2e6; border-radius: 6px; padding: 8px;
            background: #fff; width: 200px; text-align: center;
            transition: box-shadow .15s; }
.fig-card:hover { box-shadow: 0 3px 10px rgba(0,0,0,.15); }
.fig-card a { text-decoration: none; color: inherit; display: block; }
.fig-card img { max-width: 100%; height: 120px; object-fit: contain;
                display: block; margin: 0 auto 6px; }
.fig-label { font-size: 0.78em; font-weight: 600; color: #495057;
             word-break: break-word; }
.fig-type { font-size: 0.70em; color: #6c757d; margin-top: 2px; }
.html-link { color: #4c72b0; font-weight: 600; font-size: 0.85em;
             padding: 4px; }
.file-list { font-size: 0.82em; column-count: 2; column-gap: 24px; padding-left: 18px; }
.file-list li { margin-bottom: 4px; }
.config-box { background: #f1f3f5; border-radius: 6px; padding: 12px;
              font-size: 0.82em; font-family: monospace; white-space: pre-wrap;
              max-height: 200px; overflow-y: auto; }
"""


def _build_html(
    out_dir: Path,
    report_html: Path,
    all_figs: list[_FigFile],
    cfg: dict,
    timestamp: str,
) -> str:
    """Assemble full HTML document string."""
    n_datasets = 0
    if "datasets" in cfg:
        n_datasets = len(cfg["datasets"])
    else:
        # Infer from 07_clustering/ subdirectories (new layout).
        # Fallback: legacy per_dataset/ for runs produced before migration.
        clustering_d = out_dir / "07_clustering"
        per_dataset_d = out_dir / "per_dataset"
        if clustering_d.is_dir():
            n_datasets = sum(1 for x in clustering_d.iterdir() if x.is_dir())
        elif per_dataset_d.is_dir():
            n_datasets = sum(1 for x in per_dataset_d.iterdir() if x.is_dir())

    # Group by section
    by_section: dict[str, list[_FigFile]] = {}
    for f in all_figs:
        # Exclude the report itself
        if f.path == report_html:
            continue
        by_section.setdefault(f.section, []).append(f)

    section_order = ["tier1", "tier2", "tier3", "tier4", "other"]
    sections_html = "".join(
        _section_html(sec, by_section.get(sec, []), report_html)
        for sec in section_order
    )

    # Filter non-report figs for appendix
    appendix_figs = [f for f in all_figs if f.path != report_html]

    cfg_json = json.dumps(cfg, indent=2) if cfg else "(no run_config.json found)"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>PeakATail run report</title>
  <style>{_CSS}</style>
</head>
<body>
<header>
  <h1>PeakATail Run Report</h1>
  <div class="meta">
    Generated: {timestamp} &nbsp;|&nbsp;
    Output dir: {out_dir} &nbsp;|&nbsp;
    Datasets: {n_datasets}
  </div>
</header>
<main>
<section id="config-summary">
  <h2>Run configuration</h2>
  <div class="config-box">{cfg_json}</div>
</section>
{sections_html}
{_appendix_html(appendix_figs, report_html)}
</main>
</body>
</html>
"""
    return html


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_run_report(
    out_dir: Path,
    output_html: Path | None = None,
) -> Path:
    """Scan ``out_dir`` and write a composite HTML run report.

    The function is designed to be called at the end of the pipeline.  It
    never raises: all errors are caught and logged as warnings so that report
    failure does not crash the analysis.

    Args:
        out_dir: Root pipeline output directory (contains ``per_dataset/``,
            ``figures/``, ``run_config.json``, …).
        output_html: Destination for the HTML report.  Defaults to
            ``out_dir/figures/run_report.html``.

    Returns:
        Path to the written HTML file (may not exist if writing failed).
    """
    out_dir = Path(out_dir)
    if output_html is None:
        output_html = out_dir / "figures" / "run_report.html"
    output_html = Path(output_html)

    try:
        output_html.parent.mkdir(parents=True, exist_ok=True)

        all_figs = _scan_figures(out_dir)
        cfg = _read_run_config(out_dir)
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        html = _build_html(out_dir, output_html, all_figs, cfg, timestamp)
        output_html.write_text(html, encoding="utf-8")

        log.info(
            "Run report written: %s (%d figure files indexed)",
            output_html, len(all_figs),
        )
    except Exception as exc:
        log.warning("generate_run_report failed: %s", exc)

    return output_html
