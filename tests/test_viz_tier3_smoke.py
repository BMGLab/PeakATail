"""Smoke tests for Tier 3 viz strategies: pas_overlap, atlas_snap_diag, run_report.

All tests follow the same TDD contract as test_viz_strategies_smoke.py:
  - synthetic data only (no real BAM/BED files)
  - each strategy writes non-empty files to tmp_path
  - run_report assembles an HTML from a fake out_dir tree
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from ema.viz import get_viz_strategy, list_viz_strategies


# ---------------------------------------------------------------------------
# pas_overlap
# ---------------------------------------------------------------------------

@pytest.fixture
def pas_sets_3() -> dict[str, set[str]]:
    """Three datasets with partial overlap — suitable for Venn or UpSet."""
    return {
        "ds_A": {"PAS_1", "PAS_2", "PAS_3", "PAS_4"},
        "ds_B": {"PAS_2", "PAS_3", "PAS_5", "PAS_6"},
        "ds_C": {"PAS_3", "PAS_4", "PAS_6", "PAS_7"},
    }


@pytest.fixture
def pas_sets_5() -> dict[str, set[str]]:
    """Five datasets — forces UpSet path (>3)."""
    rng = np.random.default_rng(42)
    base = [f"PAS_{i}" for i in range(30)]
    return {
        f"ds_{c}": set(rng.choice(base, size=15, replace=False).tolist())
        for c in "ABCDE"
    }


@pytest.mark.parametrize("name", list_viz_strategies(plot_type="pas_overlap"))
def test_pas_overlap_3_datasets(name: str, pas_sets_3: dict, tmp_path: Path) -> None:
    strat = get_viz_strategy(name)
    paths = strat.render(pas_sets_3, tmp_path / f"pas_overlap_3_{name}")
    assert paths, f"{name} returned no paths"
    for p in paths:
        assert p.exists(), f"{p} not written"
        assert p.stat().st_size > 50, f"{p} is suspiciously small"


@pytest.mark.parametrize("name", list_viz_strategies(plot_type="pas_overlap"))
def test_pas_overlap_5_datasets(name: str, pas_sets_5: dict, tmp_path: Path) -> None:
    strat = get_viz_strategy(name)
    paths = strat.render(pas_sets_5, tmp_path / f"pas_overlap_5_{name}")
    assert paths, f"{name} returned no paths"
    for p in paths:
        assert p.exists(), f"{p} not written"
        assert p.stat().st_size > 50, f"{p} is suspiciously small"


@pytest.mark.parametrize("name", list_viz_strategies(plot_type="pas_overlap"))
def test_pas_overlap_2_datasets(name: str, tmp_path: Path) -> None:
    """Edge case: two datasets (minimum for overlap)."""
    data = {
        "ds_X": {"PAS_1", "PAS_2", "PAS_3"},
        "ds_Y": {"PAS_2", "PAS_3", "PAS_4"},
    }
    strat = get_viz_strategy(name)
    paths = strat.render(data, tmp_path / f"pas_overlap_2_{name}")
    assert paths, f"{name} returned no paths with 2 datasets"
    for p in paths:
        assert p.exists()


@pytest.mark.parametrize("name", list_viz_strategies(plot_type="pas_overlap"))
def test_pas_overlap_single_dataset_graceful(name: str, tmp_path: Path) -> None:
    """Single dataset: strategy must not crash — may return empty list."""
    strat = get_viz_strategy(name)
    # Should not raise
    paths = strat.render({"ds_A": {"PAS_1"}}, tmp_path / f"pas_overlap_1_{name}")
    # Returning empty list is acceptable for 1 dataset (no overlap possible)
    assert isinstance(paths, list)


# ---------------------------------------------------------------------------
# atlas_snap_diag
# ---------------------------------------------------------------------------

@pytest.fixture
def atlas_snap_stats() -> dict:
    """Typical atlas-snap diagnostic payload."""
    rng = np.random.default_rng(7)
    return {
        "snapped": 850,
        "unsnapped": 150,
        "snap_distances": rng.integers(0, 51, size=850).tolist(),
    }


@pytest.fixture
def atlas_snap_stats_empty_distances() -> dict:
    """Edge case: no snapped peaks (0 distances)."""
    return {
        "snapped": 0,
        "unsnapped": 500,
        "snap_distances": [],
    }


@pytest.mark.parametrize("name", list_viz_strategies(plot_type="atlas_snap_diag"))
def test_atlas_snap_diag_typical(name: str, atlas_snap_stats: dict, tmp_path: Path) -> None:
    strat = get_viz_strategy(name)
    paths = strat.render(atlas_snap_stats, tmp_path / f"atlas_snap_{name}")
    assert paths, f"{name} returned no paths"
    for p in paths:
        assert p.exists(), f"{p} not written"
        assert p.stat().st_size > 50


@pytest.mark.parametrize("name", list_viz_strategies(plot_type="atlas_snap_diag"))
def test_atlas_snap_diag_empty_distances(
    name: str, atlas_snap_stats_empty_distances: dict, tmp_path: Path
) -> None:
    """Zero distances must not crash."""
    strat = get_viz_strategy(name)
    paths = strat.render(atlas_snap_stats_empty_distances, tmp_path / f"atlas_snap_empty_{name}")
    assert isinstance(paths, list)
    for p in paths:
        assert p.exists()


# ---------------------------------------------------------------------------
# run_report  (generate_run_report function, not a VizStrategy subclass)
# ---------------------------------------------------------------------------

def _make_fake_out_dir(root: Path) -> Path:
    """Create a minimal fake pipeline output directory tree."""
    # per_dataset/sampleA/figures/
    a_figs = root / "per_dataset" / "sampleA" / "figures"
    a_figs.mkdir(parents=True)
    (a_figs / "umap_sampleA.html").write_text("<html>umap-A</html>")
    (a_figs / "umap_sampleA.png").write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * 100)
    (a_figs / "clusters_sampleA.html").write_text("<html>clusters-A</html>")

    b_figs = root / "per_dataset" / "sampleB" / "figures"
    b_figs.mkdir(parents=True)
    (b_figs / "umap_sampleB.html").write_text("<html>umap-B</html>")
    (b_figs / "clusters_sampleB.html").write_text("<html>clusters-B</html>")

    # top-level figures/ (Tier 3 outputs)
    top_figs = root / "figures"
    top_figs.mkdir(parents=True)
    (top_figs / "pas_overlap.html").write_text("<html>pas-overlap</html>")
    (top_figs / "atlas_snap.html").write_text("<html>atlas-snap</html>")

    # run_config.json
    import json
    (root / "run_config.json").write_text(json.dumps({"datasets": ["sampleA", "sampleB"]}))

    return root


def test_run_report_creates_html(tmp_path: Path) -> None:
    """generate_run_report writes a non-empty HTML file."""
    from ema.viz.run_report import generate_run_report

    out_dir = _make_fake_out_dir(tmp_path / "fake_run")
    report_path = out_dir / "figures" / "run_report.html"
    generate_run_report(out_dir, output_html=report_path)

    assert report_path.exists(), "run_report.html was not created"
    content = report_path.read_text()
    assert len(content) > 200, "run_report.html is suspiciously small"
    assert "<html" in content.lower(), "output does not look like HTML"


def test_run_report_links_per_dataset_figures(tmp_path: Path) -> None:
    """run_report HTML contains relative links to per-dataset UMAP files."""
    from ema.viz.run_report import generate_run_report

    out_dir = _make_fake_out_dir(tmp_path / "fake_run2")
    report_path = out_dir / "figures" / "run_report.html"
    generate_run_report(out_dir, output_html=report_path)

    content = report_path.read_text()
    # Should mention the umap files for both samples
    assert "umap_sampleA" in content, "sampleA UMAP not in report"
    assert "umap_sampleB" in content, "sampleB UMAP not in report"


def test_run_report_links_top_level_figures(tmp_path: Path) -> None:
    """run_report HTML references Tier 3 outputs (pas_overlap, atlas_snap)."""
    from ema.viz.run_report import generate_run_report

    out_dir = _make_fake_out_dir(tmp_path / "fake_run3")
    report_path = out_dir / "figures" / "run_report.html"
    generate_run_report(out_dir, output_html=report_path)

    content = report_path.read_text()
    assert "pas_overlap" in content, "pas_overlap not in report"
    assert "atlas_snap" in content, "atlas_snap not in report"


def test_run_report_files_appendix(tmp_path: Path) -> None:
    """run_report includes a Files appendix that lists every figure file."""
    from ema.viz.run_report import generate_run_report

    out_dir = _make_fake_out_dir(tmp_path / "fake_run4")
    report_path = out_dir / "figures" / "run_report.html"
    generate_run_report(out_dir, output_html=report_path)

    content = report_path.read_text()
    # appendix section must be present
    assert re.search(r"files|appendix|figure", content, re.IGNORECASE), (
        "No files appendix section found in report"
    )
    # every .html figure file should appear as a link somewhere
    for p in out_dir.rglob("*.html"):
        if p == report_path:
            continue
        assert p.name in content, f"Figure file {p.name} not listed in report"


def test_run_report_empty_out_dir(tmp_path: Path) -> None:
    """generate_run_report must not crash when out_dir has no figures at all."""
    from ema.viz.run_report import generate_run_report

    empty_run = tmp_path / "empty_run"
    empty_run.mkdir()
    report_path = empty_run / "figures" / "run_report.html"
    # Must not raise
    generate_run_report(empty_run, output_html=report_path)
    assert report_path.exists()


def test_run_report_missing_sections_graceful(tmp_path: Path) -> None:
    """If only some sections exist (no atlas snap), report still renders."""
    from ema.viz.run_report import generate_run_report

    partial_dir = tmp_path / "partial_run"
    a_figs = partial_dir / "per_dataset" / "sampleA" / "figures"
    a_figs.mkdir(parents=True)
    (a_figs / "umap_sampleA.html").write_text("<html>umap-A</html>")
    # No top-level figures/ directory, no atlas snap

    report_path = partial_dir / "figures" / "run_report.html"
    generate_run_report(partial_dir, output_html=report_path)
    assert report_path.exists()
    content = report_path.read_text()
    assert "umap_sampleA" in content
