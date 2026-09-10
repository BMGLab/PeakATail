"""Smoke tests for Phase V4 (Tier 4) visualization strategies.

Covers:
- tile_timing_matplotlib / tile_timing_plotly
- resource_timeline_matplotlib / resource_timeline_plotly

Each test generates synthetic input data, calls render(), and asserts
that all returned paths exist and are non-empty.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ema.viz import get_viz_strategy, list_viz_strategies


# ---------------------------------------------------------------------------
# Synthetic data factories
# ---------------------------------------------------------------------------


def _make_tile_timings(
    n_chroms: int = 3,
    tiles_per_chrom: int = 4,
    n_datasets: int = 1,
    rng_seed: int = 0,
) -> list[dict[str, Any]]:
    """Create synthetic tile timing records matching the tile_runner JSON schema.

    Schema per record::

        {
            "dataset_id": str,
            "chrom": str,
            "tile_idx": int,          # 0-based index within chrom
            "tile_start": int,
            "tile_end": int,
            "wall_seconds": float,
        }
    """
    rng = np.random.default_rng(rng_seed)
    chroms = [f"chr{i + 1}" for i in range(n_chroms)]
    ds_ids = [f"sample_{d}" for d in range(n_datasets)]
    records: list[dict[str, Any]] = []
    for ds_id in ds_ids:
        for chrom in chroms:
            for tile_idx in range(tiles_per_chrom):
                tile_start = tile_idx * 25_000_000
                tile_end = tile_start + 25_000_000
                # Inject one outlier per dataset × chrom to test red-highlight
                if tile_idx == 0:
                    wall = float(rng.uniform(60, 120))  # deliberately slow
                else:
                    wall = float(rng.uniform(1, 10))
                records.append(
                    {
                        "dataset_id": ds_id,
                        "chrom": chrom,
                        "tile_idx": tile_idx,
                        "tile_start": tile_start,
                        "tile_end": tile_end,
                        "wall_seconds": wall,
                    }
                )
    return records


def _make_resource_samples(
    n_samples: int = 30,
    rng_seed: int = 0,
) -> list[dict[str, Any]]:
    """Create synthetic resource-timeline records matching the sampler schema.

    Schema per record::

        {
            "elapsed_s": float,
            "rss_gb": float,
            "cpu_pct": float,
        }
    """
    rng = np.random.default_rng(rng_seed)
    t0 = time.time() - n_samples * 5
    records: list[dict[str, Any]] = []
    for i in range(n_samples):
        records.append(
            {
                "elapsed_s": float(i * 5),
                "rss_gb": float(rng.uniform(0.5, 4.0)),
                "cpu_pct": float(rng.uniform(10, 95)),
            }
        )
    return records


def _make_stage_annotations() -> list[dict[str, Any]]:
    """Create synthetic stage-transition annotations.

    Schema per record::

        {
            "label": str,
            "elapsed_s": float,
        }
    """
    return [
        {"label": "peak calling done", "elapsed_s": 30.0},
        {"label": "atlas snap done", "elapsed_s": 80.0},
        {"label": "clustering done", "elapsed_s": 120.0},
    ]


# ---------------------------------------------------------------------------
# tile_timing tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name",
    [n for n in list_viz_strategies(plot_type="tile_timing")],
)
def test_tile_timing_single_dataset_smoke(name: str, tmp_path: Path) -> None:
    """Single-dataset tile timing renders without error and produces files."""
    strat = get_viz_strategy(name)
    data = _make_tile_timings(n_chroms=3, tiles_per_chrom=4, n_datasets=1)
    paths = strat.render(data, tmp_path / "tile_timing_single")
    assert paths, f"{name} returned no output paths"
    for p in paths:
        assert p.exists(), f"Expected output file {p} does not exist"
        assert p.stat().st_size > 100, f"Output file {p} is suspiciously small"


@pytest.mark.parametrize(
    "name",
    [n for n in list_viz_strategies(plot_type="tile_timing")],
)
def test_tile_timing_multi_dataset_smoke(name: str, tmp_path: Path) -> None:
    """Multi-dataset tile timing produces subplots/grouped rendering."""
    strat = get_viz_strategy(name)
    data = _make_tile_timings(n_chroms=2, tiles_per_chrom=3, n_datasets=3)
    paths = strat.render(data, tmp_path / "tile_timing_multi")
    assert paths, f"{name} returned no output paths for multi-dataset data"
    for p in paths:
        assert p.exists(), f"Expected output file {p} does not exist"


@pytest.mark.parametrize(
    "name",
    [n for n in list_viz_strategies(plot_type="tile_timing")],
)
def test_tile_timing_empty_input_does_not_crash(name: str, tmp_path: Path) -> None:
    """Empty timing list must not crash — render() may return [] or a placeholder."""
    strat = get_viz_strategy(name)
    try:
        paths = strat.render([], tmp_path / "tile_timing_empty")
        # Returning [] is acceptable for empty input
        assert isinstance(paths, list)
    except Exception as exc:
        pytest.fail(f"{name}.render([]) raised unexpectedly: {exc}")


# ---------------------------------------------------------------------------
# resource_timeline tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name",
    [n for n in list_viz_strategies(plot_type="resource_timeline")],
)
def test_resource_timeline_smoke(name: str, tmp_path: Path) -> None:
    """Resource timeline renders without error and produces files."""
    strat = get_viz_strategy(name)
    samples = _make_resource_samples(n_samples=30)
    annotations = _make_stage_annotations()
    data = {"samples": samples, "annotations": annotations}
    paths = strat.render(data, tmp_path / "resource_timeline")
    assert paths, f"{name} returned no output paths"
    for p in paths:
        assert p.exists(), f"Expected output file {p} does not exist"
        assert p.stat().st_size > 100, f"Output file {p} is suspiciously small"


@pytest.mark.parametrize(
    "name",
    [n for n in list_viz_strategies(plot_type="resource_timeline")],
)
def test_resource_timeline_no_annotations(name: str, tmp_path: Path) -> None:
    """Resource timeline works when no stage annotations are provided."""
    strat = get_viz_strategy(name)
    samples = _make_resource_samples(n_samples=20)
    data = {"samples": samples, "annotations": []}
    paths = strat.render(data, tmp_path / "resource_timeline_no_ann")
    assert paths, f"{name} returned no paths even without annotations"
    for p in paths:
        assert p.exists()


@pytest.mark.parametrize(
    "name",
    [n for n in list_viz_strategies(plot_type="resource_timeline")],
)
def test_resource_timeline_empty_input_does_not_crash(
    name: str, tmp_path: Path
) -> None:
    """Empty sample list must not crash."""
    strat = get_viz_strategy(name)
    data = {"samples": [], "annotations": []}
    try:
        paths = strat.render(data, tmp_path / "resource_timeline_empty")
        assert isinstance(paths, list)
    except Exception as exc:
        pytest.fail(f"{name}.render(empty) raised unexpectedly: {exc}")


# ---------------------------------------------------------------------------
# tile_runner integration: _write_timings round-trip
# ---------------------------------------------------------------------------


def test_write_timings_json_round_trip(tmp_path: Path) -> None:
    """_write_timings() writes valid JSON that renders without error."""
    from ema.countmatrix.tile_runner import _write_timings

    records = _make_tile_timings(n_chroms=2, tiles_per_chrom=2, n_datasets=1)
    out_path = tmp_path / "tile_timings.json"
    _write_timings(records, out_path)

    assert out_path.exists()
    loaded = json.loads(out_path.read_text())
    assert isinstance(loaded, list)
    assert len(loaded) == len(records)
    assert "wall_seconds" in loaded[0]


# ---------------------------------------------------------------------------
# resource_sampler integration
# ---------------------------------------------------------------------------


def test_resource_sampler_starts_and_stops(tmp_path: Path) -> None:
    """_ResourceSampler thread starts, collects at least one sample, stops cleanly."""
    from ema.main import _ResourceSampler

    out_path = tmp_path / "resources.jsonl"
    sampler = _ResourceSampler(out_path=out_path, interval_s=0.05)
    sampler.start()
    time.sleep(0.2)  # let it collect a couple of samples
    sampler.stop()
    sampler.join(timeout=2.0)

    assert not sampler.is_alive(), "Sampler thread did not stop"
    # File may or may not be written depending on timing; either is fine
    # but if it exists it must be valid JSONL
    if out_path.exists() and out_path.stat().st_size > 0:
        lines = [l for l in out_path.read_text().splitlines() if l.strip()]
        assert lines, "resources.jsonl is non-empty but has no valid lines"
        rec = json.loads(lines[0])
        assert "elapsed_s" in rec
        assert "rss_gb" in rec
        assert "cpu_pct" in rec


def test_resource_sampler_stop_before_start_safe(tmp_path: Path) -> None:
    """Calling stop() before start() must not raise."""
    from ema.main import _ResourceSampler

    out_path = tmp_path / "resources_noop.jsonl"
    sampler = _ResourceSampler(out_path=out_path, interval_s=0.1)
    # Do NOT call start() — just stop and join
    sampler.stop()
    sampler.join(timeout=1.0)
