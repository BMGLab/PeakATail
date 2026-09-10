"""B4 regression: switch run_diff resolves pasbed via --pasbed / run layout.

Bug B4: run_diff hardcoded ``Path(h5ad_path).parent / "pasbed.bed"`` and ignored
its own ``pasbed`` argument, so under the current run layout (pasbed a sibling
dir away from the clustering h5ad) the chrom/start/end/strand columns came back
silently blank.
"""
from __future__ import annotations

from pathlib import Path

from ema.switch_test.runner import _resolve_pasbed


def test_explicit_pasbed_wins(tmp_path: Path) -> None:
    h5ad = tmp_path / "07_clustering" / "dsA" / "clusters.h5ad"
    h5ad.parent.mkdir(parents=True)
    explicit = tmp_path / "my_pasbed.bed"
    explicit.write_text("chr1\t1\t2\t1\t0\t+\n")
    assert _resolve_pasbed(str(h5ad), str(explicit)) == explicit


def test_walk_up_finds_run_layout_pasbed(tmp_path: Path) -> None:
    # pasbed lives a couple dirs up from the clustering h5ad, NOT as a sibling.
    run = tmp_path / "run_TS"
    h5ad = run / "07_clustering" / "dsA" / "clusters.h5ad"
    h5ad.parent.mkdir(parents=True)
    pasbed = run / "07_clustering" / "pasbed.bed"
    pasbed.write_text("chr1\t1\t2\t1\t0\t+\n")
    got = _resolve_pasbed(str(h5ad), None)
    assert got == pasbed.resolve()


def test_returns_none_when_absent(tmp_path: Path) -> None:
    h5ad = tmp_path / "a" / "b" / "clusters.h5ad"
    h5ad.parent.mkdir(parents=True)
    assert _resolve_pasbed(str(h5ad), None) is None


def test_nonexistent_explicit_falls_back_to_walkup(tmp_path: Path) -> None:
    run = tmp_path / "run"
    h5ad = run / "07_clustering" / "dsA" / "clusters.h5ad"
    h5ad.parent.mkdir(parents=True)
    pasbed = run / "pasbed.bed"
    pasbed.write_text("chr1\t1\t2\t1\t0\t+\n")
    # explicit path that does NOT exist → fall back to walk-up
    got = _resolve_pasbed(str(h5ad), str(tmp_path / "missing.bed"))
    assert got == pasbed.resolve()
