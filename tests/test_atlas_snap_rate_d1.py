"""D1 regression: atlas snap-rate denominator counts BOTH strands.

Bug D1: ``_render_atlas_snap_diag`` counted the snap-rate denominator from
``*.pos.bed`` only (positive strand), then clamped ``unsnapped`` with
``max(0, all_called - snapped)``. Together these pinned the reported snap
rate near 100%. The fix globs both strands (``count_called_pas``) and drops
the clamp.
"""
from __future__ import annotations

from pathlib import Path

from ema.viz.pipeline_hooks import count_called_pas, parse_atlas_mapping


def _write_bed(path: Path, n_lines: int, strand: str) -> None:
    with open(path, "w") as f:
        for i in range(n_lines):
            start = 100 + i * 10
            f.write(f"chr1\t{start}\t{start + 5}\t{i + 1}\t0\t{strand}\n")


def test_count_called_pas_counts_both_strands(tmp_path: Path) -> None:
    pc = tmp_path / "peakcalling"
    pc.mkdir()
    _write_bed(pc / "dsA.pos.bed", 7, "+")
    _write_bed(pc / "dsA.neg.bed", 5, "-")
    _write_bed(pc / "dsB.pos.bed", 3, "+")
    _write_bed(pc / "dsB.neg.bed", 4, "-")
    # Old (buggy) behaviour would have returned only the pos-strand total (10).
    assert count_called_pas(pc) == 7 + 5 + 3 + 4


def test_count_called_pas_ignores_blank_lines(tmp_path: Path) -> None:
    pc = tmp_path / "peakcalling"
    pc.mkdir()
    (pc / "d.pos.bed").write_text("chr1\t1\t5\t1\t0\t+\n\n\nchr1\t10\t15\t2\t0\t+\n")
    (pc / "d.neg.bed").write_text("\n")
    assert count_called_pas(pc) == 2


def test_count_called_pas_missing_dir_returns_zero(tmp_path: Path) -> None:
    assert count_called_pas(tmp_path / "does_not_exist") == 0


def test_unsnapped_is_not_clamped_when_denominator_complete(tmp_path: Path) -> None:
    """With both strands counted, unsnapped = all_called - snapped is honest."""
    pc = tmp_path / "peakcalling"
    pc.mkdir()
    _write_bed(pc / "d.pos.bed", 6, "+")
    _write_bed(pc / "d.neg.bed", 6, "-")  # 12 called total

    unified = tmp_path / "unified"
    unified.mkdir()
    # 4 snapped input PAS (mapping has one row per kept input PAS).
    with open(unified / "atlas_mapping.tsv", "w") as f:
        f.write("dataset_id\told_pasnumber\tnew_pas_id\tsnap_distance_bp\n")
        for i in range(4):
            f.write(f"d\t{i + 1}\t{i + 1}\t{i * 3}\n")

    snapped, distances = parse_atlas_mapping(unified / "atlas_mapping.tsv")
    all_called = count_called_pas(pc)
    assert snapped == 4
    assert all_called == 12
    assert all_called - snapped == 8  # would have been near-0 under pos-only glob
    assert distances == [0, 3, 6, 9]
