"""B3 regression: filtered_cb.tsv is always written (even when empty).

Bug B3: the per-dataset write was gated on ``if filtered_cb_path.exists()``,
which silently skipped when the run-root file was missing — leaving
``02_cb_filter/<ds>/`` with only the stats JSON, so the min_read survivor set
was unrecoverable. The write helper must always emit a file with the header.
"""
from __future__ import annotations

from pathlib import Path

from ema.config import set_directory_config
from ema.outputs import write_filtered_cb


def test_write_filtered_cb_creates_file_with_header(tmp_path: Path) -> None:
    set_directory_config(output_dir=tmp_path / "run")
    dst = write_filtered_cb(tmp_path / "run", "dsA", ["AAAA", "CCCC", "GGGG"], 1500)
    assert dst.exists()
    lines = dst.read_text().splitlines()
    assert lines[0] == "barcode\tmin_read=1500"
    assert lines[1:] == ["AAAA", "CCCC", "GGGG"]


def test_write_filtered_cb_empty_list_still_writes_header(tmp_path: Path) -> None:
    """Even with zero surviving barcodes, a header-only file must exist."""
    set_directory_config(output_dir=tmp_path / "run")
    dst = write_filtered_cb(tmp_path / "run", "dsB", [], 900)
    assert dst.exists(), "filtered_cb.tsv must be written even when empty (B3)"
    lines = dst.read_text().splitlines()
    assert lines == ["barcode\tmin_read=900"]


def test_write_filtered_cb_lands_in_cb_filter_dir(tmp_path: Path) -> None:
    set_directory_config(output_dir=tmp_path / "run")
    dst = write_filtered_cb(tmp_path / "run", "dsC", ["TTTT"], 100)
    # Canonical location: <run>/02_cb_filter/<ds>/filtered_cb.tsv
    assert dst.parent.name == "dsC"
    assert dst.parent.parent.name == "02_cb_filter"
