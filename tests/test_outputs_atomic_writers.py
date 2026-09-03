"""Every text artifact ema.outputs writes is temp-file-then-rename.

tests/test_pas_gene_artifacts_atomic_write.py covers pas_gene.tsv /
annotatedpas.bed (the two files the reported corruption was observed on).
The SAME non-atomic pattern — build the final path with a plain
``open(..., "w")`` / ``write_text()`` — was used by every other text writer
in the module, so a reader (or a second concurrent writer) could observe a
truncated stage stats JSON, filtered_cb.tsv, run_manifest.json or matrix
index sidecar in exactly the same way. These tests pin the fix by spying on
``os.replace``: the destination must still hold its COMPLETE previous
content at the moment of the rename, and the temp file must live in the
destination's own directory (same filesystem => the rename is one atomic
syscall, not a copy).
"""
from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from ema.config import set_directory_config
from ema.outputs import (
    OutputManager,
    atomic_write,
    write_annotated_matrix,
    write_filtered_cb,
)


class _ReplaceSpy:
    """Records ``(tmp_src, dst, dst_content_just_before_the_rename)``."""

    def __init__(self) -> None:
        self.calls: list[tuple[Path, Path, str | None]] = []
        self._real = os.replace

    def __call__(self, src, dst):
        d = Path(dst)
        self.calls.append((Path(src), d, d.read_text() if d.exists() else None))
        self._real(src, dst)

    def for_dst(self, dst: Path):
        return [c for c in self.calls if c[1] == Path(dst)]


def _assert_atomic(spy: _ReplaceSpy, dst: Path, old_content: str) -> None:
    calls = spy.for_dst(dst)
    assert len(calls) == 1, f"{dst} was not written via a single atomic rename: {spy.calls}"
    src, _, before = calls[0]
    assert src.parent == dst.parent, (
        f"temp file {src} is not in {dst}'s directory — os.replace() across "
        "filesystems is not atomic"
    )
    assert before == old_content, (
        f"{dst} was mutated before the rename — a reader could see a partial file"
    )


def test_save_stats_is_atomic(tmp_path):
    om = OutputManager(str(tmp_path / "run"))
    om.setup()
    dst = Path(om.path("cb_filter", "cb_filter_stats.json"))
    dst.write_text('{"old": true}')

    spy = _ReplaceSpy()
    with patch("os.replace", new=spy):
        om.save_stats("cb_filter", {"n_cells": 42})

    _assert_atomic(spy, dst, '{"old": true}')
    assert '"n_cells": 42' in dst.read_text()


def test_save_run_config_and_manifest_are_atomic(tmp_path):
    om = OutputManager(str(tmp_path / "run"))
    om.setup()
    cfg = Path(om.base_dir) / "run_config.json"
    cfg.write_text('{"old": true}')

    spy = _ReplaceSpy()
    with patch("os.replace", new=spy):
        om.save_run_config({"min_read": 1500})
        manifest_path = Path(om.write_manifest({"min_read": 1500}))

    _assert_atomic(spy, cfg, '{"old": true}')
    _assert_atomic(spy, manifest_path, None)  # first write: dst absent until the rename
    assert '"min_read": 1500' in cfg.read_text()


def test_write_filtered_cb_is_atomic(tmp_path):
    set_directory_config(output_dir=tmp_path / "run")
    from ema.config import directory_config

    dst = directory_config.filtered_cb_for("dsA")
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("barcode\tmin_read=1\nOLD\n")

    spy = _ReplaceSpy()
    with patch("os.replace", new=spy):
        write_filtered_cb(tmp_path / "run", "dsA", ["AAA", "CCC"], 10)

    _assert_atomic(spy, dst, "barcode\tmin_read=1\nOLD\n")
    assert dst.read_text() == "barcode\tmin_read=10\nAAA\nCCC\n"


def test_annotated_matrix_index_sidecars_are_atomic(tmp_path):
    sparse = pytest.importorskip("scipy.sparse")
    set_directory_config(output_dir=tmp_path / "run")
    from ema.config import directory_config

    ids = directory_config.annotated_pas_ids_for("dsA")
    cells = directory_config.annotated_cells_for("dsA")
    ids.parent.mkdir(parents=True, exist_ok=True)
    ids.write_text("pas_id\nOLD\n")
    cells.write_text("barcode\nOLD\n")

    mat = sparse.csr_matrix([[1, 0], [0, 2]])
    spy = _ReplaceSpy()
    with patch("os.replace", new=spy):
        write_annotated_matrix(tmp_path / "run", "dsA", mat, ["p1", "p2"], ["cbA", "cbB"])

    _assert_atomic(spy, ids, "pas_id\nOLD\n")
    _assert_atomic(spy, cells, "barcode\nOLD\n")
    assert ids.read_text() == "pas_id\np1\np2\n"
    assert cells.read_text() == "barcode\ncbA\ncbB\n"


def test_atomic_write_failure_leaves_the_old_file_and_no_temp_behind(tmp_path):
    dst = tmp_path / "artifact.tsv"
    dst.write_text("COMPLETE OLD CONTENT\n")

    with pytest.raises(RuntimeError, match="boom"):
        with atomic_write(dst) as f:
            f.write("half a fi")
            raise RuntimeError("boom")

    assert dst.read_text() == "COMPLETE OLD CONTENT\n"
    assert list(tmp_path.glob(f".{dst.name}.tmp*")) == []
