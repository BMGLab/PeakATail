"""Byte-identical regression test for ``ema.outputs._concat_beds``.

Perf pass: replaced the per-line `if line.strip(): out.write(line)` loop
with a single `out.writelines(...)` call per source file (same bytes
written — each surviving line, including its own original line ending, is
passed through unmodified — just fewer Python-level I/O calls).
"""
from __future__ import annotations

from pathlib import Path

from ema.outputs import _concat_beds


def _reference_concat_beds(srcs, dst: Path) -> None:
    """Verbatim copy of the pre-change implementation."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w") as out:
        for src in srcs:
            with open(src) as f:
                for line in f:
                    if line.strip():
                        out.write(line)


def test_concat_beds_matches_reference_byte_identical(tmp_path):
    src1 = tmp_path / "a.bed"
    src2 = tmp_path / "b.bed"
    # Blank lines, whitespace-only lines, and a final line with NO trailing
    # newline -- all must be preserved/dropped exactly as before.
    src1.write_text("chr1\t1\t5\t1\t0\t+\n\n   \nchr1\t10\t15\t2\t0\t+\n")
    src2.write_text("chr2\t1\t5\t1\t0\t-\nchr2\t20\t25\t2\t0\t-")  # no trailing \n

    ref_dst = tmp_path / "ref_out.bed"
    new_dst = tmp_path / "new_out.bed"
    _reference_concat_beds([src1, src2], ref_dst)
    _concat_beds([src1, src2], new_dst)

    assert ref_dst.read_bytes() == new_dst.read_bytes()
    assert new_dst.read_text() == (
        "chr1\t1\t5\t1\t0\t+\n"
        "chr1\t10\t15\t2\t0\t+\n"
        "chr2\t1\t5\t1\t0\t-\n"
        "chr2\t20\t25\t2\t0\t-"
    )


def test_concat_beds_creates_parent_dir_and_handles_empty_srcs(tmp_path):
    dst = tmp_path / "nested" / "dir" / "out.bed"
    _concat_beds([], dst)
    assert dst.exists()
    assert dst.read_text() == ""
