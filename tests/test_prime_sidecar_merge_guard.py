"""The two sidecar MERGES must never produce a silently ragged file.

``pas_support.tsv`` is built by concatenation in two places -- once per
direction in ``chrom_parallel.merge_chrom_results`` (over the per-contig
workers) and once per run in ``ema.main._write_run_support`` (over the two
strands).  Both write ONE header and then copy rows through unread.

peakAtail-prime appends columns to that file from three different flags
(``--pas-features``, ``--emit-inferred-cleavage``, ``--pas-score``), two of
which have to cross a spawn boundary to reach the producer.  If one producer
gets the setting and another does not, the merged file has a header of one
width over rows of another and every column after ``tier`` is mislabelled for
part of the rows -- which is exactly the defect commit ``bfb58d1`` fixed for
the compile-time header, one level up.

Nothing downstream reads the sidecar, so the failure is invisible until an
offline model is fitted on the wrong columns.  These tests require it to be
loud instead.
"""
from __future__ import annotations

import logging

import pytest

from ema.countmatrix.paswrite import SUPPORT_COLUMNS, support_columns


def _write_sidecar(path, header_cols, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header_cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


# ---------------------------------------------------------------------------
# 1. chrom_parallel: the per-direction header
# ---------------------------------------------------------------------------
def test_the_merged_header_is_the_childrens(tmp_path):
    from ema.countmatrix.chrom_parallel import _merged_support_header

    wide = list(SUPPORT_COLUMNS) + ["clip_positions", "clip_span",
                                    "clip_offset_mean"]
    sp = tmp_path / "c1.support.tsv"
    _write_sidecar(sp, wide, [[1, 0, 0, 0, 0, 5, 1, 2, 10, "0.00"]])
    got = _merged_support_header(
        [{"direction": False, "support": str(sp)}], False)
    assert got == "\t".join(wide) + "\n"


def test_a_child_that_missed_the_flag_is_reported_not_swallowed(tmp_path, caplog):
    """The header of worker A over the rows of worker B is the whole hazard."""
    from ema.countmatrix.chrom_parallel import _merged_support_header

    wide = list(SUPPORT_COLUMNS) + ["clip_positions", "clip_span"]
    a = tmp_path / "a.support.tsv"
    b = tmp_path / "b.support.tsv"
    _write_sidecar(a, wide, [[1, 0, 0, 0, 0, 5, 1, 2, 10]])
    _write_sidecar(b, list(SUPPORT_COLUMNS), [[2, 0, 0, 0, 0, 5, 1]])
    with caplog.at_level(logging.ERROR,
                         logger="ema.countmatrix.chrom_parallel"):
        _merged_support_header(
            [{"direction": False, "support": str(a)},
             {"direction": False, "support": str(b)}], False)
    assert any("sidecar header mismatch" in r.message for r in caplog.records), (
        "a worker that did not receive --pas-features merged silently"
    )


@pytest.mark.parametrize("features,expect_extra", [("off", 0), ("on", 3)])
def test_the_no_child_fallback_uses_this_runs_column_set(tmp_path, features,
                                                         expect_extra):
    """A direction whose contigs all came back empty must still get THIS
    run's header, or the run-root merge concatenates two widths."""
    from ema.config import variable_config
    from ema.countmatrix.chrom_parallel import _merged_support_header

    saved = getattr(variable_config, "pas_features", "off")
    variable_config.pas_features = features
    try:
        got = _merged_support_header([], False).rstrip("\n").split("\t")
    finally:
        variable_config.pas_features = saved
    assert got == list(support_columns(features == "on"))
    assert len(got) == len(SUPPORT_COLUMNS) + expect_extra


# ---------------------------------------------------------------------------
# 2. ema.main: the run-root two-strand merge
# ---------------------------------------------------------------------------
def _bed_for(sidecar_path):
    """`support_path_for` maps x.bed -> x.support.tsv; invert it."""
    return str(sidecar_path)[: -len(".support.tsv")] + ".bed"


def test_the_run_root_merge_keeps_a_consistent_header(tmp_path):
    from ema.main import _write_run_support

    wide = list(SUPPORT_COLUMNS) + ["clip_positions", "clip_span"]
    pos = tmp_path / "pos.support.tsv"
    neg = tmp_path / "neg.support.tsv"
    _write_sidecar(pos, wide, [[1, 0, 0, 0, 0, 5, 1, 2, 10]])
    _write_sidecar(neg, wide, [[2, 1, 1, 1, 1, 6, 1, 1, 0]])
    out = tmp_path / "pas_support.tsv"
    _write_run_support([_bed_for(pos), _bed_for(neg)], out)
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == wide
    assert len(lines) == 3
    for line in lines[1:]:
        assert len(line.split("\t")) == len(wide)


def test_a_strand_with_a_different_header_is_reported(tmp_path, caplog):
    from ema.main import _write_run_support

    wide = list(SUPPORT_COLUMNS) + ["clip_positions", "clip_span"]
    pos = tmp_path / "pos.support.tsv"
    neg = tmp_path / "neg.support.tsv"
    _write_sidecar(pos, wide, [[1, 0, 0, 0, 0, 5, 1, 2, 10]])
    _write_sidecar(neg, list(SUPPORT_COLUMNS), [[2, 1, 1, 1, 1, 6, 1]])
    out = tmp_path / "pas_support.tsv"
    with caplog.at_level(logging.ERROR, logger="ema.main"):
        _write_run_support([_bed_for(pos), _bed_for(neg)], out)
    assert any("RAGGED" in r.message for r in caplog.records), (
        "the two strands disagreed on the column set and nothing said so"
    )


def test_the_merge_is_still_v2_when_nobody_added_a_column(tmp_path):
    """The guard must not change the v2 bytes: seven columns in, seven out."""
    from ema.main import _write_run_support

    pos = tmp_path / "pos.support.tsv"
    neg = tmp_path / "neg.support.tsv"
    _write_sidecar(pos, list(SUPPORT_COLUMNS), [[1, 0, 0, 0, 0, 5, 2]])
    _write_sidecar(neg, list(SUPPORT_COLUMNS), [[2, 3, 3, 3, 3, 9, 1]])
    out = tmp_path / "pas_support.tsv"
    _write_run_support([_bed_for(pos), _bed_for(neg)], out)
    assert out.read_text() == (
        "\t".join(SUPPORT_COLUMNS) + "\n"
        "1\t0\t0\t0\t0\t5\t2\n"
        "2\t3\t3\t3\t3\t9\t1\n"
    )
