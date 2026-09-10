"""`peakatail reannotate` refuses an ``--out`` that another branch already holds.

Companion to tests/test_pas_gene_artifacts_atomic_write.py. That file proves
a SINGLE artifact write can't be torn; this one proves the two runs never get
to share the output dir in the first place.

Root cause: 5 branches of the corrected-rerun A2/A3 trim/cluster grid shared
a `branch_name`, so `${out_root}/reannotate/${branch_name}` resolved to ONE
`--out` for all of them. Nextflow has no dependency between grid rows, so
they ran concurrently against the same run dir and produced 3 different
`pas_gene.tsv` row counts for identical declared params. Atomic writes stop
a reader seeing half a file; only this guard stops branch B from overwriting
branch A's artifacts wholesale. It must fail fast, with a clear message and
a non-zero exit — never silently interleave.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from ema.reannotate import OUT_DIR_LOCK_NAME, ReannotateError, _claim_out_dir, reannotate_run


def _fake_base_run(tmp_path: Path) -> tuple[Path, Path]:
    """A base-run dir and a GTF that merely EXIST — the out-dir guard must
    fire before any base artifact is read, so they need no real content."""
    base = tmp_path / "base"
    base.mkdir()
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("#\n")
    return base, gtf


def test_second_claim_on_the_same_out_dir_is_refused(tmp_path):
    out = tmp_path / "branch_A"
    with _claim_out_dir(out):
        assert (out / OUT_DIR_LOCK_NAME).exists()
        with pytest.raises(ReannotateError) as ei:
            with _claim_out_dir(out):
                pass
    msg = str(ei.value)
    assert "already in use" in msg
    assert str(out) in msg


def test_claim_is_released_so_a_later_run_can_reuse_the_dir(tmp_path):
    """Sequential re-runs of the same branch must still work — the guard is
    against CONCURRENT runs, not against re-running a finished branch."""
    out = tmp_path / "branch_A"
    with _claim_out_dir(out):
        pass
    with _claim_out_dir(out):  # must not raise
        pass


def test_a_crashed_branch_leaves_no_stale_lock(tmp_path):
    """flock is dropped by the kernel when the fd closes, so an exception
    inside the branch must not wedge the dir for every future run."""
    out = tmp_path / "branch_A"
    with pytest.raises(RuntimeError, match="boom"):
        with _claim_out_dir(out):
            raise RuntimeError("boom")
    with _claim_out_dir(out):  # must not raise
        pass


def test_reannotate_run_refuses_an_out_dir_held_by_another_branch(tmp_path):
    """The real entry point: a duplicate --out fails immediately, BEFORE any
    base-run artifact is touched or written."""
    base, gtf = _fake_base_run(tmp_path)
    out = tmp_path / "reannotate" / "dup_branch_name"

    with _claim_out_dir(out):
        with pytest.raises(ReannotateError) as ei:
            reannotate_run(base_run=base, out=out, gtf=gtf)
    assert "already in use" in str(ei.value)
    # Failed fast: no branch artifacts were produced under the held dir.
    assert not (out / "04_pas_gene_assignment").exists()


def test_reannotate_run_still_reports_out_equals_base_run(tmp_path):
    """The pre-existing --out == --base-run check keeps its own message (and
    must not be shadowed by the lock, nor drop a lock file into the base)."""
    base, gtf = _fake_base_run(tmp_path)
    with pytest.raises(ReannotateError, match="must differ from --base-run"):
        reannotate_run(base_run=base, out=base, gtf=gtf)
    assert not (base / OUT_DIR_LOCK_NAME).exists()


def test_cli_exits_non_zero_on_a_duplicate_out(tmp_path):
    from click.testing import CliRunner

    from ema.cli.reannotate import reannotate

    base, gtf = _fake_base_run(tmp_path)
    out = tmp_path / "reannotate" / "dup_branch_name"

    with _claim_out_dir(out):
        result = CliRunner().invoke(
            reannotate,
            ["--base-run", str(base), "--out", str(out), "--gtf", str(gtf), "--no-log-file"],
        )
    assert result.exit_code != 0, result.output
    assert "already in use" in result.output
