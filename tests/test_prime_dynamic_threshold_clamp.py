"""``--dynamic-threshold-clamp``: the crash, now bounded, and the identity rule.

MERGE NOTE (develop -> peakAtail-prime, issues #101/#102).  This module used
to pin ``off`` == v2 == ``IndexError``.  It no longer can: the merge made the
look-back bound UNCONDITIONAL and reduced ``--dynamic-threshold-clamp`` to an
accepted no-op, so neither value can reproduce the abort.  The two tests that
asserted the crash now assert the flag's no-op-ness instead -- the same
fixture, the opposite expectation -- and the identity guarantee they rested
on is unchanged and still checked exhaustively below.

The 2026-08 parameter sweep proved (``results/paramsweep/VERDICTS.md`` §5) that
``--dynamic-threshold --lambda-fold-change 2.0`` — the parameter reference's
own documented "find more PAS" setting — aborts the caller with a bare
``IndexError`` at::

    l_end = data_array[-current_threshold]

because the dynamic estimator sets ``current_threshold`` from the local read
density and nothing bounds it by ``len(data_array)``.  The sweep's accuracy
verifier reproduced the same abort on the mouse dev slice
(``results/paramsweep/VERIFY/rerun_identity.tsv``), so the defect is
species-independent.  It is a v2 (``9dfdefb``) defect, reachable only with
``--dynamic-threshold`` (off by default).

This module pins what the merged tree does:

1. the crash fixture — a 14-read synthetic BAM — COMPLETES at the flag's
   default, because the bound no longer needs asking for;
2. ``on`` completes the identical run and still emits the peak;
3. ``on`` is byte-identical to ``off`` on a dynamic-threshold run that does
   NOT trip the bound — the bound can only change a run v2 would have
   aborted;
4. :func:`ema.countmatrix.dynamic_threshold.resolve_l_end` honours the same
   identity exhaustively at the unit level;
5. both loop implementations are covered: the monolithic loop
   (``peackcalling.py``) and the 3-stage pipeline's finder
   (``peak_pipeline.py``), which the sweep found carries the same expression.

Run against code WITHOUT the fix (the frozen v2 worktree), tests 1 and 2 fail
with ``IndexError`` / ``TypeError: unexpected keyword argument
'dynamic_threshold_clamp'`` and test 4 with ``ModuleNotFoundError``.  That the
fixture really reaches the bound is proved directly by test 4's exhaustive
unit check of the same rule, and by the fixture arithmetic spelled out below.

CRASH MECHANICS OF THE FIXTURE (all forward-strand, one chromosome):
14 reads start at ``920+i``; under v2 "fixed" geometry each end is rewritten
to ``start + seq_len``, so the ends land at ``1011..1024`` — inside one
``lambda_window`` of each other, so the background deque never prunes.  With
``default_threshold=floor_threshold=3`` the signal fires at read 3.  At read
11 the deque tops 10 entries (the estimator's minimum), and with
``lambda_window=1000, lambda_fold_change=10000`` the threshold jumps to
``int(11/1000 * 10000) = 110`` while the live window holds 10 ends:
``data_array[-110]`` → ``IndexError``.  With ``lambda_fold_change=2.0`` the
same arithmetic yields ``max(3, int(0.011*2)) = 3`` and nothing ever crashes,
which is what makes fixture 3's identity comparison meaningful.
"""
from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import pysam
import pytest

CHROM = "1"
CHROM_LEN = 200_000
SEQ_LEN = 91
N_READS = 14

# Crash arm: threshold jumps to int(11/1000 * 10000) = 110 >> window size 10.
CRASH_KWARGS = dict(
    dynamic_threshold=True,
    floor_threshold=3,
    lambda_window=1000,
    lambda_fold_change=10000.0,
)
# Benign arm: same run, threshold stays at the floor; v2 completes it.
BENIGN_KWARGS = dict(
    dynamic_threshold=True,
    floor_threshold=3,
    lambda_window=1000,
    lambda_fold_change=2.0,
)

# spawn (pipeline path) re-parses argv in the children — same convention as
# tests/test_polya_three_path_agreement.py.
_ARGV = [
    "ema",
    "--sequenceLen", str(SEQ_LEN),
    "--CellBarcodeLen", "16",
    "--BarcodeTag", "CB",
]


def _barcode(i: int) -> str:
    alphabet = "ACGT"
    out = [alphabet[(i >> (2 * k)) & 3] for k in range(8)]
    return ("".join(out) + "ACGTACGT")[:16]


def _write_crash_bam(path: Path) -> Path:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": CHROM_LEN}],
        "RG": [{"ID": "testsample", "SM": "testsample"}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i in range(N_READS):
            # The last four reads carry a genuine terminal poly(A) soft clip:
            # clip_seeded's seeder then guarantees a tier-1 call whether or
            # not a coverage peak completes, so the rescued runs have a row
            # to assert on.  The clip does not touch the crash mechanics —
            # read_check's coverage 5-tuple ignores soft clips.
            clip = 12 if i >= N_READS - 4 else 0
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "C" * 80 + "A" * clip
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 920 + i
            a.mapping_quality = 60
            a.cigartuples = [(0, 80)] + ([(4, clip)] if clip else [])
            a.query_qualities = pysam.qualitystring_to_array(
                "I" * (80 + clip))
            a.set_tag("CB", _barcode(i))
            a.set_tag("UB", f"UMI{i:06d}")
            a.set_tag("RG", "testsample")
            out.write(a)
    pysam.index(str(path))
    return path


@pytest.fixture(scope="module")
def crash_bam(tmp_path_factory) -> Path:
    d = tmp_path_factory.mktemp("dyn_clamp_bam")
    return _write_crash_bam(d / "crash.bam")


@pytest.fixture(autouse=True)
def _slice_config():
    """Point ``variable_config`` at the fixture BAM's parameters; restore after."""
    from ema.config import variable_config

    saved = {
        k: getattr(variable_config, k)
        for k in ("seqlen", "cb_len", "barcode_tag", "ignore_chro")
    }
    variable_config.seqlen = SEQ_LEN
    variable_config.cb_len = 16
    variable_config.barcode_tag = "CB"
    variable_config.ignore_chro = []
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)


def _call(bam: Path, out: Path, tag: str, **kwargs):
    from ema.countmatrix.indexing import reset_index
    from ema.countmatrix.peackcalling import peak_calling
    from ema.countmatrix.peak import Peak
    from ema.strategies import get_strategy

    reset_index()
    Peak.reset_pasnumber()
    bed = out / f"{tag}.bed"
    mtx = out / f"{tag}.mtx"
    with patch.object(sys, "argv", _ARGV):
        peak_calling(
            False,
            bedfilepath=str(bed),
            matrixpath=str(mtx),
            bamfile_dir=str(bam),
            default_threshold=3,
            merge_len=100,
            strategy=get_strategy("clip_seeded"),
            **kwargs,
        )
    return bed, mtx


# ---------------------------------------------------------------------------
# 1. the crash, reproduced — and reproduced as V2 BEHAVIOUR (default off)
# ---------------------------------------------------------------------------

def test_the_default_no_longer_aborts(crash_bam, tmp_path):
    """The sweep's crash, in 14 reads — and it does not happen any more.

    Before the merge this asserted ``pytest.raises(IndexError)``, because the
    branch defaulted the clamp ``off`` to keep v2 reproducible down to its
    abort.  Issue #101 settled that a crash is not a behaviour worth
    preserving, so the bound is on with nothing asked for."""
    bed, _ = _call(crash_bam, tmp_path, "default_bounded", **CRASH_KWARGS)
    rows = [ln for ln in bed.read_text().splitlines() if ln.strip()]
    assert rows, "the bounded run must still emit the coverage peak"


def test_explicit_off_cannot_bring_the_abort_back(crash_bam, tmp_path):
    """``--dynamic-threshold-clamp off`` is an accepted NO-OP.

    It is still parsed, still threaded through every dispatch path and still
    in ``V2_COMPAT_FLAGS`` so existing command lines and configs keep working
    — but it can no longer select the unbounded index.  If this ever raises
    again, someone re-introduced a user-reachable route to issue #101."""
    bed, _ = _call(crash_bam, tmp_path, "off_is_a_noop",
                   dynamic_threshold_clamp=False, **CRASH_KWARGS)
    rows = [ln for ln in bed.read_text().splitlines() if ln.strip()]
    assert rows, "clamp=off must behave exactly like clamp=on now"


def test_off_and_on_agree_on_the_run_that_used_to_abort(crash_bam, tmp_path):
    """The no-op claim, stated as bytes rather than as an absence of raise."""
    bed_off, mtx_off = _call(crash_bam, tmp_path, "noop_off",
                             dynamic_threshold_clamp=False, **CRASH_KWARGS)
    bed_on, mtx_on = _call(crash_bam, tmp_path, "noop_on",
                           dynamic_threshold_clamp=True, **CRASH_KWARGS)
    assert bed_off.read_bytes() == bed_on.read_bytes()
    assert mtx_off.read_bytes() == mtx_on.read_bytes()


# ---------------------------------------------------------------------------
# 2. clamp on: the identical run completes and still emits the peak
# ---------------------------------------------------------------------------

def test_clamp_on_completes_where_v2_aborts(crash_bam, tmp_path):
    bed, _ = _call(crash_bam, tmp_path, "on_rescued",
                   dynamic_threshold_clamp=True, **CRASH_KWARGS)
    rows = [l for l in bed.read_text().splitlines() if l.strip()]
    assert rows, "the rescued run must still emit the coverage peak"
    for l in rows:
        assert len(l.split("\t")) == 6, "pasbed must stay BED6-parseable"


# ---------------------------------------------------------------------------
# 3. clamp on == clamp off, byte for byte, when v2 would not have crashed
# ---------------------------------------------------------------------------

def test_clamp_is_byte_identity_on_a_run_that_does_not_crash(crash_bam, tmp_path):
    bed_off, mtx_off = _call(crash_bam, tmp_path, "benign_off",
                             dynamic_threshold_clamp=False, **BENIGN_KWARGS)
    bed_on, mtx_on = _call(crash_bam, tmp_path, "benign_on",
                           dynamic_threshold_clamp=True, **BENIGN_KWARGS)
    assert bed_off.read_bytes() == bed_on.read_bytes(), (
        "--dynamic-threshold-clamp on changed the BED of a run that would "
        "not have crashed — it must only be able to rescue an abort"
    )
    assert mtx_off.read_bytes() == mtx_on.read_bytes(), (
        "--dynamic-threshold-clamp on changed the count matrix of a run "
        "that would not have crashed"
    )


# ---------------------------------------------------------------------------
# 4. the single copy of the rule, exhaustively
# ---------------------------------------------------------------------------

def test_resolve_l_end_identity_exhaustive():
    from sortedcontainers import SortedList

    from ema.countmatrix.dynamic_threshold import resolve_l_end

    for n in range(1, 9):
        arr = SortedList(range(100, 100 + n))
        # in range: both settings return the same element (the identity
        # guarantee that makes the clamp safe to offer at all)
        for thr in range(1, n + 1):
            assert resolve_l_end(arr, thr, clamp=False) == arr[-thr]
            assert resolve_l_end(arr, thr, clamp=True) == arr[-thr]
        # out of range: off raises exactly as v2, on returns the oldest end
        for thr in (n + 1, n + 7, 10 * n):
            with pytest.raises(IndexError):
                resolve_l_end(arr, thr, clamp=False)
            assert resolve_l_end(arr, thr, clamp=True) == arr[0]
    # an empty window has no oldest end: clamping must not invent one
    empty = SortedList()
    for clamp in (False, True):
        with pytest.raises(IndexError):
            resolve_l_end(empty, 3, clamp=clamp)


def test_guard_counts_clamp_hits():
    from sortedcontainers import SortedList

    from ema.countmatrix.dynamic_threshold import DynamicThresholdGuard

    g = DynamicThresholdGuard(True)
    arr = SortedList([5, 6, 7])
    assert g.l_end(arr, 2) == 6
    assert g.n_clamped == 0
    assert g.l_end(arr, 9) == 5
    assert g.n_clamped == 1
    assert (g.worst_threshold, g.worst_len) == (9, 3)


# ---------------------------------------------------------------------------
# 5. the pipeline finder carries the same expression — cover it too
# ---------------------------------------------------------------------------

def test_pipeline_path_is_bounded_too(crash_bam, tmp_path):
    """``peak_pipeline.finder_loop`` carries the identical index (VERDICTS §5
    names ``peak_pipeline.py`` alongside the monolithic loop).  Before the
    merge the finder died in a spawned subprocess and the abort surfaced as
    ``run_pipeline``'s RuntimeError; the bound must cover this path too."""
    bed_default, _ = _call(crash_bam, tmp_path, "pipe_default",
                           use_pipeline=True, **CRASH_KWARGS)
    rows = [l for l in bed_default.read_text().splitlines() if l.strip()]
    assert rows, "the bounded pipeline run must still emit the coverage peak"

    bed_off, _ = _call(crash_bam, tmp_path, "pipe_off", use_pipeline=True,
                       dynamic_threshold_clamp=False, **CRASH_KWARGS)
    assert bed_off.read_bytes() == bed_default.read_bytes(), (
        "the pipeline finder must honour the bound with the no-op flag too"
    )
