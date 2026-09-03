"""issue #101: ``--dynamic-threshold`` must not abort with a bare ``IndexError``.

``--dynamic-threshold --lambda-fold-change 2.0`` — the parameter reference's
own documented "find more PAS" setting — used to abort the caller at::

    l_end = data_array[-current_threshold]

because the dynamic estimator sets ``current_threshold`` from the local read
density and nothing bounded it by ``len(data_array)``.  It was reproduced on
the PBMC chr19+21 dev slice and, independently, on the GSE104556 mouse1
chr18+19 slice, so the defect is species-independent.  It is reachable only
with ``--dynamic-threshold`` (off by default), which is why no published
number is affected.

What is pinned here:

1. the crash fixture — 14 reads — now COMPLETES on both loop implementations
   (``peackcalling.py``'s monolithic loop and ``peak_pipeline.py``'s finder,
   which carries the identical expression).  Run against the unfixed tree,
   these two tests fail with ``IndexError`` / ``RuntimeError: finder``;
2. a dynamic-threshold run that never trips the bound is BYTE-IDENTICAL
   before and after the fix — "before" being the unbounded expression itself,
   reached by handing the loop a ``DynamicThresholdGuard(clamp=False)``;
3. :func:`~ema.countmatrix.dynamic_threshold.resolve_l_end`'s identity
   guarantee, exhaustively: for every in-range threshold the bounded and
   unbounded branches return the same element, so the bound can only change a
   run that would otherwise have aborted;
4. ``--floor-threshold 0`` is refused rather than silently returning
   ``data_array[-0] == data_array[0]``, the OLDEST end in the window;
5. the two help-text defects the same report filed — ``--pas-gap`` (a
   multi-dataset merge parameter, inert on a single-BAM run) and
   ``--min-cells`` / ``--min-pas-per-cell`` (AnnData filters applied after
   ``pasbed.bed`` is written).

CRASH MECHANICS OF THE FIXTURE (all forward-strand, one chromosome):
14 reads start at ``920+i``; under "fixed" geometry each end is rewritten to
``start + seq_len``, so the ends land at ``1011..1024`` — inside one
``lambda_window`` of each other, so the background deque never prunes.  With
``default_threshold = floor_threshold = 3`` the signal fires at read 3.  At
read 11 the deque tops 10 entries (the estimator's minimum), and with
``lambda_window=1000, lambda_fold_change=10000`` the threshold jumps to
``int(11/1000 * 10000) = 110`` while the live window holds 10 ends:
``data_array[-110]`` → ``IndexError``.  With ``lambda_fold_change=2.0`` the
same arithmetic yields ``max(3, int(0.011*2)) = 3`` and nothing ever trips
the bound, which is what makes the byte-identity comparison meaningful.
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

# Tripping arm: threshold jumps to int(11/1000 * 10000) = 110 >> window 10.
TRIPPING_KWARGS = dict(
    dynamic_threshold=True,
    floor_threshold=3,
    lambda_window=1000,
    lambda_fold_change=10000.0,
)
# Benign arm: same run, threshold stays at the floor; the bound never fires.
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


def _write_bam(path: Path) -> Path:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": CHROM_LEN}],
        "RG": [{"ID": "testsample", "SM": "testsample"}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i in range(N_READS):
            # The last four reads carry a genuine terminal poly(A) soft clip
            # so the clip_seeded seeder guarantees a tier-1 call whether or
            # not a coverage peak completes.  The clip does not touch the
            # mechanics above — read_check's coverage 5-tuple ignores soft
            # clips.
            clip = 12 if i >= N_READS - 4 else 0
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "C" * 80 + "A" * clip
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 920 + i
            a.mapping_quality = 60
            a.cigartuples = [(0, 80)] + ([(4, clip)] if clip else [])
            a.query_qualities = pysam.qualitystring_to_array("I" * (80 + clip))
            a.set_tag("CB", _barcode(i))
            a.set_tag("UB", f"UMI{i:06d}")
            a.set_tag("RG", "testsample")
            out.write(a)
    pysam.index(str(path))
    return path


@pytest.fixture(scope="module")
def dyn_bam(tmp_path_factory) -> Path:
    d = tmp_path_factory.mktemp("dyn_bounds_bam")
    return _write_bam(d / "dyn.bam")


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
# 1. the crash, fixed — both loop implementations
# ---------------------------------------------------------------------------

def test_dynamic_threshold_completes_where_it_used_to_abort(dyn_bam, tmp_path):
    """On the unfixed tree this raises ``IndexError: list index out of range``
    from ``peackcalling.py``'s ``l_end = data_array[-current_threshold]``."""
    bed, _ = _call(dyn_bam, tmp_path, "mono", **TRIPPING_KWARGS)
    rows = [ln for ln in bed.read_text().splitlines() if ln.strip()]
    assert rows, "the bounded run must still emit the coverage peak"
    for ln in rows:
        assert len(ln.split("\t")) == 6, "pasbed must stay BED6-parseable"


def test_pipeline_finder_completes_where_it_used_to_abort(dyn_bam, tmp_path):
    """``peak_pipeline.finder_loop`` carries the identical unbounded index.
    The finder runs in a spawned subprocess, so on the unfixed tree the abort
    surfaces as ``run_pipeline``'s ``RuntimeError: ... finder ...``."""
    bed, _ = _call(dyn_bam, tmp_path, "pipe", use_pipeline=True,
                   **TRIPPING_KWARGS)
    rows = [ln for ln in bed.read_text().splitlines() if ln.strip()]
    assert rows, "the bounded pipeline run must still emit the coverage peak"


# ---------------------------------------------------------------------------
# 2. a dynamic run that never trips the bound is byte-identical before/after
# ---------------------------------------------------------------------------

def test_non_crashing_dynamic_run_is_byte_identical_before_and_after(
        dyn_bam, tmp_path):
    """"Before" is the unbounded expression itself: a guard built with
    ``clamp=False`` evaluates ``data_array[-current_threshold]`` with nothing
    wrapped around it, which is exactly the pre-fix code path."""
    import ema.countmatrix.peackcalling as pc
    from ema.countmatrix.dynamic_threshold import DynamicThresholdGuard

    with patch.object(pc, "DynamicThresholdGuard",
                      lambda *a, **k: DynamicThresholdGuard(clamp=False)):
        bed_before, mtx_before = _call(dyn_bam, tmp_path, "benign_before",
                                       **BENIGN_KWARGS)
        before = (bed_before.read_bytes(), mtx_before.read_bytes())

    bed_after, mtx_after = _call(dyn_bam, tmp_path, "benign_after",
                                 **BENIGN_KWARGS)

    assert before[0] == bed_after.read_bytes(), (
        "bounding the look-back index changed the BED of a run that never "
        "trips the bound — it must only be able to rescue an abort"
    )
    assert before[1] == mtx_after.read_bytes(), (
        "bounding the look-back index changed the count matrix of a run "
        "that never trips the bound"
    )


# ---------------------------------------------------------------------------
# 3. the single copy of the rule, exhaustively
# ---------------------------------------------------------------------------

def test_resolve_l_end_identity_exhaustive():
    from sortedcontainers import SortedList

    from ema.countmatrix.dynamic_threshold import resolve_l_end

    for n in range(1, 9):
        arr = SortedList(range(100, 100 + n))
        # in range: both settings return the same element (the identity
        # guarantee that makes the bound safe to apply by default)
        for thr in range(1, n + 1):
            assert resolve_l_end(arr, thr, clamp=False) == arr[-thr]
            assert resolve_l_end(arr, thr, clamp=True) == arr[-thr]
        # out of range: unbounded raises, bounded returns the oldest end
        for thr in (n + 1, n + 7, 10 * n):
            with pytest.raises(IndexError):
                resolve_l_end(arr, thr, clamp=False)
            assert resolve_l_end(arr, thr, clamp=True) == arr[0]
    # an empty window has no oldest end: bounding must not invent one
    empty = SortedList()
    for clamp in (False, True):
        with pytest.raises(IndexError):
            resolve_l_end(empty, 3, clamp=clamp)


def test_resolve_l_end_bounds_by_default():
    from ema.countmatrix.dynamic_threshold import (
        DYNAMIC_THRESHOLD_CLAMP_DEFAULT, resolve_l_end,
    )

    assert DYNAMIC_THRESHOLD_CLAMP_DEFAULT is True
    assert resolve_l_end([10, 11, 12], 99) == 10


def test_guard_counts_bound_hits():
    from sortedcontainers import SortedList

    from ema.countmatrix.dynamic_threshold import DynamicThresholdGuard

    g = DynamicThresholdGuard()
    arr = SortedList([5, 6, 7])
    assert g.l_end(arr, 2) == 6
    assert g.n_clamped == 0
    assert g.l_end(arr, 9) == 5
    assert g.n_clamped == 1
    assert (g.worst_threshold, g.worst_len) == (9, 3)


# ---------------------------------------------------------------------------
# 4. --floor-threshold 0 is a silent wrong answer, not a crash — refuse it
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad", [0, -1, -7])
def test_validate_floor_threshold_refuses_useless_values(bad):
    from ema.countmatrix.dynamic_threshold import validate_floor_threshold

    with pytest.raises(ValueError, match="floor-threshold"):
        validate_floor_threshold(bad)


def test_validate_floor_threshold_passes_usable_values():
    from ema.countmatrix.dynamic_threshold import validate_floor_threshold

    for good in (1, 3, 500):
        assert validate_floor_threshold(good) == good


def test_peak_calling_refuses_floor_threshold_zero(dyn_bam, tmp_path):
    """The guard sits in ``peak_calling`` so YAML and library callers — and
    every spawned tile / chromosome worker, which all re-enter through it —
    get the same refusal Click gives on the command line."""
    with pytest.raises(ValueError, match="floor-threshold"):
        _call(dyn_bam, tmp_path, "floor0", dynamic_threshold=True,
              floor_threshold=0, lambda_window=1000, lambda_fold_change=2.0)


def test_cli_rejects_floor_threshold_zero():
    from ema.cli.config_schema import RunConfig

    spec = {f.name: f.metadata["spec"]
            for f in __import__("dataclasses").fields(RunConfig)}
    click_type = spec["floor_threshold"].click_type
    assert click_type is not None, "--floor-threshold must carry a range type"
    assert getattr(click_type, "min", None) == 1


# ---------------------------------------------------------------------------
# 5. the two help-text defects filed with the same report
# ---------------------------------------------------------------------------

def _description(field_name: str) -> str:
    import dataclasses

    from ema.cli.config_schema import RunConfig

    for f in dataclasses.fields(RunConfig):
        if f.name == field_name:
            return f.metadata["spec"].description
    raise AssertionError(f"no such RunConfig field: {field_name}")


def test_pas_gap_help_text_says_it_is_a_multi_dataset_merge_parameter():
    """``pas_gap`` is consumed at exactly one place, ``merge_pas_beds`` on the
    multi-dataset unified path — never on a single-BAM run.  The old help text
    ("minimum gap between PAS within a peak") described a filter that does not
    exist."""
    lowered = _description("pas_gap").lower()
    assert "gap between pas within a peak" not in lowered, (
        "the old help text described a within-peak filter that does not exist"
    )
    assert "multi-dataset" in lowered or "multi dataset" in lowered
    assert "no effect on a single-bam run" in lowered


def test_pas_gap_is_read_only_by_the_multi_dataset_merge():
    """Ground truth for the claim above: if a second consumer ever appears,
    this test fails and the help text has to be revisited."""
    import re
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "ema"
    # A *consumer* reads the value off the legacy args namespace; every other
    # mention in the tree is a declaration surface (config_schema, wizard,
    # yaml_loader) or prose in a docstring.
    reader = re.compile(r"""(getattr\(\s*args\s*,\s*['"]pas_gap['"]"""
                        r"""|args\.pas_gap|args\[['"]pas_gap['"]\])""")
    hits = []
    for py in sorted(root.rglob("*.py")):
        for n, line in enumerate(py.read_text().splitlines(), 1):
            if reader.search(line):
                hits.append(f"{py.relative_to(root)}:{n}: {line.strip()}")
    assert len(hits) == 1, f"pas_gap gained a consumer: {hits!r}"
    assert hits[0].startswith("main.py:"), hits[0]


@pytest.mark.parametrize("field_name", ["min_cells", "min_pas_per_cell"])
def test_min_cells_help_text_says_it_filters_the_anndata(field_name):
    """Both act in ``preprocessing()``, after ``pasbed.bed`` is written, so
    reading them as call-set filters is wrong."""
    text = _description(field_name).lower()
    assert "anndata" in text
    assert "pasbed.bed" in text
