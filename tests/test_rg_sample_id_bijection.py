"""Composite CB ``"<sample_id>_<barcode>"`` must be a BIJECTION.

Regression tests for Stage 0d/0e of ``the Stage-0 fix campaign (PeakATail_wd analysis, 2026-08; see PR description and issues #67-#72)

Two independent defects made the mapping non-injective, and both merge cells
that belong to *different* samples into a single count-matrix column:

0e  ``read_check`` DISCARDED any ``RG`` containing ``_`` and substituted the
    run-level sample id (commit ``eb13529``).  ``peakatail merge`` stamps
    ``RG = dataset_id`` and ``samtools merge`` derives RG ids from file names,
    so ``sampleA_rep1.bam`` / ``sampleB_rep1.bam`` merged into one run
    collapsed onto ONE read group -> identical barcodes from the two samples
    shared a matrix column.  Sanitising instead (``"_"`` -> ``"-"``) is also
    non-injective (``a_b`` and ``a-b`` both become ``a-b``) and additionally
    breaks the ``cb.startswith(f"{ds_id}_")`` per-dataset column selector.

pre-0e ``BarcodeIndex`` split the composite on the FIRST underscore, so a
    sample id containing ``_`` was truncated and the remainder was fed to
    ``encode_cb`` as if it were a barcode -- yielding -1 for EVERY cell of the
    run, i.e. one shared column for the whole dataset.

The fix is an identity mapping plus a parser that splits on the LAST
underscore: the barcode half is a fixed-length ACGTN string and can never
contain ``_``, so ``"<sample_id>_<barcode>"`` decodes exactly for any sample
id.  See :func:`ema.countmatrix.indexing.split_cb`.
"""
from __future__ import annotations

import importlib
import sys
from unittest.mock import patch

import pytest

from ema.clustering.cross_dataset.jaccard import _strip_prefix
from ema.countmatrix.cb_encode import encode_cb
from ema.countmatrix.indexing import BarcodeIndex
from ema.datasets.collapse import _library, canonical_cb

BARCODE = "AAAACCCCGGGGTTTT"
BARCODE2 = "TTTTGGGGCCCCAAAA"


class _MockRead:
    """Minimal ``pysam.AlignedSegment`` stand-in for ``read_check``."""

    def __init__(
        self,
        cb: str = BARCODE,
        chrom: str = "chr1",
        start: int = 100,
        end: int = 250,
        is_reverse: bool = False,
        rg: str | None = None,
        is_unmapped: bool = False,
    ) -> None:
        self._cb = cb
        self._rg = rg
        self.reference_name = chrom
        self.reference_start = start
        self.reference_end = end
        self.is_reverse = is_reverse
        self.is_unmapped = is_unmapped

    def get_tag(self, tag: str) -> str:
        if tag == "CB":
            return self._cb
        if tag == "RG":
            if self._rg is None:
                raise KeyError(tag)
            return self._rg
        raise KeyError(tag)


_READ_CHECK_KWARGS: dict = dict(
    barcode="CB",
    barcode_len=16,
    seq_len=150,
    ignore_chro=set(),
)


def _read_check():
    """Import ``read_check`` with a patched argv (ema.config parses at import)."""
    with patch.object(sys, "argv", ["ema"]):
        import ema.countmatrix.read as _read_mod

        importlib.reload(_read_mod)
        return _read_mod.read_check


# ---------------------------------------------------------------------------
# 0e -- the headline test: two underscore-bearing RGs must stay distinct
# ---------------------------------------------------------------------------


def test_two_underscore_bearing_rgs_never_share_a_column() -> None:
    """`sampleA_rep1` and `sampleB_rep1` must occupy two matrix columns.

    This is the cross-sample contamination case.  Before the fix both RGs were
    replaced by the run-level sample id, so the SAME barcode seen in the two
    samples produced the SAME composite cb and therefore the SAME column.
    """
    read_check = _read_check()

    cbs = []
    for rg in ("sampleA_rep1", "sampleB_rep1"):
        read = _MockRead(cb=BARCODE, rg=rg)
        result = read_check(
            read, direction=False, **_READ_CHECK_KWARGS, sample_id="merged_run"
        )
        assert result != (0, 0, 0, 0, 0), f"read with RG {rg!r} was dropped"
        cbs.append(result[4])

    assert cbs[0] != cbs[1], (
        f"two distinct read groups produced the same composite cb: {cbs[0]!r}"
    )

    index = BarcodeIndex()
    col_a = index.get_index(cbs[0])
    col_b = index.get_index(cbs[1])
    assert col_a != col_b, (
        "sampleA_rep1 and sampleB_rep1 share matrix column "
        f"{col_a} -- their cells are silently pooled"
    )
    assert len(index.mapping) == 2

    # ... and the batch path must agree with the per-read path.
    batch = BarcodeIndex().get_indices_batch(cbs)
    assert batch[0] != batch[1]


def test_read_check_preserves_underscore_rg_verbatim() -> None:
    """The RG must survive byte-for-byte: `startswith(f"{ds_id}_")` depends on it.

    ``ema/main.py`` and ``ema/reannotate.py`` select each dataset's columns with
    ``cb.startswith(f"{ds_id}_")``.  Any rewrite of the prefix (dropping it, or
    ``"_"`` -> ``"-"``) makes that selector match nothing and the dataset is
    silently skipped with "no cells found".
    """
    read_check = _read_check()
    ds_id = "pbmc_10k_v3"

    result = read_check(
        _MockRead(cb=BARCODE, rg=ds_id), direction=False, **_READ_CHECK_KWARGS
    )
    cb = result[4]
    assert cb == f"{ds_id}_{BARCODE}", f"RG was rewritten: {cb!r}"
    assert cb.startswith(f"{ds_id}_")


def test_read_check_preserves_underscore_sample_id_fallback() -> None:
    """The no-RG fallback must not rewrite the run-level sample id either."""
    read_check = _read_check()

    result = read_check(
        _MockRead(cb=BARCODE, rg=None),
        direction=False,
        **_READ_CHECK_KWARGS,
        sample_id="pbmc_10k_v3",
    )
    assert result[4] == f"pbmc_10k_v3_{BARCODE}"


@pytest.mark.parametrize(
    "rg_a, rg_b",
    [
        ("sampleA_rep1", "sampleB_rep1"),   # the peakatail merge / samtools case
        ("lib_1", "lib-1"),                 # what "_"->"-" sanitising collides
        ("a_b_c", "a_b-c"),                 # ... and again, one level deeper
        ("pbmc_10k_v3", "pbmc"),            # truncation to the first token
    ],
)
def test_distinct_sample_ids_get_distinct_columns(rg_a: str, rg_b: str) -> None:
    """No pair of distinct sample ids may ever land in the same column."""
    index = BarcodeIndex()
    col_a = index.get_index(f"{rg_a}_{BARCODE}")
    col_b = index.get_index(f"{rg_b}_{BARCODE}")
    assert col_a != col_b, f"{rg_a!r} and {rg_b!r} collided in column {col_a}"


# ---------------------------------------------------------------------------
# the parser: split on the LAST underscore
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "sample_id",
    ["plain", "with_underscore", "pbmc_10k_v3", "GSM123-StageI", "a_b_c_d"],
)
def test_split_cb_roundtrips_any_sample_id(sample_id: str) -> None:
    """``split_cb(f"{s}_{bc}")`` must return exactly ``(s, bc)``."""
    # imported inside the test (not at module scope) so this file still
    # COLLECTS on revisions predating split_cb -- the regression then shows up
    # as these two tests failing rather than the whole module erroring out.
    from ema.countmatrix.indexing import split_cb

    assert split_cb(f"{sample_id}_{BARCODE}") == (sample_id, BARCODE)


def test_split_cb_without_separator() -> None:
    from ema.countmatrix.indexing import split_cb

    assert split_cb(BARCODE) == ("", BARCODE)


def test_barcode_index_encodes_barcode_of_underscore_sample() -> None:
    """The barcode half must reach ``encode_cb`` intact, not as -1.

    Before the fix, ``"pbmc_10k_v3_AAAACCCCGGGGTTTT"`` parsed as sample
    ``"pbmc"`` + barcode ``"10k_v3_AAAACCCCGGGGTTTT"``; the non-ACGT barcode
    encoded to -1, so every cell of the dataset shared the key
    ``("pbmc", -1)`` -- a one-column matrix.
    """
    index = BarcodeIndex()
    col1 = index.get_index(f"pbmc_10k_v3_{BARCODE}")
    col2 = index.get_index(f"pbmc_10k_v3_{BARCODE2}")
    assert col1 != col2, "two barcodes of one dataset collapsed into one column"

    keys = list(index._cb_total)
    assert {k[0] for k in keys} == {"pbmc_10k_v3"}
    assert -1 not in {k[1] for k in keys}, "barcode half was mangled before encode_cb"
    assert {k[1] for k in keys} == {encode_cb(BARCODE), encode_cb(BARCODE2)}

    # mapping must reconstruct the original strings byte-for-byte
    assert set(index.mapping) == {
        f"pbmc_10k_v3_{BARCODE}",
        f"pbmc_10k_v3_{BARCODE2}",
    }


def test_batch_and_single_index_paths_agree_on_underscore_sample() -> None:
    """``get_indices_batch`` and ``get_index`` must not diverge."""
    cbs = [f"pbmc_10k_v3_{BARCODE}", f"pbmc_10k_v3_{BARCODE2}", f"other_lib_{BARCODE}"]
    single = BarcodeIndex()
    assert [single.get_index(c) for c in cbs] == BarcodeIndex().get_indices_batch(cbs)
    assert len(single.mapping) == 3


# ---------------------------------------------------------------------------
# downstream consumers of the composite
# ---------------------------------------------------------------------------


def test_canonical_cb_strips_run_suffix_from_underscore_library() -> None:
    """``peakatail collapse`` must pool runs of a library whose id contains ``_``."""
    assert (
        canonical_cb(f"lib_A-4B7F9BA8_{BARCODE}") == f"lib_A_{BARCODE}"
    )
    assert (
        canonical_cb(f"lib_A-C0FFEE12_{BARCODE}") == f"lib_A_{BARCODE}"
    )
    # different libraries must NOT pool
    assert canonical_cb(f"lib_B-4B7F9BA8_{BARCODE}") != f"lib_A_{BARCODE}"


def test_library_prefix_keeps_underscores() -> None:
    assert _library(f"pbmc_10k_v3_{BARCODE}") == "pbmc_10k_v3"


def test_jaccard_strip_prefix_handles_underscore_sample() -> None:
    """Cross-dataset Jaccard must recover the biological barcode."""
    assert _strip_prefix(f"pbmc_10k_v3_{BARCODE}") == BARCODE
    assert _strip_prefix(f"pbmc_10k_v3_{BARCODE}-1") == f"{BARCODE}-1"
    # the synthetic-name guard still holds
    assert _strip_prefix("cellA_0") == "cellA_0"


def test_reannotate_dataset_ids_from_composite_cbs() -> None:
    """The ds-id derivation in ``ema/reannotate.py`` must not truncate ids.

    Mirrors ``reannotate.py``'s ``cb.rsplit("_", 1)[0]`` over ``unified_cbs``
    and the ``startswith(f"{ds_id}_")`` selector that consumes its output.
    """
    all_cbs = [
        f"sampleA_rep1_{BARCODE}",
        f"sampleA_rep1_{BARCODE2}",
        f"sampleB_rep1_{BARCODE}",
    ]
    unique_ds_ids = list(dict.fromkeys(cb.rsplit("_", 1)[0] for cb in all_cbs))
    assert unique_ds_ids == ["sampleA_rep1", "sampleB_rep1"]

    selected = {
        ds: [i for i, cb in enumerate(all_cbs) if cb.rsplit("_", 1)[0] == ds]
        for ds in unique_ds_ids
    }
    assert selected == {"sampleA_rep1": [0, 1], "sampleB_rep1": [2]}


def test_per_dataset_selector_is_not_prefix_ambiguous() -> None:
    """Dataset "a" must not claim dataset "a_b"'s cells.

    Guards the selector used in ``ema/main.py`` and ``ema/reannotate.py``.
    ``cb.startswith(f"{ds_id}_")`` matched both ``a_<barcode>`` and
    ``a_b_<barcode>`` once sample ids are allowed to contain ``_``.
    """
    all_cbs = [f"a_{BARCODE}", f"a_b_{BARCODE}", f"a_b_{BARCODE2}"]
    selected = {
        ds: [i for i, cb in enumerate(all_cbs) if cb.rsplit("_", 1)[0] == ds]
        for ds in ("a", "a_b")
    }
    assert selected == {"a": [0], "a_b": [1, 2]}

    # the old selector was ambiguous — kept here so the regression is explicit
    old_selector = [i for i, cb in enumerate(all_cbs) if cb.startswith("a_")]
    assert old_selector == [0, 1, 2]


# ---------------------------------------------------------------------------
# 0d -- the unmapped guard (commit 19790c8) needs a test of its own
# ---------------------------------------------------------------------------


def test_read_check_skips_unmapped_read() -> None:
    """Unmapped records carry a CB but no reference_end; they must be skipped."""
    read_check = _read_check()
    read = _MockRead(cb=BARCODE, is_unmapped=True, chrom=None, start=0, end=None)
    assert read_check(read, direction=False, **_READ_CHECK_KWARGS) == (0, 0, 0, 0, 0)


def test_read_check_skips_mapped_read_without_reference_end() -> None:
    """A CIGAR-less record (reference_end None) must not raise TypeError."""
    read_check = _read_check()
    read = _MockRead(cb=BARCODE, is_unmapped=False, end=None)
    assert read_check(read, direction=False, **_READ_CHECK_KWARGS) == (0, 0, 0, 0, 0)


def test_read_check_strips_cellranger_gem_suffix() -> None:
    """CellRanger's ``-1`` GEM-group suffix is stripped before the length check."""
    read_check = _read_check()
    result = read_check(
        _MockRead(cb=f"{BARCODE}-1", rg="sampleA_rep1"),
        direction=False,
        **_READ_CHECK_KWARGS,
    )
    assert result[4] == f"sampleA_rep1_{BARCODE}"


# ---------------------------------------------------------------------------
# same thing again, against REAL pysam.AlignedSegment objects
# ---------------------------------------------------------------------------
# 0d was a mock-drift bug: `_MockRead` did not carry `is_unmapped`, so a guard
# that works on real reads looked broken in the suite. These tests run
# read_check over a BAM written by pysam, so no attribute of a real
# AlignedSegment can silently go missing again.


def _write_bam(path, records):
    """Write a tiny BAM. *records* are (name, rg, cb, start, unmapped) tuples."""
    import pysam

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 100_000}],
        "RG": [
            {"ID": "sampleA_rep1", "SM": "sampleA"},
            {"ID": "sampleB_rep1", "SM": "sampleB"},
        ],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        h = out.header
        for name, rg, cb, start, unmapped in records:
            seg = pysam.AlignedSegment(h)
            seg.query_name = name
            seg.query_sequence = "A" * 150
            seg.query_qualities = pysam.qualitystring_to_array("I" * 150)
            if unmapped:
                seg.flag = 4
                seg.reference_id = -1
                seg.reference_start = -1
            else:
                seg.flag = 0
                seg.reference_id = 0
                seg.reference_start = start
                seg.mapping_quality = 60
                seg.cigarstring = "150M"
            seg.set_tag("CB", cb)
            seg.set_tag("RG", rg)
            out.write(seg)
    return path


def test_real_bam_two_underscore_rgs_keep_separate_columns(tmp_path) -> None:
    """End-to-end on a real BAM: merged samples must not share a column."""
    import pysam

    read_check = _read_check()
    bam = _write_bam(
        tmp_path / "merged.bam",
        [
            # CellRanger-style CB with the "-1" GEM suffix, same cell barcode
            # observed in BOTH samples — the contamination case.
            ("rA1", "sampleA_rep1", f"{BARCODE}-1", 100, False),
            ("rB1", "sampleB_rep1", f"{BARCODE}-1", 200, False),
            ("rA2", "sampleA_rep1", f"{BARCODE2}-1", 300, False),
            ("rU", "sampleA_rep1", f"{BARCODE}-1", 0, True),
        ],
    )

    index = BarcodeIndex()
    seen: dict[str, int] = {}
    n_skipped = 0
    with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as fh:
        for read in fh:
            result = read_check(read, direction=False, **_READ_CHECK_KWARGS)
            if result == (0, 0, 0, 0, 0):
                n_skipped += 1
                continue
            seen[result[4]] = index.get_index(result[4])

    assert n_skipped == 1, "the unmapped record should be the only skipped read"
    assert set(seen) == {
        f"sampleA_rep1_{BARCODE}",
        f"sampleB_rep1_{BARCODE}",
        f"sampleA_rep1_{BARCODE2}",
    }
    assert len(set(seen.values())) == 3, (
        f"expected 3 distinct matrix columns, got {seen}"
    )
    assert seen[f"sampleA_rep1_{BARCODE}"] != seen[f"sampleB_rep1_{BARCODE}"], (
        "the same barcode in two merged samples landed in ONE column"
    )
