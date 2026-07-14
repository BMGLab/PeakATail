"""A4: `ema collapse` — RG-suffix run-tag pooling.

samtools-merge renames colliding @RG ids with an 8-hex suffix, so one library's
runs get per-run cb prefixes and stop pooling. collapse strips the suffix and
SUMS the columns that then share a canonical cb, without re-running peakcalling.
"""
from __future__ import annotations

import numpy as np
import scipy.io as sio
import scipy.sparse as sp
from click.testing import CliRunner

from ema.datasets.collapse import canonical_cb, collapse_columns, collapse_run


# ---- canonical_cb (pure) --------------------------------------------------

def test_canonical_cb_strips_8hex_run_suffix():
    assert canonical_cb("GSM123-StageI-4B7F9BA8_AAACGT") == "GSM123-StageI_AAACGT"


def test_canonical_cb_pools_two_runs_of_same_library():
    a = canonical_cb("GSM1-StageI-4B7F9BA8_AAAC")
    b = canonical_cb("GSM1-StageI-9C1D2E3F_AAAC")
    assert a == b == "GSM1-StageI_AAAC"


def test_canonical_cb_keeps_distinct_libraries_distinct():
    assert canonical_cb("GSM1-StageI-4B7F9BA8_AAAC") != canonical_cb("GSM2-StageI-4B7F9BA8_AAAC")


def test_canonical_cb_no_barcode_separator_unchanged():
    assert canonical_cb("nobarcodehere") == "nobarcodehere"


def test_canonical_cb_only_touches_prefix_not_barcode():
    # A barcode ending in -<8hex> must NOT be stripped (only the prefix is).
    assert canonical_cb("GSM1-StageI_AAAC-4B7F9BA8") == "GSM1-StageI_AAAC-4B7F9BA8"


def test_canonical_cb_prefix_without_suffix_unchanged():
    assert canonical_cb("GSM1-StageI_AAAC") == "GSM1-StageI_AAAC"


# ---- collapse_columns (pure, sparse) --------------------------------------

def _matrix(cols):
    # 2 PAS x len(cols) cells; column j has counts [j+1, j+1].
    data = np.array([[j + 1 for j in range(len(cols))],
                     [j + 1 for j in range(len(cols))]], dtype=int)
    return sp.csc_matrix(data)


def test_collapse_columns_sums_shared_canonical_cb():
    cbs = ["L1-AAAAAAAA_C1", "L1-BBBBBBBB_C1", "L1-AAAAAAAA_C2"]
    m = _matrix(cbs)  # cols: C1=1, C1=2, C2=3
    collapsed, uniq = collapse_columns(m, cbs)
    assert uniq == ["L1_C1", "L1_C2"]
    dense = collapsed.toarray()
    # L1_C1 pools cols 0 (=1) + 1 (=2) = 3; L1_C2 = col 2 = 3.
    assert dense[0].tolist() == [3, 3]
    assert dense[1].tolist() == [3, 3]


def test_collapse_columns_preserves_total_counts():
    cbs = ["L1-AAAAAAAA_C1", "L1-BBBBBBBB_C1", "L2-AAAAAAAA_C1"]
    m = _matrix(cbs)
    collapsed, uniq = collapse_columns(m, cbs)
    assert collapsed.sum() == m.sum()  # collapsing sums, never drops, counts


def test_collapse_columns_distinct_libraries_stay_separate():
    cbs = ["L1-AAAAAAAA_C1", "L2-AAAAAAAA_C1"]
    _, uniq = collapse_columns(_matrix(cbs), cbs)
    assert set(uniq) == {"L1_C1", "L2_C1"}


def test_collapse_columns_shape_mismatch_raises():
    m = _matrix(["a_C1", "b_C1"])
    try:
        collapse_columns(m, ["only_one"])
        assert False, "expected ValueError"
    except ValueError:
        pass


# ---- collapse_run (on-disk orchestration) ---------------------------------

def _write_run(run_dir, cbs):
    unified = run_dir / "unified"
    unified.mkdir(parents=True)
    sio.mmwrite(str(unified / "concatenated.mtx"), _matrix(cbs).tocoo(), field="integer")
    (unified / "concatenated_cbs.tsv").write_text("\n".join(cbs) + "\n")
    (run_dir / "posbed.bed").write_text("chr1\t10\t20\tpas1\t.\t+\n")
    (run_dir / "negbed.bed").write_text("chr1\t30\t40\tpas2\t.\t-\n")


def test_collapse_run_end_to_end(tmp_path):
    cbs = ["L1-AAAAAAAA_C1", "L1-BBBBBBBB_C1", "L1-AAAAAAAA_C2"]
    in_run = tmp_path / "in"
    out_run = tmp_path / "out"
    _write_run(in_run, cbs)
    summary = collapse_run(in_run, out_run)
    assert summary["n_cells_before"] == 3
    assert summary["n_cells_after"] == 2
    assert summary["n_pooled"] == 1
    # 2 distinct run-prefixes (L1-AAAAAAAA, L1-BBBBBBBB) collapse to 1 library.
    assert summary["n_libraries_before"] == 2
    assert summary["n_libraries_after"] == 1
    # corrected artifacts written; RG-independent beds copied verbatim.
    out_cbs = (out_run / "unified" / "concatenated_cbs.tsv").read_text().split()
    assert out_cbs == ["L1_C1", "L1_C2"]
    assert (out_run / "posbed.bed").exists()
    assert (out_run / "negbed.bed").exists()
    m = sio.mmread(str(out_run / "unified" / "concatenated.mtx"))
    assert m.sum() == _matrix(cbs).sum()  # counts preserved


def test_collapse_run_missing_input_raises(tmp_path):
    try:
        collapse_run(tmp_path / "nope", tmp_path / "out")
        assert False, "expected FileNotFoundError"
    except FileNotFoundError:
        pass


def test_collapse_cli_registered_and_runs(tmp_path):
    from ema.cli import main
    cbs = ["L1-AAAAAAAA_C1", "L1-BBBBBBBB_C1"]
    in_run = tmp_path / "in"
    _write_run(in_run, cbs)
    out_run = tmp_path / "out"
    res = CliRunner().invoke(
        main, ["collapse", "--in-run", str(in_run), "--out-run", str(out_run)]
    )
    assert res.exit_code == 0, res.output
    assert (out_run / "unified" / "concatenated_cbs.tsv").exists()
