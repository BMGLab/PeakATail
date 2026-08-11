"""Byte-identical regression test for the vectorized ``ema.matrixfilter.filter_cb``.

Perf pass: filter_cb used to do two full line-by-line passes over every
matrix file (Python ``split()`` + ``int()`` per token) plus a Python
``sorted()`` and a per-row f-string write loop. That's replaced with a
pandas C-parser read (once per file, reused for both passes), a
groupby-sum, an ``isin`` filter, and a vectorized contiguous re-index —
falling back to the untouched original algorithm (``_filter_cb_legacy``)
whenever a file's tokenization is too ambiguous for a fixed-width
vectorized parse.

``_reference_filter_cb`` below is a verbatim copy of the ORIGINAL
(pre-vectorization) implementation, captured from git history before the
edit. The test builds a synthetic multi-file input (one MatrixMarket-header
file, one legacy headerless file) with some barcodes below ``min_read``,
runs BOTH implementations, and asserts every output artifact — the sparse
matrix file bytes, the filtered-cb-list file bytes, and the module-level
``filtered_cb_list`` — are identical.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import ema.matrixfilter as mf


# ---------------------------------------------------------------------------
# Reference implementation: verbatim copy of the pre-vectorization
# ema.matrixfilter.filter_cb body, parameterized instead of reading globals.
# ---------------------------------------------------------------------------
def _reference_filter_cb(matrix_paths, cb_lookup, min_read,
                          sorted_corrected_sparse_path, filter_cb_file):
    _cb_lookup = list(cb_lookup)
    filtered_cb_list_ref: list = []

    # Pass 1: sum counts per column index (1-based in file)
    cb_counts: dict[str, int] = defaultdict(int)
    for path in matrix_paths:
        with open(path, "r") as fh:
            first_data_line = True
            for line in fh:
                if not line.strip() or line.startswith("%"):
                    continue
                columns = line.split()
                if len(columns) < 3:
                    continue
                if first_data_line:
                    first_data_line = False
                    try:
                        col_idx = int(columns[1])
                    except ValueError:
                        continue
                    if col_idx > len(_cb_lookup) and len(_cb_lookup) > 0:
                        continue
                    cb_counts[columns[1]] += int(columns[2])
                    continue
                cb_counts[columns[1]] += int(columns[2])

    keep_cb = {int(cb) for cb, count in cb_counts.items() if count >= min_read}

    # Pass 2: collect data rows for kept CBs
    lines_to_keep = []
    pas_set: set[int] = set()
    for path in matrix_paths:
        with open(path, "r") as fh:
            first_data_line = True
            for line in fh:
                if not line.strip() or line.startswith("%"):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                if first_data_line:
                    first_data_line = False
                    try:
                        col_idx = int(parts[1])
                    except ValueError:
                        continue
                    if col_idx > len(_cb_lookup) and len(_cb_lookup) > 0:
                        continue
                try:
                    item = [int(x) for x in parts[:3]]
                except ValueError:
                    continue
                if item[1] in keep_cb:
                    lines_to_keep.append(item)
                    pas_set.add(item[0])

    matrix_pas_header = max(pas_set) if pas_set else 0
    matrix_cb_header = len(keep_cb)
    matrix_nzero_header = len(lines_to_keep)
    sorted_spars_matrix = sorted(lines_to_keep, key=lambda x: x[1])
    last_corrected_index, corrected_index = 0, 0

    with open(sorted_corrected_sparse_path, "w") as out:
        out.write(mf.MATRIX_MARKET_HEADER)
        out.write(f"{matrix_pas_header} {matrix_cb_header} {matrix_nzero_header}\n")
        for item in sorted_spars_matrix:
            col_index = item[1]
            if item[1] != last_corrected_index:
                corrected_index += 1
                last_corrected_index = col_index
                filtered_cb_list_ref.append(_cb_lookup[col_index - 1])
                item[1] = corrected_index
                out.write(f"{item[0]} {item[1]} {item[2]}\n")
            else:
                item[1] = corrected_index
                out.write(f"{item[0]} {item[1]} {item[2]}\n")

    with open(filter_cb_file, "w") as f:
        for item in filtered_cb_list_ref:
            f.write(f"{item}\n")

    return filtered_cb_list_ref


CB_LOOKUP = ["BC1", "BC2", "BC3", "BC4", "BC5"]


def _write_matrix_market_file(path: Path) -> None:
    """File WITH a MatrixMarket header (comment lines + oversized dim line)."""
    path.write_text(
        "%%MatrixMarket matrix coordinate integer general\n"
        "%\n"
        "999 999 6\n"   # dimension header: col[1]=999 > n_cb(5) -> detected as header
        "1 1 3\n"
        "1 2 4\n"
        "2 1 2\n"
        "2 2 10\n"
        "3 3 1\n"
        "3 5 1\n"
    )


def _write_legacy_headerless_file(path: Path) -> None:
    """File with NO header — first line is real data."""
    path.write_text(
        "4 1 1\n"
        "4 3 6\n"
        "5 4 2\n"
        "5 5 10\n"
        "\n"          # blank line — must be skipped
    )


def test_filter_cb_vectorized_matches_reference_byte_identical(tmp_path):
    file_a = tmp_path / "mm_header.mtx"
    file_b = tmp_path / "legacy_headerless.mtx"
    _write_matrix_market_file(file_a)
    _write_legacy_headerless_file(file_b)
    matrix_paths = [str(file_a), str(file_b)]
    min_read = 5

    ref_mtx = tmp_path / "ref_matrix.mtx"
    ref_cb = tmp_path / "ref_cb.tsv"
    ref_filtered_cb_list = _reference_filter_cb(
        matrix_paths, CB_LOOKUP, min_read, ref_mtx, ref_cb,
    )

    new_mtx = tmp_path / "new_matrix.mtx"
    new_cb = tmp_path / "new_cb.tsv"
    mf.filter_cb(
        input_matrix_paths=matrix_paths,
        cb_list=CB_LOOKUP,
        sorted_corrected_sparse_path=new_mtx,
        min_read=min_read,
        filter_cb_file=new_cb,
    )

    assert ref_mtx.read_bytes() == new_mtx.read_bytes()
    assert ref_cb.read_bytes() == new_cb.read_bytes()
    assert mf.filtered_cb_list == ref_filtered_cb_list

    # Pin down the actual expected content too, so a future accidental
    # semantic change to BOTH implementations still gets caught.
    assert new_mtx.read_text() == (
        "%%MatrixMarket matrix coordinate integer general\n"
        "5 4 9\n"
        "1 1 3\n"
        "2 1 2\n"
        "4 1 1\n"
        "1 2 4\n"
        "2 2 10\n"
        "3 3 1\n"
        "4 3 6\n"
        "3 4 1\n"
        "5 4 10\n"
    )
    assert new_cb.read_text() == "BC1\nBC2\nBC3\nBC5\n"


def test_filter_cb_vectorized_matches_reference_no_kept_cbs(tmp_path):
    """Every cb below min_read -> empty outputs, both implementations agree."""
    file_a = tmp_path / "a.mtx"
    file_a.write_text("1 1 1\n1 2 2\n")
    matrix_paths = [str(file_a)]
    min_read = 100

    ref_mtx = tmp_path / "ref_matrix.mtx"
    ref_cb = tmp_path / "ref_cb.tsv"
    _reference_filter_cb(matrix_paths, CB_LOOKUP, min_read, ref_mtx, ref_cb)

    new_mtx = tmp_path / "new_matrix.mtx"
    new_cb = tmp_path / "new_cb.tsv"
    mf.filter_cb(
        input_matrix_paths=matrix_paths,
        cb_list=CB_LOOKUP,
        sorted_corrected_sparse_path=new_mtx,
        min_read=min_read,
        filter_cb_file=new_cb,
    )

    assert ref_mtx.read_bytes() == new_mtx.read_bytes()
    assert ref_cb.read_bytes() == new_cb.read_bytes()
    assert new_mtx.read_text() == "%%MatrixMarket matrix coordinate integer general\n0 0 0\n"
    assert new_cb.read_text() == ""


def test_filter_cb_vectorized_falls_back_on_ragged_row(tmp_path):
    """A data row with >3 whitespace tokens can't be vectorized safely —
    the fast path must detect this and transparently fall back to the
    legacy row-by-row implementation, still producing correct output
    (only the first 3 tokens of the ragged row are used, matching the
    original ``parts[:3]`` slicing).
    """
    file_a = tmp_path / "ragged.mtx"
    file_a.write_text(
        "1 1 3\n"
        "1 2 4 extra_token\n"  # ragged: 4 whitespace tokens
        "2 1 2\n"
    )
    matrix_paths = [str(file_a)]
    min_read = 1

    ref_mtx = tmp_path / "ref_matrix.mtx"
    ref_cb = tmp_path / "ref_cb.tsv"
    _reference_filter_cb(matrix_paths, CB_LOOKUP, min_read, ref_mtx, ref_cb)

    new_mtx = tmp_path / "new_matrix.mtx"
    new_cb = tmp_path / "new_cb.tsv"
    mf.filter_cb(
        input_matrix_paths=matrix_paths,
        cb_list=CB_LOOKUP,
        sorted_corrected_sparse_path=new_mtx,
        min_read=min_read,
        filter_cb_file=new_cb,
    )

    assert ref_mtx.read_bytes() == new_mtx.read_bytes()
    assert ref_cb.read_bytes() == new_cb.read_bytes()
