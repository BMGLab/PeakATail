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


# ---------------------------------------------------------------------------
# Streaming (low-memory) reader: same outputs, chunk boundaries and all.
#
# filter_cb no longer keeps a pandas frame with a per-row ``cb_str`` object
# column alive (~95 B/non-zero, an estimated 20-27 GB on the PBMC 10k
# matrices); it parses each file once in chunks into three integer arrays.
# These tests pin the byte-identity of the result against the reference
# implementation on randomised and adversarial inputs, and pin the
# deferral cases that keep the pathological semantics exact.
# ---------------------------------------------------------------------------

import random

import numpy as np
import pytest


def _run_both(tmp_path, matrix_paths, cb_lookup, min_read):
    ref_mtx, ref_cb = tmp_path / "ref.mtx", tmp_path / "ref_cb.tsv"
    ref_list = _reference_filter_cb(matrix_paths, cb_lookup, min_read, ref_mtx, ref_cb)
    new_mtx, new_cb = tmp_path / "new.mtx", tmp_path / "new_cb.tsv"
    mf.filter_cb(
        input_matrix_paths=matrix_paths, cb_list=cb_lookup,
        sorted_corrected_sparse_path=new_mtx, min_read=min_read,
        filter_cb_file=new_cb,
    )
    assert ref_mtx.read_bytes() == new_mtx.read_bytes()
    assert ref_cb.read_bytes() == new_cb.read_bytes()
    assert mf.filtered_cb_list == ref_list
    return new_mtx, new_cb


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
@pytest.mark.parametrize("min_read", [0, 1, 40])
def test_streaming_matches_reference_on_random_matrices(tmp_path, seed, min_read):
    rng = random.Random(seed)
    n_cb = 60
    lookup = [f"sample_BC{i:04d}" for i in range(n_cb)]
    paths = []
    for f in range(2):
        rows = []
        if f == 0:  # MatrixMarket dimension header (col[1] > n_cb)
            rows.append("%%MatrixMarket matrix coordinate integer general")
            rows.append("%")
            rows.append(f"5000 {n_cb + 500} 0")
        for _ in range(700):
            rows.append(f"{rng.randint(1, 400)} {rng.randint(1, n_cb)} {rng.randint(1, 9)}")
        rows.append("")  # blank line
        p = tmp_path / f"m{f}.mtx"
        p.write_text("\n".join(rows) + "\n")
        paths.append(str(p))
    _run_both(tmp_path, paths, lookup, min_read)


def test_streaming_crosses_chunk_boundaries(tmp_path, monkeypatch):
    """Header detection, grouping and the stable sort must not depend on
    where a chunk ends."""
    monkeypatch.setattr(mf, "_CHUNK_ROWS", 7)
    rng = random.Random(11)
    lookup = [f"BC{i}" for i in range(20)]
    rows = ["%%MatrixMarket matrix coordinate integer general", "999 999 0"]
    rows += [f"{rng.randint(1, 50)} {rng.randint(1, 20)} {rng.randint(1, 5)}" for _ in range(101)]
    p = tmp_path / "chunky.mtx"
    p.write_text("\n".join(rows) + "\n")
    _run_both(tmp_path, [str(p)], lookup, min_read=6)


def test_streaming_defers_non_canonical_barcode_token(tmp_path):
    """``007`` and ``7`` are two Pass-1 buckets in the reference algorithm;
    the integer grouping cannot express that, so the call defers."""
    p = tmp_path / "leading_zero.mtx"
    p.write_text("1 1 3\n1 01 4\n2 1 2\n")
    lookup = ["BC1", "BC2", "BC3"]
    with pytest.raises(mf._NonCanonicalBarcodeToken):
        mf._read_matrices_streaming([str(p)], len(lookup))
    _run_both(tmp_path, [str(p)], lookup, min_read=4)


def test_streaming_defers_float_valued_token(tmp_path):
    """A float count crashes the reference implementation in Pass 1
    (``int("4.5")``); the fast path must defer to it rather than silently
    truncating, so the SAME error surfaces."""
    p = tmp_path / "float.mtx"
    p.write_text("1 1 3\n1 2 4.5\n2 1 2\n")
    with pytest.raises(mf._RaggedMatrixData):
        mf._read_matrices_streaming([str(p)], 3)
    with pytest.raises(ValueError):
        _reference_filter_cb([str(p)], ["BC1", "BC2", "BC3"], 1,
                             tmp_path / "r.mtx", tmp_path / "r_cb.tsv")
    with pytest.raises(ValueError):
        mf.filter_cb(input_matrix_paths=[str(p)], cb_list=["BC1", "BC2", "BC3"],
                     sorted_corrected_sparse_path=tmp_path / "n.mtx",
                     min_read=1, filter_cb_file=tmp_path / "n_cb.tsv")


def test_streaming_defers_absurd_column_index(tmp_path):
    """A column index far outside the barcode list would need a multi-GB
    accumulator; the call defers to the reference implementation, which
    raises the same IndexError it always did."""
    p = tmp_path / "absurd.mtx"
    # first row is consumed as the dimension header, the second is data
    p.write_text("9 9 9\n1 999999999 2\n1 1 5\n")
    with pytest.raises(IndexError):
        _reference_filter_cb([str(p)], ["BC1", "BC2"], 1,
                             tmp_path / "r.mtx", tmp_path / "r_cb.tsv")
    with pytest.raises(IndexError):
        mf.filter_cb(input_matrix_paths=[str(p)], cb_list=["BC1", "BC2"],
                     sorted_corrected_sparse_path=tmp_path / "n.mtx",
                     min_read=1, filter_cb_file=tmp_path / "n_cb.tsv")


def test_streaming_handles_empty_and_comment_only_files(tmp_path):
    empty = tmp_path / "empty.mtx"
    empty.write_text("")
    comments = tmp_path / "comments.mtx"
    comments.write_text("%%MatrixMarket matrix coordinate integer general\n%\n")
    data = tmp_path / "data.mtx"
    data.write_text("1 1 9\n2 1 1\n")
    _run_both(tmp_path, [str(empty), str(comments), str(data)], ["BC1", "BC2"], min_read=5)


def test_streaming_drops_short_rows_like_the_reference(tmp_path):
    p = tmp_path / "short.mtx"
    p.write_text("1 1 5\n2 2\n3 1 6\n")
    _run_both(tmp_path, [str(p)], ["BC1", "BC2"], min_read=1)


def test_streaming_reader_returns_downcast_integer_arrays(tmp_path):
    p = tmp_path / "ints.mtx"
    p.write_text("7 7 7\n1 2 3\n4 5 6\n")   # first row consumed as header
    pas, cb, count = mf._read_matrices_streaming([str(p)], 3)
    assert pas.dtype == np.int32 and cb.dtype == np.int32 and count.dtype == np.int32
    assert pas.tolist() == [1, 4] and cb.tolist() == [2, 5] and count.tolist() == [3, 6]
