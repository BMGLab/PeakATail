"""B1 regression: strand-safe count routing in merge + concat.

Bug B1: within one dataset, pos and neg PAS are both numbered from 1, so the
count-routing key ``(dataset_id, pasnumber)`` collided across strands and neg
counts were routed onto the pos PAS row (or vice versa). The fix threads a
strand through the key: ``(dataset_id, strand, pasnumber)``.

These tests exercise the pure-Python key helpers plus concat_matrices without
bedtools (they feed a hand-written mapping TSV).
"""
from __future__ import annotations

from pathlib import Path

import scipy.io as sio

from ema.datasets.pas_merge import (
    _decode_pas_key,
    _encode_pas_key,
    concat_matrices,
)


def test_encode_decode_roundtrip() -> None:
    assert _encode_pas_key("dsA", "+", "5") == "dsA::+::5"
    assert _encode_pas_key("dsA", None, "5") == "dsA::5"
    assert _decode_pas_key("dsA::+::5") == ("dsA", "+", "5")
    assert _decode_pas_key("dsA::5") == ("dsA", "", "5")
    # dataset id containing '::' stays intact (rsplit from the right)
    assert _decode_pas_key("d::x::-::9") == ("d::x", "-", "9")


def _write_mtx(path: Path, triples, n_rows, n_cols) -> None:
    with open(path, "w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{n_rows} {n_cols} {len(triples)}\n")
        for r, c, v in triples:
            f.write(f"{r} {c} {v}\n")


def _write_cb(path: Path, cbs) -> None:
    path.write_text("\n".join(cbs) + "\n")


def test_concat_routes_pos_and_neg_pas_n_to_distinct_rows(tmp_path: Path) -> None:
    """dsA pos PAS #1 and dsA neg PAS #1 must land on DIFFERENT unified rows."""
    # Mapping (strand-aware): pos PAS 1 -> unified 10 ; neg PAS 1 -> unified 20.
    mapping = tmp_path / "map.tsv"
    mapping.write_text(
        "dataset_id\tstrand\told_pasnumber\tnew_pas_id\n"
        "dsA\t+\t1\t10\n"
        "dsA\t-\t1\t20\n"
    )

    # Two per-dataset matrices for the SAME dataset id, one per strand. Each has
    # its own row #1 with a distinctive count.
    pos_mtx = tmp_path / "pos.mtx"
    neg_mtx = tmp_path / "neg.mtx"
    _write_mtx(pos_mtx, [(1, 1, 7)], n_rows=1, n_cols=1)   # pos count 7
    _write_mtx(neg_mtx, [(1, 1, 99)], n_rows=1, n_cols=1)  # neg count 99
    pos_cb = tmp_path / "pos.cb"
    neg_cb = tmp_path / "neg.cb"
    _write_cb(pos_cb, ["CELL_POS"])
    _write_cb(neg_cb, ["CELL_NEG"])

    out_mtx = tmp_path / "out.mtx"
    out_cb = tmp_path / "out.cb"
    concat_matrices(
        mtx_paths=[pos_mtx, neg_mtx],
        dataset_ids=["dsA", "dsA"],
        cb_paths=[pos_cb, neg_cb],
        mapping_path=mapping,
        output_mtx=out_mtx,
        output_cb=out_cb,
        strands=["+", "-"],
    )

    m = sio.mmread(str(out_mtx)).tocsr()
    # unified rows sorted: 10 -> row 0, 20 -> row 1 ; columns: pos cell 0, neg cell 1
    # pos count 7 must be at unified row for id 10, neg cell column.
    assert m.shape == (2, 2)
    dense = m.toarray()
    # row0 (id10=pos) has the 7 in the pos cell col0; row1 (id20=neg) has 99 in col1
    assert dense[0, 0] == 7
    assert dense[1, 1] == 99
    # Crucially neither strand's count leaked onto the other's row.
    assert dense[0, 1] == 0
    assert dense[1, 0] == 0


def test_concat_without_strands_collapses_key(tmp_path: Path) -> None:
    """Legacy path (no strands, no strand column) still routes by (ds, pas)."""
    mapping = tmp_path / "map.tsv"
    mapping.write_text(
        "dataset_id\told_pasnumber\tnew_pas_id\n"
        "dsA\t1\t10\n"
    )
    mtx = tmp_path / "a.mtx"
    _write_mtx(mtx, [(1, 1, 5)], n_rows=1, n_cols=1)
    cb = tmp_path / "a.cb"
    _write_cb(cb, ["CELL"])
    out_mtx = tmp_path / "out.mtx"
    out_cb = tmp_path / "out.cb"
    concat_matrices(
        mtx_paths=[mtx], dataset_ids=["dsA"], cb_paths=[cb],
        mapping_path=mapping, output_mtx=out_mtx, output_cb=out_cb,
    )
    m = sio.mmread(str(out_mtx)).tocsr()
    assert m.toarray()[0, 0] == 5
