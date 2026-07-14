"""Contract tests for ema.provenance -- the E3 append-only drop ledgers."""
import csv
from pathlib import Path

import pytest

from ema.provenance import (
    CELL_LEDGER_COLUMNS,
    PAS_LEDGER_COLUMNS,
    ProvenanceLedger,
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        rows = [dict(zip(header, row)) for row in reader]
    return rows, header


def test_pas_and_cell_ledgers_survivors_and_drops(tmp_path):
    run_dir = tmp_path / "run1"
    ledger = ProvenanceLedger(run_dir)

    # Surviving PAS.
    ledger.record_pas(
        orig_pas_key="chr1:100-200:+",
        chrom="chr1",
        start=100,
        end=200,
        strand="+",
        unified_pas_id="PAS_1",
        snap_distance_bp=5,
        gene_id="ENSG1",
        gene_distance_bp=10,
        tier="high",
        last_stage="matrix_concat",
    )
    # Dropped PAS.
    ledger.record_pas(
        orig_pas_key="chr1:300-400:-",
        chrom="chr1",
        start=300,
        end=400,
        strand="-",
        dropped_at="atlas_snap",
        drop_reason="no nearby atlas PAS",
    )
    # Another surviving PAS.
    ledger.record_pas(orig_pas_key="chr2:1-2:+", unified_pas_id="PAS_2")

    # Surviving cell.
    ledger.record_cell(
        barcode="AAAA",
        dataset_id="ds1",
        total_reads=1000,
        n_pas=42,
        cluster="0",
    )
    # Dropped cell.
    ledger.record_cell(
        barcode="TTTT",
        dataset_id="ds1",
        total_reads=3,
        dropped_at="cb_filter",
        drop_reason="below min read count",
    )

    ledger.flush()

    pas_path = run_dir / "provenance" / "pas_ledger.tsv"
    cell_path = run_dir / "provenance" / "cell_ledger.tsv"
    assert pas_path.exists()
    assert cell_path.exists()

    pas_rows, pas_header = _read_tsv(pas_path)
    cell_rows, cell_header = _read_tsv(cell_path)

    assert pas_header == PAS_LEDGER_COLUMNS
    assert cell_header == CELL_LEDGER_COLUMNS

    assert len(pas_rows) == 3
    assert len(cell_rows) == 2

    # Surviving semantics: dropped_at == "" and drop_reason == "".
    surviving_pas = [r for r in pas_rows if r["dropped_at"] == ""]
    assert len(surviving_pas) == 2
    for r in surviving_pas:
        assert r["drop_reason"] == ""

    dropped_pas = [r for r in pas_rows if r["dropped_at"] != ""]
    assert len(dropped_pas) == 1
    assert dropped_pas[0]["dropped_at"] == "atlas_snap"
    assert dropped_pas[0]["drop_reason"] == "no nearby atlas PAS"

    surviving_cells = [r for r in cell_rows if r["dropped_at"] == ""]
    assert len(surviving_cells) == 1
    assert surviving_cells[0]["barcode"] == "AAAA"

    # count_surviving helper / property must match the on-disk survivor count.
    assert ledger.count_surviving("pas") == len(surviving_pas) == 2
    assert ledger.count_surviving("cell") == len(surviving_cells) == 1
    assert ledger.surviving_pas_count == 2


def test_none_to_empty_string_and_sanitization(tmp_path):
    run_dir = tmp_path / "run2"
    ledger = ProvenanceLedger(run_dir)

    ledger.record_pas(
        orig_pas_key="chr1:1-2:+",
        chrom=None,
        start=None,
        end=None,
        strand=None,
        unified_pas_id=None,
        snap_distance_bp=None,
        gene_id="gene\twith\ttabs",
        gene_distance_bp=None,
        tier=None,
        last_stage=None,
    )
    ledger.record_cell(
        barcode="bc\nwith\nnewlines",
        dataset_id="ds\rwith\rcr",
        total_reads=None,
        n_pas=None,
        cluster=None,
    )
    ledger.flush()

    pas_rows, _ = _read_tsv(run_dir / "provenance" / "pas_ledger.tsv")
    cell_rows, _ = _read_tsv(run_dir / "provenance" / "cell_ledger.tsv")

    row = pas_rows[0]
    for col in ("chrom", "start", "end", "strand", "unified_pas_id",
                "snap_distance_bp", "gene_distance_bp", "tier", "last_stage",
                "dropped_at", "drop_reason"):
        assert row[col] == ""
    # No literal tabs leaked into the value (the TSV parses cleanly at all,
    # which already proves this, but check no stray tab char remains).
    assert "\t" not in row["gene_id"]
    assert row["gene_id"] == "gene with tabs"

    crow = cell_rows[0]
    assert "\n" not in crow["barcode"]
    assert crow["barcode"] == "bc with newlines"
    assert "\r" not in crow["dataset_id"]
    assert crow["dataset_id"] == "ds with cr"
    assert crow["total_reads"] == ""
    assert crow["n_pas"] == ""
    assert crow["cluster"] == ""


def test_record_drop_requires_dropped_at_and_reason(tmp_path):
    run_dir = tmp_path / "run3"
    ledger = ProvenanceLedger(run_dir)

    ledger.record_drop(
        "pas",
        orig_pas_key="chr1:1-2:+",
        dropped_at="pas_gene",
        drop_reason="no gene within max distance",
    )
    ledger.record_drop(
        "cell",
        barcode="GGGG",
        dataset_id="ds1",
        dropped_at="preprocess",
        drop_reason="doublet",
    )

    with pytest.raises(ValueError):
        ledger.record_drop("pas", orig_pas_key="x", dropped_at="", drop_reason="reason")

    with pytest.raises(ValueError):
        ledger.record_drop("pas", orig_pas_key="x", dropped_at="stage", drop_reason="")

    with pytest.raises(ValueError):
        ledger.record_drop("bogus", dropped_at="stage", drop_reason="reason")

    ledger.flush()
    pas_rows, _ = _read_tsv(run_dir / "provenance" / "pas_ledger.tsv")
    cell_rows, _ = _read_tsv(run_dir / "provenance" / "cell_ledger.tsv")
    assert len(pas_rows) == 1
    assert pas_rows[0]["dropped_at"] == "pas_gene"
    assert pas_rows[0]["drop_reason"] == "no gene within max distance"
    assert len(cell_rows) == 1
    assert cell_rows[0]["dropped_at"] == "preprocess"
    assert cell_rows[0]["drop_reason"] == "doublet"


def test_ledgers_land_under_provenance_subdir(tmp_path):
    run_dir = tmp_path / "some_run"
    ledger = ProvenanceLedger(run_dir)
    ledger.record_pas(orig_pas_key="k")
    ledger.record_cell(barcode="b", dataset_id="d")
    ledger.write()

    provenance_dir = run_dir / "provenance"
    assert provenance_dir.is_dir()
    assert sorted(p.name for p in provenance_dir.iterdir()) == [
        "cell_ledger.tsv",
        "pas_ledger.tsv",
    ]


def test_flush_is_idempotent_and_rewritable(tmp_path):
    run_dir = tmp_path / "run4"
    ledger = ProvenanceLedger(run_dir)
    ledger.record_pas(orig_pas_key="k1")
    ledger.flush()
    ledger.flush()  # calling twice should not duplicate rows
    rows, _ = _read_tsv(run_dir / "provenance" / "pas_ledger.tsv")
    assert len(rows) == 1

    ledger.record_pas(orig_pas_key="k2", dropped_at="atlas_snap", drop_reason="x")
    ledger.flush()
    rows, _ = _read_tsv(run_dir / "provenance" / "pas_ledger.tsv")
    assert len(rows) == 2
