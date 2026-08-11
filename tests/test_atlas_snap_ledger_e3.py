"""E3 wiring: atlas snap records the drop funnel into the provenance ledger.

The atlas snap is the biggest silent drop point (~59-63% of called PAS). With a
ProvenanceLedger passed in, snap_beds_to_atlas records every input PAS as either
a survivor (dropped_at="") with its snap distance, or a drop (dropped_at=
"atlas_snap") with a reason.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(
    not (shutil.which("bedtools") and shutil.which("sort")),
    reason="requires bedtools + sort on PATH",
)

from ema.datasets.atlas_snap import snap_beds_to_atlas
from ema.provenance import ProvenanceLedger, PAS_LEDGER_COLUMNS


def _read_ledger(run_dir: Path) -> list[dict]:
    text = (run_dir / "provenance" / "pas_ledger.tsv").read_text().splitlines()
    header = text[0].split("\t")
    assert header == PAS_LEDGER_COLUMNS
    return [dict(zip(header, line.split("\t"))) for line in text[1:]]


def test_ledger_records_survivor_and_drop(tmp_path: Path) -> None:
    run = tmp_path / "run"
    run.mkdir()
    unified = run / "unified"

    # Atlas site at chr1:1000 (+). One near PAS (summit within 50bp) survives;
    # one far PAS (summit ~5000bp away) is dropped.
    atlas = tmp_path / "atlas.bed"
    atlas.write_text("chr1\t1000\t1001\tA1\t0\t+\n")

    bed = tmp_path / "dsA.bed"
    bed.write_text(
        "chr1\t980\t1001\t1\t0\t+\n"     # summit end-1=1000 -> dist 0, survives
        "chr1\t6000\t6021\t2\t0\t+\n"    # summit 6020 -> ~5000bp, dropped
    )

    ledger = ProvenanceLedger(run)
    # mode="filter": this test asserts the pre-D9 drop behaviour (the far
    # PAS is dropped, not annotated-and-kept). See test_atlas_snap_annotate_d9.py
    # for the new default ("annotate") mode's ledger semantics.
    snap_beds_to_atlas(
        [bed], ["dsA"], atlas_bed=atlas, output_dir=unified,
        distance=50, strands=["+"], ledger=ledger, mode="filter",
    )
    ledger.flush()

    rows = _read_ledger(run)
    assert len(rows) == 2
    survivors = [r for r in rows if r["dropped_at"] == ""]
    drops = [r for r in rows if r["dropped_at"] == "atlas_snap"]
    assert len(survivors) == 1
    assert len(drops) == 1
    assert survivors[0]["orig_pas_key"] == "dsA::+::1"
    assert survivors[0]["snap_distance_bp"] == "0"
    assert drops[0]["orig_pas_key"] == "dsA::+::2"
    assert "atlas_distance" in drops[0]["drop_reason"]
    assert ledger.count_surviving("pas") == 1
