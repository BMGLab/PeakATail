"""D9: atlas-snap and internal-priming become ANNOTATE-not-drop by default.

The scientist is hunting ALTERNATIVE polyadenylation, so a PAS that doesn't
match the reference atlas (or looks internally-primed) is signal, not noise.
This module tests the atlas-snap half of that change:

  * ``snap_beds_to_atlas(..., mode="annotate")`` (the new default) keeps
    EVERY input PAS -- matched PAS snap onto the atlas summit, unmatched PAS
    keep their own coordinates -- and writes ``atlas_status.tsv`` /
    ``atlas_stats.json`` sidecars with the per-PAS / run-level match stats.
  * ``mode="filter"`` restores the pre-D9 drop behaviour (see
    test_atlas_snap_summit.py / test_atlas_snap_ledger_e3.py for those
    regression tests, updated to pass ``mode="filter"`` explicitly now that
    it is no longer the default).
  * The three new PAS-ledger columns (``atlas_match``, ``atlas_distance_bp``,
    ``internal_priming``) are additive and empty-by-default.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(
    not (shutil.which("bedtools") and shutil.which("sort")),
    reason="requires bedtools + sort on PATH",
)

from ema.datasets.atlas_snap import snap_beds_to_atlas
from ema.provenance import PAS_LEDGER_COLUMNS, ProvenanceLedger


def _write_bed(path: Path, rows: list[tuple]) -> None:
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")


def _setup(tmp_path: Path):
    atlas_bed = tmp_path / "atlas.bed"
    _write_bed(atlas_bed, [("chr1", 1000, 1001, "atlas1", "0", "+")])

    input_bed = tmp_path / "sample.bed"
    _write_bed(input_bed, [
        # summit end-1 = 999 -> dist 1 from atlas1 (1000) -> MATCHES within 50bp
        ("chr1", 990, 1000, "1", "0", "+"),
        # summit end-1 = 6019 -> ~5000bp away -> no atlas match
        ("chr1", 6000, 6020, "2", "0", "+"),
    ])
    return atlas_bed, input_bed


# ---------------------------------------------------------------------------
# 1. annotate mode (default): unmatched PAS KEPT, not dropped.
# ---------------------------------------------------------------------------
def test_annotate_mode_keeps_unmatched_pas(tmp_path: Path) -> None:
    atlas_bed, input_bed = _setup(tmp_path)
    out_dir = tmp_path / "out"

    # mode omitted -> default is "annotate".
    snapped_bed, mapping_path = snap_beds_to_atlas(
        [input_bed], ["dsA"], atlas_bed=atlas_bed, output_dir=out_dir,
        distance=50, strands=["+"],
    )

    mapping_rows = mapping_path.read_text().splitlines()[1:]
    assert len(mapping_rows) == 2, "both matched AND unmatched PAS must appear in the mapping"

    snapped_rows = snapped_bed.read_text().splitlines()
    assert len(snapped_rows) == 2, "both PAS must be written to the output BED"


def test_filter_mode_drops_unmatched_pas(tmp_path: Path) -> None:
    atlas_bed, input_bed = _setup(tmp_path)
    out_dir = tmp_path / "out"

    _snapped, mapping_path = snap_beds_to_atlas(
        [input_bed], ["dsA"], atlas_bed=atlas_bed, output_dir=out_dir,
        distance=50, strands=["+"], mode="filter",
    )

    mapping_rows = mapping_path.read_text().splitlines()[1:]
    assert len(mapping_rows) == 1, "filter mode must drop the unmatched PAS"


# ---------------------------------------------------------------------------
# 2. Stats sidecars: atlas_status.tsv + atlas_stats.json.
# ---------------------------------------------------------------------------
def test_annotate_mode_writes_status_and_stats_sidecars(tmp_path: Path) -> None:
    atlas_bed, input_bed = _setup(tmp_path)
    out_dir = tmp_path / "out"

    snap_beds_to_atlas(
        [input_bed], ["dsA"], atlas_bed=atlas_bed, output_dir=out_dir,
        distance=50, strands=["+"],
    )

    stats = json.loads((out_dir / "atlas_stats.json").read_text())
    assert stats["mode"] == "annotate"
    assert stats["n_atlas_matched"] == 1
    assert stats["n_atlas_unmatched"] == 1
    assert stats["atlas_match_rate"] == 0.5

    status_lines = (out_dir / "atlas_status.tsv").read_text().splitlines()
    header = status_lines[0].split("\t")
    assert header == ["new_pas_id", "atlas_match", "atlas_distance_bp"]
    rows = [dict(zip(header, line.split("\t"))) for line in status_lines[1:]]
    assert len(rows) == 2
    matched = [r for r in rows if r["atlas_match"] == "True"]
    unmatched = [r for r in rows if r["atlas_match"] == "False"]
    assert len(matched) == 1 and len(unmatched) == 1
    assert matched[0]["atlas_distance_bp"] == "1"
    # Unmatched PAS #2 was too far, but bedtools DID find a closest feature,
    # so the distance is a real (large) int, not "" (that's reserved for "no
    # atlas feature on this contig/strand at all").
    assert unmatched[0]["atlas_distance_bp"] not in ("", None)


# ---------------------------------------------------------------------------
# 3. Ledger: atlas_match/atlas_distance_bp populated for EVERY PAS in
#    annotate mode (both survive, dropped_at="").
# ---------------------------------------------------------------------------
def test_annotate_mode_ledger_both_survive_with_atlas_status(tmp_path: Path) -> None:
    atlas_bed, input_bed = _setup(tmp_path)
    out_dir = tmp_path / "out"
    run_dir = tmp_path / "run"

    ledger = ProvenanceLedger(run_dir)
    snap_beds_to_atlas(
        [input_bed], ["dsA"], atlas_bed=atlas_bed, output_dir=out_dir,
        distance=50, strands=["+"], ledger=ledger,
    )
    ledger.flush()

    text = (run_dir / "provenance" / "pas_ledger.tsv").read_text().splitlines()
    header = text[0].split("\t")
    assert header == PAS_LEDGER_COLUMNS
    rows = [dict(zip(header, line.split("\t"))) for line in text[1:]]
    assert len(rows) == 2
    # annotate mode: nothing dropped at the atlas-snap stage.
    assert all(r["dropped_at"] == "" for r in rows)

    by_key = {r["orig_pas_key"]: r for r in rows}
    assert by_key["dsA::+::1"]["atlas_match"] == "True"
    assert by_key["dsA::+::1"]["atlas_distance_bp"] == "1"
    assert by_key["dsA::+::2"]["atlas_match"] == "False"
    # internal_priming column always empty here -- ip runs at a later stage.
    assert by_key["dsA::+::1"]["internal_priming"] == ""
    assert by_key["dsA::+::2"]["internal_priming"] == ""


# ---------------------------------------------------------------------------
# 4. Default-off: when atlas isn't configured at all, none of this code runs
#    (covered structurally in ema/main.py -- directory_config.atlas gates
#    the whole snap_beds_to_atlas call). This test locks the CONTRACT that
#    an un-passed ledger produces no atlas-snap side effects and no crash.
# ---------------------------------------------------------------------------
def test_no_ledger_no_crash_annotate_mode(tmp_path: Path) -> None:
    atlas_bed, input_bed = _setup(tmp_path)
    out_dir = tmp_path / "out"
    # ledger=None (default) -- must not raise, must still write sidecars.
    snap_beds_to_atlas(
        [input_bed], ["dsA"], atlas_bed=atlas_bed, output_dir=out_dir,
        distance=50, strands=["+"],
    )
    assert (out_dir / "atlas_stats.json").exists()
