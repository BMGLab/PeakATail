"""Atlas-as-pure-overlay: atlas_mode="annotate" no longer routes through
``snap_beds_to_atlas`` (which decides PAS *identity* -- matched PAS get
snapped onto atlas coordinates, unmatched PAS get their own 1:1 id per INPUT
row, before any cross-dataset merge). Run per-dataset before a merge, that
approach fragments a genuinely novel PAS seen in two datasets at the same
locus into two "unmatched" features instead of merging them into one.

The fix: ``ema/main.py``'s atlas_mode="annotate" (default) path now builds
the unified PAS set EXACTLY like a no-atlas run --
``ema.datasets.pas_merge.merge_pas_beds`` proximity-merges ALL called PAS
across datasets, atlas plays NO role in merging or id assignment -- and
``ema.datasets.atlas_annotate.annotate_pas_against_atlas`` overlays
atlas_match/atlas_distance_bp onto that ALREADY-FINAL unified set afterwards.

This module tests that overlay function directly, and the
merge-then-overlay COMPOSITION that ``ema/main.py`` now performs for
atlas_mode="annotate" (mirroring how test_atlas_snap_annotate_d9.py tests
``snap_beds_to_atlas`` directly for the old routing). ``mode="filter"``
(``snap_beds_to_atlas``, unchanged) keeps its own dedicated regression tests
in test_atlas_snap_annotate_d9.py / test_atlas_snap_summit.py /
test_atlas_snap_ledger_e3.py.
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

from ema.datasets.atlas_annotate import annotate_pas_against_atlas
from ema.datasets.atlas_snap import snap_beds_to_atlas
from ema.datasets.pas_merge import merge_pas_beds


def _write_bed(path: Path, rows: list[tuple]) -> None:
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")


def _read_status(status_path: Path) -> list[dict]:
    lines = status_path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


# ---------------------------------------------------------------------------
# 1. annotate_pas_against_atlas is a pure overlay: never drops a row.
# ---------------------------------------------------------------------------
def test_overlay_never_drops_any_unified_pas(tmp_path: Path) -> None:
    atlas_bed = tmp_path / "atlas.bed"
    _write_bed(atlas_bed, [("chr1", 1000, 1001, "atlas1", "0", "+")])

    unified_bed = tmp_path / "unified.bed"
    _write_bed(unified_bed, [
        ("chr1", 990, 1000, "1", "0", "+"),   # summit 999 -> dist 1 -> matches
        ("chr1", 6000, 6020, "2", "0", "+"),  # summit 6019 -> far -> no match
        ("chr2", 100, 110, "3", "0", "-"),    # different contig -> no atlas feature at all
    ])

    out_dir = tmp_path / "unified_out"
    status_path, stats_path = annotate_pas_against_atlas(
        unified_pasbed_path=unified_bed,
        atlas_bed_path=atlas_bed,
        distance=50,
        output_dir=out_dir,
    )

    rows = _read_status(status_path)
    assert {r["unified_pas_id"] for r in rows} == {"1", "2", "3"}
    assert len(rows) == 3  # NOTHING dropped

    by_id = {r["unified_pas_id"]: r for r in rows}
    assert by_id["1"]["atlas_match"] == "True"
    assert by_id["1"]["atlas_distance_bp"] == "1"
    assert by_id["2"]["atlas_match"] == "False"
    assert by_id["2"]["atlas_distance_bp"] not in ("", None)
    # chr2 has no atlas feature at all (different contig/strand) -> real
    # distance is unknowable, "" is reserved for exactly this case.
    assert by_id["3"]["atlas_match"] == "False"
    assert by_id["3"]["atlas_distance_bp"] == ""

    # Input unified bed itself is untouched (pure overlay, no rewrite).
    assert unified_bed.read_text().splitlines() == [
        "chr1\t990\t1000\t1\t0\t+",
        "chr1\t6000\t6020\t2\t0\t+",
        "chr2\t100\t110\t3\t0\t-",
    ]


def test_overlay_stats_sidecar(tmp_path: Path) -> None:
    atlas_bed = tmp_path / "atlas.bed"
    _write_bed(atlas_bed, [("chr1", 1000, 1001, "atlas1", "0", "+")])

    unified_bed = tmp_path / "unified.bed"
    _write_bed(unified_bed, [
        ("chr1", 990, 1000, "1", "0", "+"),
        ("chr1", 6000, 6020, "2", "0", "+"),
    ])

    out_dir = tmp_path / "unified_out"
    _status_path, stats_path = annotate_pas_against_atlas(
        unified_pasbed_path=unified_bed,
        atlas_bed_path=atlas_bed,
        distance=50,
        output_dir=out_dir,
    )
    stats = json.loads(stats_path.read_text())
    assert stats == {
        "mode": "annotate",
        "n_atlas_matched": 1,
        "n_atlas_unmatched": 1,
        "atlas_match_rate": 0.5,
    }


# ---------------------------------------------------------------------------
# 2. The merge-then-overlay COMPOSITION ema/main.py now performs for
#    atlas_mode="annotate": unified PAS ids are IDENTICAL to a no-atlas
#    merge_pas_beds run of the same input -- atlas never changes the merge.
# ---------------------------------------------------------------------------
def test_annotate_composition_ids_match_no_atlas_merge(tmp_path: Path) -> None:
    bed_a = tmp_path / "dsA.bed"
    _write_bed(bed_a, [
        ("chr1", 990, 1000, "1", "0", "+"),
        ("chr1", 6000, 6020, "2", "0", "+"),
    ])
    bed_b = tmp_path / "dsB.bed"
    _write_bed(bed_b, [
        ("chr1", 991, 1001, "1", "0", "+"),  # close to dsA's PAS #1 -> should merge
    ])

    # No-atlas baseline.
    no_atlas_dir = tmp_path / "no_atlas"
    no_atlas_bed, _mapping = merge_pas_beds(
        [bed_a, bed_b], ["dsA", "dsB"], output_dir=no_atlas_dir,
        gap=100, strands=["+", "+"],
    )
    no_atlas_ids = sorted(
        line.split("\t")[3] for line in no_atlas_bed.read_text().splitlines() if line
    )
    no_atlas_coords = {
        line.split("\t")[3]: tuple(line.split("\t")[:3])
        for line in no_atlas_bed.read_text().splitlines() if line
    }

    # atlas_mode="annotate" composition: merge (atlas-blind) + overlay.
    atlas_bed = tmp_path / "atlas.bed"
    _write_bed(atlas_bed, [("chr1", 1000, 1001, "atlas1", "0", "+")])

    annotate_dir = tmp_path / "annotate"
    annotate_bed, _mapping2 = merge_pas_beds(
        [bed_a, bed_b], ["dsA", "dsB"], output_dir=annotate_dir,
        gap=100, strands=["+", "+"],
    )
    status_path, _stats_path = annotate_pas_against_atlas(
        unified_pasbed_path=annotate_bed,
        atlas_bed_path=atlas_bed,
        distance=50,
        output_dir=annotate_dir,
    )
    annotate_ids = sorted(
        line.split("\t")[3] for line in annotate_bed.read_text().splitlines() if line
    )
    annotate_coords = {
        line.split("\t")[3]: tuple(line.split("\t")[:3])
        for line in annotate_bed.read_text().splitlines() if line
    }

    # Same ids, same coordinates -- the atlas overlay changed NOTHING about
    # the merge.
    assert annotate_ids == no_atlas_ids
    assert annotate_coords == no_atlas_coords

    # And the overlay covers exactly those ids.
    status_ids = {r["unified_pas_id"] for r in _read_status(status_path)}
    assert status_ids == set(annotate_ids)


# ---------------------------------------------------------------------------
# 3. Anti-fragmentation guarantee: a genuinely NOVEL PAS (no atlas match)
#    seen in TWO datasets at the same locus becomes ONE unified feature, not
#    two -- exactly what snap_beds_to_atlas's per-dataset id-minting would
#    have fragmented.
# ---------------------------------------------------------------------------
def test_novel_pas_in_two_datasets_merges_to_one_feature(tmp_path: Path) -> None:
    # Atlas has NOTHING near this locus -> both datasets' PAS here are novel.
    atlas_bed = tmp_path / "atlas.bed"
    _write_bed(atlas_bed, [("chr5", 50000, 50001, "atlas_far", "0", "+")])

    bed_a = tmp_path / "dsA.bed"
    _write_bed(bed_a, [("chr1", 2000, 2010, "1", "0", "+")])  # summit 2009
    bed_b = tmp_path / "dsB.bed"
    _write_bed(bed_b, [("chr1", 2003, 2013, "1", "0", "+")])  # summit 2012, 3bp away -> merges

    out_dir = tmp_path / "out"
    unified_bed, mapping_path = merge_pas_beds(
        [bed_a, bed_b], ["dsA", "dsB"], output_dir=out_dir,
        gap=100, strands=["+", "+"],
    )
    status_path, stats_path = annotate_pas_against_atlas(
        unified_pasbed_path=unified_bed,
        atlas_bed_path=atlas_bed,
        distance=50,
        output_dir=out_dir,
    )

    unified_rows = [l for l in unified_bed.read_text().splitlines() if l]
    assert len(unified_rows) == 1, "novel PAS from two datasets at one locus must merge to ONE feature"

    mapping_rows = mapping_path.read_text().splitlines()[1:]
    assert len(mapping_rows) == 2, "both original PAS still map onto that one unified feature"
    unified_id = unified_rows[0].split("\t")[3]
    assert all(row.split("\t")[-1] == unified_id for row in mapping_rows)

    status_rows = _read_status(status_path)
    assert len(status_rows) == 1
    assert status_rows[0]["unified_pas_id"] == unified_id
    assert status_rows[0]["atlas_match"] == "False"  # novel PAS, correctly unmatched

    stats = json.loads(stats_path.read_text())
    assert stats["n_atlas_matched"] == 0
    assert stats["n_atlas_unmatched"] == 1


def test_novel_pas_in_two_datasets_would_fragment_under_old_snap_routing(tmp_path: Path) -> None:
    """Contrast case: confirms the OLD per-dataset snap_beds_to_atlas routing
    really did fragment this exact scenario, motivating the fix above.
    """
    atlas_bed = tmp_path / "atlas.bed"
    _write_bed(atlas_bed, [("chr5", 50000, 50001, "atlas_far", "0", "+")])

    bed_a = tmp_path / "dsA.bed"
    _write_bed(bed_a, [("chr1", 2000, 2010, "1", "0", "+")])
    bed_b = tmp_path / "dsB.bed"
    _write_bed(bed_b, [("chr1", 2003, 2013, "1", "0", "+")])

    out_dir = tmp_path / "out"
    snapped_bed, _mapping = snap_beds_to_atlas(
        [bed_a, bed_b], ["dsA", "dsB"], atlas_bed=atlas_bed, output_dir=out_dir,
        distance=50, strands=["+", "+"], mode="annotate",
    )
    rows = [l for l in snapped_bed.read_text().splitlines() if l]
    assert len(rows) == 2, "old routing fragments the same novel PAS into two features"


# ---------------------------------------------------------------------------
# 4. filter mode is untouched by this change (still goes through
#    snap_beds_to_atlas, still snaps+drops).
# ---------------------------------------------------------------------------
def test_filter_mode_still_snaps_and_drops(tmp_path: Path) -> None:
    atlas_bed = tmp_path / "atlas.bed"
    _write_bed(atlas_bed, [("chr1", 1000, 1001, "atlas1", "0", "+")])

    input_bed = tmp_path / "sample.bed"
    _write_bed(input_bed, [
        ("chr1", 990, 1000, "1", "0", "+"),   # matches -> kept
        ("chr1", 6000, 6020, "2", "0", "+"),  # far -> dropped
    ])

    out_dir = tmp_path / "out"
    snapped_bed, mapping_path = snap_beds_to_atlas(
        [input_bed], ["dsA"], atlas_bed=atlas_bed, output_dir=out_dir,
        distance=50, strands=["+"], mode="filter",
    )
    rows = [l for l in snapped_bed.read_text().splitlines() if l]
    assert len(rows) == 1  # the far PAS was dropped
    mapping_rows = mapping_path.read_text().splitlines()[1:]
    assert len(mapping_rows) == 1
