"""D2 regression test: atlas snap distance must reflect the PAS 3' summit,
not the raw peak interval.

Before the fix, ``snap_beds_to_atlas`` ran ``bedtools closest -s -d`` on the
full called-peak interval. A wide peak that merely *overlaps* a nearby atlas
entry is reported as distance 0 even though its actual 3' end (the
polyadenylation site) sits far away. This test builds a synthetic wide peak
whose interval overlaps an atlas site but whose 3' summit does not, and
asserts the ``snap_distance_bp`` written to ``atlas_mapping.tsv`` is the
summit distance.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from ema.datasets.atlas_snap import snap_beds_to_atlas

pytestmark = pytest.mark.skipif(
    shutil.which("bedtools") is None or shutil.which("sort") is None,
    reason="bedtools/sort not on PATH",
)


def _write_bed(path: Path, rows: list[tuple[str, int, int, str, str, str]]) -> None:
    with open(path, "w") as f:
        for chrom, start, end, name, score, strand in rows:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")


def test_plus_strand_summit_distance_not_interval_distance(tmp_path: Path) -> None:
    # Atlas site: chr1:1000-1001 (+)
    atlas_bed = tmp_path / "atlas.bed"
    _write_bed(atlas_bed, [("chr1", 1000, 1001, "atlas1", "0", "+")])

    # Wide "+" peak: interval chr1:900-1010 overlaps the atlas site (interval
    # distance == 0), but its 3' summit is end-1 == 1009 (1bp interval
    # [1009, 1010)), which bedtools reports as 9bp from the atlas entry.
    input_bed = tmp_path / "sample.bed"
    _write_bed(input_bed, [("chr1", 900, 1010, "pas1", "0", "+")])

    out_dir = tmp_path / "out"
    _snapped, mapping_path = snap_beds_to_atlas(
        bed_paths=[input_bed],
        dataset_ids=["sample"],
        atlas_bed=atlas_bed,
        output_dir=out_dir,
        distance=50,
        mode="filter",
    )

    rows = _read_mapping(mapping_path)
    assert len(rows) == 1
    dataset_id, old_pasnumber, new_pas_id, snap_distance_bp = rows[0]
    assert dataset_id == "sample"
    assert old_pasnumber == "pas1"
    # The overlapping-interval distance would be 0; the correct summit
    # distance is 9. Assert we get the summit value, not the interval one.
    assert snap_distance_bp != 0
    assert snap_distance_bp == 9


def test_minus_strand_summit_distance_not_interval_distance(tmp_path: Path) -> None:
    # Atlas site: chr1:5000-5001 (-)
    atlas_bed = tmp_path / "atlas.bed"
    _write_bed(atlas_bed, [("chr1", 5000, 5001, "atlas2", "0", "-")])

    # Wide "-" peak: interval chr1:4990-5100 overlaps the atlas site (interval
    # distance == 0), but its 3' summit (on "-", summit == start == 4990,
    # 1bp interval [4990, 4991)) is 10bp from the atlas entry per bedtools.
    input_bed = tmp_path / "sample.bed"
    _write_bed(input_bed, [("chr1", 4990, 5100, "pas2", "0", "-")])

    out_dir = tmp_path / "out"
    _snapped, mapping_path = snap_beds_to_atlas(
        bed_paths=[input_bed],
        dataset_ids=["sample"],
        atlas_bed=atlas_bed,
        output_dir=out_dir,
        distance=50,
        mode="filter",
    )

    rows = _read_mapping(mapping_path)
    assert len(rows) == 1
    dataset_id, old_pasnumber, new_pas_id, snap_distance_bp = rows[0]
    assert dataset_id == "sample"
    assert old_pasnumber == "pas2"
    assert snap_distance_bp != 0
    assert snap_distance_bp == 10


def test_wide_peak_correctly_dropped_when_only_summit_is_far(tmp_path: Path) -> None:
    """A wide peak whose interval overlaps the atlas but whose summit is
    farther than the distance cutoff must be DROPPED — proving the filter
    now acts on the summit distance, not the (misleadingly close) interval.
    """
    atlas_bed = tmp_path / "atlas.bed"
    _write_bed(atlas_bed, [("chr1", 1000, 1001, "atlas1", "0", "+")])

    input_bed = tmp_path / "sample.bed"
    # Same wide peak as above: interval overlaps atlas (dist 0) but summit
    # is 9bp away.
    _write_bed(input_bed, [("chr1", 900, 1010, "pas1", "0", "+")])

    out_dir = tmp_path / "out"
    _snapped, mapping_path = snap_beds_to_atlas(
        bed_paths=[input_bed],
        dataset_ids=["sample"],
        atlas_bed=atlas_bed,
        output_dir=out_dir,
        distance=5,  # summit distance (9) > 5, so this must be dropped
        mode="filter",
    )

    rows = _read_mapping(mapping_path)
    assert rows == []


def _read_mapping(mapping_path: Path) -> list[tuple[str, str, str, int]]:
    """Header-aware read → (dataset_id, old_pasnumber, new_pas_id, distance).

    The mapping gained a ``strand`` column (bug B1), so column positions are
    resolved by header name rather than fixed index.
    """
    rows: list[tuple[str, str, str, int]] = []
    with open(mapping_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        assert {"dataset_id", "old_pasnumber", "new_pas_id", "snap_distance_bp"} <= set(header)
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            rows.append(
                (
                    parts[col["dataset_id"]],
                    parts[col["old_pasnumber"]],
                    parts[col["new_pas_id"]],
                    int(parts[col["snap_distance_bp"]]),
                )
            )
    return rows
