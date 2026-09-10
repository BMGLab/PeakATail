"""E1: content-addressed stable pas_uid = chrom:pos:strand.

pas_uid is emitted as a sidecar (pas_uid.tsv) alongside the merged/atlas BED so
cross-run joins can key on a stable id instead of the run-local integer.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from ema.datasets.pas_merge import merge_pas_beds, pas_uid_of


def test_pas_uid_of_is_strand_aware_3prime() -> None:
    # + strand 3' end = end-1 (last base)
    assert pas_uid_of("chr1", 100, 110, "+") == "chr1:109:+"
    # - strand 3' end = start (first base)
    assert pas_uid_of("chr2", 200, 260, "-") == "chr2:200:-"


def test_pas_uid_of_nonint_fallback() -> None:
    assert pas_uid_of("chrX", "a", "b", "+") == "chrX:b:+"


@pytest.mark.skipif(
    not (shutil.which("bedtools") and shutil.which("sort")),
    reason="requires bedtools + sort on PATH",
)
def test_merge_writes_pas_uid_sidecar(tmp_path: Path) -> None:
    bed = tmp_path / "dsA.bed"
    bed.write_text(
        "chr1\t100\t110\t1\t0\t+\n"
        "chr1\t500\t560\t2\t0\t-\n"
    )
    out = tmp_path / "unified"
    merged, mapping = merge_pas_beds(
        [bed], ["dsA"], output_dir=out, gap=10, strands=["+"],
    )
    uid_path = out / "pas_uid.tsv"
    assert uid_path.exists()
    lines = uid_path.read_text().splitlines()
    assert lines[0] == "new_pas_id\tpas_uid"
    uids = dict(l.split("\t") for l in lines[1:])
    # Two merged PAS, each with a content-addressed strand-aware 3' uid.
    assert set(uids.values()) == {"chr1:109:+", "chr1:500:-"}
