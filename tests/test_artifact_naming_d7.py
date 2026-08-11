"""D7 regression: an artifact's path must describe its actual scope.

Bug D7 (the "artifacts that describe something other than their name" class):
a multi-sample run collapsed per-dataset outputs onto a single shared path, so
one dataset's ``pasbed.bed`` silently overwrote the next and downstream tooling
reported one dataset's numbers as if they were the whole run. The fix (with
B2/B5) scopes every per-dataset artifact under ``<peak_calling>/<dataset_id>/``
via ``directory_config.pasbed_for(ds)`` / ``posbed_for(ds)`` / ``negbed_for(ds)``.

This test pins the invariant: given two datasets, ``write_per_dataset_beds``
must write TWO distinct pasbeds at two distinct dataset-scoped paths, each
containing only that dataset's peaks — never a single overwritten file.
"""
from __future__ import annotations

from pathlib import Path

import ema.config as cfg
from ema.outputs import write_per_dataset_beds


def _write_bed(path: Path, rows: list[tuple[str, int, int, int, int, str]]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(f"{c}\t{s}\t{e}\t{i}\t{sc}\t{st}" for c, s, e, i, sc, st in rows) + "\n")
    return str(path)


def test_per_dataset_pasbeds_are_dataset_scoped_not_overwritten(tmp_path: Path) -> None:
    run = tmp_path / "run_TS"
    run.mkdir(parents=True)
    cfg.set_directory_config(output_dir=run)

    # Two datasets, each with a distinct, recognisable peak.
    src = tmp_path / "src"
    ds1_pos = _write_bed(src / "ds1.pos.bed", [("chr1", 100, 101, 11, 0, "+")])
    ds1_neg = _write_bed(src / "ds1.neg.bed", [("chr1", 200, 201, 12, 0, "-")])
    ds2_pos = _write_bed(src / "ds2.pos.bed", [("chr2", 300, 301, 21, 0, "+")])
    ds2_neg = _write_bed(src / "ds2.neg.bed", [("chr2", 400, 401, 22, 0, "-")])

    pasbeds = write_per_dataset_beds(
        run,
        pos_beds_by_ds={"ds1": [ds1_pos], "ds2": [ds2_pos]},
        neg_beds_by_ds={"ds1": [ds1_neg], "ds2": [ds2_neg]},
    )

    # Two DISTINCT dataset-scoped paths — the name carries the dataset id.
    assert set(pasbeds) == {"ds1", "ds2"}
    p1, p2 = Path(pasbeds["ds1"]), Path(pasbeds["ds2"])
    assert p1 != p2
    assert "ds1" in p1.parts and "ds2" in p2.parts
    assert p1.exists() and p2.exists()

    # Each pasbed holds only its own dataset's peaks (no cross-contamination,
    # no overwrite): ds1 -> chr1 ids 11/12, ds2 -> chr2 ids 21/22.
    t1, t2 = p1.read_text(), p2.read_text()
    assert "chr1" in t1 and "chr2" not in t1
    assert "chr2" in t2 and "chr1" not in t2
    assert "\t11\t" in t1 and "\t12\t" in t1
    assert "\t21\t" in t2 and "\t22\t" in t2


def test_directory_config_paths_are_per_dataset(tmp_path: Path) -> None:
    """The path scheme itself must differ by dataset id (name == scope)."""
    cfg.set_directory_config(output_dir=tmp_path / "run")
    assert cfg.directory_config.pasbed_for("A") != cfg.directory_config.pasbed_for("B")
    assert "A" in cfg.directory_config.pasbed_for("A").parts
    assert "B" in cfg.directory_config.pasbed_for("B").parts
