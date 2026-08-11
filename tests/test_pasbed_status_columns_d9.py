"""D9: annotatedpas.bed carries atlas_match/atlas_distance_bp/internal_priming
as trailing columns (after gene_id), and ema.data.Run.pasbed reads both the
old (BED6+gene_id) and new (BED6+gene_id+3 status cols) shapes correctly.
"""
from __future__ import annotations

from pathlib import Path

from ema.config import set_directory_config
from ema.outputs import write_pas_gene_artifacts


def _write_pasbed(path: Path, rows: list[tuple]) -> None:
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")


def test_write_pas_gene_artifacts_appends_status_columns(tmp_path):
    run_dir = tmp_path / "run"
    set_directory_config(output_dir=run_dir)
    from ema.config import directory_config

    pasbed = directory_config.pasbed_for("dsA")
    pasbed.parent.mkdir(parents=True, exist_ok=True)
    _write_pasbed(pasbed, [
        ("chr1", 10, 11, "1", 0, "+"),
        ("chr1", 20, 21, "2", 0, "-"),
    ])

    write_pas_gene_artifacts(
        run_dir, "dsA",
        pas_ids=["1", "2"], gene_ids=["ENSG1", "ENSG2"],
        atlas_of={"1": (True, 5), "2": (False, "")},
        ip_of={"1": False, "2": True},
    )

    annot = directory_config.annotatedpas_for("dsA")
    lines = annot.read_text().splitlines()
    assert lines[0].split("\t") == ["chr1", "10", "11", "1", "0", "+", "ENSG1", "True", "5", "False"]
    assert lines[1].split("\t") == ["chr1", "20", "21", "2", "0", "-", "ENSG2", "False", "", "True"]


def test_write_pas_gene_artifacts_defaults_to_empty_status_when_not_supplied(tmp_path):
    """No atlas_of/ip_of passed -> status columns are still appended (additive
    contract) but "" (atlas/ip didn't run for this PAS)."""
    run_dir = tmp_path / "run"
    set_directory_config(output_dir=run_dir)
    from ema.config import directory_config

    pasbed = directory_config.pasbed_for("dsB")
    pasbed.parent.mkdir(parents=True, exist_ok=True)
    _write_pasbed(pasbed, [("chr1", 10, 11, "1", 0, "+")])

    write_pas_gene_artifacts(run_dir, "dsB", pas_ids=["1"], gene_ids=["ENSG1"])

    annot = directory_config.annotatedpas_for("dsB")
    line = annot.read_text().splitlines()[0]
    assert line.split("\t") == ["chr1", "10", "11", "1", "0", "+", "ENSG1", "", "", ""]


def test_run_pasbed_reads_extended_annotatedpas_bed(tmp_path):
    """ema.data.Run.pasbed must not silently corrupt columns when the file
    has MORE than 6 columns (a plain names=_PASBED_COLUMNS read mis-shifts
    everything -- see ema/data/run.py::_read_pasbed_like). Written as
    pasbed.bed (the primary lookup) rather than the annotatedpas.bed
    fallback, since Run._resolve_artifact raises when the primary artifact
    is entirely absent from disk -- a pre-existing behaviour independent of
    this reader change, and out of scope here; this test isolates the
    reader itself.
    """
    from ema.data.run import Run

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "run_manifest.json").write_text("{}")
    (run_dir / "pasbed.bed").write_text(
        "chr1\t10\t11\t1\t0\t+\tENSG1\tTrue\t5\tFalse\n"
        "chr1\t20\t21\t2\t0\t-\tENSG2\tFalse\t\tTrue\n"
    )

    run = Run.from_dir(run_dir)
    df = run.pasbed

    assert list(df.columns) == [
        "chrom", "start", "end", "pas_id", "score", "strand",
        "gene_id", "atlas_match", "atlas_distance_bp", "internal_priming",
    ]
    assert list(df["chrom"]) == ["chr1", "chr1"]
    assert list(df["pas_id"]) == [1, 2]
    assert list(df["gene_id"]) == ["ENSG1", "ENSG2"]
    assert list(df["atlas_match"]) == [True, False]


def test_run_pasbed_reads_legacy_bed6_unchanged(tmp_path):
    """Plain BED6 pasbed.bed (no gene_id / status columns at all) still
    parses exactly as before."""
    from ema.data.run import Run

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "run_manifest.json").write_text("{}")
    (run_dir / "pasbed.bed").write_text("chr1\t10\t11\t1\t0\t+\n")

    run = Run.from_dir(run_dir)
    df = run.pasbed
    assert list(df.columns) == ["chrom", "start", "end", "pas_id", "score", "strand"]
    assert df.iloc[0]["pas_id"] == 1
