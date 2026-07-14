"""B0 regression: run_config.json records the RESOLVED config, not argparse defaults.

Bug B0: ``save_run_config(vars(args))`` serialized the argparse namespace, which
holds defaults for everything supplied via YAML / ``set_directory_config`` — so a
run that snapped to an atlas recorded ``atlas: null``. ``build_resolved_run_config``
reads the resolved ``directory_config`` singleton instead.
"""
from __future__ import annotations

import json
from pathlib import Path

from ema.config import set_directory_config
from ema.outputs import OutputManager, build_resolved_run_config


def test_resolved_config_captures_atlas_and_inputs(tmp_path: Path) -> None:
    atlas = tmp_path / "atlas.bed"
    atlas.write_text("chr1\t100\t101\tA1\t0\t+\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text("# gtf\n")

    set_directory_config(
        output_dir=tmp_path / "run",
        atlas=str(atlas),
        atlas_distance=42,
        gtf_dir=str(gtf),
        datasets=[{"id": "dsA", "bams": ["a.bam"]}],
    )

    resolved = build_resolved_run_config()

    # The atlas path must be present and NOT null (the B0 symptom).
    assert resolved["directories"]["atlas"] == str(atlas)
    assert resolved["directories"]["atlas_distance"] == 42
    assert resolved["directories"]["gtf_dir"] == str(gtf)
    assert resolved["directories"]["datasets"] == [{"id": "dsA", "bams": ["a.bam"]}]
    # Resolved filter/variable blocks exist.
    assert "min_read" in resolved["filters"]
    assert "seqlen" in resolved["variables"]


def test_save_run_config_writes_json_with_resolved_atlas(tmp_path: Path) -> None:
    atlas = tmp_path / "atlas.bed"
    atlas.write_text("chr1\t100\t101\tA1\t0\t+\n")
    out = tmp_path / "run2"
    set_directory_config(output_dir=out, atlas=str(atlas), atlas_distance=7)
    out.mkdir(parents=True, exist_ok=True)

    mgr = OutputManager(base_dir=str(out))
    mgr.save_run_config(build_resolved_run_config())

    data = json.loads((out / "run_config.json").read_text())
    assert "timestamp" in data
    assert data["directories"]["atlas"] == str(atlas)
    assert data["directories"]["atlas_distance"] == 7
