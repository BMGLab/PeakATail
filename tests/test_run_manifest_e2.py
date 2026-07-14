"""E2: run_manifest.json — the contract artifact the hub indexes."""
from __future__ import annotations

import json
from pathlib import Path

from ema.outputs import OutputManager


def test_manifest_has_contract_version_and_id_grammar(tmp_path: Path) -> None:
    mgr = OutputManager(base_dir=str(tmp_path))
    path = mgr.write_manifest({"directories": {"atlas": "a.bed"}})
    data = json.loads(Path(path).read_text())
    assert data["contract_version"] == "0.1.0"
    assert data["resolved_config"]["directories"]["atlas"] == "a.bed"
    assert "pas_uid" in data["id_grammar"]
    assert data["id_grammar"]["cell_uid"] == "{dataset_id}:{barcode}"
    assert data["stratum_to_label"] == {}


def test_register_and_stratum_to_label(tmp_path: Path) -> None:
    mgr = OutputManager(base_dir=str(tmp_path))
    art = tmp_path / "unified" / "concatenated.mtx"
    art.parent.mkdir(parents=True)
    art.write_text("%%MatrixMarket\n")
    mgr.register_artifact(str(art), stage="merge", fmt="mtx", schema="counts@1", n_rows=5)
    path = mgr.write_manifest(
        {"filters": {"min_read": 1500}},
        stratum_to_label={"CD8_T_cells_effector_memory_termi": "CD8 T cells effector memory terminally differentiated"},
    )
    data = json.loads(Path(path).read_text())
    paths = {a["path"] for a in data["artifacts"]}
    assert "unified/concatenated.mtx" in paths
    reg = next(a for a in data["artifacts"] if a["path"] == "unified/concatenated.mtx")
    assert reg["n_rows"] == 5 and reg["schema"] == "counts@1"
    # D10: full label recoverable from truncated stratum dir name.
    assert data["stratum_to_label"]["CD8_T_cells_effector_memory_termi"].startswith("CD8 T cells")


def test_auto_discovers_standard_artifacts(tmp_path: Path) -> None:
    mgr = OutputManager(base_dir=str(tmp_path))
    (tmp_path / "unified").mkdir()
    (tmp_path / "unified" / "pas_uid.tsv").write_text("new_pas_id\tpas_uid\n1\tchr1:9:+\n")
    (tmp_path / "provenance").mkdir()
    (tmp_path / "provenance" / "pas_ledger.tsv").write_text("orig_pas_key\n")
    path = mgr.write_manifest({})
    data = json.loads(Path(path).read_text())
    paths = {a["path"] for a in data["artifacts"]}
    assert "unified/pas_uid.tsv" in paths
    assert "provenance/pas_ledger.tsv" in paths
