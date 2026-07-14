"""E2: run_manifest.json — the contract artifact the hub indexes.

Conforms to peakatail_contract.models.RunManifest (frozen v0.1.0).
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from ema.outputs import OutputManager

# The frozen contract package lives in the sibling hub repo. If it (and pydantic)
# can be imported, we additionally assert the manifest parses cleanly there.
_CONTRACT_SRC = Path("/home/user/D/peakatail-hub/packages/contract/src")


def _contract_available() -> bool:
    if not _CONTRACT_SRC.exists():
        return False
    import sys
    if str(_CONTRACT_SRC) not in sys.path:
        sys.path.insert(0, str(_CONTRACT_SRC))
    return importlib.util.find_spec("peakatail_contract") is not None and \
        importlib.util.find_spec("pydantic") is not None


def test_manifest_has_run_identity_and_id_grammar(tmp_path: Path) -> None:
    run = tmp_path / "run_TS"
    run.mkdir()
    mgr = OutputManager(base_dir=str(run))
    path = mgr.write_manifest({"directories": {"atlas": "a.bed", "datasets": []}})
    data = json.loads(Path(path).read_text())
    assert data["run_id"] == "run_TS"
    assert data["root"] == str(run.resolve())
    assert data["contract_version"] == "0.1.0"
    assert data["resolved_config"]["directories"]["atlas"] == "a.bed"
    assert data["id_grammar"]["cell_uid"] == "{dataset_id}:{barcode}"
    assert data["entity_counts"]["n_datasets"] == 0


def test_register_artifact_and_datasets(tmp_path: Path) -> None:
    mgr = OutputManager(base_dir=str(tmp_path))
    art = tmp_path / "unified" / "concatenated.mtx"
    art.parent.mkdir(parents=True)
    art.write_text("%%MatrixMarket\n")
    mgr.register_artifact(
        str(art), stage="merge", fmt="mtx", schema_name="counts",
        entity_counts={"n_pas": 5},
    )
    path = mgr.write_manifest(
        {"directories": {"datasets": [{"id": "dsA", "bams": ["a.bam"]}]}},
        stratum_to_label={"CD8_trunc": "CD8 T cells effector memory"},
    )
    data = json.loads(Path(path).read_text())
    reg = next(a for a in data["artifacts"] if a["path"] == "unified/concatenated.mtx")
    assert reg["schema_name"] == "counts"
    assert reg["schema_version"] == "0.1.0"
    assert reg["entity_counts"] == {"n_pas": 5}
    # per-artifact content hash (staleness detection) — sha256 of the file bytes.
    import hashlib
    expect = "sha256:" + hashlib.sha256(art.read_bytes()).hexdigest()
    assert reg["content_hash"] == expect
    assert data["datasets"] == [{"dataset_id": "dsA", "bam_paths": ["a.bam"], "label": None}]
    assert data["stratum_to_label"]["CD8_trunc"].startswith("CD8 T cells")


def test_content_hash_changes_when_artifact_changes(tmp_path: Path) -> None:
    """A changed artifact yields a different content_hash even if the manifest
    shape is identical — the staleness signal the hub indexer needs."""
    def _hash_for(text: str) -> str:
        mgr = OutputManager(base_dir=str(tmp_path))
        art = tmp_path / "provenance" / "pas_ledger.tsv"
        art.parent.mkdir(parents=True, exist_ok=True)
        art.write_text(text)
        mgr.register_artifact(str(art), stage="provenance", fmt="tsv",
                              schema_name="PasLedgerRow")
        return mgr._artifacts[0]["content_hash"]

    h1 = _hash_for("orig_pas_key\na\n")
    h2 = _hash_for("orig_pas_key\na\nb\n")
    assert h1 and h2 and h1 != h2


def test_auto_discovered_artifacts_carry_content_hash(tmp_path: Path) -> None:
    mgr = OutputManager(base_dir=str(tmp_path))
    (tmp_path / "provenance").mkdir()
    (tmp_path / "provenance" / "pas_ledger.tsv").write_text("orig_pas_key\n1\n")
    data = json.loads(Path(mgr.write_manifest({})).read_text())
    led = next(a for a in data["artifacts"] if a["path"] == "provenance/pas_ledger.tsv")
    assert led["content_hash"].startswith("sha256:")


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
    # format values must be contract Format enum members
    assert {a["format"] for a in data["artifacts"]} <= {"bed", "mtx", "tsv", "h5ad", "parquet", "json"}


@pytest.mark.skipif(not _contract_available(), reason="peakatail_contract/pydantic not importable")
def test_manifest_parses_against_frozen_contract(tmp_path: Path) -> None:
    from peakatail_contract.models import RunManifest  # type: ignore

    run = tmp_path / "run_TS"
    run.mkdir()
    (run / "unified").mkdir()
    (run / "unified" / "pas_uid.tsv").write_text("new_pas_id\tpas_uid\n1\tchr1:9:+\n")
    mgr = OutputManager(base_dir=str(run))
    mgr.register_artifact(
        str(run / "unified" / "pas_uid.tsv"),
        stage="unified", fmt="tsv", schema_name="pas_uid",
    )
    path = mgr.write_manifest(
        {"directories": {"datasets": [{"id": "dsA", "bams": ["a.bam"]}]}},
    )
    data = json.loads(Path(path).read_text())
    # Must construct without raising — extras (timestamp, id_grammar) are ignored.
    manifest = RunManifest.model_validate(data)
    assert manifest.run_id == "run_TS"
    assert manifest.datasets[0].dataset_id == "dsA"
    assert any(a.format.value == "tsv" for a in manifest.artifacts)
