"""E4: ema.data.Run — the manifest-driven, pandas-based run loader.

Builds a synthetic run directory (mirroring the on-disk shapes documented in
``ema/outputs.py::OutputManager`` and ``ema/provenance.py``) under
``tmp_path`` and exercises every accessor, the manifest-vs-conventional-path
fallback in ``_resolve_artifact``, and the parquet-vs-tsv fallback for the E5
long tables (pyarrow is not installed in this environment, so the tsv branch
of ``write_long_table`` is exercised for real, not simulated).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ema.data import Run, RunReadError
from ema.provenance import CELL_LEDGER_COLUMNS, PAS_LEDGER_COLUMNS
from ema.switch_test.long_output import FINDING_LONG_COLUMNS, LENGTH_LONG_COLUMNS, write_long_table

ad = pytest.importorskip("anndata")


# --------------------------------------------------------------------------- #
# fixture builder
# --------------------------------------------------------------------------- #

def _write_tsv(path: Path, columns: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=columns)
    df.to_csv(path, sep="\t", index=False)


def _write_h5ad(path: Path, leiden: list[str], canonical: list[int], n_vars: int = 4) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = len(leiden)
    X = np.zeros((n, n_vars), dtype=np.float32)
    obs = pd.DataFrame({
        "leiden": pd.Categorical(leiden),
        "canonical_cluster": canonical,
    })
    obs.index = [f"cell{i}" for i in range(n)]
    var = pd.DataFrame(index=[f"pas{i}" for i in range(n_vars)])
    ad.AnnData(X=X, obs=obs, var=var).write(path)


def _build_run(tmp_path: Path, *, register_clusters_artifact: bool = False) -> Path:
    root = tmp_path / "run_TS"
    root.mkdir()

    # -- provenance ledgers (conventional root-level, per ema/provenance.py) --
    _write_tsv(
        root / "provenance" / "pas_ledger.tsv",
        PAS_LEDGER_COLUMNS,
        [
            ["p1", "chr1", 100, 110, "+", "u1", "", "ENSG1", 0, "CORE", "pas_gene", "", ""],
            ["p2", "chr2", 200, 260, "-", "u2", "", "ENSG2", 5, "CORE", "pas_gene", "", ""],
            ["p3", "chr3", 300, 310, "+", "u3", "", "", "", "", "atlas_snap", "atlas_snap", "no_gene"],
        ],
    )
    _write_tsv(
        root / "provenance" / "cell_ledger.tsv",
        CELL_LEDGER_COLUMNS,
        [
            ["AAAA", "dsA", 100, 5, "", "", "0"],
            ["CCCC", "dsA", 50, 2, "cb_filter", "low_reads", ""],
            ["GGGG", "dsB", 80, 4, "", "", "1"],
        ],
    )
    # E3 reconcile outputs: flat concatenated by_dataset/{pas,cell}_ledger.tsv
    # (real engine layout — see ema/main.py's reconcile_dataset_ledgers call).
    _write_tsv(
        root / "provenance" / "by_dataset" / "cell_ledger.tsv",
        CELL_LEDGER_COLUMNS,
        [
            ["AAAA", "dsA", 100, 5, "", "", "0"],
            ["GGGG", "dsB", 80, 4, "", "", "1"],
        ],
    )
    (root / "provenance" / "reconcile_summary.json").write_text(json.dumps({
        "datasets": [
            {"ds_id": "dsA", "surviving_pas": 1, "n_vars": 4, "invariant_ok": True},
            {"ds_id": "dsB", "surviving_pas": 1, "n_vars": 2, "invariant_ok": True},
        ],
        "all_ok": True,
        "n_datasets": 2,
    }))

    # -- pasbed (BED6, no header) --
    bed_path = root / "pasbed.bed"
    bed_path.write_text(
        "chr1\t100\t110\tpas1\t0\t+\n"
        "chr2\t200\t260\tpas2\t0\t-\n"
    )

    # -- E5 long tables: written via the real writer so the tsv fallback path
    # (pyarrow not installed here) is exercised exactly as production does. --
    findings_df = pd.DataFrame(
        [["f1", "chr1:109:+", "ENSG1", "0", "1", None, "fisher", "switch_diff:fisher",
          "flat", None, 0.9, 0.5, 0.01, 0.1, 1.0, 100, 40, 60, 40, 30, 10, "significance"]],
        columns=FINDING_LONG_COLUMNS,
    )
    write_long_table(findings_df, str(root / "findings_long"))

    length_df = pd.DataFrame(
        [["fisher", "ENSG1", "T1", "dsA:AAAA", "0", 0.4, "chr1:109:+", 1, "flat", "significance"]],
        columns=LENGTH_LONG_COLUMNS,
    )
    write_long_table(length_df, str(root / "length_long"))

    # -- clusters.h5ad (dsA only, by default — dsB added later for the
    # ambiguity test) --
    _write_h5ad(
        root / "07_clustering" / "dsA" / "clusters.h5ad",
        leiden=["0", "0", "1"],
        canonical=[2, 2, 3],
        n_vars=4,
    )

    artifacts = [
        {
            "path": "provenance/pas_ledger.tsv", "stage": "provenance", "format": "tsv",
            "schema_name": "PasLedgerRow", "schema_version": "0.1.0",
            "entity_counts": {}, "content_hash": None,
        },
    ]
    if register_clusters_artifact:
        artifacts.append({
            "path": "07_clustering/dsA/clusters.h5ad", "stage": "clustering", "format": "h5ad",
            "schema_name": "clusters.h5ad", "schema_version": "0.1.0",
            "entity_counts": {}, "content_hash": None,
        })

    manifest = {
        "run_id": "run_TS",
        "root": str(root),
        "contract_version": "0.1.0",
        "datasets": [
            {"dataset_id": "dsA", "bam_paths": ["a.bam"], "label": None},
            {"dataset_id": "dsB", "bam_paths": ["b.bam"], "label": None},
        ],
        # CellLedgerRow / pasbed.bed / FindingRow / LengthRow deliberately
        # NOT registered here -- exercises the conventional-path fallback.
        "artifacts": artifacts,
        "resolved_config": {
            "directories": {"output_dir": str(root), "atlas": None, "datasets": []},
            "variables": {"seqlen": 10, "cb_len": 16},
            "filters": {"min_read": 5},
            "args": {"strategy": "fisher"},
        },
        "stratum_to_label": {},
        "entity_counts": {"n_datasets": 2},
    }
    (root / "run_manifest.json").write_text(json.dumps(manifest, indent=2))
    return root


# --------------------------------------------------------------------------- #
# construction
# --------------------------------------------------------------------------- #

def test_from_dir_accepts_run_dir_or_manifest_path(tmp_path):
    root = _build_run(tmp_path)
    run_a = Run.from_dir(root)
    run_b = Run.from_dir(root / "run_manifest.json")
    assert run_a.run_id == run_b.run_id == "run_TS"
    assert run_a.root == root.resolve()


def test_missing_manifest_raises_run_read_error(tmp_path):
    empty = tmp_path / "not_a_run"
    empty.mkdir()
    with pytest.raises(RunReadError):
        Run.from_dir(empty)


# --------------------------------------------------------------------------- #
# config / identity
# --------------------------------------------------------------------------- #

def test_config_subdicts(tmp_path):
    run = Run.from_dir(_build_run(tmp_path))
    assert run.variables["seqlen"] == 10
    assert run.filters["min_read"] == 5
    assert run.args["strategy"] == "fisher"
    assert run.directories["output_dir"] == str(run.root)
    assert run.datasets == ["dsA", "dsB"]


def test_config_subdicts_default_to_empty_when_absent(tmp_path):
    # Directly constructed (not via from_dir) with a minimal manifest, to
    # check the "empty dict if absent" fallback independent of file I/O.
    run = Run(root=tmp_path, manifest={"run_id": "bare"})
    assert run.args == {}
    assert run.filters == {}
    assert run.variables == {}
    assert run.directories == {}
    assert run.datasets == []


# --------------------------------------------------------------------------- #
# _resolve_artifact: manifest-first, conventional-fallback
# --------------------------------------------------------------------------- #

def test_resolve_artifact_uses_manifest_entry_when_present(tmp_path):
    root = _build_run(tmp_path)
    run = Run.from_dir(root)
    resolved = run._resolve_artifact("PasLedgerRow", "provenance/pas_ledger.tsv")
    assert resolved == (root / "provenance" / "pas_ledger.tsv").resolve()


def test_resolve_artifact_falls_back_to_conventional_when_manifest_omits(tmp_path):
    root = _build_run(tmp_path)
    run = Run.from_dir(root)
    # CellLedgerRow is NOT registered in the manifest for this fixture.
    resolved = run._resolve_artifact("CellLedgerRow", "provenance/cell_ledger.tsv")
    assert resolved == (root / "provenance" / "cell_ledger.tsv").resolve()


def test_resolve_artifact_raises_when_neither_exists(tmp_path):
    root = _build_run(tmp_path)
    run = Run.from_dir(root)
    with pytest.raises(RunReadError):
        run._resolve_artifact("NoSuchSchema", "nowhere/at/all.tsv")


# --------------------------------------------------------------------------- #
# ledgers / pasbed / long tables
# --------------------------------------------------------------------------- #

def test_pas_ledger_and_cell_ledger_shapes(tmp_path):
    run = Run.from_dir(_build_run(tmp_path))
    pas = run.pas_ledger
    assert list(pas.columns) == PAS_LEDGER_COLUMNS
    assert len(pas) == 3

    cell = run.cell_ledger
    assert list(cell.columns) == CELL_LEDGER_COLUMNS
    assert len(cell) == 3
    assert set(cell["dataset_id"]) == {"dsA", "dsB"}


def test_ledgers_fall_back_to_by_dataset(tmp_path):
    """Real multi-sample no-atlas runs write the reconciled ledgers ONLY under
    provenance/by_dataset/ (no run-level provenance/pas_ledger.tsv, because no
    atlas-snap drop ledger exists) — the loader must resolve them there.

    Regression for a bug caught on a real Laughney run: .pas_ledger raised
    RunReadError because it only knew the run-level path.
    """
    root = tmp_path / "run_byds"
    (root / "provenance" / "by_dataset").mkdir(parents=True)
    # ONLY the by_dataset ledgers exist; no run-level provenance/*.tsv.
    pd.DataFrame([{c: "" for c in PAS_LEDGER_COLUMNS} for _ in range(2)]).to_csv(
        root / "provenance" / "by_dataset" / "pas_ledger.tsv", sep="\t", index=False
    )
    pd.DataFrame([{c: "" for c in CELL_LEDGER_COLUMNS} for _ in range(2)]).to_csv(
        root / "provenance" / "by_dataset" / "cell_ledger.tsv", sep="\t", index=False
    )
    # Minimal manifest with NO ledger artifacts registered (mirrors the real
    # run: auto-discovery is what surfaces them).
    (root / "run_manifest.json").write_text(json.dumps({
        "run_id": "run_byds", "root": str(root), "contract_version": "0.1.0",
        "artifacts": [], "resolved_config": {}, "datasets": [],
    }))

    run = Run.from_dir(root)
    assert list(run.pas_ledger.columns) == PAS_LEDGER_COLUMNS
    assert len(run.pas_ledger) == 2
    assert list(run.cell_ledger.columns) == CELL_LEDGER_COLUMNS
    assert len(run.cell_ledger) == 2


def test_pasbed_columns_and_rows(tmp_path):
    run = Run.from_dir(_build_run(tmp_path))
    bed = run.pasbed
    assert list(bed.columns) == ["chrom", "start", "end", "pas_id", "score", "strand"]
    assert len(bed) == 2
    assert bed.iloc[0]["pas_id"] == "pas1"
    assert bed.iloc[1]["strand"] == "-"


def test_findings_and_length_use_tsv_fallback(tmp_path):
    root = _build_run(tmp_path)
    # `write_long_table` writes .parquet when a parquet engine is installed and
    # falls back to a sibling .tsv otherwise. Assert whichever the environment
    # actually produced — the point of the test is that the loader resolves
    # EITHER extension, so pin the resolver to whichever branch this env took.
    try:
        import pyarrow  # noqa: F401
        _has_parquet = True
    except ImportError:
        try:
            import fastparquet  # noqa: F401
            _has_parquet = True
        except ImportError:
            _has_parquet = False
    if _has_parquet:
        assert (root / "findings_long.parquet").exists()
        assert (root / "length_long.parquet").exists()
    else:
        assert (root / "findings_long.tsv").exists()
        assert not (root / "findings_long.parquet").exists()
        assert (root / "length_long.tsv").exists()

    run = Run.from_dir(root)
    findings = run.findings
    assert list(findings.columns) == FINDING_LONG_COLUMNS
    assert len(findings) == 1
    assert findings.iloc[0]["gene_id"] == "ENSG1"

    length = run.length
    assert list(length.columns) == LENGTH_LONG_COLUMNS
    assert len(length) == 1
    assert length.iloc[0]["gene_id"] == "ENSG1"


# --------------------------------------------------------------------------- #
# genes
# --------------------------------------------------------------------------- #

def test_genes_dedupes_and_drops_empties(tmp_path):
    run = Run.from_dir(_build_run(tmp_path))
    assert run.genes == ["ENSG1", "ENSG2"]


# --------------------------------------------------------------------------- #
# provenance namespace
# --------------------------------------------------------------------------- #

def test_provenance_namespace(tmp_path):
    run = Run.from_dir(_build_run(tmp_path))
    prov = run.provenance
    assert len(prov.pas) == 3
    assert len(prov.cell) == 3
    assert prov.reconcile_summary["all_ok"] is True
    assert prov.reconcile_summary["n_datasets"] == 2
    # by_dataset: split from the flat concatenated by_dataset/cell_ledger.tsv
    # (real engine layout — pas_ledger.tsv has no dataset_id column so it is
    # intentionally not split; see Run.provenance docstring).
    assert set(prov.by_dataset.keys()) == {"dsA", "dsB"}
    assert len(prov.by_dataset["dsA"]) == 1
    assert prov.by_dataset["dsA"].iloc[0]["barcode"] == "AAAA"


def test_provenance_by_dataset_empty_when_dir_absent(tmp_path):
    root = _build_run(tmp_path)
    import shutil
    shutil.rmtree(root / "provenance" / "by_dataset")
    run = Run.from_dir(root)
    assert run.provenance.by_dataset == {}
    assert run.provenance.reconcile_summary is not None  # unaffected


# --------------------------------------------------------------------------- #
# clusters.h5ad
# --------------------------------------------------------------------------- #

def test_clusters_single_dataset_via_conventional_glob(tmp_path):
    root = _build_run(tmp_path)  # only dsA has a clusters.h5ad on disk
    run = Run.from_dir(root)
    adata = run.clusters()
    assert adata.n_vars == 4
    assert run.n_vars() == 4
    # dataset_id-scoped call resolves the same file.
    assert run.clusters("dsA").n_vars == 4


def test_clusters_ambiguous_without_dataset_id_raises(tmp_path):
    root = _build_run(tmp_path)
    # Fabricate a second h5ad for dsB -- now the run has two.
    _write_h5ad(
        root / "07_clustering" / "dsB" / "clusters.h5ad",
        leiden=["0", "1"],
        canonical=[5, 2],
        n_vars=2,
    )
    run = Run.from_dir(root)  # fresh instance: no stale cache from other tests
    with pytest.raises(RunReadError):
        run.clusters()
    # Scoped lookups remain unambiguous.
    assert run.clusters("dsA").n_vars == 4
    assert run.clusters("dsB").n_vars == 2
    assert run.n_vars("dsB") == 2


def test_clusters_via_registered_manifest_artifact(tmp_path):
    root = _build_run(tmp_path, register_clusters_artifact=True)
    run = Run.from_dir(root)
    assert run.clusters("dsA").n_vars == 4


# --------------------------------------------------------------------------- #
# canonical_cluster
# --------------------------------------------------------------------------- #

def test_canonical_cluster_mapped_and_unknown(tmp_path):
    run = Run.from_dir(_build_run(tmp_path))
    assert run.canonical_cluster("dsA", "0") == "2"
    assert run.canonical_cluster("dsA", "1") == "3"
    assert run.canonical_cluster("dsA", "99") is None


def test_canonical_cluster_none_when_h5ad_lacks_columns(tmp_path):
    root = _build_run(tmp_path)
    # Overwrite dsA's h5ad with one that has no canonical_cluster column.
    path = root / "07_clustering" / "dsA" / "clusters.h5ad"
    obs = pd.DataFrame({"leiden": pd.Categorical(["0", "1"])})
    obs.index = ["c0", "c1"]
    ad.AnnData(X=np.zeros((2, 3), dtype=np.float32), obs=obs).write(path)
    run = Run.from_dir(root)
    assert run.canonical_cluster("dsA", "0") is None
