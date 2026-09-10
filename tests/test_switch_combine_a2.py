"""A2: `peakatail switch combine` — stitch stage-labelled h5ads for grouped analysis.

Promotes the stage-combining glue (previously only in
scripts/b3_stage_celltype_switch.py) into a package feature: stamp obs[group_key]
with a per-input stage label, optionally split by a cell-type obs column, drop
groups below a min-cell floor, and keep only splits spanning >= 2 groups.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

import anndata as ad

from ema.switch_test.combine import combine_labeled, combine_to_dir, slug


def _adata(n, celltypes, gene_ids=("g1", "g2", "g3")):
    X = np.arange(n * len(gene_ids)).reshape(n, len(gene_ids)).astype(float)
    obs = pd.DataFrame({"celltype": celltypes}, index=[f"c{i}" for i in range(n)])
    var = pd.DataFrame({"gene_id": list(gene_ids)}, index=list(gene_ids))
    return ad.AnnData(X=X, obs=obs, var=var)


# ---- slug -----------------------------------------------------------------

def test_slug_sanitizes():
    assert slug("T cells / CD8+") == "T_cells_CD8"
    assert slug("") == "group"


# ---- combine_labeled (pure) ----------------------------------------------

def test_combine_no_split_single_all_group():
    a = _adata(3, ["x"] * 3)
    b = _adata(2, ["x"] * 2)
    out = combine_labeled([a, b], ["Normal", "Tumor"], group_key="stage")
    assert set(out) == {"__all__"}
    comb = out["__all__"]
    assert comb.n_obs == 5
    assert set(comb.obs["stage"]) == {"Normal", "Tumor"}


def test_combine_drops_single_group_split():
    # Only one stage present → cannot contrast → omitted.
    a = _adata(3, ["x"] * 3)
    out = combine_labeled([a], ["Normal"], group_key="stage")
    assert out == {}


def test_combine_split_by_celltype():
    a = _adata(4, ["T", "T", "B", "B"])
    b = _adata(4, ["T", "T", "B", "B"])
    out = combine_labeled([a, b], ["Normal", "Tumor"],
                          group_key="stage", split_key="celltype")
    assert set(out) == {"T", "B"}
    for ct, comb in out.items():
        assert (comb.obs["celltype"].astype(str) == ct).all()
        assert set(comb.obs["stage"]) == {"Normal", "Tumor"}


def test_combine_min_cells_filters_group():
    # Tumor contributes only 1 cell → below min_cells=2 → dropped → only Normal
    # remains → <2 groups → split omitted.
    a = _adata(4, ["T"] * 4)     # Normal: 4 T cells
    b = _adata(1, ["T"])         # Tumor: 1 T cell
    out = combine_labeled([a, b], ["Normal", "Tumor"],
                          group_key="stage", split_key="celltype", min_cells=2)
    assert out == {}


def test_combine_preserves_gene_id_var():
    a = _adata(2, ["x"] * 2)
    b = _adata(2, ["x"] * 2)
    out = combine_labeled([a, b], ["N", "T"])
    assert "gene_id" in out["__all__"].var.columns


def test_combine_length_mismatch_raises():
    with pytest.raises(ValueError):
        combine_labeled([_adata(2, ["x"] * 2)], ["N", "T"])


def test_combine_empty_raises():
    with pytest.raises(ValueError):
        combine_labeled([], [])


# ---- combine_to_dir + CLI -------------------------------------------------

def test_combine_to_dir_writes_h5ads(tmp_path):
    pa = tmp_path / "normal.h5ad"
    pb = tmp_path / "tumor.h5ad"
    _adata(3, ["T", "T", "B"]).write_h5ad(pa)
    _adata(3, ["T", "B", "B"]).write_h5ad(pb)
    manifest = combine_to_dir(
        [("Normal", str(pa)), ("Tumor", str(pb))], tmp_path / "out",
        group_key="stage", split_key="celltype", min_cells=1,
    )
    assert set(manifest) == {"T", "B"}
    for info in manifest.values():
        assert (tmp_path / "out" / f"{info['slug']}.h5ad").exists()
        assert set(info["groups"]) == {"Normal", "Tumor"}


def test_combine_to_dir_missing_input_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        combine_to_dir([("N", str(tmp_path / "nope.h5ad")),
                        ("T", str(tmp_path / "nope2.h5ad"))], tmp_path / "out")


def test_switch_combine_cli(tmp_path):
    from ema.cli import main
    pa = tmp_path / "normal.h5ad"
    pb = tmp_path / "tumor.h5ad"
    _adata(3, ["T"] * 3).write_h5ad(pa)
    _adata(3, ["T"] * 3).write_h5ad(pb)
    out = tmp_path / "combined"
    res = CliRunner().invoke(main, [
        "switch", "combine",
        "-i", f"Normal={pa}", "-i", f"Tumor={pb}",
        "--split-key", "celltype", "-o", str(out),
    ])
    assert res.exit_code == 0, res.output
    assert (out / "combine_manifest.json").exists()
    assert (out / "T.h5ad").exists()


def test_switch_combine_cli_bad_input_format(tmp_path):
    from ema.cli import main
    res = CliRunner().invoke(main, [
        "switch", "combine", "-i", "noequalssign", "-i", "also=bad_but_ok",
    ])
    assert res.exit_code != 0
