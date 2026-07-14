"""B7 regression: repair NaN var['gene_id'] from the pas_gene map on var_names.

Bug B7: after an across-dataset concat, adata.var['gene_id'] is ~59% NaN (outer
join drops it for PAS absent in some datasets), so joins on the COLUMN drop most
genes. The var_names (pas_id) INDEX is intact, so gene_id is repaired from the
per-dataset pas_gene mapping keyed on var_names.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ema.switch_test.runner import _repair_gene_id, _resolve_pas_gene_map


def test_repair_fills_nan_from_map() -> None:
    s = pd.Series(["ENSG1", np.nan, "", "nan"], index=["1", "2", "3", "4"])
    pas_gene = {"2": "ENSG2", "3": "ENSG3", "4": "ENSG4"}
    repaired, n = _repair_gene_id(s, pas_gene)
    assert n == 3
    assert list(repaired) == ["ENSG1", "ENSG2", "ENSG3", "ENSG4"]


def test_repair_noop_without_map() -> None:
    s = pd.Series([np.nan, "ENSG1"], index=["1", "2"])
    repaired, n = _repair_gene_id(s, {})
    assert n == 0
    assert repaired.isna().iloc[0]


def test_repair_keeps_existing_values() -> None:
    s = pd.Series(["ENSG_keep", np.nan], index=["1", "2"])
    repaired, n = _repair_gene_id(s, {"1": "ENSG_other", "2": "ENSG2"})
    assert n == 1
    assert repaired["1"] == "ENSG_keep"  # existing value untouched
    assert repaired["2"] == "ENSG2"


def test_resolve_pas_gene_map_from_annotatedpas_bed(tmp_path: Path) -> None:
    run = tmp_path / "run"
    h5ad = run / "07_clustering" / "dsA" / "clusters.h5ad"
    h5ad.parent.mkdir(parents=True)
    # annotatedpas.bed a couple dirs up: chrom start end pas_id score strand gene_id
    annot = run / "07_clustering" / "annotatedpas.bed"
    annot.write_text(
        "chr1\t10\t11\t1\t0\t+\tENSG1\n"
        "chr1\t20\t21\t2\t0\t-\tENSG2\n"
    )
    m = _resolve_pas_gene_map(str(h5ad))
    assert m == {"1": "ENSG1", "2": "ENSG2"}


def test_resolve_pas_gene_map_from_pas_gene_tsv(tmp_path: Path) -> None:
    run = tmp_path / "run"
    h5ad = run / "04" / "dsA" / "clusters.h5ad"
    h5ad.parent.mkdir(parents=True)
    pg = run / "04" / "pas_gene.tsv"
    pg.write_text("pas_id\tgene_id\n1\tENSG1\n2\tENSG2\n")
    m = _resolve_pas_gene_map(str(h5ad))
    assert m == {"1": "ENSG1", "2": "ENSG2"}


def test_resolve_returns_empty_when_absent(tmp_path: Path) -> None:
    h5ad = tmp_path / "a" / "b" / "clusters.h5ad"
    h5ad.parent.mkdir(parents=True)
    assert _resolve_pas_gene_map(str(h5ad)) == {}
