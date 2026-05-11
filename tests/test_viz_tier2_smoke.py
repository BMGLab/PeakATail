"""Parametrized smoke tests for Tier 2 viz strategies.

Each strategy receives synthetic data matching the real runner output shapes
and must write at least one non-empty file without raising.
"""
from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from ema.viz import get_viz_strategy, list_viz_strategies


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(42)


@pytest.fixture
def small_adata(rng: np.random.Generator) -> ad.AnnData:
    n_cells, n_pas = 60, 80
    X = sp.csr_matrix(rng.poisson(0.5, (n_cells, n_pas)))
    adata = ad.AnnData(
        X=X,
        obs={"leiden": [str(i % 4) for i in range(n_cells)]},
        obsm={"X_umap": rng.standard_normal((n_cells, 2))},
    )
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"PAS_{i}" for i in range(n_pas)]
    return adata


@pytest.fixture
def volcano_df(rng: np.random.Generator) -> pd.DataFrame:
    """Minimal DataFrame matching run_diff per-pair TSV shape."""
    n = 200
    return pd.DataFrame(
        {
            "pas_id": [f"PAS_{i}" for i in range(n)],
            "log2fc": rng.normal(0, 2, n),
            "qvalue": rng.uniform(0, 1, n),
        }
    )


@pytest.fixture
def diff_agreement_data(rng: np.random.Generator) -> dict[str, set]:
    """dict[strategy_name -> set[sig_pas_ids]]."""
    n_pas = 300
    all_pas = {f"PAS_{i}" for i in range(n_pas)}
    return {
        "fisher": {f"PAS_{i}" for i in rng.choice(n_pas, 60, replace=False)},
        "nbinom": {f"PAS_{i}" for i in rng.choice(n_pas, 50, replace=False)},
        "wilcoxon": {f"PAS_{i}" for i in rng.choice(n_pas, 70, replace=False)},
    }


@pytest.fixture
def pdui_df(small_adata: ad.AnnData, rng: np.random.Generator) -> ad.AnnData:
    """AnnData with .obs['leiden'] + .obs['pdui_score'] for PDUI distribution plot."""
    adata = small_adata.copy()
    adata.obs["pdui_score"] = rng.uniform(0, 1, adata.n_obs)
    return adata


@pytest.fixture
def length_shifts_df(rng: np.random.Generator) -> pd.DataFrame:
    """DataFrame indexed by gene_id with per-cluster-pair columns."""
    n_genes = 80
    return pd.DataFrame(
        rng.normal(0, 0.3, (n_genes, 3)),
        index=[f"GENE_{i}" for i in range(n_genes)],
        columns=["0_vs_1", "0_vs_2", "1_vs_2"],
    )


@pytest.fixture
def match_df(rng: np.random.Generator) -> pd.DataFrame:
    """DataFrame matching ClusterMatchStrategy.match() return shape."""
    datasets = ["sampleA", "sampleB", "sampleC"]
    clusters = ["0", "1", "2", "3"]
    rows = []
    can_id = 1
    for ds in datasets:
        for cl in clusters:
            confidence = float(rng.uniform(0.5, 1.0))
            matched = [
                [other_ds, cl]
                for other_ds in datasets
                if other_ds != ds
            ]
            rows.append(
                {
                    "dataset_id": ds,
                    "original_cluster": cl,
                    "canonical_cluster": can_id,
                    "match_confidence": confidence,
                    "matched_to": json.dumps(matched),
                }
            )
        can_id += 1
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _assert_paths_written(paths: list[Path], strategy_name: str) -> None:
    assert paths, f"{strategy_name} returned no paths"
    for p in paths:
        assert p.exists(), f"{strategy_name}: file not created: {p}"
        assert p.stat().st_size > 50, f"{strategy_name}: file suspiciously small: {p}"


# ---------------------------------------------------------------------------
# Volcano
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", list_viz_strategies(plot_type="volcano"))
def test_volcano_smoke(name: str, volcano_df: pd.DataFrame, tmp_path: Path) -> None:
    strat = get_viz_strategy(name)
    paths = strat.render(volcano_df, tmp_path / f"volcano_{name}")
    _assert_paths_written(paths, name)


# ---------------------------------------------------------------------------
# Diff agreement
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", list_viz_strategies(plot_type="diff_agreement"))
def test_diff_agreement_smoke(
    name: str, diff_agreement_data: dict[str, set], tmp_path: Path
) -> None:
    strat = get_viz_strategy(name)
    paths = strat.render(diff_agreement_data, tmp_path / f"diff_agreement_{name}")
    _assert_paths_written(paths, name)


# ---------------------------------------------------------------------------
# PDUI distribution
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", list_viz_strategies(plot_type="pdui_distribution"))
def test_pdui_distribution_smoke(
    name: str, pdui_df: ad.AnnData, tmp_path: Path
) -> None:
    strat = get_viz_strategy(name)
    paths = strat.render((pdui_df, "pdui_score"), tmp_path / f"pdui_{name}")
    _assert_paths_written(paths, name)


# ---------------------------------------------------------------------------
# Length shifts
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", list_viz_strategies(plot_type="length_shifts"))
def test_length_shifts_smoke(
    name: str, length_shifts_df: pd.DataFrame, tmp_path: Path
) -> None:
    strat = get_viz_strategy(name)
    paths = strat.render(length_shifts_df, tmp_path / f"length_shifts_{name}")
    _assert_paths_written(paths, name)


# ---------------------------------------------------------------------------
# Cluster match sankey
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", list_viz_strategies(plot_type="cluster_match_sankey"))
def test_cluster_match_sankey_smoke(
    name: str, match_df: pd.DataFrame, tmp_path: Path
) -> None:
    strat = get_viz_strategy(name)
    paths = strat.render(match_df, tmp_path / f"sankey_{name}")
    _assert_paths_written(paths, name)


# ---------------------------------------------------------------------------
# Match confidence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", list_viz_strategies(plot_type="match_confidence"))
def test_match_confidence_smoke(
    name: str, match_df: pd.DataFrame, tmp_path: Path
) -> None:
    strat = get_viz_strategy(name)
    paths = strat.render(match_df, tmp_path / f"confidence_{name}")
    _assert_paths_written(paths, name)
