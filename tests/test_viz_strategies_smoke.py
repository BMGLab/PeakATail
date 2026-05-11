"""Each registered strategy renders without error on synthetic data + writes expected files."""
from pathlib import Path
import numpy as np
import anndata as ad
import scipy.sparse as sp
import pytest

from ema.viz import list_viz_strategies, get_viz_strategy


@pytest.fixture
def small_adata():
    rng = np.random.default_rng(0)
    n_cells, n_pas = 50, 100
    X = sp.csr_matrix(rng.poisson(0.5, (n_cells, n_pas)))
    obs_idx = [f"cell_{i}" for i in range(n_cells)]
    var_idx = [f"PAS_{i}" for i in range(n_pas)]
    adata = ad.AnnData(X=X, obs={"leiden": [str(i % 3) for i in range(n_cells)]},
                       var={}, obsm={"X_umap": rng.standard_normal((n_cells, 2))})
    adata.obs_names = obs_idx
    adata.var_names = var_idx
    return adata


@pytest.mark.parametrize("name", [n for n in list_viz_strategies(plot_type="umap")])
def test_umap_strategies_smoke(name, small_adata, tmp_path):
    strat = get_viz_strategy(name)
    paths = strat.render((small_adata, "smoke"), tmp_path / "umap_smoke")
    assert paths, f"{name} returned no paths"
    for p in paths:
        assert p.exists() and p.stat().st_size > 100, f"{p} not written or empty"


@pytest.mark.parametrize("name", [n for n in list_viz_strategies(plot_type="cluster_sizes")])
def test_cluster_sizes_smoke(name, small_adata, tmp_path):
    strat = get_viz_strategy(name)
    paths = strat.render((small_adata, "smoke"), tmp_path / "cluster_sizes_smoke")
    assert paths
    for p in paths:
        assert p.exists() and p.stat().st_size > 100


def test_peak_qc_smoke(tmp_path):
    """peak_qc accepts a dict of arrays (peak counts per chrom, etc.)."""
    rng = np.random.default_rng(0)
    data = {
        "peaks_per_chrom": {"chr1": 1000, "chr2": 800, "chr3": 600, "chrX": 400},
        "per_cell_pas": rng.integers(50, 200, 100),
        "peak_widths": rng.integers(50, 500, 1000),
        "per_cell_reads": rng.integers(2000, 10000, 100),
    }
    for name in list_viz_strategies(plot_type="peak_qc"):
        strat = get_viz_strategy(name)
        paths = strat.render(data, tmp_path / f"peak_qc_{name}")
        assert paths
        for p in paths:
            assert p.exists()
