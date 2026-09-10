"""D4 regression: fisher counts CELLS not reads (de-pseudoreplication).

Bug D4: the 2x2 table summed READS, so significance scaled with sequencing
depth (reads within a cell are correlated). count_mode='cells' (default) builds
the table from per-cell detection among gene-expressing cells.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ema.switch_test.strategies.fisher import FisherStrategy


def _matrix(cells_per_group: int, depth: int):
    """Two clusters, one gene with 2 PAS. In BOTH clusters every cell expresses
    both PAS at the SAME per-cell rate — so per-cell there is NO differential
    usage. ``depth`` scales reads-per-cell; a reads-based test would grow more
    significant with depth (pseudoreplication), a per-cell test should not.
    """
    rng = np.random.RandomState(0)
    rows = []
    idx = []
    labels = []
    for grp in ("c1", "c2"):
        for i in range(cells_per_group):
            # identical expected usage in both groups: PAS1 ~= PAS2
            p1 = depth
            p2 = depth
            rows.append([p1, p2])
            idx.append(f"{grp}_{i}")
            labels.append(grp)
    cm = pd.DataFrame(rows, index=idx, columns=["1", "2"])
    return cm, pd.Series(labels, index=idx)


def test_cells_mode_not_significant_when_no_percell_difference():
    cm, labels = _matrix(cells_per_group=40, depth=100)
    pas_gene = {"1": "G", "2": "G"}
    df = FisherStrategy().test(
        cm, labels, cluster1="c1", cluster2="c2",
        min_cells_per_group=5, pas_gene_map=pas_gene, count_mode="cells",
    )
    # No per-cell usage difference → not significant.
    assert (df["qvalue"] >= 0.05).all()
    # delta_proportion is per-cell and ~0 here.
    assert df["delta_proportion"].abs().max() < 1e-9


def test_reads_mode_still_available_and_returns_columns():
    cm, labels = _matrix(cells_per_group=40, depth=100)
    pas_gene = {"1": "G", "2": "G"}
    df = FisherStrategy().test(
        cm, labels, cluster1="c1", cluster2="c2",
        min_cells_per_group=5, pas_gene_map=pas_gene, count_mode="reads",
    )
    assert {"pvalue", "qvalue", "delta_proportion", "log2fc"} <= set(df.columns)


def test_invalid_count_mode_raises():
    cm, labels = _matrix(10, 10)
    try:
        FisherStrategy().test(cm, labels, cluster1="c1", cluster2="c2",
                              pas_gene_map={"1": "G", "2": "G"}, count_mode="bogus")
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_cells_mode_detects_real_percell_switch():
    """Cluster c1: all cells use PAS1; cluster c2: all cells use PAS2. A real
    per-cell switch → significant with a large delta_proportion."""
    rows, idx, labels = [], [], []
    for i in range(40):
        rows.append([10, 0]); idx.append(f"c1_{i}"); labels.append("c1")
    for i in range(40):
        rows.append([0, 10]); idx.append(f"c2_{i}"); labels.append("c2")
    cm = pd.DataFrame(rows, index=idx, columns=["1", "2"])
    labels = pd.Series(labels, index=idx)
    df = FisherStrategy().test(
        cm, labels, cluster1="c1", cluster2="c2",
        min_cells_per_group=5, pas_gene_map={"1": "G", "2": "G"}, count_mode="cells",
    )
    assert (df["qvalue"] < 0.05).any()
    assert df["delta_proportion"].abs().max() > 0.9
