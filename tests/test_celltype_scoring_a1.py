"""A1 core: marker-signature cell-typing + PAS-cluster concordance.

The pure, dependency-light core the eventual `ema celltype` subcommand will
call. GEX I/O + scanpy scoring/clustering that PRODUCE the score matrix are
real-data-gated and out of scope here (see scripts/gex_celltyping.py); this
locks the assignment + concordance math.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from ema.celltype.scoring import assign_celltypes, concordance, zscore_across_cells


def test_zscore_across_cells_centers_and_scales():
    df = pd.DataFrame({"T": [1.0, 2.0, 3.0], "B": [10.0, 10.0, 10.0]})
    z = zscore_across_cells(df)
    assert abs(z["T"].mean()) < 1e-9
    assert abs(z["T"].std(ddof=0) - 1.0) < 1e-9
    # constant column → left at 0 (no divide-by-zero).
    assert (z["B"] == 0.0).all()


def _scores(n_per, layout):
    """Build a cells x signatures score DataFrame + cluster labels.

    layout: list of (cluster, dominant_sig). Each cluster gets n_per cells
    strongly enriched for its dominant signature.
    """
    sigs = ["T", "B", "Myeloid"]
    rows, idx, clusters = [], [], []
    k = 0
    for clu, dom in layout:
        for _ in range(n_per):
            vec = [0.1, 0.1, 0.1]
            vec[sigs.index(dom)] = 5.0
            rows.append(vec)
            idx.append(f"cell{k}"); k += 1
            clusters.append(clu)
    return pd.DataFrame(rows, index=idx, columns=sigs), pd.Series(clusters, index=idx)


def test_assign_celltypes_m1_m2_agree_on_clean_layout():
    df, labels = _scores(20, [("0", "T"), ("1", "B"), ("2", "Myeloid")])
    res = assign_celltypes(df, labels)
    assert res["m1_by_cluster"] == {"0": "T", "1": "B", "2": "Myeloid"}
    assert res["m2_by_cluster"] == {"0": "T", "1": "B", "2": "Myeloid"}
    assert res["agreement"] == 1.0
    assert res["n_clusters"] == 3


def test_assign_celltypes_cell_argmax_present():
    df, labels = _scores(5, [("0", "T"), ("1", "B")])
    res = assign_celltypes(df, labels)
    assert set(res["cell_argmax"].values()) == {"T", "B"}
    assert len(res["cell_argmax"]) == 10


def test_assign_celltypes_disagreement_lowers_agreement():
    # Build a T<->B symmetric dataset so the z-transform is IDENTICAL for both
    # signatures (same global mean/std) — this isolates a genuine mean-vs-
    # majority split within a cluster from any magnitude/scale artifact.
    #   mix : 6 cells (T=2,B=0) + 5 cells (T=0,B=3)
    #   mix2: mirror image (6 cells B=2 + 5 cells T=3) so global T,B match.
    # In "mix": per-cell argmax majority = T (6 vs 5) but cluster-mean tilts to
    # B (mean raw 1.36 > 1.09) → M1 != M2. mix2 disagrees symmetrically.
    sigs = ["T", "B", "Myeloid"]
    rows, idx, clusters = [], [], []
    for i in range(6):
        rows.append([2.0, 0.0, 0.0]); idx.append(f"mix_t{i}"); clusters.append("mix")
    for i in range(5):
        rows.append([0.0, 3.0, 0.0]); idx.append(f"mix_b{i}"); clusters.append("mix")
    for i in range(6):
        rows.append([0.0, 2.0, 0.0]); idx.append(f"m2_b{i}"); clusters.append("mix2")
    for i in range(5):
        rows.append([3.0, 0.0, 0.0]); idx.append(f"m2_t{i}"); clusters.append("mix2")
    df = pd.DataFrame(rows, index=idx, columns=sigs)
    labels = pd.Series(clusters, index=idx)
    res = assign_celltypes(df, labels)
    assert res["n_clusters"] == 2
    assert res["m2_by_cluster"]["mix"] == "T"     # majority of cells
    assert res["m1_by_cluster"]["mix"] == "B"     # cluster-mean argmax
    assert res["agreement"] == 0.0                # both clusters disagree


def test_assign_celltypes_empty_signatures_raises():
    df = pd.DataFrame(index=["a", "b"])
    with pytest.raises(ValueError):
        assign_celltypes(df, ["0", "0"])


def test_assign_celltypes_length_mismatch_raises():
    df, _ = _scores(3, [("0", "T")])
    with pytest.raises(ValueError):
        assign_celltypes(df, ["0", "0"])  # 2 labels, 3 cells


# ---- concordance ----------------------------------------------------------

def test_concordance_identical_labelings_ari_one():
    a = ["x", "x", "y", "y", "z"]
    r = concordance(a, list(a))
    assert r["ARI"] == 1.0 and r["AMI"] == 1.0 and r["n_cells"] == 5


def test_concordance_permuted_labels_still_one():
    # ARI/AMI are invariant to label renaming (co-membership only).
    a = ["T", "T", "B", "B"]
    b = ["0", "0", "1", "1"]
    assert concordance(a, b)["ARI"] == 1.0


def test_concordance_random_is_low():
    rng = np.random.RandomState(0)
    a = rng.randint(0, 4, size=200).tolist()
    b = rng.randint(0, 4, size=200).tolist()
    assert concordance(a, b)["ARI"] < 0.1


def test_concordance_length_mismatch_raises():
    with pytest.raises(ValueError):
        concordance(["a"], ["a", "b"])


def test_concordance_empty_raises():
    with pytest.raises(ValueError):
        concordance([], [])
