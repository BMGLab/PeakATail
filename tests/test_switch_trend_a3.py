"""A3: `peakatail switch trend` — ordered-stage APA trend aggregation.

Generalizes the Laughney-specific stage-PDUI-trend glue (previously only in
scripts/b3_stage_celltype_switch.py) into a package feature: given a per-cluster
PDUI long table and a caller-supplied stage order, report slope + Spearman
monotonicity + direction, overall and per gene.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from ema.switch_test.trend import (
    build_stage_rank,
    distal_proportion_trend,
    pdui_trend,
    pdui_trend_by_gene,
    trend_from_stage_means,
)


ORDER = ["Normal", "StageI", "StageII", "StageIII"]


# ---- build_stage_rank -----------------------------------------------------

def test_build_stage_rank_orders_from_zero():
    assert build_stage_rank(ORDER) == {"Normal": 0, "StageI": 1, "StageII": 2, "StageIII": 3}


def test_build_stage_rank_ignores_blanks_and_dupes():
    assert build_stage_rank(["A", "", " ", "A", "B"]) == {"A": 0, "B": 1}


def test_build_stage_rank_empty_raises():
    with pytest.raises(ValueError):
        build_stage_rank(["", "  "])


# ---- trend_from_stage_means ----------------------------------------------

def test_trend_decreasing_slope_and_direction():
    r = trend_from_stage_means({0: 0.9, 1: 0.6, 2: 0.3, 3: 0.0})
    assert r["n_stages"] == 4
    assert r["slope"] < 0
    assert r["direction"] == "decreasing"
    assert math.isclose(r["spearman"], -1.0, abs_tol=1e-9)


def test_trend_increasing():
    r = trend_from_stage_means({0: 0.1, 1: 0.4, 2: 0.8})
    assert r["slope"] > 0 and r["direction"] == "increasing"


def test_trend_flat():
    r = trend_from_stage_means({0: 0.5, 1: 0.5, 2: 0.5})
    assert r["direction"] == "flat"
    assert math.isnan(r["spearman"])  # constant → undefined correlation


def test_trend_single_stage_undetermined():
    r = trend_from_stage_means({2: 0.5})
    assert r["n_stages"] == 1 and r["direction"] == "undetermined"
    assert math.isnan(r["slope"])


def test_trend_ignores_nan_means():
    r = trend_from_stage_means({0: 0.9, 1: float("nan"), 2: 0.3})
    assert r["n_stages"] == 2 and r["slope"] < 0


# ---- pdui_trend (overall) -------------------------------------------------

def _long_table():
    # PDUI (distal fraction) falls with stage → 3'UTR shortening. Multiple
    # cells per stage; means are 0.9, 0.6, 0.3, 0.0. Unknown stage is ignored.
    rows = []
    for stage, base in [("Normal", 0.9), ("StageI", 0.6), ("StageII", 0.3), ("StageIII", 0.0)]:
        for d in (-0.05, 0.0, 0.05):
            rows.append({"gene_id": "G1", "cluster": stage, "pdui": base + d})
    rows.append({"gene_id": "G1", "cluster": "UnknownStage", "pdui": 0.99})
    return pd.DataFrame(rows)


def test_pdui_trend_overall_shortening():
    r = pdui_trend(_long_table(), ORDER)
    assert r["direction"] == "decreasing"       # PDUI decreasing == 3'UTR shortening
    assert set(r["mean_by_stage"]) == {"Normal", "StageI", "StageII", "StageIII"}
    assert "UnknownStage" not in r["mean_by_stage"]  # out-of-order stage dropped
    assert r["n_stages"] == 4


def test_pdui_trend_missing_column_raises():
    df = pd.DataFrame({"cluster": ["Normal"], "nope": [0.1]})
    with pytest.raises(ValueError):
        pdui_trend(df, ORDER)


def test_pdui_trend_coerces_non_numeric_values():
    df = pd.DataFrame({"cluster": ["Normal", "StageI", "StageI"],
                       "pdui": ["0.9", "0.3", "bad"]})
    r = pdui_trend(df, ORDER)
    # "bad" dropped; Normal=0.9, StageI=0.3 → decreasing.
    assert r["direction"] == "decreasing" and r["n_stages"] == 2


# ---- pdui_trend_by_gene ---------------------------------------------------

def test_pdui_trend_by_gene_separates_genes():
    rows = []
    # G_short shortens, G_long lengthens.
    for stage, s, l in [("Normal", 0.9, 0.1), ("StageII", 0.5, 0.5), ("StageIII", 0.1, 0.9)]:
        rows.append({"gene_id": "G_short", "cluster": stage, "pdui": s})
        rows.append({"gene_id": "G_long", "cluster": stage, "pdui": l})
    out = pdui_trend_by_gene(pd.DataFrame(rows), ORDER)
    assert out.loc["G_short", "direction"] == "decreasing"
    assert out.loc["G_long", "direction"] == "increasing"
    # sorted by slope ascending → most-shortening gene first.
    assert out.index[0] == "G_short"


def test_pdui_trend_by_gene_requires_gene_col():
    df = pd.DataFrame({"cluster": ["Normal", "StageI"], "pdui": [0.9, 0.3]})
    with pytest.raises(ValueError):
        pdui_trend_by_gene(df, ORDER)


# ---- CLI ------------------------------------------------------------------

def test_switch_trend_cli(tmp_path):
    from ema.cli import main
    pdui = tmp_path / "pdui.tsv"
    _long_table().to_csv(pdui, sep="\t", index=False)
    out = tmp_path / "trendout"
    res = CliRunner().invoke(main, [
        "switch", "trend", "--pdui", str(pdui),
        "--stage-order", "Normal,StageI,StageII,StageIII",
        "-o", str(out),
    ])
    assert res.exit_code == 0, res.output
    assert (out / "length_trend.json").exists()
    assert (out / "length_trend_by_gene.tsv").exists()
    import json
    summary = json.loads((out / "length_trend.json").read_text())
    assert summary["direction"] == "decreasing"


# ---- distal_proportion_trend ----------------------------------------------
#
# The corrected trend for the `proportion` strategy: restrict to multi-PAS
# genes (single-PAS genes have proportion == 1.0 trivially and would flood
# the trend with a constant), take only the distal-PAS (max rank) row per
# (gene, cell), then fit the same slope/spearman as pdui_trend.

def _proportion_long_table():
    """2 genes x 4 stages x few cells; ranks are 1-based (per_gene synthesis
    convention -- see rank_pas_by_genomic_position).

    G_MULTI (2 PAS): distal-PAS proportion falls Normal->StageIII
        (0.9, 0.6, 0.3, 0.0) -- real multi-PAS APA signal.
    G_SINGLE (1 PAS): proportion == 1.0 in every stage/cell (trivial; must
        be excluded by min_pas=2).
    """
    rows = []
    stage_vals = [("Normal", 0.9), ("StageI", 0.6), ("StageII", 0.3), ("StageIII", 0.0)]
    cell_i = 0
    for stage, distal_val in stage_vals:
        for d in (-0.02, 0.0, 0.02):
            cell_i += 1
            cell = f"cell{cell_i}"
            proximal_val = 1.0 - (distal_val + d)
            rows.append({"gene_id": "G_MULTI", "cell": cell, "cluster": stage,
                         "rank": 1, "proportion": proximal_val})
            rows.append({"gene_id": "G_MULTI", "cell": cell, "cluster": stage,
                         "rank": 2, "proportion": distal_val + d})
            rows.append({"gene_id": "G_SINGLE", "cell": cell, "cluster": stage,
                         "rank": 1, "proportion": 1.0})
    return pd.DataFrame(rows)


def test_distal_proportion_trend_excludes_single_pas_genes():
    df = _proportion_long_table()
    r = distal_proportion_trend(df, ORDER, stage_col="cluster")
    # Only G_MULTI's distal (rank=2) rows should contribute.
    assert r["n_genes"] == 1
    assert r["direction"] == "decreasing"
    assert r["n_stages"] == 4


def test_distal_proportion_trend_all_single_pas_is_empty_result():
    df = _proportion_long_table()
    only_single = df[df["gene_id"] == "G_SINGLE"]
    r = distal_proportion_trend(only_single, ORDER, stage_col="cluster")
    assert r["n_genes"] == 0
    assert r["n_stages"] == 0
    assert math.isnan(r["slope"])
    assert math.isnan(r["spearman"])


def test_distal_proportion_trend_picks_max_rank_row():
    """Only the max-rank (distal) row per (gene, cell) should feed the
    trend -- not the proximal row, which moves in the opposite direction."""
    df = _proportion_long_table()
    r = distal_proportion_trend(df, ORDER, stage_col="cluster", min_pas=2)
    # distal values fall 0.9 -> 0.0 so slope must be negative; if the
    # proximal row leaked in (rises 0.1 -> 1.0) the sign would flip / the
    # mean would be wrong.
    assert r["mean_by_stage"]["Normal"] > r["mean_by_stage"]["StageIII"]
    assert abs(r["mean_by_stage"]["Normal"] - 0.9) < 1e-6
    assert abs(r["mean_by_stage"]["StageIII"] - 0.0) < 1e-6


def test_distal_proportion_trend_n_stages_lt_3_spearman_nan():
    """A 2-point trend has a spearman of +/-1 by construction -- meaningless.
    distal_proportion_trend nulls it out (report-analyst convention) while
    still reporting n_stages and a real (non-NaN) slope."""
    df = _proportion_long_table()
    two_stage = df[df["cluster"].isin(["Normal", "StageI"])]
    r = distal_proportion_trend(two_stage, ORDER, stage_col="cluster")
    assert r["n_stages"] == 2
    assert math.isnan(r["spearman"])
    assert r["slope"] < 0


def test_distal_proportion_trend_missing_columns_raises():
    df = pd.DataFrame({"cluster": ["Normal"], "proportion": [0.5]})
    with pytest.raises(ValueError):
        distal_proportion_trend(df, ORDER, stage_col="cluster")


def test_distal_proportion_trend_empty_input():
    df = pd.DataFrame(columns=["gene_id", "rank", "cell", "cluster", "proportion"])
    r = distal_proportion_trend(df, ORDER, stage_col="cluster")
    assert r["n_genes"] == 0 and r["n_stages"] == 0
