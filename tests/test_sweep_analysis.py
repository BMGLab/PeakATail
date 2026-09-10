"""Unit tests for the sweep-level cumulative analysis layer.

All tests build small synthetic sweep directory trees under ``tmp_path`` --
no real sweep data is required. Focused on the aggregation helpers in
``ema.benchmark.sweep_analysis``, not on plotting internals (figures are
smoke-tested only: they must not raise and must produce a file).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ema.benchmark import sweep_analysis as sa


# ---------------------------------------------------------------------------
# Fixtures / builders
# ---------------------------------------------------------------------------


def _write_bed(path: Path, n_lines: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for i in range(n_lines):
            f.write(f"1\t{100 * i}\t{100 * i + 50}\tpas{i}\tGENE{i}\t+\n")


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def _bench_json(n_predicted: int, n_reference: int = 1000) -> dict:
    return {
        "n_predicted": n_predicted,
        "n_reference": n_reference,
        "cutoffs": {
            "50": {"precision": 0.90, "recall": 0.10, "f1": 0.18},
            "100": {"precision": 0.92, "recall": 0.12, "f1": 0.21},
            "500": {"precision": 0.95, "recall": 0.20, "f1": 0.33},
            "1000": {"precision": 0.97, "recall": 0.30, "f1": 0.46},
        },
    }


def _make_h5ad(path: Path, n_cells: int, n_clusters: int) -> None:
    ad = pytest.importorskip("anndata")
    import scipy.sparse as sp

    rng = np.random.default_rng(0)
    n_vars = 10
    X = sp.csr_matrix(rng.integers(0, 5, size=(n_cells, n_vars)))
    obs = pd.DataFrame(
        {"leiden": [str(i % n_clusters) for i in range(n_cells)]},
        index=[f"cell{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"pas{i}" for i in range(n_vars)])
    a = ad.AnnData(X=X, obs=obs, var=var)
    path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(path)


def _build_grid_run(root: Path, run_id: str, n_pas: int, n_predicted: int,
                     atlas_match_rate: float = 0.8) -> None:
    d = root / "runs" / "grid" / run_id
    _write_bed(d / "pasbed.bed", n_pas)
    _write_json(d / "benchmark_vs_polyasite_v3.json", _bench_json(n_predicted))
    _write_json(d / "unified" / "atlas_stats.json", {
        "mode": "annotate",
        "n_atlas_matched": int(n_pas * atlas_match_rate),
        "n_atlas_unmatched": int(n_pas * (1 - atlas_match_rate)),
        "atlas_match_rate": atlas_match_rate,
    })


def _build_reannotate_branch(root: Path, branch: str, trim: dict | None,
                              clustering: dict | None, n_pas: int,
                              datasets: dict[str, tuple[int, int]]) -> None:
    """datasets: {dataset_id: (n_cells, n_clusters)}"""
    d = root / "runs" / "reannotate" / branch
    _write_bed(d / "pasbed.bed", n_pas)
    manifest = {}
    if trim is not None:
        manifest["trim"] = trim
    if clustering is not None:
        manifest["clustering"] = clustering
    _write_json(d / "branch_manifest.json", manifest)
    for ds, (n_cells, n_clusters) in datasets.items():
        _make_h5ad(d / "07_clustering" / ds / "clusters.h5ad", n_cells, n_clusters)


def _build_switch_celltype(root: Path, celltype: str, slope: float, spearman: float,
                            direction: str, mean_by_stage: dict, contrasts: dict[str, list[float]],
                            gene_trends: list[tuple[str, float, float, str]]) -> None:
    """contrasts: {"A_vs_B": [qvalue, qvalue, ...]}; gene_trends: [(gene_id, slope, spearman, direction)]"""
    switch_dir = root / "runs" / "B1_cohort_full" / "B3_switch"
    _write_json(switch_dir / "trend" / celltype / "length_trend.json", {
        "n_stages": len(mean_by_stage),
        "slope": slope,
        "spearman": spearman,
        "direction": direction,
        "mean_by_stage": mean_by_stage,
        "value_col": "pdui",
    })

    gene_path = switch_dir / "trend" / celltype / "length_trend_by_gene.tsv"
    gene_path.parent.mkdir(parents=True, exist_ok=True)
    gdf = pd.DataFrame(gene_trends, columns=["gene_id", "slope", "spearman", "direction"])
    gdf.insert(1, "n_stages", 3)
    gdf.to_csv(gene_path, sep="\t", index=False)

    fisher_dir = switch_dir / "diff" / celltype / "fisher" / "differential"
    fisher_dir.mkdir(parents=True, exist_ok=True)
    for contrast, qvalues in contrasts.items():
        df = pd.DataFrame({
            "pas_id": range(len(qvalues)),
            "gene_id": [f"G{i}" for i in range(len(qvalues))],
            "qvalue": qvalues,
        })
        df.to_csv(fisher_dir / f"fisher_{contrast}.tsv", sep="\t", index=False)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def test_count_bed_lines(tmp_path):
    bed = tmp_path / "x.bed"
    _write_bed(bed, 7)
    assert sa.count_bed_lines(bed) == 7


def test_count_bed_lines_missing(tmp_path):
    assert sa.count_bed_lines(tmp_path / "nope.bed") is None


def test_fmt_handles_none_and_nan():
    assert sa._fmt(None) == "n/a"
    assert sa._fmt(float("nan")) == "n/a"
    assert sa._fmt(1234, ",") == "1,234"
    assert sa._fmt(0.5, ".1%") == "50.0%"


# ---------------------------------------------------------------------------
# 1. Strategy comparison
# ---------------------------------------------------------------------------


def test_harvest_strategy_comparison_ranks_by_widest_f1(tmp_path):
    _build_grid_run(tmp_path, "lg_annotate", n_pas=100, n_predicted=100)
    _build_grid_run(tmp_path, "si_annotate", n_pas=200, n_predicted=200)
    # si_annotate gets a higher f1@1000 by overriding its bench json.
    _write_json(
        tmp_path / "runs" / "grid" / "si_annotate" / "benchmark_vs_polyasite_v3.json",
        {**_bench_json(200), "cutoffs": {**_bench_json(200)["cutoffs"],
                                          "1000": {"precision": 0.99, "recall": 0.99, "f1": 0.99}}},
    )

    df = sa.harvest_strategy_comparison(tmp_path)
    assert set(df["run"]) == {"lg_annotate", "si_annotate"}
    assert df.iloc[0]["run"] == "si_annotate"  # ranked first: higher f1@1000
    assert df.iloc[0]["rank"] == 1
    assert df.loc[df["run"] == "lg_annotate", "n_pas"].iloc[0] == 100
    assert df.loc[df["run"] == "lg_annotate", "atlas_match_rate"].iloc[0] == 0.8


def test_harvest_strategy_comparison_missing_run_skipped(tmp_path):
    _build_grid_run(tmp_path, "lg_annotate", n_pas=50, n_predicted=50)
    df = sa.harvest_strategy_comparison(tmp_path)
    assert list(df["run"]) == ["lg_annotate"]


def test_harvest_strategy_comparison_empty_sweep(tmp_path):
    df = sa.harvest_strategy_comparison(tmp_path)
    assert df.empty


# ---------------------------------------------------------------------------
# 2. Internal-priming axis
# ---------------------------------------------------------------------------


def test_harvest_ip_filter_axis_computes_drop_fraction():
    strategy_df = pd.DataFrame([
        {"run": "lg_ip_off", "peak_strategy": "lambda_gradient", "ip_mode": "off", "n_pas": 1000},
        {"run": "lg_ip_filter", "peak_strategy": "lambda_gradient", "ip_mode": "filter", "n_pas": 900},
        {"run": "lg_annotate", "peak_strategy": "lambda_gradient", "ip_mode": "annotate", "n_pas": 1000},
        {"run": "lp_annotate", "peak_strategy": "lambda_poisson", "ip_mode": "annotate", "n_pas": 800},
    ])
    ip_df = sa.harvest_ip_filter_axis(strategy_df)
    # Only lambda_gradient rows participate in this axis.
    assert set(ip_df["run"]) == {"lg_ip_off", "lg_ip_filter", "lg_annotate"}
    filt_row = ip_df[ip_df["ip_mode"] == "filter"].iloc[0]
    assert filt_row["n_pas_dropped_vs_off"] == 100
    assert filt_row["pct_flagged_internal_priming"] == pytest.approx(10.0)
    off_row = ip_df[ip_df["ip_mode"] == "off"].iloc[0]
    assert off_row["n_pas_dropped_vs_off"] == 0


def test_harvest_ip_filter_axis_no_off_run_returns_unmodified():
    strategy_df = pd.DataFrame([
        {"run": "lg_annotate", "peak_strategy": "lambda_gradient", "ip_mode": "annotate", "n_pas": 1000},
    ])
    ip_df = sa.harvest_ip_filter_axis(strategy_df)
    assert "n_pas_dropped_vs_off" not in ip_df.columns


def test_harvest_ip_filter_axis_empty_input():
    assert sa.harvest_ip_filter_axis(pd.DataFrame()).empty


# ---------------------------------------------------------------------------
# 3. Trim / clustering comparison
# ---------------------------------------------------------------------------


def test_harvest_trim_comparison_reads_manifest_and_h5ad(tmp_path):
    _build_reannotate_branch(
        tmp_path, "A2_trim_default",
        trim={"max_gene_distance": 5000, "utr_multiplier": 2.0, "include_extended": False},
        clustering={"method": "leiden_tfidf", "resolution": 1.0, "n_neighbors": 30},
        n_pas=500,
        datasets={"ds1": (30, 3), "ds2": (20, 2)},
    )
    df = sa.harvest_trim_comparison(tmp_path)
    assert len(df) == 1
    row = df.iloc[0]
    assert row["branch"] == "A2_trim_default"
    assert row["max_gene_distance"] == 5000
    assert row["n_pas"] == 500
    assert row["n_datasets"] == 2
    assert row["n_cells_total"] == 50
    assert row["n_clusters_mean"] == pytest.approx(2.5)


def test_harvest_trim_comparison_falls_back_when_manifest_missing(tmp_path):
    d = tmp_path / "runs" / "reannotate" / "A2_trim_d1000_ext"
    _write_bed(d / "pasbed.bed", 42)
    # no branch_manifest.json written -> falls back to TRIM_META
    df = sa.harvest_trim_comparison(tmp_path)
    assert df.iloc[0]["max_gene_distance"] == 1000
    assert bool(df.iloc[0]["include_extended"]) is True


def test_harvest_clustering_comparison(tmp_path):
    _build_reannotate_branch(
        tmp_path, "A3_res0.5", trim=None,
        clustering={"method": "leiden_tfidf", "resolution": 0.5, "n_neighbors": 30},
        n_pas=100, datasets={"ds1": (10, 4)},
    )
    df = sa.harvest_clustering_comparison(tmp_path)
    assert len(df) == 1
    assert df.iloc[0]["resolution"] == 0.5
    assert df.iloc[0]["n_clusters_mean"] == 4


def test_gex_concordance_summary(tmp_path):
    conc_path = tmp_path / "runs" / "B1_cohort_full" / "B2_gex_celltyping" / "concordance.csv"
    conc_path.parent.mkdir(parents=True)
    pd.DataFrame({
        "gsm": ["GSM1", "GSM2", "GSM3"],
        "ARI_celltype_vs_pas": [0.5, 0.6, 0.7],
        "AMI_celltype_vs_pas": [0.4, 0.5, 0.6],
    }).to_csv(conc_path, index=False)

    df = sa.harvest_gex_concordance(tmp_path)
    assert len(df) == 3
    summary = sa.summarize_gex_concordance(df)
    assert summary["n_gsm"] == 3
    assert summary["ARI_celltype_vs_pas"]["mean"] == pytest.approx(0.6)


def test_gex_concordance_missing_file(tmp_path):
    assert sa.harvest_gex_concordance(tmp_path).empty
    assert sa.summarize_gex_concordance(pd.DataFrame()) == {}


# ---------------------------------------------------------------------------
# 4. Switch cumulative
# ---------------------------------------------------------------------------


def _default_sweep_with_switch(tmp_path) -> Path:
    _build_switch_celltype(
        tmp_path, "CT_SHORTEN",
        slope=-0.01, spearman=-1.0, direction="decreasing",
        mean_by_stage={"Normal": 0.20, "StageI": 0.15, "IVprimary": 0.10, "Met": 0.05},
        contrasts={"IVprimary_vs_Normal": [0.001, 0.2, 0.5], "Met_vs_Normal": [0.001, 0.001]},
        gene_trends=[("GENEA", -0.3, -0.9, "decreasing"), ("GENEB", -0.1, -0.5, "decreasing")],
    )
    _build_switch_celltype(
        tmp_path, "CT_LENGTHEN",
        slope=0.02, spearman=1.0, direction="increasing",
        mean_by_stage={"Normal": 0.05, "StageI": 0.10, "IVprimary": 0.15, "Met": 0.20},
        contrasts={"IVprimary_vs_Normal": [0.001, 0.001, 0.5]},
        gene_trends=[("GENEA", 0.25, 0.85, "increasing")],
    )
    _build_switch_celltype(
        tmp_path, "CT_SHORTEN2",
        slope=-0.005, spearman=-0.9, direction="decreasing",
        mean_by_stage={"Normal": 0.30, "Met": 0.10},
        contrasts={"Met_vs_Normal": [0.001, 0.9]},
        gene_trends=[("GENEA", -0.28, -0.95, "decreasing")],
    )
    return tmp_path


def test_harvest_switch_trends_and_shortening_summary(tmp_path):
    root = _default_sweep_with_switch(tmp_path)
    trend_df = sa.harvest_switch_trends(root)
    assert set(trend_df["celltype"]) == {"CT_SHORTEN", "CT_LENGTHEN", "CT_SHORTEN2"}
    assert list(trend_df["celltype"])[0] in ("CT_SHORTEN", "CT_SHORTEN2")  # sorted by slope ascending

    summary = sa.summarize_shortening_lengthening(trend_df)
    assert summary["n_celltypes"] == 3
    assert summary["n_shortening"] == 2
    assert summary["n_lengthening"] == 1
    assert summary["n_flat_or_other"] == 0


def test_summarize_shortening_lengthening_empty():
    summary = sa.summarize_shortening_lengthening(pd.DataFrame())
    assert summary == {"n_celltypes": 0, "n_shortening": 0, "n_lengthening": 0, "n_flat_or_other": 0}


def test_celltype_stage_pdui_matrix(tmp_path):
    root = _default_sweep_with_switch(tmp_path)
    trend_df = sa.harvest_switch_trends(root)
    matrix = sa.celltype_stage_pdui_matrix(trend_df)
    assert "Normal" in matrix.columns and "Met" in matrix.columns
    assert matrix.loc["CT_SHORTEN", "Normal"] == pytest.approx(0.20)
    assert pd.isna(matrix.loc["CT_SHORTEN2", "StageI"])


def test_harvest_fisher_hit_counts_and_contrast_summary(tmp_path):
    root = _default_sweep_with_switch(tmp_path)
    fisher_df = sa.harvest_fisher_hit_counts(root, fdr=0.05)
    row = fisher_df[(fisher_df["celltype"] == "CT_SHORTEN") & (fisher_df["contrast"] == "IVprimary_vs_Normal")].iloc[0]
    assert row["n_tests"] == 3
    assert row["n_sig"] == 1  # only qvalue=0.001 < 0.05

    contrast_summary = sa.summarize_fisher_hits(fisher_df)
    met_vs_normal = contrast_summary[contrast_summary["contrast"] == "Met_vs_Normal"].iloc[0]
    # CT_SHORTEN (2 sig of [0.001, 0.001]) + CT_SHORTEN2 (1 sig of [0.001, 0.9]) = 3 sig, 2 celltypes
    assert met_vs_normal["n_celltypes"] == 2
    assert met_vs_normal["total_sig"] == 3
    assert met_vs_normal["total_tests"] == 4


def test_harvest_fisher_hit_counts_no_switch_data(tmp_path):
    # grid/reannotate runs have no B3_switch at all -- must not raise.
    _build_grid_run(tmp_path, "lg_annotate", n_pas=10, n_predicted=10)
    df = sa.harvest_fisher_hit_counts(tmp_path)
    assert df.empty
    assert sa.summarize_fisher_hits(df).empty


def test_harvest_recurrent_switch_genes_min_celltypes(tmp_path):
    root = _default_sweep_with_switch(tmp_path)
    # GENEA trends (|spearman|>=0.8) in all 3 celltypes; GENEB only in 1.
    recurrent = sa.harvest_recurrent_switch_genes(root, min_celltypes=3, spearman_thresh=0.8)
    assert list(recurrent["gene_id"]) == ["GENEA"]
    assert recurrent.iloc[0]["n_celltypes_trending"] == 3
    assert recurrent.iloc[0]["dominant_direction"] == "decreasing"
    assert recurrent.iloc[0]["n_consistent"] == 2  # 2/3 calls were "decreasing"


def test_harvest_recurrent_switch_genes_threshold_excludes_weak_trends(tmp_path):
    root = _default_sweep_with_switch(tmp_path)
    recurrent = sa.harvest_recurrent_switch_genes(root, min_celltypes=2, spearman_thresh=0.99)
    # only spearman==-1.0 / 1.0 entries qualify -> GENEA appears in just 2 celltypes now
    assert set(recurrent["gene_id"]) <= {"GENEA"}


def test_harvest_recurrent_switch_genes_no_data_returns_empty_with_columns(tmp_path):
    df = sa.harvest_recurrent_switch_genes(tmp_path)
    assert df.empty
    assert list(df.columns) == [
        "gene_id", "n_celltypes_trending", "dominant_direction",
        "n_consistent", "mean_slope", "mean_abs_spearman",
    ]


# ---------------------------------------------------------------------------
# Figures (smoke only)
# ---------------------------------------------------------------------------


def test_plot_strategy_precision_recall_writes_png_and_svg(tmp_path):
    df = pd.DataFrame([
        {"run": "a", "n_pas": 100, "recall_1000": 0.3},
        {"run": "b", "n_pas": 200, "recall_1000": 0.5},
    ])
    out = tmp_path / "fig" / "strategy"
    sa.plot_strategy_precision_recall(df, out)
    assert out.with_suffix(".png").exists()
    assert out.with_suffix(".svg").exists()


def test_plot_functions_skip_gracefully_on_empty_input(tmp_path):
    out = tmp_path / "fig" / "x"
    # None of these should raise on empty input.
    sa.plot_strategy_precision_recall(pd.DataFrame(), out)
    sa.plot_trim_cluster_sensitivity(pd.DataFrame(), pd.DataFrame(), out)
    sa.plot_celltype_stage_heatmap(pd.DataFrame(), out)
    sa.plot_shortening_lengthening_summary(pd.DataFrame(), out)
    sa.plot_recurrent_genes(pd.DataFrame(), out)
    assert not out.with_suffix(".png").exists()


# ---------------------------------------------------------------------------
# End-to-end
# ---------------------------------------------------------------------------


def test_analyze_sweep_end_to_end(tmp_path):
    root = tmp_path / "sweep"
    _build_grid_run(root, "lg_annotate", n_pas=100, n_predicted=100)
    _build_grid_run(root, "lg_ip_filter", n_pas=90, n_predicted=90)
    _build_grid_run(root, "lg_ip_off", n_pas=100, n_predicted=100)
    _build_reannotate_branch(
        root, "A2_trim_default",
        trim={"max_gene_distance": 5000, "utr_multiplier": 2.0, "include_extended": False},
        clustering={"method": "leiden_tfidf", "resolution": 1.0, "n_neighbors": 30},
        n_pas=500, datasets={"ds1": (30, 3)},
    )
    _default_sweep_with_switch(root)

    out_dir = tmp_path / "out"
    summary = sa.analyze_sweep(root, out_dir, min_recurrent_celltypes=2, spearman_thresh=0.8)

    assert (out_dir / "cumulative_summary.json").exists()
    assert (out_dir / "cumulative_summary.md").exists()
    assert (out_dir / "tables" / "strategy_comparison.csv").exists()
    assert (out_dir / "tables" / "switch_trend_by_celltype.csv").exists()
    assert (out_dir / "tables" / "recurrent_switch_genes.csv").exists()
    assert (out_dir / "figures" / "strategy_precision_recall.png").exists()
    assert (out_dir / "figures" / "shortening_lengthening_summary.svg").exists()

    assert summary["switch"]["shortening_vs_lengthening"]["n_shortening"] == 2
    md_text = (out_dir / "cumulative_summary.md").read_text()
    assert "Cumulative Sweep Analysis" in md_text
    assert "lg_annotate" in md_text

    # Reloadable JSON with no leftover numpy types.
    reloaded = json.loads((out_dir / "cumulative_summary.json").read_text())
    assert reloaded["switch"]["n_celltypes_analyzed"] == 3


# ---------------------------------------------------------------------------
# FILTER-EFFECT scenario comparison (independent named run dirs, not a
# sweep_root grid/reannotate tier)
# ---------------------------------------------------------------------------


def _build_scenario_run(run_dir: Path, n_pas: int, atlas_match_rate: float | None = None,
                         peak_filters: dict | None = None) -> None:
    """A minimal `peakatail run` output dir: annotatedpas.bed + optional filter stats."""
    _write_bed(run_dir / "annotatedpas.bed", n_pas)
    if atlas_match_rate is not None:
        _write_json(run_dir / "atlas_stats.json", {
            "atlas_match_rate": atlas_match_rate,
            "n_atlas_matched": int(n_pas * atlas_match_rate),
            "n_atlas_unmatched": int(n_pas * (1 - atlas_match_rate)),
        })
    if peak_filters is not None:
        _write_json(run_dir / "04_pas_gene_assignment" / "peak_filters_stats.json", peak_filters)


def _build_scenario_switch(run_dir: Path, celltype: str, slope: float, direction: str,
                            contrasts: dict) -> None:
    """Switch trend + fisher outputs directly under <run_dir>/B3_switch (no B1_cohort_full)."""
    switch_dir = run_dir / "B3_switch"
    _write_json(switch_dir / "trend" / celltype / "length_trend.json", {
        "n_stages": 4, "slope": slope, "spearman": -0.9, "direction": direction,
        "mean_by_stage": {"Normal": 0.5, "StageI": 0.45, "IVprimary": 0.4, "Met": 0.35},
    })
    fisher_dir = switch_dir / "diff" / celltype / "fisher" / "differential"
    fisher_dir.mkdir(parents=True, exist_ok=True)
    for contrast, qvalues in contrasts.items():
        pd.DataFrame({"qvalue": qvalues}).to_csv(
            fisher_dir / f"fisher_{contrast}.tsv", sep="\t", index=False
        )


def _default_scenario_dirs(root: Path) -> dict:
    scenarios = {
        "baseline": dict(n_pas=1000, atlas_match_rate=0.80, slope=-0.010, direction="decreasing"),
        "atlas_filter": dict(n_pas=950, atlas_match_rate=0.82, slope=-0.012, direction="decreasing"),
        "annot_filter_3utr": dict(n_pas=900, atlas_match_rate=0.79, slope=-0.008, direction="decreasing"),
        "ip_filter": dict(n_pas=850, atlas_match_rate=0.78, slope=0.004, direction="increasing"),
    }
    scenario_dirs = {}
    for name, cfg in scenarios.items():
        d = root / name
        _build_scenario_run(
            d, n_pas=cfg["n_pas"], atlas_match_rate=cfg["atlas_match_rate"],
            peak_filters={
                "ip_flag_rate": 0.08, "n_ip_flagged": 40,
                "pos": {"total": 500, "filtered": 10,
                        "annotation": {"total": 500, "filtered": 5},
                        "internal_priming": {"filtered": 10}},
                "neg": {"total": 500, "filtered": 8,
                        "annotation": {"total": 500, "filtered": 4},
                        "internal_priming": {"filtered": 8}},
            },
        )
        _build_scenario_switch(
            d, "tcell", slope=cfg["slope"], direction=cfg["direction"],
            contrasts={"Normal_vs_StageI": [0.01, 0.5, 0.2, 0.001]},
        )
        conc_path = d / "B2_gex_celltyping" / "concordance.csv"
        conc_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"GSM": ["a", "b"], "ARI_leiden": [0.4, 0.5], "AMI_leiden": [0.6, 0.55]}).to_csv(
            conc_path, index=False
        )
        scenario_dirs[name] = d
    return scenario_dirs


def test_harvest_filter_effect_pas_counts_computes_delta_vs_baseline(tmp_path):
    scenario_dirs = _default_scenario_dirs(tmp_path)
    df = sa.harvest_filter_effect_pas_counts(scenario_dirs)

    row = df.set_index("scenario")
    assert row.loc["baseline", "n_pas"] == 1000
    assert row.loc["baseline", "n_pas_dropped_vs_baseline"] == 0
    assert row.loc["atlas_filter", "n_pas"] == 950
    assert row.loc["atlas_filter", "n_pas_dropped_vs_baseline"] == 50
    assert row.loc["ip_filter", "n_pas_dropped_vs_baseline"] == 150
    # peak_filters_stats.json fields surfaced
    assert row.loc["baseline", "n_ip_flagged"] == 40
    assert row.loc["baseline", "annot_filtered_pos"] == 5


def test_harvest_filter_effect_pas_counts_missing_baseline_no_delta_column(tmp_path):
    scenario_dirs = _default_scenario_dirs(tmp_path)
    del scenario_dirs["baseline"]
    df = sa.harvest_filter_effect_pas_counts(scenario_dirs)
    assert "n_pas_dropped_vs_baseline" not in df.columns


def test_harvest_filter_effect_clustering(tmp_path):
    scenario_dirs = _default_scenario_dirs(tmp_path)
    for name, d in scenario_dirs.items():
        _make_h5ad(d / "07_clustering" / "GSMfoo" / "clusters.h5ad", n_cells=20, n_clusters=4)
    df = sa.harvest_filter_effect_clustering(scenario_dirs)
    assert set(df["scenario"]) == set(scenario_dirs)
    assert (df["n_clusters_mean"] == 4).all()
    assert (df["ARI_leiden_mean"] == 0.45).all()


def test_harvest_filter_effect_fisher_and_trend(tmp_path):
    scenario_dirs = _default_scenario_dirs(tmp_path)
    fisher_df = sa.harvest_filter_effect_fisher(scenario_dirs, fdr=0.05)
    assert set(fisher_df["scenario"]) == set(scenario_dirs)
    # 2 of the 4 qvalues (0.01, 0.001) are < 0.05
    assert (fisher_df["n_sig"] == 2).all()

    trend_df = sa.harvest_filter_effect_trend(scenario_dirs)
    assert set(trend_df["scenario"]) == set(scenario_dirs)
    ip_row = trend_df.set_index("scenario").loc["ip_filter"]
    assert ip_row["direction"] == "increasing"

    fisher_summary = sa.summarize_filter_effect_fisher(fisher_df)
    assert set(fisher_summary["scenario"]) == set(scenario_dirs)
    assert (fisher_summary["total_sig"] == 2).all()


def test_analyze_filter_effect_end_to_end(tmp_path):
    scenario_dirs = _default_scenario_dirs(tmp_path)
    for name, d in scenario_dirs.items():
        _make_h5ad(d / "07_clustering" / "GSMfoo" / "clusters.h5ad", n_cells=20, n_clusters=4)

    out_dir = tmp_path / "analysis"
    summary = sa.analyze_filter_effect(scenario_dirs, out_dir, fdr=0.05)

    assert (out_dir / "filter_effect_summary.json").exists()
    assert (out_dir / "filter_effect_summary.md").exists()
    assert (out_dir / "tables" / "filter_effect_pas_counts.csv").exists()
    assert (out_dir / "tables" / "filter_effect_clustering.csv").exists()
    assert (out_dir / "tables" / "filter_effect_fisher_summary.csv").exists()
    assert (out_dir / "tables" / "filter_effect_trend_summary.csv").exists()
    assert (out_dir / "figures" / "filter_effect_summary.png").exists()
    assert (out_dir / "figures" / "filter_effect_summary.svg").exists()

    assert summary["baseline"] == "baseline"
    assert summary["scenarios"] == list(scenario_dirs)

    md_text = (out_dir / "filter_effect_summary.md").read_text()
    assert "atlas_filter" in md_text
    assert "ip_filter" in md_text

    # Reloadable JSON with no leftover numpy/Path types.
    reloaded = json.loads((out_dir / "filter_effect_summary.json").read_text())
    assert reloaded["pas_counts"][0]["scenario"] == "baseline"


def test_analyze_filter_effect_missing_scenario_dir_does_not_crash(tmp_path):
    scenario_dirs = {"baseline": tmp_path / "baseline", "atlas_filter": tmp_path / "nonexistent"}
    (tmp_path / "baseline").mkdir()
    out_dir = tmp_path / "analysis"
    summary = sa.analyze_filter_effect(scenario_dirs, out_dir)
    assert (out_dir / "filter_effect_summary.json").exists()
    assert summary["scenarios"] == ["baseline", "atlas_filter"]
