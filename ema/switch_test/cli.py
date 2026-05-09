"""CLI entry point for the APA switch test.

Runs differential APA + PDUI on user-selected cluster pairs, on a
marker-restricted PAS subset (top-N per cluster, union across clusters).
This is the heavy-compute companion to the main `ema` pipeline — kept as a
separate command so the user controls scope (which cluster pairs, how many
markers, which strategies) instead of running everything unconditionally.
"""

from __future__ import annotations
import argparse
from itertools import combinations
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.sparse as sp

from joblib import parallel_backend

from ema.switch_test.strategies import get_diff_strategy, list_diff_strategies
from ema.quantification.strategies import get_pdui_strategy, list_pdui_strategies
from ema.quantification.marker_selector import (
    select_marker_pas,
    restrict_count_matrix,
    save_markers,
)
from ema.utils import ResourceManager


def parse_cell_combinations(s: str) -> list[tuple[str, str]]:
    """'0,1;0,2;1,2' -> [('0','1'), ('0','2'), ('1','2')]"""
    out = []
    for pair in s.split(";"):
        parts = pair.split(",")
        if len(parts) != 2:
            raise ValueError(f"Each pair must be 'A,B', got: {pair}")
        out.append((parts[0].strip(), parts[1].strip()))
    return out


def build_count_dfs(adata: ad.AnnData):
    """Return (pdui_df, diff_df, cell_index, pas_index)."""
    X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    try:
        pas_index = [int(v) for v in adata.var_names]
    except (TypeError, ValueError):
        pas_index = list(adata.var_names)
    cell_index = list(adata.obs_names)
    pdui_df = pd.DataFrame(X.T, index=pas_index, columns=cell_index)
    diff_df = pd.DataFrame(X, index=cell_index, columns=pas_index)
    return pdui_df, diff_df, cell_index, pas_index


def cli():
    parser = argparse.ArgumentParser(prog="ema_switch")

    # Input
    parser.add_argument(
        "--h5ad", type=str, required=True,
        help="Path to clustered AnnData h5ad (output of `ema` main pipeline). "
             "For multi-sample, point at one dataset's h5ad: "
             "emaout/per_dataset/{ds_id}/clusters.h5ad",
    )
    parser.add_argument(
        "--pasbed", type=str, default=None,
        help="Optional PAS BED path for isoform-aware PDUI. "
             "Defaults to emaout/pasbed.bed",
    )
    parser.add_argument(
        "--gtf", type=str, default=None,
        help="Optional GTF for isoform-aware PDUI. Required if --pdui-method "
             "includes per_isoform aggregation.",
    )
    parser.add_argument(
        "--output-dir", dest="output_dir", type=str, required=True,
        help="Where to write differential/* and pdui_*.tsv files.",
    )

    # Cluster pair selection
    parser.add_argument(
        "--cell_combinations", "--cluster-pairs", dest="cell_combinations",
        type=str, default=None,
        help="Cluster pairs to compare. Format: '0,1;0,2;1,2'. "
             "If omitted: all C(K,2) pairs (only with --diff-method=fisher; "
             "NB methods require explicit selection).",
    )
    parser.add_argument(
        "--cluster-key", dest="cluster_key", type=str, default="leiden",
        help="adata.obs column with cluster labels (default: leiden)",
    )

    # Marker-based subsetting
    parser.add_argument(
        "--marker-top-n", dest="marker_top_n", type=int, default=200,
        help="Top-N marker PAS per cluster to retain for differential testing "
             "(default: 200). Set to 0 to disable subsetting (use all PAS).",
    )
    parser.add_argument(
        "--marker-method", dest="marker_method", type=str, default="wilcoxon",
        choices=["wilcoxon", "t-test", "logreg"],
        help="Marker ranking method (default: wilcoxon)",
    )

    # Differential strategy
    parser.add_argument(
        "--diff-method", dest="diff_method", type=str, default="fisher",
        choices=list_diff_strategies(),
        help=f"Differential APA test (default: fisher). Available: "
             f"{', '.join(list_diff_strategies())}",
    )
    parser.add_argument(
        "--fdr-threshold", dest="fdr_threshold", type=float, default=0.05,
        help="FDR q-value threshold for significance (default: 0.05)",
    )

    # PDUI strategies (comma-separated)
    parser.add_argument(
        "--pdui-method", dest="pdui_method", type=str, default="classic",
        help=f"Comma-separated PDUI methods. Available: "
             f"{', '.join(list_pdui_strategies())}. Default: classic. "
             f"Use 'none' to skip PDUI.",
    )
    parser.add_argument(
        "--pdui-isoform-agg", dest="pdui_isoform_agg", type=str,
        default="per_gene", choices=["per_gene", "per_isoform"],
        help="PDUI aggregation level (default: per_gene)",
    )
    parser.add_argument(
        "--pdui-isoform-collapse", dest="pdui_isoform_collapse", type=str,
        default="none", choices=["none", "mean", "majority"],
        help="When per_isoform: how to summarize isoforms (default: none)",
    )

    # Resource control
    parser.add_argument(
        "--max-jobs", dest="max_jobs", type=int, default=None,
        help="Override n_jobs for parallel work. Default: auto from ResourceManager.",
    )
    parser.add_argument(
        "--per-worker-mb", dest="per_worker_mb", type=int, default=300,
        help="Estimated peak RAM per parallel worker (default: 300 MB).",
    )

    args = parser.parse_args()

    # Resource decisions
    rm = ResourceManager()
    n_jobs = args.max_jobs if args.max_jobs else rm.get_n_jobs(args.per_worker_mb)
    print(f"[ema_switch] resources: {rm.report()}")
    print(f"[ema_switch] using n_jobs={n_jobs}")

    # Load
    print(f"[ema_switch] loading {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Marker selection
    if args.marker_top_n > 0:
        print(f"[ema_switch] selecting top {args.marker_top_n} markers per cluster "
              f"({args.marker_method})")
        markers = select_marker_pas(
            adata,
            cluster_key=args.cluster_key,
            top_n_per_cluster=args.marker_top_n,
            method=args.marker_method,
        )
        save_markers(markers, output_dir / "markers.tsv")
        print(f"[ema_switch] {len(markers)} unique marker PAS selected "
              f"(union across {adata.obs[args.cluster_key].nunique()} clusters)")
    else:
        markers = None
        print("[ema_switch] marker subsetting disabled — using all PAS")

    # Build oriented count matrices
    pdui_df_full, diff_df_full, cell_idx, pas_idx = build_count_dfs(adata)
    if markers is not None:
        marker_set = set(markers)
        pdui_df = restrict_count_matrix(pdui_df_full, list(marker_set), axis="rows")
        diff_df = restrict_count_matrix(diff_df_full, list(marker_set), axis="cols")
        print(f"[ema_switch] count matrix restricted: "
              f"{pdui_df.shape[0]} PAS x {pdui_df.shape[1]} cells")
    else:
        pdui_df, diff_df = pdui_df_full, diff_df_full

    cluster_labels = pd.Series(adata.obs[args.cluster_key].values, index=adata.obs_names)
    unique_clusters = sorted(cluster_labels.unique().astype(str).tolist())

    # PDUI
    pdui_methods = [m.strip() for m in args.pdui_method.split(",") if m.strip() and m.strip() != "none"]
    if pdui_methods:
        # PDUI strategies need pas_isoform_map for gene grouping (even per_gene aggregation
        # uses the gene_id from the map). Auto-default --pasbed/--gtf if not given.
        pasbed = args.pasbed or "emaout/pasbed.bed"
        gtf = args.gtf or "data/Homo_sapiens.GRCh38.99.gtf"
        if not Path(pasbed).exists() or not Path(gtf).exists():
            print(f"[ema_switch] PDUI requires --gtf and --pasbed (or default paths). "
                  f"pasbed={pasbed} (exists={Path(pasbed).exists()}), "
                  f"gtf={gtf} (exists={Path(gtf).exists()}). Skipping PDUI.")
            pdui_methods = []
        else:
            from ema.annotate.gtf2isoform_utr import parse_isoform_utrs
            from ema.quantification.pas_to_isoform import map_pas_to_isoforms
            isoform_utrs = parse_isoform_utrs(Path(gtf))
            pas_isoform_map = map_pas_to_isoforms(Path(pasbed), isoform_utrs)
            print(f"[ema_switch] isoform map: {len(pas_isoform_map)} PAS mapped")

        for method in pdui_methods:
            strat = get_pdui_strategy(method)
            with parallel_backend("loky", n_jobs=n_jobs):
                df = strat.compute(
                    count_matrix=pdui_df,
                    pas_isoform_map=pas_isoform_map,
                    aggregation=args.pdui_isoform_agg,
                    isoform_collapse=args.pdui_isoform_collapse,
                )
            out_path = output_dir / f"pdui_{method}.tsv"
            df.to_csv(out_path, sep="\t", index=False)
            print(f"[ema_switch] PDUI ({method}): {len(df)} rows -> {out_path}")
    else:
        print("[ema_switch] PDUI skipped")

    # Differential APA
    diff_strat = get_diff_strategy(args.diff_method)
    diff_dir = output_dir / "differential"
    diff_dir.mkdir(exist_ok=True)

    if diff_strat.supports_multi_condition:
        print(f"[ema_switch] running {args.diff_method} omnibus across "
              f"{len(unique_clusters)} clusters")
        df = diff_strat.test(
            count_matrix=diff_df,
            cluster_labels=cluster_labels,
            n_jobs=n_jobs,
        )
        out_path = diff_dir / f"{args.diff_method}_omnibus.tsv"
        df.to_csv(out_path, sep="\t")
        sig = (df["qvalue"] < args.fdr_threshold).sum() if "qvalue" in df.columns else 0
        print(f"[ema_switch] omnibus: {sig} significant PAS (q<{args.fdr_threshold}) -> {out_path}")
    else:
        if args.cell_combinations:
            pairs = parse_cell_combinations(args.cell_combinations)
        else:
            pairs = list(combinations(unique_clusters, 2))
            if args.diff_method != "fisher":
                print(f"[ema_switch] WARNING: --diff-method={args.diff_method} on "
                      f"all {len(pairs)} cluster pairs may be slow. "
                      f"Use --cell_combinations to subset.")
        print(f"[ema_switch] running {args.diff_method} on {len(pairs)} cluster pairs")
        total_sig = 0
        for c1, c2 in pairs:
            df = diff_strat.test(
                count_matrix=diff_df,
                cluster_labels=cluster_labels,
                cluster1=c1,
                cluster2=c2,
                n_jobs=n_jobs,
            )
            out_path = diff_dir / f"{args.diff_method}_{c1}_vs_{c2}.tsv"
            df.to_csv(out_path, sep="\t")
            if "qvalue" in df.columns:
                total_sig += (df["qvalue"] < args.fdr_threshold).sum()
        print(f"[ema_switch] {len(pairs)} pairs: {total_sig} total significant PAS "
              f"(q<{args.fdr_threshold}) -> {diff_dir}")

    print("[ema_switch] done.")


if __name__ == "__main__":
    cli()
