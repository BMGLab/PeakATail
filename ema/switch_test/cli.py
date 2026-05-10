"""CLI entry point for the APA switch test.

Runs differential APA + PDUI on user-selected cluster pairs, on a
marker-restricted PAS subset (top-N per cluster, union across clusters).
This is the heavy-compute companion to the main `ema` pipeline — kept as a
separate command so the user controls scope (which cluster pairs, how many
markers, which strategies) instead of running everything unconditionally.
"""

from __future__ import annotations
import argparse
import functools
import logging
import multiprocessing
from itertools import combinations
from pathlib import Path

log = logging.getLogger(__name__)

import anndata as ad
import pandas as pd
import scipy.sparse as sp

from joblib import parallel_backend

from ema.switch_test.strategies import get_diff_strategy, list_diff_strategies
from ema.switch_test.pair_runner import run_one_pair
from ema.quantification.strategies import get_pdui_strategy, list_pdui_strategies
from ema.quantification.marker_selector import (
    select_marker_pas,
    restrict_count_matrix,
    save_markers,
)
import warnings

from ema.utils import ResourceManager, get_resource_manager, reset_resource_manager


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


def _dispatch_pair(
    pair: tuple[str, str],
    *,
    strategy_name: str,
    diff_df: pd.DataFrame,
    cluster_labels: pd.Series,
    n_jobs_inner: int,
) -> tuple[str, str, pd.DataFrame]:
    """Top-level wrapper for ``run_one_pair`` suitable for ``Pool.imap_unordered``.

    ``functools.partial`` cannot carry keyword-only args in all Python versions,
    so this thin wrapper unpacks the positional ``pair`` argument and forwards
    the rest as keyword args to :func:`run_one_pair`.

    This function must remain at module level (not a closure) so that the
    ``spawn`` pool can pickle it reliably.

    Args:
        pair: ``(cluster1_label, cluster2_label)`` tuple from the pair list.
        strategy_name: Registry key for the differential strategy.
        diff_df: Full count matrix (cells × PAS).
        cluster_labels: Cell-to-cluster assignment Series.
        n_jobs_inner: Inner worker budget from ``ResourceManager.split_jobs()``.

    Returns:
        Forwarded ``(c1, c2, result_df)`` from :func:`run_one_pair`.
    """
    c1, c2 = pair
    return run_one_pair(
        strategy_name,
        diff_df,
        cluster_labels,
        c1,
        c2,
        n_jobs_inner=n_jobs_inner,
    )


def run_diff(
    h5ad_paths: list[str],
    pasbed: str | None,
    gtf: str | None,
    output_dir: str,
    cluster_pairs: str | None,
    cluster_key: str,
    marker_top_n: int,
    marker_method: str,
    strategy: str,
    fdr: float,
    threads: int | None,
    per_worker_mb: int,
) -> None:
    """Library-level entry point for differential APA testing.

    Runs differential APA (Fisher / NB regression) across cluster pairs on
    one or more h5ad files.  Body is the post-argparse logic of the legacy
    ``cli()`` diff pathway, promoted so ``ema/cli/switch_diff.py`` can call it
    without going through argparse.

    Args:
        h5ad_paths: One or more paths to clustered AnnData h5ad files.
        pasbed: Optional PAS BED path for isoform-aware PDUI.
        gtf: Optional GTF path for isoform-aware PDUI.
        output_dir: Directory to write differential/* and pdui_*.tsv files.
        cluster_pairs: ``'c1,c2;c3,c4'``-style pair string, or None for all.
        cluster_key: adata.obs column with cluster labels.
        marker_top_n: Top-N marker PAS per cluster (0 = disabled).
        marker_method: Marker ranking method (wilcoxon / t-test / logreg).
        strategy: Registered differential APA strategy name.
        fdr: FDR q-value threshold for significance.
        threads: Max parallel workers ceiling (or None for auto).
        per_worker_mb: Estimated peak RAM per parallel worker (MB).
    """
    reset_resource_manager()
    rm = ResourceManager(
        user_max_threads=threads,
        user_per_worker_mb=per_worker_mb if per_worker_mb != 300 else None,
    )
    import ema.utils as _ema_utils
    _ema_utils._RM_INSTANCE = rm

    n_jobs = rm.get_n_jobs(per_worker_mb=per_worker_mb, stage="run_diff")
    log.info("run_diff: resources: %s", rm.report())
    log.info("run_diff: n_jobs=%d", n_jobs)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    diff_dir = out_dir / "differential"
    diff_dir.mkdir(exist_ok=True)

    diff_strat = get_diff_strategy(strategy)

    for h5ad_path in h5ad_paths:
        log.info("run_diff: loading %s", h5ad_path)
        adata = ad.read_h5ad(h5ad_path)

        # Marker selection
        if marker_top_n > 0:
            log.info("run_diff: selecting top %d markers per cluster (%s)", marker_top_n, marker_method)
            markers = select_marker_pas(
                adata,
                cluster_key=cluster_key,
                top_n_per_cluster=marker_top_n,
                method=marker_method,
            )
            save_markers(markers, out_dir / "markers.tsv")
            log.info(
                "run_diff: %d unique marker PAS selected (union across %d clusters)",
                len(markers), adata.obs[cluster_key].nunique(),
            )
        else:
            markers = None
            log.info("run_diff: marker subsetting disabled — using all PAS")

        _, diff_df_full, _, _ = build_count_dfs(adata)
        if markers is not None:
            marker_set = set(markers)
            diff_df = restrict_count_matrix(diff_df_full, list(marker_set), axis="cols")
            log.info("run_diff: count matrix restricted to %d PAS", diff_df.shape[1])
        else:
            diff_df = diff_df_full

        cluster_labels = pd.Series(adata.obs[cluster_key].values, index=adata.obs_names)
        unique_clusters = sorted(cluster_labels.unique().astype(str).tolist())

        if diff_strat.supports_multi_condition:
            log.info("run_diff: running %s omnibus across %d clusters", strategy, len(unique_clusters))
            df = diff_strat.test(
                count_matrix=diff_df,
                cluster_labels=cluster_labels,
                n_jobs=n_jobs,
            )
            out_path = diff_dir / f"{strategy}_omnibus.tsv"
            df.to_csv(out_path, sep="\t")
            sig = (df["qvalue"] < fdr).sum() if "qvalue" in df.columns else 0
            log.info("run_diff: omnibus: %d significant PAS (q<%s) -> %s", sig, fdr, out_path)
        else:
            if cluster_pairs:
                pairs = parse_cell_combinations(cluster_pairs)
            else:
                pairs = list(combinations(unique_clusters, 2))
                if strategy != "fisher":
                    log.warning(
                        "run_diff: --strategy=%s on all %d pairs may be slow. "
                        "Use --cluster-pairs to subset.",
                        strategy, len(pairs),
                    )
            log.info("run_diff: running %s on %d cluster pairs", strategy, len(pairs))

            n_outer, n_inner = rm.split_jobs(n_outer=len(pairs), per_inner_mb=300)
            log.info(
                "run_diff: parallelism: %d pairs -> (%d outer, %d inner per pair)",
                len(pairs), n_outer, n_inner,
            )

            pair_results: dict[tuple[str, str], pd.DataFrame] = {}

            if n_outer <= 1 or len(pairs) <= 1:
                for c1, c2 in pairs:
                    _, _, df = run_one_pair(
                        strategy, diff_df, cluster_labels, c1, c2, n_jobs_inner=n_inner,
                    )
                    pair_results[(c1, c2)] = df
            else:
                worker_fn = functools.partial(
                    _dispatch_pair,
                    strategy_name=strategy,
                    diff_df=diff_df,
                    cluster_labels=cluster_labels,
                    n_jobs_inner=n_inner,
                )
                ctx = multiprocessing.get_context("spawn")
                with ctx.Pool(n_outer) as pool:
                    for c1, c2, df in pool.imap_unordered(worker_fn, pairs, chunksize=1):
                        pair_results[(c1, c2)] = df

            total_sig = 0
            for c1, c2 in pairs:
                df = pair_results[(c1, c2)]
                out_path = diff_dir / f"{strategy}_{c1}_vs_{c2}.tsv"
                df.to_csv(out_path, sep="\t")
                if "qvalue" in df.columns:
                    total_sig += (df["qvalue"] < fdr).sum()

            log.info(
                "run_diff: %d pairs: %d total significant PAS (q<%s) -> %s",
                len(pairs), total_sig, fdr, diff_dir,
            )

    log.info("run_diff: done.")


def run_length(
    h5ad_paths: list[str],
    gtf: str | None,
    output_dir: str,
    cluster_pairs: str | None,
    cluster_key: str,
    strategy: str,
    isoform_agg: str,
    isoform_collapse: str,
    threads: int | None,
) -> None:
    """Library-level entry point for 3'UTR length / PDUI quantification.

    Computes per-cluster PDUI scores using the requested strategy on one or
    more h5ad files.  Body is the PDUI pathway of the legacy ``cli()`` function,
    promoted so ``ema/cli/switch_length.py`` can call it without argparse.

    Args:
        h5ad_paths: One or more paths to clustered AnnData h5ad files.
        gtf: Optional GTF path required for isoform-level aggregation.
        output_dir: Directory to write pdui_*.tsv files.
        cluster_pairs: Reserved for future filtering; currently unused.
        cluster_key: adata.obs column with cluster labels.
        strategy: PDUI strategy name (classic / proportion / shannon).
        isoform_agg: Aggregation level — ``"gene"`` or ``"isoform"``.
        isoform_collapse: How to summarise isoforms — ``"none"`` / ``"mean"`` /
            ``"majority"``.
        threads: Max parallel workers ceiling (or None for auto).
    """
    reset_resource_manager()
    rm = ResourceManager(user_max_threads=threads)
    import ema.utils as _ema_utils
    _ema_utils._RM_INSTANCE = rm

    n_jobs = rm.get_n_jobs(per_worker_mb=300, stage="run_length")
    log.info("run_length: resources: %s", rm.report())

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdui_methods = [strategy] if strategy and strategy != "none" else []

    for h5ad_path in h5ad_paths:
        log.info("run_length: loading %s", h5ad_path)
        adata = ad.read_h5ad(h5ad_path)

        pdui_df_full, _, _, _ = build_count_dfs(adata)

        if not pdui_methods:
            log.info("run_length: PDUI skipped (no strategy)")
            continue

        # Build isoform map if GTF provided
        pas_isoform_map = {}
        if gtf and Path(gtf).exists():
            from ema.annotate.gtf2isoform_utr import parse_isoform_utrs
            from ema.quantification.pas_to_isoform import map_pas_to_isoforms
            isoform_utrs = parse_isoform_utrs(Path(gtf))
            # Use emaout/pasbed.bed as default PAS BED (consistent with old cli)
            pasbed_default = Path("emaout/pasbed.bed")
            if pasbed_default.exists():
                pas_isoform_map = map_pas_to_isoforms(pasbed_default, isoform_utrs)
                log.info("run_length: isoform map: %d PAS mapped", len(pas_isoform_map))
        else:
            if isoform_agg == "isoform":
                log.warning("run_length: --isoform-agg=isoform requires --gtf; falling back to gene")
                isoform_agg = "gene"

        for method in pdui_methods:
            strat = get_pdui_strategy(method)
            with parallel_backend("loky", n_jobs=n_jobs):
                df = strat.compute(
                    count_matrix=pdui_df_full,
                    pas_isoform_map=pas_isoform_map,
                    aggregation=isoform_agg,
                    isoform_collapse=isoform_collapse,
                )
            out_path = out_dir / f"pdui_{method}.tsv"
            df.to_csv(out_path, sep="\t", index=False)
            log.info("run_length: PDUI (%s): %d rows -> %s", method, len(df), out_path)

    log.info("run_length: done.")


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
        "--threads", dest="threads", type=int, default=None,
        help="Max parallel workers (overrides auto-detected default). "
             "Respected by ResourceManager as an absolute ceiling for all "
             "parallel stages.",
    )
    parser.add_argument(
        "--max-jobs", dest="max_jobs", type=int, default=None,
        help="[DEPRECATED] Use --threads instead. Kept for backwards "
             "compatibility; forwards to --threads.",
    )
    parser.add_argument(
        "--per-worker-mb", dest="per_worker_mb", type=int, default=300,
        help="Estimated peak RAM per parallel worker (default: 300 MB).",
    )

    args = parser.parse_args()

    # Handle deprecated --max-jobs: forward to --threads with a warning.
    if args.max_jobs is not None:
        warnings.warn(
            "--max-jobs is deprecated; use --threads instead.",
            DeprecationWarning,
            stacklevel=1,
        )
        if args.threads is None:
            args.threads = args.max_jobs

    # Initialise the process-wide ResourceManager singleton with the user
    # ceiling so every downstream call to get_resource_manager() picks it up.
    reset_resource_manager()
    rm = ResourceManager(
        user_max_threads=args.threads,
        user_per_worker_mb=args.per_worker_mb if args.per_worker_mb != 300 else None,
    )
    # Override the lazy singleton with our already-constructed instance.
    import ema.utils as _ema_utils
    _ema_utils._RM_INSTANCE = rm

    n_jobs = rm.get_n_jobs(per_worker_mb=args.per_worker_mb, stage="ema_switch")
    log.info("resources: %s", rm.report())
    log.info("using n_jobs=%d", n_jobs)

    # Load
    log.info("loading %s", args.h5ad)
    adata = ad.read_h5ad(args.h5ad)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Marker selection
    if args.marker_top_n > 0:
        log.info("selecting top %d markers per cluster (%s)", args.marker_top_n, args.marker_method)
        markers = select_marker_pas(
            adata,
            cluster_key=args.cluster_key,
            top_n_per_cluster=args.marker_top_n,
            method=args.marker_method,
        )
        save_markers(markers, output_dir / "markers.tsv")
        log.info(
            "%d unique marker PAS selected (union across %d clusters)",
            len(markers), adata.obs[args.cluster_key].nunique(),
        )
    else:
        markers = None
        log.info("marker subsetting disabled — using all PAS")

    # Build oriented count matrices
    pdui_df_full, diff_df_full, cell_idx, pas_idx = build_count_dfs(adata)
    if markers is not None:
        marker_set = set(markers)
        pdui_df = restrict_count_matrix(pdui_df_full, list(marker_set), axis="rows")
        diff_df = restrict_count_matrix(diff_df_full, list(marker_set), axis="cols")
        log.info("count matrix restricted: %d PAS x %d cells", pdui_df.shape[0], pdui_df.shape[1])
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
            log.warning(
                "PDUI requires --gtf and --pasbed (or default paths). "
                "pasbed=%s (exists=%s), gtf=%s (exists=%s). Skipping PDUI.",
                pasbed, Path(pasbed).exists(), gtf, Path(gtf).exists(),
            )
            pdui_methods = []
        else:
            from ema.annotate.gtf2isoform_utr import parse_isoform_utrs
            from ema.quantification.pas_to_isoform import map_pas_to_isoforms
            isoform_utrs = parse_isoform_utrs(Path(gtf))
            pas_isoform_map = map_pas_to_isoforms(Path(pasbed), isoform_utrs)
            log.info("isoform map: %d PAS mapped", len(pas_isoform_map))

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
            log.info("PDUI (%s): %d rows -> %s", method, len(df), out_path)
    else:
        log.info("PDUI skipped")

    # Differential APA
    diff_strat = get_diff_strategy(args.diff_method)
    diff_dir = output_dir / "differential"
    diff_dir.mkdir(exist_ok=True)

    if diff_strat.supports_multi_condition:
        log.info("running %s omnibus across %d clusters", args.diff_method, len(unique_clusters))
        df = diff_strat.test(
            count_matrix=diff_df,
            cluster_labels=cluster_labels,
            n_jobs=n_jobs,
        )
        out_path = diff_dir / f"{args.diff_method}_omnibus.tsv"
        df.to_csv(out_path, sep="\t")
        sig = (df["qvalue"] < args.fdr_threshold).sum() if "qvalue" in df.columns else 0
        log.info("omnibus: %d significant PAS (q<%s) -> %s", sig, args.fdr_threshold, out_path)
    else:
        if args.cell_combinations:
            pairs = parse_cell_combinations(args.cell_combinations)
        else:
            pairs = list(combinations(unique_clusters, 2))
            if args.diff_method != "fisher":
                log.warning(
                    "--diff-method=%s on all %d cluster pairs may be slow. "
                    "Use --cell_combinations to subset.",
                    args.diff_method, len(pairs),
                )
        log.info("running %s on %d cluster pairs", args.diff_method, len(pairs))

        # Split workers between pair-level (outer) and PAS-level (inner) parallelism
        # to avoid CPU over-subscription when both stages run concurrently.
        n_outer, n_inner = rm.split_jobs(n_outer=len(pairs), per_inner_mb=300)
        log.info(
            "parallelism: %d pairs -> split_jobs gave (%d outer, %d inner per pair)",
            len(pairs), n_outer, n_inner,
        )

        pair_results: dict[tuple[str, str], pd.DataFrame] = {}

        if n_outer <= 1 or len(pairs) <= 1:
            # Sequential path: avoids spawn overhead for small jobs (K=2 clusters
            # → 1 pair) or when --threads 1 forces serial execution.
            for c1, c2 in pairs:
                _, _, df = run_one_pair(
                    args.diff_method,
                    diff_df,
                    cluster_labels,
                    c1,
                    c2,
                    n_jobs_inner=n_inner,
                )
                pair_results[(c1, c2)] = df
        else:
            # Parallel path: spawn-based Pool over cluster pairs.
            # functools.partial creates a pickle-safe partial with fixed args;
            # each call receives (c1, c2) from the pool.
            worker_fn = functools.partial(
                _dispatch_pair,
                strategy_name=args.diff_method,
                diff_df=diff_df,
                cluster_labels=cluster_labels,
                n_jobs_inner=n_inner,
            )
            ctx = multiprocessing.get_context("spawn")
            with ctx.Pool(n_outer) as pool:
                for c1, c2, df in pool.imap_unordered(
                    worker_fn, pairs, chunksize=1
                ):
                    pair_results[(c1, c2)] = df

        # Write results in deterministic order (matches the pair list order).
        total_sig = 0
        for c1, c2 in pairs:
            df = pair_results[(c1, c2)]
            out_path = diff_dir / f"{args.diff_method}_{c1}_vs_{c2}.tsv"
            df.to_csv(out_path, sep="\t")
            if "qvalue" in df.columns:
                total_sig += (df["qvalue"] < args.fdr_threshold).sum()

        log.info(
            "%d pairs: %d total significant PAS (q<%s) -> %s",
            len(pairs), total_sig, args.fdr_threshold, diff_dir,
        )

    log.info("done.")


if __name__ == "__main__":
    cli()
