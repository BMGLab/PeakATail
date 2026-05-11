"""Runner functions for the APA switch test (differential APA + PDUI).

Provides library-level ``run_diff`` and ``run_length`` entry points consumed
by ``ema/cli/switch_diff.py`` and ``ema/cli/switch_length.py``.  The old
argparse ``cli()`` has been removed; use ``ema switch diff`` / ``ema switch
length`` instead.
"""

from __future__ import annotations
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
) -> dict[tuple[str, str], pd.DataFrame]:
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

    Returns:
        Dict mapping ``(c1, c2)`` pairs to their result DataFrames (all h5ads
        accumulated). For omnibus strategies the key is ``("omnibus", "")`` and
        the value is the omnibus DataFrame. Empty dict if no pairs were run.
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

    # Accumulate all pair results across h5ads (returned to caller for viz).
    all_pair_results: dict[tuple[str, str], pd.DataFrame] = {}

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
            all_pair_results[("omnibus", "")] = df
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
                all_pair_results[(c1, c2)] = df

            log.info(
                "run_diff: %d pairs: %d total significant PAS (q<%s) -> %s",
                len(pairs), total_sig, fdr, diff_dir,
            )

    log.info("run_diff: done.")
    return all_pair_results


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
) -> tuple[pd.DataFrame | None, "ad.AnnData | None"]:
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
        isoform_agg: Aggregation level — ``"per_gene"`` or ``"per_isoform"``.
            Legacy values ``"gene"`` / ``"isoform"`` are translated to the
            ``per_*`` form so existing scripts/YAMLs do not silently invoke
            the wrong branch.
        isoform_collapse: How to summarise isoforms — ``"none"`` / ``"mean"`` /
            ``"majority"``.
        threads: Max parallel workers ceiling (or None for auto).

    Returns:
        Tuple ``(pdui_df, adata)`` from the last h5ad processed, or
        ``(None, None)`` if nothing was computed.  ``pdui_df`` is the PDUI
        result DataFrame; ``adata`` is the loaded AnnData (for cluster obs).
    """
    # Backward-compat: translate legacy {"gene","isoform"} tokens to the
    # canonical {"per_gene","per_isoform"} the strategy code branches on.
    # Without this, isoform_agg="gene" silently runs the per_isoform path.
    _LEGACY_AGG = {"gene": "per_gene", "isoform": "per_isoform"}
    if isoform_agg in _LEGACY_AGG:
        log.warning(
            "isoform_agg=%r is a legacy alias; using %r. Please update callers.",
            isoform_agg, _LEGACY_AGG[isoform_agg],
        )
        isoform_agg = _LEGACY_AGG[isoform_agg]
    if isoform_agg not in ("per_gene", "per_isoform"):
        raise ValueError(
            f"isoform_agg={isoform_agg!r} invalid; expected per_gene or per_isoform"
        )
    if isoform_collapse not in ("none", "mean", "majority"):
        raise ValueError(
            f"isoform_collapse={isoform_collapse!r} invalid; expected none|mean|majority"
        )

    reset_resource_manager()
    rm = ResourceManager(user_max_threads=threads)
    import ema.utils as _ema_utils
    _ema_utils._RM_INSTANCE = rm

    n_jobs = rm.get_n_jobs(per_worker_mb=300, stage="run_length")
    log.info("run_length: resources: %s", rm.report())

    # Loud warning for the documented no-op flag rather than silent acceptance.
    if cluster_pairs is not None and cluster_pairs != "":
        log.warning(
            "ema switch length: --cluster-pairs is currently unused for length "
            "analysis (your value %r will NOT filter the output).",
            cluster_pairs,
        )

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdui_methods = [strategy] if strategy and strategy != "none" else []

    # Track the last computed PDUI df + adata for viz (returned to caller).
    _last_pdui_df: pd.DataFrame | None = None
    _last_adata: "ad.AnnData | None" = None

    for h5ad_path in h5ad_paths:
        log.info("run_length: loading %s", h5ad_path)
        adata = ad.read_h5ad(h5ad_path)
        _last_adata = adata

        pdui_df_full, _, _, _ = build_count_dfs(adata)

        if not pdui_methods:
            log.info("run_length: PDUI skipped (no strategy)")
            continue

        # Build pas_isoform_map. Per_isoform needs the GTF (transcript-level
        # UTRs); per_gene only needs the PAS->gene assignment which now lives
        # directly on adata.var (written by ema.matrixfilter.preprocessing when
        # gene_ids are forwarded from annotate).
        pas_isoform_map: dict[int, list[tuple[str, str, int, int, int]]] = {}
        if gtf and Path(gtf).exists() and isoform_agg == "per_isoform":
            # Isoform-aware path requires both GTF and a PAS BED with strand info.
            # The pasbed path is derived from the run dir (set by the new
            # auto-routing helper) so the deprecated `emaout/pasbed.bed`
            # global lookup is gone.
            from ema.annotate.gtf2isoform_utr import parse_isoform_utrs
            from ema.quantification.pas_to_isoform import map_pas_to_isoforms
            # Cache the parsed isoform UTR map next to the source h5ad's run
            # dir (or fall back to the user-level ~/.cache).  Without this the
            # human GTF is re-parsed (~1 GB on disk, ~3-5 GB in memory) on
            # every `ema switch length --isoform-agg per_isoform` invocation.
            # Cap workers via the same ResourceManager n_jobs we already
            # computed for this run so we don't spawn one-per-chromosome.
            _gtf_cache = Path(h5ad_path).resolve().parent
            for _ in range(4):
                if (_gtf_cache / "gtf_cache").exists() or (_gtf_cache / "run_config.json").exists():
                    _gtf_cache = _gtf_cache / "gtf_cache"
                    break
                _gtf_cache = _gtf_cache.parent
            else:
                _gtf_cache = Path.home() / ".cache" / "peakatail" / "gtf"
            _gtf_cache.mkdir(parents=True, exist_ok=True)
            isoform_utrs = parse_isoform_utrs(
                Path(gtf),
                cache_dir=_gtf_cache,
                n_workers=min(n_jobs, 4),
            )
            # Find pasbed.bed: walk up from the h5ad path, then fall back to
            # concatenating the per-dataset BEDs the pipeline already writes
            # under <run>/peakcalling/.  The deprecated hardcoded
            # `emaout/pasbed.bed` lookup is gone.
            pasbed_candidates: list[Path] = []
            h5ad_run_root = Path(h5ad_path).resolve().parent
            search_root = h5ad_run_root
            for _ in range(4):  # walk up at most 4 dirs to the run root
                pasbed_candidates.append(search_root / "pasbed.bed")
                search_root = search_root.parent
            chosen_pasbed = next((p for p in pasbed_candidates if p.exists()), None)

            if chosen_pasbed is None:
                # Build pasbed.bed on the fly from peakcalling/*.{pos,neg}.bed.
                peak_dir = h5ad_run_root
                for _ in range(4):
                    if (peak_dir / "peakcalling").is_dir():
                        peak_dir = peak_dir / "peakcalling"
                        break
                    peak_dir = peak_dir.parent
                if peak_dir.name == "peakcalling":
                    beds = sorted(peak_dir.glob("*.pos.bed")) + sorted(peak_dir.glob("*.neg.bed"))
                    if beds:
                        combined = peak_dir.parent / "pasbed.bed"
                        with open(combined, "w") as out:
                            for b in beds:
                                with open(b) as f:
                                    for line in f:
                                        if line.strip():
                                            out.write(line)
                        chosen_pasbed = combined
                        log.info(
                            "run_length: built pasbed.bed by concatenating "
                            "%d peak BEDs -> %s",
                            len(beds), combined,
                        )
            if chosen_pasbed is not None:
                pas_isoform_map = map_pas_to_isoforms(chosen_pasbed, isoform_utrs)
                log.info(
                    "run_length: isoform map: %d PAS mapped from %s",
                    len(pas_isoform_map), chosen_pasbed,
                )
            else:
                log.warning(
                    "run_length: --isoform-agg=per_isoform requires a pasbed; "
                    "none found near %s. Falling back to per_gene.",
                    Path(h5ad_path).resolve(),
                )
                isoform_agg = "per_gene"
        elif isoform_agg == "per_isoform" and (not gtf or not Path(gtf).exists()):
            log.warning(
                "run_length: --isoform-agg=per_isoform requires --gtf; "
                "falling back to per_gene"
            )
            isoform_agg = "per_gene"

        # Per_gene path: synthesise a minimal pas_isoform_map from adata.var.
        # Each PAS maps to a single (gene, "_gene_", rank=1) entry — enough
        # for the strategy's per_gene aggregation, which collapses isoforms.
        if isoform_agg == "per_gene" and not pas_isoform_map:
            if "gene_id" in adata.var.columns:
                gene_id_col = adata.var["gene_id"]
                for pas_id_str, gene_id in gene_id_col.items():
                    if not gene_id or (isinstance(gene_id, float) and pd.isna(gene_id)):
                        continue
                    try:
                        pas_id_int = int(pas_id_str)
                    except (TypeError, ValueError):
                        continue
                    pas_isoform_map[pas_id_int] = [
                        (str(gene_id), "_gene_", 0, 1, 1)
                    ]
                log.info(
                    "run_length: per_gene map built from adata.var['gene_id']: "
                    "%d PAS mapped",
                    len(pas_isoform_map),
                )
            else:
                log.warning(
                    "run_length: per_gene requested but adata.var has no "
                    "'gene_id' column (older h5ad?). PDUI will be empty."
                )

        # Use the **threading** joblib backend for PDUI strategies.  Two
        # reasons, both observed on the full-BAM regression run:
        #   1. The default loky backend pickles a copy of ``pdui_df_full``
        #      (a 1051 x 22419 dense DataFrame ~ 190 MB) into each worker.
        #      With n_jobs=4 the parent peaked at ~13 GB RSS and crashed the
        #      host; threads share memory so this stays at ~600 MB.
        #   2. loky workers also reinitialise OpenBLAS thread pools and
        #      deadlocked even with the BLAS env caps in ema/__init__.py.
        # Inner pandas/numpy ops release the GIL, so threading still gets
        # real parallelism for the gene-level loops.  threadpool_limits(1)
        # belt-and-braces caps any BLAS calls inside the threads.
        from contextlib import contextmanager
        try:
            from threadpoolctl import threadpool_limits as _threadpool_limits
        except ImportError:
            @contextmanager
            def _threadpool_limits(limits=1):  # type: ignore[misc]
                yield
        for method in pdui_methods:
            strat = get_pdui_strategy(method)
            with _threadpool_limits(limits=1):
                with parallel_backend("threading", n_jobs=n_jobs):
                    df = strat.compute(
                        count_matrix=pdui_df_full,
                        pas_isoform_map=pas_isoform_map,
                        aggregation=isoform_agg,
                        isoform_collapse=isoform_collapse,
                    )
            out_path = out_dir / f"pdui_{method}.tsv"
            df.to_csv(out_path, sep="\t", index=False)
            log.info("run_length: PDUI (%s): %d rows -> %s", method, len(df), out_path)
            _last_pdui_df = df

    log.info("run_length: done.")
    return _last_pdui_df, _last_adata


# NOTE: The old argparse cli() has been removed.
# Use `ema switch diff` / `ema switch length` instead.
