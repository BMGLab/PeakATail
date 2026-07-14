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
    min_cells_per_group: int = 10,
    pas_gene_map: dict[str, str] | None = None,
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
        min_cells_per_group: Minimum cells per group for a PAS to enter testing.

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
        min_cells_per_group=min_cells_per_group,
        pas_gene_map=pas_gene_map,
    )


def _resolve_pasbed(h5ad_path: str, explicit: str | None = None) -> Path | None:
    """Resolve the pasbed.bed for a clustered h5ad (B4).

    Resolution order:
      1. ``explicit`` (the user's ``--pasbed``) if it exists.
      2. Walk up from the h5ad's directory (up to 4 levels) looking for a
         ``pasbed.bed`` — the run layout puts it a sibling dir away from the
         clustering h5ad, so a naive ``h5ad.parent/pasbed.bed`` misses it.
    Returns the resolved ``Path`` or ``None`` if nothing is found.
    """
    if explicit:
        p = Path(explicit)
        if p.exists():
            return p
    search_root = Path(h5ad_path).resolve().parent
    for _ in range(4):
        candidate = search_root / "pasbed.bed"
        if candidate.exists():
            return candidate
        search_root = search_root.parent
    return None


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
    min_cells_per_group: int = 10,
    progress_manager=None,
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
        min_cells_per_group: Minimum cells (with nonzero counts for NB strategies)
            in each cluster for a PAS to enter differential testing. Default 10.
        progress_manager: Optional :class:`~ema.progress.ProgressManager`.  When
            supplied, a ``"Cluster-pair testing"`` stage is registered and
            advanced once per pair completed (both serial and parallel paths).

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
                min_cells_per_group=min_cells_per_group,
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

            # Register the cluster-pair testing progress stage.  Always show
            # the bar so the user gets uniform feedback — even a 1/1 bar is
            # better than no visible activity.
            _pair_client = None
            if progress_manager is not None:
                _pair_stage = progress_manager.add_stage(
                    "Cluster-pair testing", total=len(pairs)
                )
                _pair_client = progress_manager.client(_pair_stage)

            # Build pas_gene_map for the WITHIN-GENE Fisher framing.  The
            # fisher strategy uses this to group PAS by gene; nb_pairwise /
            # nb_multi accept-and-ignore via **_ignored.  When the h5ad has
            # no gene_id column the map is left as None and Fisher logs a
            # warning + falls back to the (less APA-correct) global path.
            pas_gene_map: dict[str, str] | None = None
            if "gene_id" in adata.var.columns:
                pas_gene_map = {
                    str(k): str(v)
                    for k, v in adata.var["gene_id"].dropna().items()
                    if str(v).strip()
                }
                log.info(
                    "run_diff: pas_gene_map: %d PAS -> %d unique genes "
                    "(within-gene Fisher comparison enabled)",
                    len(pas_gene_map), len(set(pas_gene_map.values())),
                )
            else:
                log.warning(
                    "run_diff: adata.var has no 'gene_id' column; Fisher "
                    "will fall back to cross-gene comparison."
                )

            pair_results: dict[tuple[str, str], pd.DataFrame] = {}

            if n_outer <= 1 or len(pairs) <= 1:
                for c1, c2 in pairs:
                    _, _, df = run_one_pair(
                        strategy, diff_df, cluster_labels, c1, c2,
                        n_jobs_inner=n_inner,
                        min_cells_per_group=min_cells_per_group,
                        pas_gene_map=pas_gene_map,
                    )
                    pair_results[(c1, c2)] = df
                    if _pair_client is not None:
                        _pair_client.advance(1)
            else:
                worker_fn = functools.partial(
                    _dispatch_pair,
                    strategy_name=strategy,
                    diff_df=diff_df,
                    cluster_labels=cluster_labels,
                    n_jobs_inner=n_inner,
                    min_cells_per_group=min_cells_per_group,
                    pas_gene_map=pas_gene_map,
                )
                ctx = multiprocessing.get_context("spawn")
                with ctx.Pool(n_outer) as pool:
                    for c1, c2, df in pool.imap_unordered(worker_fn, pairs, chunksize=1):
                        pair_results[(c1, c2)] = df
                        if _pair_client is not None:
                            _pair_client.advance(1)

            # Build annotation lookup once per h5ad, shared across all pairs.
            # gene_id from adata.var (index = pas_id as str).
            _gene_id_map: pd.Series | None = None
            if "gene_id" in adata.var.columns:
                _gene_id_map = adata.var["gene_id"].copy()
            else:
                log.warning(
                    "run_diff: adata.var has no 'gene_id' column; "
                    "gene_id will be blank in differential TSVs."
                )

            # chrom/start/end/strand from pasbed.bed.
            # B4: honour the explicit --pasbed argument first, then resolve via
            # the run layout (walk up from the h5ad). The old code hardcoded
            # ``h5ad.parent/pasbed.bed`` and ignored the passed ``pasbed`` arg —
            # under the current layout the pasbed is a sibling *dir* away, so the
            # coordinate columns came back silently blank.
            _pasbed_cols: pd.DataFrame | None = None
            _pasbed_path = _resolve_pasbed(h5ad_path, pasbed)
            if _pasbed_path is not None and _pasbed_path.exists():
                try:
                    _bed = pd.read_csv(
                        _pasbed_path,
                        sep="\t",
                        header=None,
                        usecols=[0, 1, 2, 3, 5],
                        names=["chrom", "start", "end", "pas_id", "strand"],
                        dtype=str,
                    )
                    _bed["pas_id"] = _bed["pas_id"].astype(str)
                    _pasbed_cols = _bed.set_index("pas_id")[["chrom", "start", "end", "strand"]]
                    log.info(
                        "run_diff: loaded pasbed for coordinate annotation: %s (%d rows)",
                        _pasbed_path, len(_pasbed_cols),
                    )
                except Exception as _e:
                    log.warning("run_diff: failed to parse pasbed.bed at %s: %s", _pasbed_path, _e)
            else:
                log.debug("run_diff: no pasbed.bed found at %s; coords will be absent", _pasbed_path)

            def _augment_diff_df(
                df_raw: pd.DataFrame,
                c1_label: str,
                c2_label: str,
            ) -> pd.DataFrame:
                """Left-join gene_id and coordinates onto a per-pair result DataFrame.

                Augmented column order: pas_id, gene_id, chrom, start, end,
                strand, cluster1, cluster2, <original stat columns>.
                If a lookup is unavailable the corresponding columns are filled
                with empty strings so the TSV structure stays consistent.
                """
                aug = df_raw.copy()
                # Normalise index name so joins work regardless of whether the
                # strategy set it or left it as a positional range index.
                if aug.index.name != "pas_id":
                    aug.index.name = "pas_id"
                aug = aug.reset_index()  # pas_id becomes a regular column
                aug["pas_id"] = aug["pas_id"].astype(str)

                # --- gene_id join ---
                if _gene_id_map is not None:
                    gene_series = _gene_id_map.rename_axis("pas_id").reset_index()
                    gene_series["pas_id"] = gene_series["pas_id"].astype(str)
                    aug = aug.merge(gene_series, on="pas_id", how="left")
                else:
                    aug["gene_id"] = ""

                # --- coordinate join ---
                if _pasbed_cols is not None:
                    coord_df = _pasbed_cols.rename_axis("pas_id").reset_index()
                    aug = aug.merge(coord_df, on="pas_id", how="left")
                    aug[["chrom", "start", "end", "strand"]] = (
                        aug[["chrom", "start", "end", "strand"]].fillna("")
                    )
                else:
                    aug["chrom"] = ""
                    aug["start"] = ""
                    aug["end"] = ""
                    aug["strand"] = ""

                # --- self-describing cluster columns ---
                aug.insert(0, "cluster2", c2_label)
                aug.insert(0, "cluster1", c1_label)

                # Reorder: pas_id, gene_id, chrom, start, end, strand, cluster1, cluster2, <stats>
                stat_cols = [
                    c for c in aug.columns
                    if c not in {"pas_id", "gene_id", "chrom", "start", "end",
                                 "strand", "cluster1", "cluster2"}
                ]
                aug = aug[
                    ["pas_id", "gene_id", "chrom", "start", "end", "strand",
                     "cluster1", "cluster2"] + stat_cols
                ]
                return aug

            total_sig = 0
            for c1, c2 in pairs:
                df = pair_results[(c1, c2)]
                df_out = _augment_diff_df(df, c1, c2)
                out_path = diff_dir / f"{strategy}_{c1}_vs_{c2}.tsv"
                df_out.to_csv(out_path, sep="\t", index=False)
                if "qvalue" in df.columns:
                    total_sig += (df["qvalue"] < fdr).sum()
                # Store the augmented df (includes gene_id, chrom, strand) so
                # downstream consumers such as the viz gene_track auto-top-N
                # ranking can access gene annotations without re-parsing disk.
                all_pair_results[(c1, c2)] = df_out

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
    pseudocount: float = 0.0,
    progress_client=None,
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
        pseudocount: Added to each per-cell count before PDUI / entropy
            computation.  Default 0.0 preserves original behaviour.  Set to
            e.g. 1.0 to eliminate NaN on zero-count cells.
        progress_client: Optional :class:`~ema.progress.ProgressClient`.  When
            supplied, ``advance(1)`` is called after each h5ad is processed so
            the CLI progress bar ticks forward.

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
                raise FileNotFoundError(
                    f"pasbed.bed not found near {Path(h5ad_path).resolve().parent}. "
                    f"Either re-run the pipeline (it should produce per_dataset/<ds>/pasbed.bed) "
                    f"or pass --pasbed explicitly."
                )
            pas_isoform_map = map_pas_to_isoforms(chosen_pasbed, isoform_utrs)
            log.info(
                "run_length: isoform map: %d PAS mapped from %s",
                len(pas_isoform_map), chosen_pasbed,
            )
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
        # Look up pasbed near the h5ad input so we can decorate each row
        # with genomic coordinates.  Same walk-up logic as the per_isoform
        # branch above, but applied unconditionally — even per_gene runs
        # benefit from coords-per-PAS.
        _pasbed_cols: pd.DataFrame | None = None
        try:
            pb_candidates: list[Path] = []
            sr = Path(h5ad_path).resolve().parent
            for _ in range(4):
                pb_candidates.append(sr / "pasbed.bed")
                sr = sr.parent
            _pb = next((p for p in pb_candidates if p.exists()), None)
            if _pb is not None:
                _pasbed_cols = pd.read_csv(
                    _pb, sep="\t", header=None,
                    names=["chrom", "start", "end", "pas_id", "score", "strand"],
                    usecols=["chrom", "start", "end", "pas_id", "strand"],
                    dtype={"pas_id": str, "chrom": str, "strand": str,
                           "start": "Int64", "end": "Int64"},
                ).set_index("pas_id")[["chrom", "start", "end", "strand"]]
        except Exception as _e:
            log.warning("run_length: could not read pasbed coords: %s", _e)
            _pasbed_cols = None

        # Cell -> cluster map from the AnnData (typically "leiden", but
        # honour whatever the user passed via --cluster-key).
        _cluster_series: pd.Series | None = None
        if cluster_key and cluster_key in adata.obs.columns:
            _cluster_series = adata.obs[cluster_key].astype(str)

        def _augment_pdui_df(df_raw: pd.DataFrame) -> pd.DataFrame:
            """Left-join cluster + per-PAS coordinates onto a strategy output.

            Adds a `cluster` column (cell -> cluster lookup) and, for every
            PAS-id-like column present (`pas_id`, `proximal_pas_id`,
            `distal_pas_id`), four coordinate columns: ``<col>_chrom``,
            ``<col>_start``, ``<col>_end``, ``<col>_strand``.  Falls back
            to empty strings when a lookup is unavailable so the TSV
            structure stays consistent across runs.
            """
            aug = df_raw.copy()

            # Cluster column
            if _cluster_series is not None and "cell" in aug.columns:
                aug["cluster"] = aug["cell"].map(_cluster_series).fillna("")
            else:
                aug["cluster"] = ""

            # Coordinate joins per PAS-id column.  Use a join on a temporary
            # str-cast column so we don't disturb the original dtype.
            pas_cols = [c for c in ("pas_id", "proximal_pas_id", "distal_pas_id")
                        if c in aug.columns]
            for col in pas_cols:
                prefix = "" if col == "pas_id" else col.rsplit("_pas_id", 1)[0] + "_"
                if _pasbed_cols is None:
                    for sub in ("chrom", "start", "end", "strand"):
                        aug[f"{prefix}{sub}"] = ""
                    continue
                key = aug[col].astype(str)
                joined = key.map(_pasbed_cols["chrom"]).rename(f"{prefix}chrom")
                aug[f"{prefix}chrom"] = joined.fillna("")
                aug[f"{prefix}start"] = key.map(_pasbed_cols["start"]).astype("Int64")
                aug[f"{prefix}end"] = key.map(_pasbed_cols["end"]).astype("Int64")
                aug[f"{prefix}strand"] = key.map(_pasbed_cols["strand"]).fillna("")
            return aug

        for method in pdui_methods:
            strat = get_pdui_strategy(method)
            with _threadpool_limits(limits=1):
                with parallel_backend("threading", n_jobs=n_jobs):
                    df = strat.compute(
                        count_matrix=pdui_df_full,
                        pas_isoform_map=pas_isoform_map,
                        aggregation=isoform_agg,
                        isoform_collapse=isoform_collapse,
                        pseudocount=pseudocount,
                    )
            # Filename comes from the strategy class (not f"pdui_{method}"):
            # proportion and shannon produce non-PDUI quantities so prefixing
            # them with "pdui_" was misleading.  See PDUIStrategy.output_filename.
            out_path = out_dir / getattr(strat, "output_filename", f"pdui_{method}.tsv")
            df_out = _augment_pdui_df(df)
            df_out.to_csv(out_path, sep="\t", index=False)
            log.info(
                "run_length: %s: %d rows -> %s (augmented with cluster + coords)",
                method, len(df_out), out_path,
            )
            _last_pdui_df = df

        # Advance the progress bar once per h5ad (all methods for this h5ad done).
        if progress_client is not None:
            try:
                progress_client.advance(1)
            except Exception:
                pass

    log.info("run_length: done.")
    return _last_pdui_df, _last_adata


# NOTE: The old argparse cli() has been removed.
# Use `ema switch diff` / `ema switch length` instead.
