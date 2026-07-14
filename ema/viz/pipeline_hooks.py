"""Per-CLI-command viz orchestrators.

The strategies in ``ema.viz`` are pure renderers — given the right ``data``
shape they produce files.  This module is the **only** glue between a CLI
command and those strategies.  Different commands produce different sets of
artifacts and therefore call different orchestrator entry points:

================  =================================================
CLI command       Entry point
================  =================================================
``ema run``       :func:`render_run_outputs`
``ema switch diff``     :func:`render_switch_diff_outputs`
``ema switch length``   :func:`render_switch_length_outputs`
``ema switch match``    :func:`render_switch_match_outputs`
================  =================================================

Each entry point reads its inputs (in-memory results passed in by the command,
plus on-disk artifacts written by the pipeline body), assembles the per-strategy
payloads, calls :func:`ema.viz.render_all`, and logs how many files each
render wrote.

All renderers swallow exceptions with a warning so a viz failure never crashes
a successful run.
"""
from __future__ import annotations

import json
import logging
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

from ema.viz import render_all

log = logging.getLogger(__name__)

DEFAULT_ENGINES: list[str] = ["matplotlib"]

# Number of top genes rendered as gene-track figures by the auto-top-N block in
# each switch orchestrator.  No CLI flag yet — Phase 2 if requested.
_AUTO_TOP_N_GENES: int = 5


# ===========================================================================
# Pure data assemblers (shared across commands; safe to unit-test in isolation)
# ===========================================================================

def build_peak_qc_data(bed_paths: Iterable[str | Path], adata: Any) -> dict[str, Any]:
    """Assemble the 4-panel input dict consumed by the ``peak_qc`` strategies."""
    import scipy.sparse as _sp

    peaks_per_chrom: dict[str, int] = {}
    peak_widths: list[int] = []
    for bed in bed_paths:
        try:
            with open(bed) as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 3:
                        continue
                    try:
                        start = int(parts[1])
                        end = int(parts[2])
                    except ValueError:
                        continue
                    chrom = parts[0]
                    peaks_per_chrom[chrom] = peaks_per_chrom.get(chrom, 0) + 1
                    peak_widths.append(end - start)
        except OSError as exc:
            log.warning("peak_qc: could not read %s: %s", bed, exc)

    per_cell_pas: list[int] = []
    per_cell_reads: list[float] = []
    if adata is not None and adata.n_obs > 0 and adata.n_vars > 0:
        X = adata.X
        if _sp.issparse(X):
            per_cell_pas = (X != 0).sum(axis=1).A1.tolist()
            per_cell_reads = X.sum(axis=1).A1.tolist()
        else:
            import numpy as _np
            X_arr = _np.asarray(X)
            per_cell_pas = (X_arr != 0).sum(axis=1).tolist()
            per_cell_reads = X_arr.sum(axis=1).tolist()

    return {
        "peaks_per_chrom": peaks_per_chrom,
        "per_cell_pas": per_cell_pas,
        "peak_widths": peak_widths,
        "per_cell_reads": per_cell_reads,
    }


def _find_pasbed_near(anchor: Path, max_walk: int = 4) -> Path | None:
    """Walk up from *anchor* (a directory or file) to find ``pasbed.bed``.

    Mirrors the same walk-up logic used by the runner.

    Args:
        anchor: Start path (file or directory).
        max_walk: Maximum number of parent directories to ascend.

    Returns:
        First ``pasbed.bed`` found, or ``None``.
    """
    start = anchor if anchor.is_dir() else anchor.parent
    search = start
    for _ in range(max_walk + 1):
        candidate = search / "pasbed.bed"
        if candidate.exists():
            return candidate
        if search == search.parent:
            break
        search = search.parent
    return None


def load_resource_samples(jsonl_path: Path) -> list[dict[str, Any]]:
    """Load resource samples from the JSONL written by ``_ResourceSampler``."""
    samples: list[dict[str, Any]] = []
    if not jsonl_path.exists():
        return samples
    for line in jsonl_path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            samples.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return samples


def load_tile_timings(json_path: Path) -> list[dict[str, Any]]:
    """Load per-tile timing records from the JSON written by ``run_all_jobs``."""
    if not json_path.exists():
        return []
    try:
        data = json.loads(json_path.read_text())
        return data if isinstance(data, list) else []
    except (json.JSONDecodeError, OSError) as exc:
        log.warning("tile_timings.json unreadable: %s", exc)
        return []


def parse_atlas_mapping(mapping_path: Path) -> tuple[int, list[int]]:
    """Read ``atlas_mapping.tsv``, returning ``(snapped_count, distances_bp)``.

    Distances are read from the optional 4th column added by
    :func:`ema.datasets.atlas_snap.snap_beds_to_atlas`.
    """
    snapped = 0
    distances: list[int] = []
    if not mapping_path.exists():
        return snapped, distances
    with open(mapping_path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                snapped += 1
            if len(parts) >= 4:
                try:
                    distances.append(int(parts[3]))
                except ValueError:
                    pass
    return snapped, distances


def count_called_pas(peakcalling_dir: Path) -> int:
    """Count all called PAS across BOTH strands under ``peakcalling_dir``.

    The atlas snap-rate denominator is the number of called PAS on both the
    positive (``*.pos.bed``) and negative (``*.neg.bed``) strands. Counting a
    single strand understates the denominator and inflates the reported snap
    rate (bug D1). Returns 0 if the directory is unreadable.
    """
    try:
        return sum(
            1
            for pattern in ("*.pos.bed", "*.neg.bed")
            for b in Path(peakcalling_dir).glob(pattern)
            for line in open(b)
            if line.strip()
        )
    except OSError:
        return 0


# ===========================================================================
# `ema run` — the main pipeline (single- or multi-sample)
# ===========================================================================

def render_run_outputs(
    *,
    output_dir: Path,
    engines: list[str] | None,
    is_single_sample: bool,
    unique_ds_ids: list[str],
    bed_paths: Iterable[str | Path],
    atlas_enabled: bool,
    clustering_dir: Path | None = None,
    single_sample_h5ad: Path | None = None,
    # Legacy alias kept for external callers that still pass per_dataset_dir.
    per_dataset_dir: Path | None = None,
) -> None:
    """One-call viz for ``ema run``.  Runs after the pipeline body completes.

    Reads only on-disk artifacts the pipeline produced:

    - ``07_clustering/<ds>/clusters.h5ad``                  (multi-sample mode)
    - ``single_sample_h5ad``                               (single-sample mode)
    - ``output_dir/peakcalling/*.{pos,neg}.bed``           (peak BEDs)
    - ``output_dir/unified/atlas_mapping.tsv``             (when --atlas)
    - ``output_dir/resources.jsonl``                       (resource sampler)
    - ``output_dir/tile_timings.json``                     (tile timings)

    Renders, in order:
      1. UMAP + cluster_sizes + peak_qc per dataset
      2. pas_overlap + atlas_snap_diag (multi-sample only)
      3. resource_timeline + tile_timing (always)
    """
    eng = engines if engines is not None else DEFAULT_ENGINES
    if not eng:
        return

    output_dir = Path(output_dir)
    beds_list = [str(b) for b in bed_paths]
    top_figs = output_dir / "figures"
    top_figs.mkdir(parents=True, exist_ok=True)

    # Resolve the clustering directory (new layout: 07_clustering/<ds>/).
    # Accept the legacy per_dataset_dir kwarg for external callers.
    _clustering_dir: Path | None = clustering_dir or per_dataset_dir

    # ----- 1. per-dataset UMAP/cluster_sizes/peak_qc -----------------------
    per_ds_adata: dict[str, Any] = {}
    pas_sets: dict[str, set[str]] = {}

    if is_single_sample:
        if single_sample_h5ad is not None and single_sample_h5ad.exists():
            try:
                import anndata as ad
                ad_obj = ad.read_h5ad(single_sample_h5ad)
                ds_id = "default"
                per_ds_adata[ds_id] = ad_obj
                _render_per_dataset_tier1(ad_obj, ds_id, beds_list, top_figs, eng)
            except Exception as exc:
                log.warning("single-sample viz failed: %s", exc)
        else:
            log.info("single-sample viz: clusters h5ad not on disk; skipping")
    else:
        if _clustering_dir is None:
            log.warning("multi-sample viz: clustering_dir not provided; skipping")
        else:
            for ds in unique_ds_ids:
                h5 = _clustering_dir / ds / "clusters.h5ad"
                if not h5.exists():
                    continue
                try:
                    import anndata as ad
                    ad_obj = ad.read_h5ad(h5)
                    per_ds_adata[ds] = ad_obj
                    pas_sets[ds] = set(ad_obj.var_names)
                    figs_ds = _clustering_dir / ds / "figures"
                    ds_beds = [b for b in beds_list if Path(b).name.startswith(f"{ds}_")]
                    _render_per_dataset_tier1(ad_obj, ds, ds_beds, figs_ds, eng)
                except Exception as exc:
                    log.warning("per-dataset viz failed for %r: %s", ds, exc)

    # ----- 2. cross-dataset (multi-sample only) ----------------------------
    if not is_single_sample:
        _render_pas_overlap(pas_sets, top_figs, eng)
        if atlas_enabled:
            _render_atlas_snap_diag(output_dir, top_figs, eng)

    # ----- 3. resource_timeline + tile_timing ------------------------------
    _render_resource_timeline(output_dir, top_figs, eng)
    _render_tile_timing(output_dir, top_figs, eng)

    # ----- 4. run_report: aggregates all figures into a single HTML --------
    try:
        from ema.viz.run_report import generate_run_report
        generate_run_report(
            output_dir,
            output_html=top_figs / "run_report.html",
        )
    except Exception as exc:
        log.warning("run_report generation failed: %s", exc)


def _render_per_dataset_tier1(
    adata: Any,
    ds_id: str,
    bed_paths: list[str],
    figures_dir: Path,
    engines: list[str],
) -> None:
    """UMAP + cluster_sizes + peak_qc for one dataset."""
    try:
        figures_dir.mkdir(parents=True, exist_ok=True)
        u = render_all("umap", (adata, ds_id), figures_dir / f"umap_{ds_id}", engines=engines)
        c = render_all(
            "cluster_sizes", (adata, ds_id),
            figures_dir / f"clusters_{ds_id}", engines=engines,
        )
        log.info(
            "viz[%s]: umap=%d, cluster_sizes=%d file(s) at %s",
            ds_id, len(u), len(c), figures_dir,
        )
        if bed_paths:
            qc = build_peak_qc_data(bed_paths, adata)
            q = render_all(
                "peak_qc", qc,
                figures_dir / f"peak_qc_{ds_id}", engines=engines,
            )
            log.info(
                "viz[%s]: peak_qc=%d file(s) (%d chroms, %d peaks)",
                ds_id, len(q),
                len(qc["peaks_per_chrom"]), len(qc["peak_widths"]),
            )
    except Exception as exc:
        log.warning("per-dataset Tier 1 viz failed for %r: %s", ds_id, exc)


def _render_pas_overlap(
    pas_sets: dict[str, set[str]],
    top_figs: Path,
    engines: list[str],
) -> None:
    if len(pas_sets) < 2:
        log.info(
            "pas_overlap skipped: need >=2 datasets with h5ad, got %d",
            len(pas_sets),
        )
        return
    try:
        w = render_all("pas_overlap", pas_sets, top_figs / "pas_overlap", engines=engines)
        log.info("viz: pas_overlap=%d file(s)", len(w))
    except Exception as exc:
        log.warning("pas_overlap viz failed: %s", exc)


def _render_atlas_snap_diag(
    output_dir: Path,
    top_figs: Path,
    engines: list[str],
) -> None:
    try:
        snapped, distances = parse_atlas_mapping(output_dir / "unified" / "atlas_mapping.tsv")
        # D1 fix: the snap-rate denominator must count ALL called PAS on BOTH
        # strands. The old glob only read ``*.pos.bed`` (positive strand), which
        # understated ``all_called`` — and the subsequent ``max(0, ...)`` clamp
        # then hid the resulting ``snapped > all_called`` inconsistency, pinning
        # the reported snap rate near 100%. Glob both strands and drop the clamp.
        all_called = count_called_pas(output_dir / "peakcalling")
        # Invariant: every snapped PAS is a called PAS, so snapped <= all_called.
        # If this fails the denominator is incomplete (missing BED files) — surface
        # it loudly rather than silently clamping unsnapped to 0.
        unsnapped = all_called - snapped
        if unsnapped < 0:
            log.warning(
                "atlas_snap_diag: snapped=%d exceeds all_called=%d — snap-rate "
                "denominator is incomplete (missing peakcalling BEDs?); reporting "
                "unsnapped=0 but the rate is unreliable.",
                snapped, all_called,
            )
            unsnapped = 0
        stats = {
            "snapped": snapped,
            "unsnapped": unsnapped,
            "snap_distances": distances,
        }
        w = render_all("atlas_snap_diag", stats, top_figs / "atlas_snap", engines=engines)
        log.info(
            "viz: atlas_snap_diag=%d file(s) (snapped=%d unsnapped=%d distances=%d)",
            len(w), stats["snapped"], stats["unsnapped"], len(distances),
        )
    except Exception as exc:
        log.warning("atlas_snap_diag viz failed: %s", exc)


def _render_resource_timeline(
    output_dir: Path,
    top_figs: Path,
    engines: list[str],
) -> None:
    try:
        samples = load_resource_samples(output_dir / "resources.jsonl")
        if not samples:
            log.info("viz: resources.jsonl missing/empty, skipping resource_timeline")
            return
        w = render_all(
            "resource_timeline",
            {"samples": samples, "annotations": []},
            top_figs / "resource_timeline",
            engines=engines,
        )
        log.info(
            "viz: resource_timeline=%d file(s) (%d samples)",
            len(w), len(samples),
        )
    except Exception as exc:
        log.warning("resource_timeline viz failed: %s", exc)


def _render_tile_timing(
    output_dir: Path,
    top_figs: Path,
    engines: list[str],
) -> None:
    try:
        tiles = load_tile_timings(output_dir / "tile_timings.json")
        if not tiles:
            log.info("viz: tile_timings.json missing/empty, skipping tile_timing")
            return
        w = render_all("tile_timing", tiles, top_figs / "tile_timing", engines=engines)
        log.info("viz: tile_timing=%d file(s) (%d tiles)", len(w), len(tiles))
    except Exception as exc:
        log.warning("tile_timing viz failed: %s", exc)


# ===========================================================================
# `ema switch diff` — volcano + diff_agreement
# ===========================================================================

def render_switch_diff_outputs(
    *,
    out_dir: Path,
    pair_results: dict[tuple[str, str | None], "Any"],
    fdr: float,
    engines: list[str] | None,
    log2fc_thresh: float = 1.0,
    h5ad_paths: list[str] | None = None,
    cluster_key: str = "leiden",
    pasbed_path: str | None = None,
    progress_manager: Any = None,
) -> None:
    """Render figures specific to ``ema switch diff``.

    ``pair_results`` maps ``(c1, c2)`` (with ``c2`` possibly ``None`` for
    omnibus tests) to the per-pair stats DataFrame.  ``fdr`` and
    ``log2fc_thresh`` flow from the CLI/YAML so the volcano cutoffs
    match the user's settings.

    Args:
        out_dir: Output directory for this diff run.
        pair_results: Per-pair result DataFrames keyed by ``(c1, c2)``.
        fdr: FDR threshold used for volcano cutoff lines.
        engines: Rendering engine names (``["matplotlib"]``, etc.).
        log2fc_thresh: log2FC threshold for volcano cutoff lines.
        h5ad_paths: Optional list of h5ad paths (from the CLI ``--h5ad``
            option).  Required for the auto-top-N gene_track block; if
            omitted or empty the block is skipped with a warning.
        cluster_key: AnnData obs column holding cluster labels.
        pasbed_path: Optional explicit path to ``pasbed.bed``.  When
            supplied and the file exists, it is preferred over the
            walk-up heuristic used internally.
    """
    if not engines or not pair_results:
        return
    try:
        figs_dir = Path(out_dir) / "figures"
        figs_dir.mkdir(parents=True, exist_ok=True)

        # --- volcano: one per pair ---
        # Pass cluster IDs + sample sizes (from n_cells_cluster1/cluster2 columns
        # added in commit 31fc463) so each volcano is self-describing in its
        # title and meta sidecar.  Researchers seeing volcano_0_vs_10.png no
        # longer have to grep the log to learn which clusters and how many
        # cells produced it.
        vol_count = 0
        for (c1, c2), df in pair_results.items():
            stem = f"volcano_{c1}_vs_{c2}" if c2 else "volcano_omnibus"
            payload: dict = {
                "df": df, "fdr": fdr, "log2fc_thresh": log2fc_thresh,
                "cluster1": c1, "cluster2": c2,
            }
            # Sample sizes are constant within a pair, so the first row's
            # value is representative.
            for _col, _key in (
                ("n_cells_cluster1", "n_cells_cluster1"),
                ("n_cells_cluster2", "n_cells_cluster2"),
            ):
                if _col in df.columns and len(df) > 0:
                    try:
                        payload[_key] = int(df[_col].iloc[0])
                    except Exception:
                        pass
            # source_tsv path lets the user open the matching differential file.
            if c2:
                payload["source_tsv"] = f"../differential/{stem.replace('volcano_', '')}.tsv"
            w = render_all(
                "volcano", payload,
                figs_dir / stem,
                engines=engines,
            )
            vol_count += len(w)
        log.info("ema switch diff: volcano=%d file(s)", vol_count)

        # --- diff_agreement: needs >=2 pair sig-sets ---
        sig_sets: dict[str, set[str]] = {}
        for (c1, c2), df in pair_results.items():
            key = f"{c1}_vs_{c2}" if c2 else "omnibus"
            if "qvalue" not in df.columns:
                continue
            if "pas_id" in df.columns:
                sig_sets[key] = set(df.loc[df["qvalue"] < fdr, "pas_id"].astype(str))
            else:
                sig_sets[key] = set(df.index[df["qvalue"] < fdr].astype(str))
        # diff_agreement uses upsetplot which enumerates 2^N intersections;
        # cap at 8 sig-sets to keep render time bounded.  Above that the chart
        # is also visually unreadable so a top-K trim would be lossy too.
        _MAX_DIFF_AGREEMENT_SETS = 8
        if 2 <= len(sig_sets) <= _MAX_DIFF_AGREEMENT_SETS:
            w = render_all(
                "diff_agreement", sig_sets,
                figs_dir / "diff_agreement", engines=engines,
            )
            log.info(
                "ema switch diff: diff_agreement=%d file(s) (%d sig-sets)",
                len(w), len(sig_sets),
            )
        elif len(sig_sets) > _MAX_DIFF_AGREEMENT_SETS:
            log.info(
                "ema switch diff: diff_agreement skipped (got %d sig-sets, "
                "max %d for upsetplot — too many for a readable chart)",
                len(sig_sets), _MAX_DIFF_AGREEMENT_SETS,
            )
        else:
            log.info(
                "ema switch diff: diff_agreement skipped (need >=2 pairs, got %d)",
                len(sig_sets),
            )

        # --- auto-top-N gene_track: rank genes by volcano score, render tracks ---
        try:
            _render_diff_gene_tracks(
                pair_results=pair_results,
                h5ad_paths=h5ad_paths or [],
                cluster_key=cluster_key,
                figs_dir=figs_dir,
                engines=engines or [],
                pasbed_path=pasbed_path,
                progress_manager=progress_manager,
            )
        except Exception as _gt_exc:
            log.warning(
                "ema switch diff: gene_track auto top-N failed: %s", _gt_exc
            )

        # Walk the figures dir once and write a researcher-readable manifest.
        from ema.viz._meta import write_figures_index
        idx = write_figures_index(figs_dir, command="ema switch diff")
        if idx is not None:
            log.info("ema switch diff: figures index -> %s", idx)
    except Exception as exc:
        log.warning("ema switch diff: viz rendering failed: %s", exc)


# ---------------------------------------------------------------------------
# Gene-track auto-top-N helpers (shared by diff and length orchestrators)
# ---------------------------------------------------------------------------

def _render_diff_gene_tracks(
    *,
    pair_results: dict[tuple[str, str | None], Any],
    h5ad_paths: list[str],
    cluster_key: str,
    figs_dir: Path,
    engines: list[str],
    pasbed_path: str | None = None,
    progress_manager: Any = None,
) -> None:
    """Rank the top-N genes by volcano score and render one gene_track per gene.

    Args:
        pair_results: Per-pair result DataFrames (same object passed to the
            orchestrator).
        h5ad_paths: Paths to the h5ad files loaded by ``run_diff``.
        cluster_key: AnnData obs column with cluster labels.
        figs_dir: Figures output directory.
        engines: Rendering engine names.
        pasbed_path: Optional explicit path to ``pasbed.bed``.  Preferred
            over the walk-up heuristic when supplied and the file exists.
    """
    if not h5ad_paths:
        log.warning(
            "ema switch diff: gene_track auto top-%d skipped — no h5ad paths "
            "provided (pass h5ad_paths= to render_switch_diff_outputs)",
            _AUTO_TOP_N_GENES,
        )
        return

    from ema.viz._gene_track_helpers import rank_top_genes, build_gene_panel

    top_genes = rank_top_genes(pair_results, n=_AUTO_TOP_N_GENES)
    if not top_genes:
        log.info(
            "ema switch diff: gene_track auto top-%d skipped — no gene_id "
            "column found in pair_results",
            _AUTO_TOP_N_GENES,
        )
        return

    # Load the last (or first available) h5ad and its pasbed.
    adata = _load_last_adata(h5ad_paths)
    if adata is None:
        log.warning(
            "ema switch diff: gene_track auto top-%d skipped — could not load "
            "any h5ad",
            _AUTO_TOP_N_GENES,
        )
        return

    # Prefer explicit pasbed_path; fall back to walk-up heuristic.
    resolved_pasbed: Path | None = None
    if pasbed_path and Path(pasbed_path).exists():
        resolved_pasbed = Path(pasbed_path)
    else:
        resolved_pasbed = _find_pasbed_near(Path(h5ad_paths[-1]))
    if resolved_pasbed is None:
        log.warning(
            "ema switch diff: gene_track auto top-%d skipped — pasbed.bed not "
            "found near %s",
            _AUTO_TOP_N_GENES, Path(h5ad_paths[-1]).parent,
        )
        return
    import pandas as _pd
    pasbed_df = _pd.read_csv(
        resolved_pasbed,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "pas_id", "score", "strand"],
    )

    _gt_client = None
    if progress_manager is not None and len(top_genes) > 0:
        _gt_stage = progress_manager.add_stage(
            f"Gene tracks (top-{len(top_genes)})", total=len(top_genes)
        )
        _gt_client = progress_manager.client(_gt_stage)

    total_files = 0
    rendered_genes: list[str] = []
    for gene_id in top_genes:
        try:
            panel = build_gene_panel(
                gene_id=gene_id,
                adata=adata,
                pasbed=pasbed_df,
                cluster_key=cluster_key,
            )
            if panel is None:
                if _gt_client is not None:
                    _gt_client.advance(1)
                continue
            paths = render_all(
                "gene_track",
                panel,
                figs_dir / f"gene_{gene_id}",
                engines=engines,
            )
            total_files += len(paths)
            rendered_genes.append(gene_id)
        except Exception as _e:
            log.warning(
                "ema switch diff: gene_track for %r failed: %s", gene_id, _e
            )
        finally:
            if _gt_client is not None:
                _gt_client.advance(1)

    log.info(
        "ema switch diff: gene_track auto top-%d = %d file(s) (%d gene(s): %s)",
        _AUTO_TOP_N_GENES,
        total_files,
        len(rendered_genes),
        ", ".join(rendered_genes),
    )


def _render_length_gene_tracks(
    *,
    adata: Any,
    pdui_df: Any,
    cluster_key: str,
    figs_dir: Path,
    engines: list[str],
    out_dir: Path | None = None,
) -> None:
    """Rank top-N genes by within-cluster score variance and render gene tracks.

    For each gene, compute the mean score per cluster then take the variance
    across clusters.  Top-N genes by variance are rendered.

    Args:
        adata: AnnData (already loaded by the length orchestrator).
        pdui_df: Per-cell score DataFrame (from the PDUI/entropy/proportion
            strategy).  May be ``None`` — skips with warning.
        cluster_key: AnnData obs column with cluster labels.
        figs_dir: Figures output directory.
        engines: Rendering engine names.
        out_dir: The switch-length output directory.  Used as an additional
            anchor for the pasbed walk-up when ``adata`` carries no filename.
    """
    if pdui_df is None or pdui_df.empty:
        log.info(
            "ema switch length: gene_track auto top-%d skipped — pdui_df empty",
            _AUTO_TOP_N_GENES,
        )
        return

    import numpy as _np
    import pandas as _pd
    from ema.viz._gene_track_helpers import build_gene_panel

    # Determine which score column to use for variance-based ranking.
    score_col: str | None = None
    for _col in ("pdui", "entropy", "proportion"):
        if _col in pdui_df.columns:
            score_col = _col
            break

    gene_col = "gene_id" if "gene_id" in pdui_df.columns else None

    if score_col is None or gene_col is None or "cell" not in pdui_df.columns:
        log.info(
            "ema switch length: gene_track auto top-%d skipped — pdui_df "
            "missing gene_id, cell, or a score column",
            _AUTO_TOP_N_GENES,
        )
        return

    # Map cell -> cluster from adata.
    if cluster_key not in adata.obs.columns:
        log.warning(
            "ema switch length: gene_track auto top-%d skipped — cluster_key "
            "%r not in adata.obs",
            _AUTO_TOP_N_GENES, cluster_key,
        )
        return

    cell_to_cluster = adata.obs[cluster_key].astype(str).to_dict()
    valid = pdui_df[pdui_df[score_col].notna()].copy()
    valid["__cluster__"] = valid["cell"].map(cell_to_cluster)
    valid = valid[valid["__cluster__"].notna()]

    if valid.empty:
        log.info(
            "ema switch length: gene_track auto top-%d skipped — no non-NaN "
            "score rows after cluster mapping",
            _AUTO_TOP_N_GENES,
        )
        return

    # Mean score per (gene, cluster); variance across clusters.
    per_gene_cluster = (
        valid.groupby([gene_col, "__cluster__"])[score_col]
        .mean()
        .unstack(fill_value=_np.nan)
    )
    # Variance across cluster columns per gene.
    gene_var = per_gene_cluster.var(axis=1, skipna=True)
    top_genes = (
        gene_var.dropna()
        .sort_values(ascending=False)
        .head(_AUTO_TOP_N_GENES)
        .index.tolist()
    )

    if not top_genes:
        log.info(
            "ema switch length: gene_track auto top-%d skipped — no genes with "
            "sufficient cluster coverage",
            _AUTO_TOP_N_GENES,
        )
        return

    # Walk up from adata filename to find pasbed.  adata may not carry its
    # source path — use the h5ad paths registered in obs_names provenance if
    # available; otherwise walk from out_dir / figs_dir.
    pasbed_path: Path | None = None
    h5ad_hint: str | None = getattr(adata, "filename", None)
    if h5ad_hint:
        pasbed_path = _find_pasbed_near(Path(str(h5ad_hint)))
    if pasbed_path is None and out_dir is not None:
        pasbed_path = _find_pasbed_near(out_dir, max_walk=6)
    if pasbed_path is None:
        pasbed_path = _find_pasbed_near(figs_dir, max_walk=6)

    if pasbed_path is None:
        log.warning(
            "ema switch length: gene_track auto top-%d skipped — pasbed.bed "
            "not found",
            _AUTO_TOP_N_GENES,
        )
        return

    pasbed_df = _pd.read_csv(
        pasbed_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "pas_id", "score", "strand"],
    )

    total_files = 0
    rendered_genes: list[str] = []
    for gene_id in top_genes:
        try:
            panel = build_gene_panel(
                gene_id=gene_id,
                adata=adata,
                pasbed=pasbed_df,
                cluster_key=cluster_key,
            )
            if panel is None:
                continue
            paths = render_all(
                "gene_track",
                panel,
                figs_dir / f"gene_{gene_id}",
                engines=engines,
            )
            total_files += len(paths)
            rendered_genes.append(gene_id)
        except Exception as _e:
            log.warning(
                "ema switch length: gene_track for %r failed: %s", gene_id, _e
            )

    log.info(
        "ema switch length: gene_track auto top-%d = %d file(s) (%d gene(s): %s)",
        _AUTO_TOP_N_GENES,
        total_files,
        len(rendered_genes),
        ", ".join(rendered_genes),
    )


def _load_last_adata(h5ad_paths: list[str]) -> Any:
    """Load the last readable h5ad from the list.

    Iterates in reverse so the most-recently-processed dataset is preferred.

    Args:
        h5ad_paths: Paths to h5ad files.

    Returns:
        AnnData object, or ``None`` if none could be loaded.
    """
    try:
        import anndata as _ad
    except ImportError:
        log.warning("gene_track: anndata not installed; skipping")
        return None
    for p in reversed(h5ad_paths):
        try:
            return _ad.read_h5ad(p)
        except Exception as _e:
            log.warning("gene_track: could not load %s: %s", p, _e)
    return None


# ===========================================================================
# `ema switch length` — pdui_distribution + length_shifts
# ===========================================================================

def render_switch_length_outputs(
    *,
    out_dir: Path,
    last_adata: Any,
    pdui_df: Any,
    cluster_key: str,
    engines: list[str] | None,
) -> None:
    """Render figures specific to ``ema switch length``."""
    if not engines or last_adata is None:
        return
    try:
        import numpy as _np
        import pandas as _pd

        figs_dir = Path(out_dir) / "figures"
        figs_dir.mkdir(parents=True, exist_ok=True)

        # --- per-strategy primary figure ---
        #
        # Each PDUI strategy produces a different quantity:
        #   * classic    → 'pdui' column          → per-cluster PDUI violin
        #   * shannon    → 'entropy' column       → per-cluster entropy violin
        #   * proportion → 'proportion' column    → PAS x cluster heatmap
        #
        # Picking the figure by which column is present (rather than by
        # strategy name) keeps the orchestrator decoupled from the strategy
        # registry — a new strategy that emits one of these columns gets the
        # right viz for free.  The previous unconditional ``pdui_distribution``
        # branch silently produced a pas_detection_rate fallback figure
        # mislabelled as PDUI when the strategy wasn't classic.
        cols = set(pdui_df.columns) if pdui_df is not None else set()

        if "pdui" in cols and "cell" in cols:
            score_key = "mean_pdui"
            per_cell = (
                pdui_df.dropna(subset=["pdui"])
                       .groupby("cell")["pdui"].mean()
            )
            last_adata.obs[score_key] = (
                last_adata.obs_names.to_series().map(per_cell)
            )
            n_cells_with_score = int(last_adata.obs[score_key].notna().sum())
            log.info(
                "ema switch length: per-cell PDUI computed (%d/%d cells have a "
                "non-NaN mean PDUI; mean=%.3f, std=%.3f)",
                n_cells_with_score, last_adata.n_obs,
                float(last_adata.obs[score_key].mean(skipna=True)),
                float(last_adata.obs[score_key].std(skipna=True)),
            )
            w = render_all(
                "pdui_distribution", (last_adata, score_key),
                figs_dir / "pdui_distribution", engines=engines,
            )
            log.info("ema switch length: pdui_distribution=%d file(s)", len(w))

        elif "entropy" in cols and "cell" in cols:
            score_key = "mean_entropy"
            per_cell = (
                pdui_df.dropna(subset=["entropy"])
                       .groupby("cell")["entropy"].mean()
            )
            last_adata.obs[score_key] = (
                last_adata.obs_names.to_series().map(per_cell)
            )
            log.info(
                "ema switch length: per-cell entropy computed (%d/%d cells have "
                "a non-NaN mean entropy)",
                int(last_adata.obs[score_key].notna().sum()), last_adata.n_obs,
            )
            w = render_all(
                "entropy_distribution", (last_adata, score_key),
                figs_dir / "entropy_distribution", engines=engines,
            )
            log.info("ema switch length: entropy_distribution=%d file(s)", len(w))

        elif "proportion" in cols and "pas_id" in cols and "cell" in cols:
            w = render_all(
                "proportion_heatmap",
                {
                    "pdui_df": pdui_df,
                    "adata": last_adata,
                    "cluster_key": cluster_key,
                    "top_n": 50,
                },
                figs_dir / "proportion_heatmap", engines=engines,
            )
            log.info("ema switch length: proportion_heatmap=%d file(s)", len(w))

        else:
            # No strategy column present — log loudly, do NOT silently emit a
            # pas_detection_rate figure that looks like a real result.
            log.warning(
                "ema switch length: pdui_df missing a known score column "
                "(pdui/entropy/proportion); skipping primary figure. cols=%s",
                sorted(cols),
            )

        # --- length_shifts: ΔPDUI per (gene × cluster pair) ---
        #
        # Previous implementation materialised the full 4.3M-row pdui_df
        # (one row per (gene, cell)) as a Python list, then computed 66
        # cluster-pair differences using the SAME per-cell vector for every
        # cluster (a bug: shifts would always be zero), spawning a 2+ GB
        # in-memory DataFrame.  On the full BAM this OOM'd the host.
        #
        # The correct computation aggregates PDUI per (gene, cluster) using
        # the cell→cluster mapping from last_adata.obs[cluster_key], then
        # forms pair-wise differences on the small gene × cluster table
        # (typically a few thousand × ~10 = <1 MB).
        can_length_shifts = (
            pdui_df is not None
            and cluster_key in last_adata.obs.columns
            and "pdui" in pdui_df.columns
            and "cell" in pdui_df.columns
            and len(pdui_df) > 0
        )
        if can_length_shifts:
            clusters = sorted(
                last_adata.obs[cluster_key].unique(),
                key=lambda x: int(x) if str(x).isdigit() else x,
            )
            gene_col = "gene_id" if "gene_id" in pdui_df.columns else pdui_df.columns[0]
            pdui_col = "pdui"

            # Map each cell to its cluster (avoid copying pdui_df).
            cell_to_cluster = last_adata.obs[cluster_key].astype(str).to_dict()
            # Pre-filter to rows with non-NaN pdui — typically <10% of the long
            # frame (8.7% in our regression run) — keeping the groupby cheap.
            valid = pdui_df[pdui_df[pdui_col].notna()].copy()
            if valid.empty:
                log.info("ema switch length: length_shifts skipped (no non-NaN PDUI rows)")
            else:
                valid["__cluster__"] = valid["cell"].map(cell_to_cluster)
                valid = valid[valid["__cluster__"].notna()]

                # mean PDUI per (gene, cluster) — small table: n_multi_pas_genes × n_clusters.
                per_gene_cluster = (
                    valid.groupby([gene_col, "__cluster__"])[pdui_col].mean().unstack(fill_value=_np.nan)
                )
                # Ensure every cluster column exists (some may have no data).
                for cl in clusters:
                    if cl not in per_gene_cluster.columns:
                        per_gene_cluster[cl] = _np.nan
                per_gene_cluster = per_gene_cluster[list(clusters)]

                shifts_data: dict[str, _pd.Series] = {}
                for c1, c2 in combinations(clusters, 2):
                    shifts_data[f"{c1}_vs_{c2}"] = per_gene_cluster[c1] - per_gene_cluster[c2]
                if shifts_data:
                    shifts_df = _pd.DataFrame(shifts_data, index=per_gene_cluster.index)
                    w = render_all(
                        "length_shifts", shifts_df,
                        figs_dir / "length_shifts", engines=engines,
                    )
                    log.info(
                        "ema switch length: length_shifts=%d file(s) "
                        "(%d genes × %d cluster pairs)",
                        len(w), len(shifts_df), len(shifts_data),
                    )

        # --- auto-top-N gene_track: rank genes by score variance across clusters ---
        try:
            _render_length_gene_tracks(
                adata=last_adata,
                pdui_df=pdui_df,
                cluster_key=cluster_key,
                figs_dir=figs_dir,
                engines=engines or [],
                out_dir=Path(out_dir),
            )
        except Exception as _gt_exc:
            log.warning(
                "ema switch length: gene_track auto top-N failed: %s", _gt_exc
            )

        from ema.viz._meta import write_figures_index
        idx = write_figures_index(figs_dir, command="ema switch length")
        if idx is not None:
            log.info("ema switch length: figures index -> %s", idx)
    except Exception as exc:
        log.warning("ema switch length: viz rendering failed: %s", exc)


# ===========================================================================
# `ema switch match` — cluster_match_sankey + match_confidence
# ===========================================================================

def render_switch_match_outputs(
    *,
    out_dir: Path,
    df: Any,
    engines: list[str] | None,
) -> None:
    """Render figures specific to ``ema switch match``."""
    if not engines or df is None or getattr(df, "empty", True):
        return
    try:
        figs_dir = Path(out_dir) / "figures"
        figs_dir.mkdir(parents=True, exist_ok=True)
        s = render_all(
            "cluster_match_sankey", df,
            figs_dir / "cluster_match_sankey", engines=engines,
        )
        m = render_all(
            "match_confidence", df,
            figs_dir / "match_confidence", engines=engines,
        )
        log.info(
            "ema switch match: sankey=%d, match_confidence=%d file(s) at %s",
            len(s), len(m), figs_dir,
        )

        from ema.viz._meta import write_figures_index
        idx = write_figures_index(figs_dir, command="ema switch match")
        if idx is not None:
            log.info("ema switch match: figures index -> %s", idx)
    except Exception as exc:
        log.warning("ema switch match: viz rendering failed: %s", exc)
