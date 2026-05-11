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

DEFAULT_ENGINES: list[str] = ["matplotlib", "plotly"]


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
    per_dataset_dir: Path | None = None,
    single_sample_h5ad: Path | None = None,
) -> None:
    """One-call viz for ``ema run``.  Runs after the pipeline body completes.

    Reads only on-disk artifacts the pipeline produced:

    - ``output_dir/per_dataset/<ds>/clusters.h5ad``        (multi-sample mode)
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
        if per_dataset_dir is None:
            log.warning("multi-sample viz: per_dataset_dir not provided; skipping")
        else:
            for ds in unique_ds_ids:
                h5 = per_dataset_dir / ds / "clusters.h5ad"
                if not h5.exists():
                    continue
                try:
                    import anndata as ad
                    ad_obj = ad.read_h5ad(h5)
                    per_ds_adata[ds] = ad_obj
                    pas_sets[ds] = set(ad_obj.var_names)
                    figs_ds = per_dataset_dir / ds / "figures"
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
        try:
            all_called = sum(
                1
                for b in (output_dir / "peakcalling").glob("*.pos.bed")
                for line in open(b)
                if line.strip()
            )
        except OSError:
            all_called = 0
        stats = {
            "snapped": snapped,
            "unsnapped": max(0, all_called - snapped),
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
) -> None:
    """Render figures specific to ``ema switch diff``.

    ``pair_results`` maps ``(c1, c2)`` (with ``c2`` possibly ``None`` for
    omnibus tests) to the per-pair stats DataFrame.
    """
    if not engines or not pair_results:
        return
    try:
        figs_dir = Path(out_dir) / "figures"
        figs_dir.mkdir(parents=True, exist_ok=True)

        # --- volcano: one per pair ---
        # Pass FDR through so the volcano cutoff matches the user's --fdr.
        # log2fc_thresh stays at the renderer default until it's surfaced via
        # RunConfig (Wave 2).
        vol_count = 0
        for (c1, c2), df in pair_results.items():
            stem = f"volcano_{c1}_vs_{c2}" if c2 else "volcano_omnibus"
            w = render_all(
                "volcano",
                {"df": df, "fdr": fdr},
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
    except Exception as exc:
        log.warning("ema switch diff: viz rendering failed: %s", exc)


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
        import scipy.sparse as _sp

        figs_dir = Path(out_dir) / "figures"
        figs_dir.mkdir(parents=True, exist_ok=True)

        # --- pdui_distribution ---
        #
        # Previous implementation assumed ``pdui_df`` was indexed by PAS (one
        # row per PAS aligned to adata.var) and tried to multiply X * pas_pdui
        # to get a per-cell weighted score.  But the PDUI strategies produce
        # LONG format — one row per (gene, transcript, cell) — so the length
        # check ``len(pdui_df) == X.shape[1]`` was *never* true and the
        # fallback ``np.zeros(...)`` made every per-cell score 0.  That's the
        # "all PDUI distributions are zero" bug the user reported.
        #
        # The correct aggregation: mean PDUI per cell across all (gene,
        # transcript) entries, joined onto adata.obs by barcode.  NaN values
        # in the long frame are skipped by groupby().mean(), so cells with
        # no detectable multi-PAS genes correctly land at NaN (not 0).
        if pdui_df is not None and "pdui" in pdui_df.columns and "cell" in pdui_df.columns:
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
        elif cluster_key in last_adata.obs.columns:
            X2 = last_adata.X.toarray() if _sp.issparse(last_adata.X) else last_adata.X
            last_adata.obs["pas_detection_rate"] = (X2 > 0).mean(axis=1)
            w = render_all(
                "pdui_distribution", (last_adata, "pas_detection_rate"),
                figs_dir / "pdui_distribution", engines=engines,
            )
            log.info(
                "ema switch length: pdui_distribution (fallback) = %d file(s)", len(w),
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
        if pdui_df is None or cluster_key not in last_adata.obs.columns:
            return
        clusters = sorted(
            last_adata.obs[cluster_key].unique(),
            key=lambda x: int(x) if str(x).isdigit() else x,
        )
        gene_col = "gene_id" if "gene_id" in pdui_df.columns else pdui_df.columns[0]
        pdui_col = "pdui" if "pdui" in pdui_df.columns else None
        if not pdui_col or len(pdui_df) == 0 or "cell" not in pdui_df.columns:
            return

        # Map each cell to its cluster (avoid copying pdui_df).
        cell_to_cluster = last_adata.obs[cluster_key].astype(str).to_dict()
        # Pre-filter to rows with non-NaN pdui — typically <10% of the long
        # frame (8.7% in our regression run) — keeping the groupby cheap.
        valid = pdui_df[pdui_df[pdui_col].notna()].copy()
        if valid.empty:
            log.info("ema switch length: length_shifts skipped (no non-NaN PDUI rows)")
            return
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
    except Exception as exc:
        log.warning("ema switch match: viz rendering failed: %s", exc)
