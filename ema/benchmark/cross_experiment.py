"""Cross-experiment analyses for the Laughney sweep.

Where ``sweep_analysis.py`` *harvests* what each run already wrote, this module
asks the questions a reviewer asks *across* runs — the ones that decide whether
an observed difference between two branches is real:

1. **Atlas null control** (``atlas_null_control``) — the PolyASite benchmark
   reports precision ~1.0 for every peak strategy. That number is only
   meaningful against a null. This shifts/permutes the called PAS and recomputes
   precision, giving "excess precision over chance".
2. **PAS set overlap** (``pas_set_overlap``) — are two strategies calling the
   *same* sites, or different sites that both happen to land near the atlas?
3. **Cluster concordance** (``cluster_concordance``) — n_clusters / ARI / AMI /
   silhouette for every branch against (a) a reference branch and (b) the GEX
   cell-type labels, so a parameter's effect on clustering is one number.
4. **Cluster -> cell type contingency** (``cluster_celltype_contingency``) —
   the Sankey input, plus how much the mapping *shifts* between branches.
5. **Switch depth confound** (``switch_depth_confound``) — 3'UTR trend slopes
   next to the sequencing depth / detection rate / PDUI-zero fraction they are
   confounded with, and the number of stages the slope was actually fit on.
6. **Fisher power check** (``fisher_power_check``) — significant-hit counts
   against cell counts, the pseudoreplication signature.
7. **Recurrent vs private switches** (``recurrent_switches``).

Everything reads only from finished run directories, so it can be re-run on a
partially complete sweep (missing branches are reported, not fatal).
"""

from __future__ import annotations

import argparse
import json
import logging
import warnings
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Number of stages the ordered-stage trend needs before a Spearman rho is
# anything other than an artefact of having two points.
MIN_STAGES_FOR_TREND = 3


# ---------------------------------------------------------------------------
# small IO helpers
# ---------------------------------------------------------------------------


def _read_json(path: Path) -> dict:
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return {}


def _read_bed_positions(path: Path) -> pd.DataFrame:
    """chrom / start / end / midpoint from a BED-like file (cols 0,1,2)."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={0: str, 1: np.int64, 2: np.int64},
        comment="#",
    )
    df["chrom"] = df["chrom"].astype(str).str.replace("^chr", "", regex=True)
    df["mid"] = ((df["start"] + df["end"]) // 2).astype(np.int64)
    return df


def _load_atlas_index(atlas_bed: Path) -> dict[str, np.ndarray]:
    """chrom -> sorted array of atlas site midpoints."""
    logger.info("loading atlas %s", atlas_bed)
    atlas = pd.read_csv(
        atlas_bed,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={0: str, 1: np.int64, 2: np.int64},
        comment="#",
    )
    atlas["chrom"] = atlas["chrom"].astype(str).str.replace("^chr", "", regex=True)
    mids = ((atlas["start"] + atlas["end"]) // 2).to_numpy(np.int64)
    out: dict[str, np.ndarray] = {}
    for chrom, idx in atlas.groupby("chrom", sort=False).indices.items():
        arr = np.sort(mids[idx])
        out[str(chrom)] = arr
    logger.info("atlas indexed: %d chromosomes, %d sites", len(out), len(atlas))
    return out


def _nearest_distance(
    index: dict[str, np.ndarray], chroms: np.ndarray, positions: np.ndarray
) -> np.ndarray:
    """Distance from each (chrom, position) to the nearest indexed site."""
    dist = np.full(positions.shape, np.inf, dtype=float)
    for chrom in np.unique(chroms):
        ref = index.get(str(chrom))
        mask = chroms == chrom
        if ref is None or ref.size == 0:
            continue
        query = positions[mask]
        pos = np.searchsorted(ref, query)
        left = ref[np.clip(pos - 1, 0, ref.size - 1)]
        right = ref[np.clip(pos, 0, ref.size - 1)]
        dist[mask] = np.minimum(np.abs(query - left), np.abs(query - right))
    return dist


# ---------------------------------------------------------------------------
# 1. atlas null control
# ---------------------------------------------------------------------------


def atlas_null_control(
    pasbed: Path,
    atlas_bed: Path,
    cutoffs: Sequence[int] = (10, 25, 50, 100, 200),
    n_shuffles: int = 5,
    shift_min: int = 5_000,
    shift_max: int = 50_000,
    seed: int = 0,
    atlas_index: Optional[dict[str, np.ndarray]] = None,
) -> pd.DataFrame:
    """Observed vs null atlas precision for one run's called PAS.

    Two nulls, both preserving the per-chromosome PAS count:

    ``shift``
        each called PAS displaced by +/- U(shift_min, shift_max) bp on its own
        chromosome. Controls for the fact that atlas sites are concentrated in
        genes/3'UTRs — a shifted PAS stays in roughly the same neighbourhood.
    ``uniform``
        positions drawn uniformly across the chromosome's observed span.
        The weaker, genome-wide "better than anywhere" control.

    ``excess_precision = observed - null_shift`` is the number that actually
    says whether peak calling put PAS where polyadenylation happens.
    """
    index = atlas_index if atlas_index is not None else _load_atlas_index(atlas_bed)
    pas = _read_bed_positions(pasbed)
    pas = pas[pas["chrom"].isin(index.keys())]
    if pas.empty:
        return pd.DataFrame()

    chroms = pas["chrom"].to_numpy()
    mids = pas["mid"].to_numpy(np.int64)
    obs_dist = _nearest_distance(index, chroms, mids)

    rng = np.random.default_rng(seed)
    # per-chromosome span, used by the uniform null
    span = pas.groupby("chrom")["mid"].agg(["min", "max"])

    null_shift: list[np.ndarray] = []
    null_unif: list[np.ndarray] = []
    for _ in range(max(1, n_shuffles)):
        offset = rng.integers(shift_min, shift_max, size=mids.size)
        offset *= rng.choice([-1, 1], size=mids.size)
        shifted = np.maximum(mids + offset, 0)
        null_shift.append(_nearest_distance(index, chroms, shifted))

        unif = np.empty_like(mids)
        for chrom in np.unique(chroms):
            mask = chroms == chrom
            lo, hi = int(span.loc[chrom, "min"]), int(span.loc[chrom, "max"])
            hi = max(hi, lo + 1)
            unif[mask] = rng.integers(lo, hi, size=int(mask.sum()))
        null_unif.append(_nearest_distance(index, chroms, unif))

    rows = []
    for cutoff in cutoffs:
        obs = float(np.mean(obs_dist <= cutoff))
        shift = float(np.mean([np.mean(d <= cutoff) for d in null_shift]))
        unif = float(np.mean([np.mean(d <= cutoff) for d in null_unif]))
        rows.append(
            {
                "cutoff_bp": cutoff,
                "n_pas": int(mids.size),
                "precision_observed": round(obs, 4),
                "precision_null_shift": round(shift, 4),
                "precision_null_uniform": round(unif, 4),
                "excess_over_shift": round(obs - shift, 4),
                "excess_over_uniform": round(obs - unif, 4),
                # how much of the *available* headroom the caller actually used
                "excess_ratio_shift": (
                    round((obs - shift) / (1.0 - shift), 4) if shift < 1 else np.nan
                ),
            }
        )
    out = pd.DataFrame(rows)
    out["median_distance_bp"] = float(np.median(obs_dist))
    out["median_distance_null_shift_bp"] = float(
        np.median(np.concatenate(null_shift))
    )
    return out


def atlas_null_control_multi(
    run_dirs: dict[str, Path],
    atlas_bed: Path,
    pas_filename: str = "pasbed.bed",
    **kwargs,
) -> pd.DataFrame:
    """``atlas_null_control`` over several runs, sharing one atlas index."""
    index = _load_atlas_index(Path(atlas_bed))
    frames = []
    for label, run_dir in run_dirs.items():
        path = Path(run_dir) / pas_filename
        if not path.exists():
            logger.warning("missing %s", path)
            continue
        df = atlas_null_control(path, atlas_bed, atlas_index=index, **kwargs)
        if df.empty:
            continue
        df.insert(0, "run", label)
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# ---------------------------------------------------------------------------
# 2. PAS yield + cross-strategy overlap
# ---------------------------------------------------------------------------


def pas_yield_table(run_dirs: dict[str, Path]) -> pd.DataFrame:
    """PAS count, peak width, gene coverage and multi-PAS genes per run."""
    rows = []
    for label, run_dir in run_dirs.items():
        run_dir = Path(run_dir)
        row: dict = {"run": label}
        pasbed = run_dir / "pasbed.bed"
        if pasbed.exists():
            pas = _read_bed_positions(pasbed)
            width = (pas["end"] - pas["start"]).to_numpy()
            row.update(
                n_pas=int(len(pas)),
                median_peak_width=float(np.median(width)),
                mean_peak_width=round(float(np.mean(width)), 1),
                p90_peak_width=float(np.percentile(width, 90)),
            )
        pas_gene = run_dir / "pas_gene.tsv"
        if pas_gene.exists():
            pg = pd.read_csv(
                pas_gene, sep="\t", header=None, names=["pas_id", "gene_id"]
            )
            per_gene = pg.groupby("gene_id")["pas_id"].nunique()
            row.update(
                n_pas_assigned=int(pg["pas_id"].nunique()),
                n_genes=int(per_gene.size),
                n_multi_pas_genes=int((per_gene >= 2).sum()),
                frac_multi_pas_genes=(
                    round(float((per_gene >= 2).mean()), 4) if per_gene.size else np.nan
                ),
                median_pas_per_gene=float(per_gene.median()) if per_gene.size else np.nan,
            )
        bench = _read_json(run_dir / "benchmark_vs_polyasite_v3.json")
        cut = (bench.get("cutoffs") or {}).get("50") or {}
        if cut:
            row["atlas_precision_50"] = cut.get("precision")
            row["atlas_recall_50"] = cut.get("recall")
        rows.append(row)
    return pd.DataFrame(rows)


def pas_set_overlap(run_dirs: dict[str, Path], tolerance: int = 50) -> pd.DataFrame:
    """Pairwise PAS-set agreement (Jaccard within ``tolerance`` bp).

    Answers: do two peak strategies find the *same* polyadenylation sites, or
    merely equally atlas-adjacent ones?
    """
    sets: dict[str, pd.DataFrame] = {}
    for label, run_dir in run_dirs.items():
        path = Path(run_dir) / "pasbed.bed"
        if path.exists():
            sets[label] = _read_bed_positions(path)
    labels = sorted(sets)
    rows = []
    for i, a in enumerate(labels):
        for b in labels[i + 1 :]:
            da, db = sets[a], sets[b]
            idx_b = {
                c: np.sort(g["mid"].to_numpy(np.int64))
                for c, g in db.groupby("chrom", sort=False)
            }
            idx_a = {
                c: np.sort(g["mid"].to_numpy(np.int64))
                for c, g in da.groupby("chrom", sort=False)
            }
            d_ab = _nearest_distance(idx_b, da["chrom"].to_numpy(), da["mid"].to_numpy())
            d_ba = _nearest_distance(idx_a, db["chrom"].to_numpy(), db["mid"].to_numpy())
            matched_a = int(np.sum(d_ab <= tolerance))
            matched_b = int(np.sum(d_ba <= tolerance))
            union = len(da) + len(db) - matched_a
            rows.append(
                {
                    "run_a": a,
                    "run_b": b,
                    "n_a": len(da),
                    "n_b": len(db),
                    "matched_a_in_b": matched_a,
                    "matched_b_in_a": matched_b,
                    "frac_a_recovered": round(matched_a / max(len(da), 1), 4),
                    "frac_b_recovered": round(matched_b / max(len(db), 1), 4),
                    "jaccard": round(matched_a / max(union, 1), 4),
                    "tolerance_bp": tolerance,
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 3-4. clustering concordance + cluster -> cell type mapping
# ---------------------------------------------------------------------------


def _load_clusters(h5ad_path: Path) -> Optional[pd.DataFrame]:
    """obs frame with a normalised ``leiden`` column, plus X_lsi if present."""
    try:
        import anndata as ad
    except ImportError:  # pragma: no cover
        logger.error("anndata not available")
        return None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata = ad.read_h5ad(h5ad_path)
    except (OSError, KeyError, ValueError) as exc:
        logger.warning("cannot read %s: %s", h5ad_path, exc)
        return None
    obs = adata.obs.copy()
    col = next(
        (c for c in ("leiden", "louvain", "cluster", "clusters") if c in obs.columns),
        None,
    )
    if col is None:
        return None
    obs["leiden"] = obs[col].astype(str)
    frame = obs[["leiden"]].copy()
    for extra in ("celltype", "n_genes", "total_counts", "canonical_cluster"):
        if extra in obs.columns:
            frame[extra] = obs[extra].values
    if "X_lsi" in adata.obsm:
        frame.attrs["X_lsi"] = np.asarray(adata.obsm["X_lsi"])
    if "X_umap" in adata.obsm:
        frame.attrs["X_umap"] = np.asarray(adata.obsm["X_umap"])
    return frame


def load_celltype_labels(gex_dir: Path) -> dict[str, pd.Series]:
    """dataset_id -> Series(cell -> celltype) from the B2 *_pas_labeled.h5ad."""
    labels: dict[str, pd.Series] = {}
    gex_dir = Path(gex_dir)
    if not gex_dir.exists():
        return labels
    for path in sorted(gex_dir.glob("*_pas_labeled.h5ad")):
        dataset = path.name.replace("_pas_labeled.h5ad", "")
        frame = _load_clusters(path)
        if frame is None or "celltype" not in frame.columns:
            continue
        labels[dataset] = frame["celltype"].astype(str)
    return labels


def cluster_concordance(
    branch_dirs: dict[str, Path],
    reference: str,
    celltype_labels: Optional[dict[str, pd.Series]] = None,
    clustering_subdir: str = "07_clustering",
    max_cells_for_silhouette: int = 4000,
) -> pd.DataFrame:
    """Per branch x dataset: n_clusters, ARI/AMI vs reference and vs cell type.

    This is the single table that makes "did this parameter change anything?"
    answerable: if ARI-vs-reference is ~1.0 the branch is a relabelling, and if
    ARI-vs-celltype does not move, the change is biologically inert.
    """
    from sklearn.metrics import (
        adjusted_mutual_info_score,
        adjusted_rand_score,
        silhouette_score,
    )

    celltype_labels = celltype_labels or {}
    cache: dict[tuple[str, str], Optional[pd.DataFrame]] = {}

    def get(branch: str, dataset: str) -> Optional[pd.DataFrame]:
        key = (branch, dataset)
        if key not in cache:
            path = Path(branch_dirs[branch]) / clustering_subdir / dataset / "clusters.h5ad"
            cache[key] = _load_clusters(path) if path.exists() else None
        return cache[key]

    ref_dir = Path(branch_dirs[reference]) / clustering_subdir
    datasets = sorted(p.name for p in ref_dir.iterdir() if p.is_dir()) if ref_dir.exists() else []

    rows = []
    for branch in sorted(branch_dirs):
        for dataset in datasets:
            frame = get(branch, dataset)
            if frame is None:
                rows.append({"branch": branch, "dataset": dataset, "status": "missing"})
                continue
            row: dict = {
                "branch": branch,
                "dataset": dataset,
                "status": "ok",
                "n_cells": int(len(frame)),
                "n_clusters": int(frame["leiden"].nunique()),
                "largest_cluster_frac": round(
                    float(frame["leiden"].value_counts(normalize=True).iloc[0]), 4
                ),
            }
            # ARI/AMI against the reference branch on shared cells
            ref_frame = get(reference, dataset)
            if ref_frame is not None and branch != reference:
                shared = frame.index.intersection(ref_frame.index)
                if len(shared) >= 10:
                    a = frame.loc[shared, "leiden"].to_numpy()
                    b = ref_frame.loc[shared, "leiden"].to_numpy()
                    row["n_shared_with_ref"] = int(len(shared))
                    row["ARI_vs_reference"] = round(adjusted_rand_score(b, a), 4)
                    row["AMI_vs_reference"] = round(
                        adjusted_mutual_info_score(b, a), 4
                    )
            elif branch == reference:
                row["ARI_vs_reference"] = 1.0
                row["AMI_vs_reference"] = 1.0

            # ARI/AMI against GEX cell types
            ct = celltype_labels.get(dataset)
            if ct is not None:
                shared = frame.index.intersection(ct.index)
                if len(shared) >= 10:
                    a = frame.loc[shared, "leiden"].to_numpy()
                    b = ct.loc[shared].to_numpy()
                    row["n_cells_with_celltype"] = int(len(shared))
                    row["n_celltypes"] = int(pd.unique(b).size)
                    row["ARI_vs_celltype"] = round(adjusted_rand_score(b, a), 4)
                    row["AMI_vs_celltype"] = round(
                        adjusted_mutual_info_score(b, a), 4
                    )

            # internal structure: silhouette of the clustering in LSI space
            lsi = frame.attrs.get("X_lsi")
            if lsi is not None and row["n_clusters"] > 1:
                take = np.arange(len(frame))
                if len(take) > max_cells_for_silhouette:
                    take = np.random.default_rng(0).choice(
                        len(frame), max_cells_for_silhouette, replace=False
                    )
                try:
                    row["silhouette_lsi"] = round(
                        float(
                            silhouette_score(
                                lsi[take], frame["leiden"].to_numpy()[take]
                            )
                        ),
                        4,
                    )
                except ValueError:
                    pass
            rows.append(row)
    return pd.DataFrame(rows)


def summarize_cluster_concordance(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse the per-dataset concordance table to one row per branch."""
    ok = df[df.get("status", "ok") == "ok"]
    if ok.empty:
        return pd.DataFrame()
    agg = {
        "n_clusters": "mean",
        "ARI_vs_reference": "mean",
        "AMI_vs_reference": "mean",
        "ARI_vs_celltype": "mean",
        "AMI_vs_celltype": "mean",
        "silhouette_lsi": "mean",
        "largest_cluster_frac": "mean",
    }
    agg = {k: v for k, v in agg.items() if k in ok.columns}
    out = ok.groupby("branch").agg(agg).round(4)
    out["n_datasets"] = ok.groupby("branch").size()
    return out.reset_index()


def cluster_celltype_contingency(
    branch_dir: Path,
    dataset: str,
    celltype_labels: dict[str, pd.Series],
    clustering_subdir: str = "07_clustering",
) -> pd.DataFrame:
    """Long-form cluster -> cell type flow table (the Sankey input)."""
    path = Path(branch_dir) / clustering_subdir / dataset / "clusters.h5ad"
    frame = _load_clusters(path) if path.exists() else None
    ct = celltype_labels.get(dataset)
    if frame is None or ct is None:
        return pd.DataFrame()
    shared = frame.index.intersection(ct.index)
    if not len(shared):
        return pd.DataFrame()
    joined = pd.DataFrame(
        {"cluster": frame.loc[shared, "leiden"].to_numpy(), "celltype": ct.loc[shared].to_numpy()}
    )
    out = (
        joined.groupby(["cluster", "celltype"]).size().reset_index(name="n_cells")
    )
    total = out.groupby("cluster")["n_cells"].transform("sum")
    out["frac_of_cluster"] = (out["n_cells"] / total).round(4)
    out.insert(0, "dataset", dataset)
    return out.sort_values(["cluster", "n_cells"], ascending=[True, False])


def celltype_assignment_shift(
    branch_dirs: dict[str, Path],
    reference: str,
    celltype_labels: dict[str, pd.Series],
    clustering_subdir: str = "07_clustering",
) -> pd.DataFrame:
    """How much each branch changes the *cluster -> cell type* call per cell.

    For every cell, the cell type is the majority cell type of its cluster.
    ``frac_cells_reassigned`` is then a directly interpretable "this parameter
    changed the cell-type call for N% of cells".
    """

    def majority_map(branch: str, dataset: str) -> Optional[pd.Series]:
        path = Path(branch_dirs[branch]) / clustering_subdir / dataset / "clusters.h5ad"
        frame = _load_clusters(path) if path.exists() else None
        ct = celltype_labels.get(dataset)
        if frame is None or ct is None:
            return None
        shared = frame.index.intersection(ct.index)
        if not len(shared):
            return None
        joined = pd.DataFrame(
            {
                "cluster": frame.loc[shared, "leiden"].to_numpy(),
                "celltype": ct.loc[shared].to_numpy(),
            },
            index=shared,
        )
        winner = joined.groupby("cluster")["celltype"].agg(
            lambda s: s.value_counts().idxmax()
        )
        return joined["cluster"].map(winner)

    ref_dir = Path(branch_dirs[reference]) / clustering_subdir
    datasets = sorted(p.name for p in ref_dir.iterdir() if p.is_dir()) if ref_dir.exists() else []
    rows = []
    for branch in sorted(branch_dirs):
        if branch == reference:
            continue
        for dataset in datasets:
            a, b = majority_map(branch, dataset), majority_map(reference, dataset)
            if a is None or b is None:
                continue
            shared = a.index.intersection(b.index)
            if not len(shared):
                continue
            changed = (a.loc[shared].to_numpy() != b.loc[shared].to_numpy()).mean()
            rows.append(
                {
                    "branch": branch,
                    "dataset": dataset,
                    "n_cells": int(len(shared)),
                    "frac_cells_reassigned": round(float(changed), 4),
                }
            )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 5. switch biology: the depth confound
# ---------------------------------------------------------------------------


def switch_depth_confound(switch_dir: Path, stage_order: Sequence[str]) -> pd.DataFrame:
    """Per cell type x stage: PDUI, its zero-fraction, cell count and depth.

    The reported "global 3'UTR shortening" is a slope through per-stage mean
    PDUI. PDUI is 0 whenever a gene's distal PAS has no reads in that cell, so
    a stage with fewer reads per cell has mechanically lower mean PDUI. This
    table puts the slope and the thing that could fake it side by side.
    """
    switch_dir = Path(switch_dir)
    length_dir = switch_dir / "length"
    if not length_dir.exists():
        return pd.DataFrame()
    rows = []
    for ct_dir in sorted(p for p in length_dir.iterdir() if p.is_dir()):
        pdui_path = ct_dir / "classic" / "pdui_classic.tsv"
        if not pdui_path.exists():
            continue
        try:
            pdui = pd.read_csv(
                pdui_path,
                sep="\t",
                usecols=[
                    "gene_id",
                    "cell",
                    "pdui",
                    "proximal_reads",
                    "distal_reads",
                    "total_reads",
                    "cluster",
                ],
            )
        except (OSError, ValueError) as exc:
            logger.warning("cannot read %s: %s", pdui_path, exc)
            continue
        for stage, grp in pdui.groupby("cluster"):
            informative = grp["total_reads"] > 0
            rows.append(
                {
                    "celltype": ct_dir.name,
                    "stage": str(stage),
                    "n_cells": int(grp["cell"].nunique()),
                    "n_genes": int(grp["gene_id"].nunique()),
                    "n_gene_cell_pairs": int(len(grp)),
                    "mean_pdui": round(float(grp["pdui"].mean()), 5),
                    # the honest version: PDUI only where the gene was seen
                    "mean_pdui_informative": (
                        round(float(grp.loc[informative, "pdui"].mean()), 5)
                        if informative.any()
                        else np.nan
                    ),
                    "frac_pdui_zero": round(float((grp["pdui"] == 0).mean()), 5),
                    "frac_uninformative": round(float((~informative).mean()), 5),
                    "mean_total_reads": round(float(grp["total_reads"].mean()), 4),
                    "median_total_reads_informative": (
                        float(grp.loc[informative, "total_reads"].median())
                        if informative.any()
                        else np.nan
                    ),
                }
            )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    order = {s: i for i, s in enumerate(stage_order)}
    df["stage_index"] = df["stage"].map(order)
    return df.sort_values(["celltype", "stage_index"]).reset_index(drop=True)


def recompute_trends(depth_df: pd.DataFrame, value_col: str = "mean_pdui") -> pd.DataFrame:
    """Refit the stage trend, recording how many stages it was fit on.

    The shipped ``length_trend.json`` reports a Spearman rho for cell types with
    only two stages, where rho is always exactly +/-1. This refits and flags
    that, and additionally fits the *informative-only* PDUI so the slope can be
    compared with and without the detection-rate confound.
    """
    from scipy import stats

    rows = []
    for celltype, grp in depth_df.dropna(subset=["stage_index"]).groupby("celltype"):
        grp = grp.sort_values("stage_index")
        row: dict = {
            "celltype": celltype,
            "n_stages": int(len(grp)),
            "stages": ",".join(grp["stage"]),
            "total_cells": int(grp["n_cells"].sum()),
            "interpretable": bool(len(grp) >= MIN_STAGES_FOR_TREND),
        }
        for col, tag in ((value_col, ""), ("mean_pdui_informative", "_informative")):
            if col not in grp.columns:
                continue
            y = grp[col].to_numpy(float)
            x = grp["stage_index"].to_numpy(float)
            ok = np.isfinite(y)
            if ok.sum() < 2:
                continue
            fit = stats.linregress(x[ok], y[ok])
            row[f"slope{tag}"] = round(float(fit.slope), 6)
            row[f"pvalue{tag}"] = (
                round(float(fit.pvalue), 6) if ok.sum() > 2 else np.nan
            )
            if ok.sum() >= MIN_STAGES_FOR_TREND:
                rho = stats.spearmanr(x[ok], y[ok])
                row[f"spearman{tag}"] = round(float(rho.statistic), 4)
            else:
                row[f"spearman{tag}"] = np.nan  # undefined with 2 points
        # does the slope just track sequencing depth?
        depth = grp["mean_total_reads"].to_numpy(float)
        x = grp["stage_index"].to_numpy(float)
        if len(grp) >= 2 and np.isfinite(depth).all():
            row["depth_slope"] = round(float(stats.linregress(x, depth).slope), 6)
        if len(grp) >= MIN_STAGES_FOR_TREND:
            row["corr_pdui_depth"] = round(
                float(np.corrcoef(grp[value_col].to_numpy(float), depth)[0, 1]), 4
            )
        rows.append(row)
    return pd.DataFrame(rows)


def depth_confound_summary(trend_df: pd.DataFrame) -> dict:
    """Headline numbers for the depth-confound control."""
    if trend_df.empty:
        return {}
    usable = trend_df[trend_df["interpretable"]]
    out = {
        "n_celltypes": int(len(trend_df)),
        "n_celltypes_with_ge3_stages": int(len(usable)),
        "n_celltypes_with_2_stages_only": int((trend_df["n_stages"] == 2).sum()),
        "pct_trends_uninterpretable": round(
            100.0 * float((trend_df["n_stages"] < MIN_STAGES_FOR_TREND).mean()), 1
        ),
    }
    if "slope" in trend_df.columns:
        out["median_slope"] = round(float(trend_df["slope"].median()), 6)
        out["frac_decreasing"] = round(float((trend_df["slope"] < 0).mean()), 4)
    if "slope_informative" in trend_df.columns:
        both = trend_df.dropna(subset=["slope", "slope_informative"])
        if not both.empty:
            out["median_slope_informative"] = round(
                float(both["slope_informative"].median()), 6
            )
            out["frac_decreasing_informative"] = round(
                float((both["slope_informative"] < 0).mean()), 4
            )
            out["sign_flips_when_depth_controlled"] = int(
                (np.sign(both["slope"]) != np.sign(both["slope_informative"])).sum()
            )
    if "corr_pdui_depth" in trend_df.columns:
        corr = trend_df["corr_pdui_depth"].dropna()
        if not corr.empty:
            out["median_corr_pdui_vs_depth"] = round(float(corr.median()), 4)
    return out


# ---------------------------------------------------------------------------
# 6. fisher power / pseudoreplication check
# ---------------------------------------------------------------------------


def fisher_power_check(switch_dir: Path, fdr: float = 0.05) -> pd.DataFrame:
    """Significant-hit counts per (cell type, contrast) with the cell counts.

    Fisher's exact test here is applied to *read* counts pooled over cells, so
    every extra cell adds pseudo-replicated reads and inflates significance.
    The signature is a strong positive relation between n_cells and n_significant
    — which this table lets you plot directly.
    """
    switch_dir = Path(switch_dir)
    diff_dir = switch_dir / "diff"
    if not diff_dir.exists():
        return pd.DataFrame()
    rows = []
    for ct_dir in sorted(p for p in diff_dir.iterdir() if p.is_dir()):
        for tsv in sorted((ct_dir / "fisher" / "differential").glob("fisher_*.tsv")):
            try:
                df = pd.read_csv(
                    tsv,
                    sep="\t",
                    usecols=lambda c: c
                    in {
                        "qvalue",
                        "n_cells",
                        "n_cells_cluster1",
                        "n_cells_cluster2",
                        "gene_id",
                        "delta_proportion",
                    },
                )
            except (OSError, ValueError) as exc:
                logger.warning("cannot read %s: %s", tsv, exc)
                continue
            if df.empty or "qvalue" not in df.columns:
                continue
            sig = df["qvalue"] <= fdr
            n1 = int(df["n_cells_cluster1"].max()) if "n_cells_cluster1" in df else np.nan
            n2 = int(df["n_cells_cluster2"].max()) if "n_cells_cluster2" in df else np.nan
            rows.append(
                {
                    "celltype": ct_dir.name,
                    "contrast": tsv.stem.replace("fisher_", ""),
                    "n_tests": int(len(df)),
                    "n_significant": int(sig.sum()),
                    "frac_significant": round(float(sig.mean()), 4),
                    "n_cells_cluster1": n1,
                    "n_cells_cluster2": n2,
                    "n_cells_total": (n1 + n2) if np.isfinite(n1) and np.isfinite(n2) else np.nan,
                    "n_genes_significant": (
                        int(df.loc[sig, "gene_id"].nunique())
                        if "gene_id" in df.columns
                        else np.nan
                    ),
                    "median_abs_delta_sig": (
                        round(float(df.loc[sig, "delta_proportion"].abs().median()), 4)
                        if "delta_proportion" in df.columns and sig.any()
                        else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


def pseudoreplication_signature(fisher_df: pd.DataFrame) -> dict:
    """Correlation between cell count and significance — the smoking gun."""
    from scipy import stats

    df = fisher_df.dropna(subset=["n_cells_total"])
    df = df[df["n_tests"] > 0]
    if len(df) < 4:
        return {"n_points": int(len(df))}
    out = {"n_points": int(len(df))}
    for target in ("n_significant", "frac_significant"):
        rho = stats.spearmanr(df["n_cells_total"], df[target])
        out[f"spearman_ncells_vs_{target}"] = round(float(rho.statistic), 4)
        out[f"pvalue_ncells_vs_{target}"] = float(rho.pvalue)
    out["median_frac_significant"] = round(float(df["frac_significant"].median()), 4)
    out["max_frac_significant"] = round(float(df["frac_significant"].max()), 4)
    return out


# ---------------------------------------------------------------------------
# 7. recurrent vs private switches
# ---------------------------------------------------------------------------


def recurrent_switches(
    switch_dir: Path,
    fdr: float = 0.05,
    min_abs_delta: float = 0.1,
    contrast: Optional[str] = None,
) -> pd.DataFrame:
    """Genes with a significant APA switch, counted across cell types.

    A gene called in one cell type only is as likely to be a power artefact as
    biology; genes recurring across many independent cell types are the ones
    worth a figure.
    """
    switch_dir = Path(switch_dir)
    diff_dir = switch_dir / "diff"
    if not diff_dir.exists():
        return pd.DataFrame()
    records = []
    for ct_dir in sorted(p for p in diff_dir.iterdir() if p.is_dir()):
        for tsv in sorted((ct_dir / "fisher" / "differential").glob("fisher_*.tsv")):
            name = tsv.stem.replace("fisher_", "")
            if contrast and name != contrast:
                continue
            try:
                df = pd.read_csv(
                    tsv, sep="\t", usecols=["gene_id", "qvalue", "delta_proportion"]
                )
            except (OSError, ValueError):
                continue
            hit = df[
                (df["qvalue"] <= fdr)
                & (df["delta_proportion"].abs() >= min_abs_delta)
            ]
            for gene, grp in hit.groupby("gene_id"):
                records.append(
                    {
                        "gene_id": gene,
                        "celltype": ct_dir.name,
                        "contrast": name,
                        "max_abs_delta": float(grp["delta_proportion"].abs().max()),
                        "direction": (
                            "shorter"
                            if grp.loc[grp["delta_proportion"].abs().idxmax(), "delta_proportion"] > 0
                            else "longer"
                        ),
                    }
                )
    if not records:
        return pd.DataFrame()
    long = pd.DataFrame(records)
    out = (
        long.groupby("gene_id")
        .agg(
            n_celltypes=("celltype", "nunique"),
            n_contrasts=("contrast", "nunique"),
            celltypes=("celltype", lambda s: ",".join(sorted(set(s))[:6])),
            max_abs_delta=("max_abs_delta", "max"),
            n_shorter=("direction", lambda s: int((s == "shorter").sum())),
            n_longer=("direction", lambda s: int((s == "longer").sum())),
        )
        .reset_index()
    )
    out["consistent_direction"] = (out["n_shorter"] == 0) | (out["n_longer"] == 0)
    return out.sort_values(
        ["n_celltypes", "max_abs_delta"], ascending=False
    ).reset_index(drop=True)


# ---------------------------------------------------------------------------
# 8. filter-effect deltas (phase 3)
# ---------------------------------------------------------------------------


def filter_effect_deltas(
    scenario_dirs: dict[str, Path],
    baseline: str = "baseline",
    celltype_labels: Optional[dict[str, pd.Series]] = None,
    stage_order: Sequence[str] = ("Normal", "StageI", "IVprimary", "Met"),
) -> dict[str, pd.DataFrame]:
    """Every downstream consequence of masking one PAS filter, vs baseline.

    Returns tables keyed ``pas``, ``clustering``, ``celltype_shift``, ``trend``,
    ``fisher`` — the same five downstream steps, so one scenario's whole blast
    radius reads off a single row.
    """
    present = {k: Path(v) for k, v in scenario_dirs.items() if Path(v).exists()}
    out: dict[str, pd.DataFrame] = {}
    if not present:
        return out

    out["pas"] = pas_yield_table(present)
    if baseline in present and "n_pas" in out["pas"].columns:
        base = out["pas"].set_index("run")
        if baseline in base.index:
            ref = base.loc[baseline, "n_pas"]
            out["pas"]["delta_n_pas"] = out["pas"]["n_pas"] - ref
            out["pas"]["frac_pas_retained"] = (out["pas"]["n_pas"] / ref).round(4)

    if baseline in present:
        out["clustering"] = cluster_concordance(
            present, reference=baseline, celltype_labels=celltype_labels
        )
        if celltype_labels:
            out["celltype_shift"] = celltype_assignment_shift(
                present, reference=baseline, celltype_labels=celltype_labels
            )

    trend_frames, fisher_frames = [], []
    for label, run_dir in present.items():
        switch_dir = run_dir / "B3_switch"
        if not switch_dir.exists():
            continue
        depth = switch_depth_confound(switch_dir, stage_order)
        if not depth.empty:
            trends = recompute_trends(depth)
            trends.insert(0, "scenario", label)
            trend_frames.append(trends)
        fisher = fisher_power_check(switch_dir)
        if not fisher.empty:
            fisher.insert(0, "scenario", label)
            fisher_frames.append(fisher)
    if trend_frames:
        out["trend"] = pd.concat(trend_frames, ignore_index=True)
    if fisher_frames:
        out["fisher"] = pd.concat(fisher_frames, ignore_index=True)
    return out


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------


def _write(df: pd.DataFrame, path: Path, label: str) -> None:
    if df is None or df.empty:
        logger.warning("no rows for %s", label)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    logger.info("wrote %s (%d rows)", path, len(df))


def run_all(
    sweep_root: Path,
    out_dir: Path,
    atlas_bed: Optional[Path] = None,
    filter_effect_root: Optional[Path] = None,
    stage_order: Sequence[str] = ("Normal", "StageI", "IVprimary", "Met"),
    n_shuffles: int = 3,
) -> dict:
    """Run every cross-experiment analysis available for this sweep."""
    sweep_root = Path(sweep_root)
    runs = sweep_root / "runs"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary: dict = {"sweep_root": str(sweep_root)}

    cohort = runs / "B1_cohort_full"
    grid_dirs = {
        p.name: p for p in sorted((runs / "grid").glob("*")) if p.is_dir()
    }
    reann_dirs = {
        p.name: p for p in sorted((runs / "reannotate").glob("*")) if p.is_dir()
    }

    # --- cell type labels (shared by every clustering comparison) -----------
    celltypes = load_celltype_labels(cohort / "B2_gex_celltyping")
    summary["n_datasets_with_celltypes"] = len(celltypes)

    # --- strategy axis ------------------------------------------------------
    if grid_dirs:
        _write(pas_yield_table(grid_dirs), out_dir / "strategy_pas_yield.tsv", "yield")
        _write(
            pas_set_overlap(grid_dirs), out_dir / "strategy_pas_overlap.tsv", "overlap"
        )
        if atlas_bed and Path(atlas_bed).exists():
            null = atlas_null_control_multi(
                grid_dirs, Path(atlas_bed), n_shuffles=n_shuffles
            )
            _write(null, out_dir / "strategy_atlas_null_control.tsv", "atlas null")
            if not null.empty:
                at50 = null[null["cutoff_bp"] == 50]
                summary["atlas_null_at_50bp"] = {
                    r["run"]: {
                        "observed": r["precision_observed"],
                        "null_shift": r["precision_null_shift"],
                        "excess": r["excess_over_shift"],
                    }
                    for _, r in at50.iterrows()
                }
        strategy_clusters = cluster_concordance(
            {**grid_dirs}, reference="lg_annotate" if "lg_annotate" in grid_dirs else sorted(grid_dirs)[0],
            celltype_labels=celltypes,
        )
        _write(strategy_clusters, out_dir / "strategy_clustering.tsv", "strategy clusters")
        _write(
            summarize_cluster_concordance(strategy_clusters),
            out_dir / "strategy_clustering_summary.tsv",
            "strategy cluster summary",
        )

    # --- trim + clustering axes (reannotate branches vs the cohort) ---------
    if reann_dirs:
        branch_dirs = {"BASE_cohort": cohort, **reann_dirs}
        _write(
            pas_yield_table(branch_dirs), out_dir / "branch_pas_yield.tsv", "branch yield"
        )
        conc = cluster_concordance(
            branch_dirs, reference="BASE_cohort", celltype_labels=celltypes
        )
        _write(conc, out_dir / "branch_clustering.tsv", "branch clusters")
        _write(
            summarize_cluster_concordance(conc),
            out_dir / "branch_clustering_summary.tsv",
            "branch cluster summary",
        )
        shift = celltype_assignment_shift(
            branch_dirs, reference="BASE_cohort", celltype_labels=celltypes
        )
        _write(shift, out_dir / "branch_celltype_shift.tsv", "celltype shift")
        if not shift.empty:
            summary["celltype_reassignment_by_branch"] = (
                shift.groupby("branch")["frac_cells_reassigned"]
                .mean()
                .round(4)
                .to_dict()
            )

    # --- switch biology -----------------------------------------------------
    switch_dir = cohort / "B3_switch"
    if switch_dir.exists():
        depth = switch_depth_confound(switch_dir, stage_order)
        _write(depth, out_dir / "switch_depth_by_stage.tsv", "depth")
        if not depth.empty:
            trends = recompute_trends(depth)
            _write(trends, out_dir / "switch_trend_recomputed.tsv", "trends")
            summary["depth_confound"] = depth_confound_summary(trends)
        fisher = fisher_power_check(switch_dir)
        _write(fisher, out_dir / "switch_fisher_power.tsv", "fisher")
        if not fisher.empty:
            summary["pseudoreplication"] = pseudoreplication_signature(fisher)
        rec = recurrent_switches(switch_dir)
        _write(rec, out_dir / "switch_recurrent_genes.tsv", "recurrent")
        if not rec.empty:
            summary["recurrent_switches"] = {
                "n_genes_total": int(len(rec)),
                "n_genes_in_1_celltype": int((rec["n_celltypes"] == 1).sum()),
                "n_genes_in_ge3_celltypes": int((rec["n_celltypes"] >= 3).sum()),
                "n_genes_in_ge5_celltypes": int((rec["n_celltypes"] >= 5).sum()),
                "frac_private": round(float((rec["n_celltypes"] == 1).mean()), 4),
                "top_genes": rec.head(15)[
                    ["gene_id", "n_celltypes", "max_abs_delta", "consistent_direction"]
                ].to_dict("records"),
            }

        # cluster -> cell type flows for the headline cohort (Sankey input)
        flows = []
        for dataset in sorted(celltypes):
            flow = cluster_celltype_contingency(cohort, dataset, celltypes)
            if not flow.empty:
                flows.append(flow)
        if flows:
            _write(
                pd.concat(flows, ignore_index=True),
                out_dir / "cohort_cluster_celltype_flows.tsv",
                "flows",
            )

    # --- filter effect (phase 3, may be incomplete) -------------------------
    if filter_effect_root:
        fe_runs = Path(filter_effect_root) / "runs"
        scenarios = {
            p.name: p for p in sorted(fe_runs.glob("*")) if p.is_dir()
        } if fe_runs.exists() else {}
        summary["filter_effect_scenarios_present"] = sorted(scenarios)
        if scenarios:
            tables = filter_effect_deltas(
                scenarios, celltype_labels=celltypes, stage_order=stage_order
            )
            for name, df in tables.items():
                _write(df, out_dir / f"filter_effect_{name}.tsv", f"fe {name}")

    with open(out_dir / "cross_experiment_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, default=str)
    logger.info("summary -> %s", out_dir / "cross_experiment_summary.json")
    return summary


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--atlas", type=Path, default=None)
    parser.add_argument("--filter-effect-root", type=Path, default=None)
    parser.add_argument("--stage-order", default="Normal,StageI,IVprimary,Met")
    parser.add_argument("--n-shuffles", type=int, default=3)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    run_all(
        sweep_root=args.sweep_root,
        out_dir=args.out,
        atlas_bed=args.atlas,
        filter_effect_root=args.filter_effect_root,
        stage_order=tuple(args.stage_order.split(",")),
        n_shuffles=args.n_shuffles,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
