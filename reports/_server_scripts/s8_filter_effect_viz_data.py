#!/usr/bin/env python
"""s8 — export UMAP + cluster/cell-type data for the four filter-effect scenarios.

Same contract as ``s6_clustering_viz_data.py`` but over the Phase-3 scenarios
rather than the sweep's clustering variants. Reads only real run artifacts:

  <scenario>/07_clustering/<dataset>/clusters.h5ad        obsm['X_umap'], obs['leiden']
  <scenario>/B2_gex_celltyping/<dataset>_gex_labeled.h5ad obs['celltype']

Cell-type calls are GEX-derived and therefore identical across scenarios; they
are joined **by barcode**, and the join rate is recorded per (scenario, dataset)
so a figure can never silently present a poorly-matched panel as complete.

Two outputs, both small enough to keep in the repo:

1. ``fe_cluster_celltype_counts.tsv`` — (scenario, dataset, leiden, celltype)
   contingency for all four scenarios and all six datasets. Drives every Sankey.
2. ``fe_umap_cells.tsv.gz`` — per-cell UMAP coordinates for all four scenarios
   on a shared reference dataset, so the panels are comparable like for like.

Usage::

    python s8_filter_effect_viz_data.py --out <dir> [--reference <GSM>]
"""
from __future__ import annotations

import argparse
import gzip
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import anndata as ad
import numpy as np
import pandas as pd

FE = Path("/mnt/ssd2/Laugney_Aligned/peakatail_experiments/FILTER_EFFECT_2026-08")
SCENARIOS = ["baseline", "atlas_filter", "annot_filter_3utr", "ip_filter"]

#: Unlabelled cells are kept and drawn grey, never dropped — dropping them makes
#: the cluster -> cell-type mapping look far cleaner than it is.
UNLABELLED = "unlabelled"


def _barcode(adata) -> np.ndarray:
    if "barcode" in adata.obs.columns:
        return adata.obs["barcode"].astype(str).to_numpy()
    return adata.obs_names.astype(str).to_numpy()


def celltype_map(scenario_dir: Path, dataset: str) -> dict[str, str]:
    """barcode -> celltype from the GEX-labelled h5ad, if present."""
    for pat in (f"{dataset}_gex_labeled.h5ad", f"{dataset}_pas_labeled.h5ad"):
        p = scenario_dir / "B2_gex_celltyping" / pat
        if not p.exists():
            continue
        try:
            a = ad.read_h5ad(p)
        except (OSError, KeyError, ValueError):
            continue
        if "celltype" not in a.obs.columns:
            continue
        return dict(zip(_barcode(a), a.obs["celltype"].astype(str)))
    return {}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--root", type=Path, default=FE)
    ap.add_argument("--reference", default=None,
                    help="dataset for the per-cell UMAP export (default: first common)")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    runs = a.root / "runs"

    # datasets present in every scenario, so panels are comparable
    per_scen = {}
    for s in SCENARIOS:
        d = runs / s / "07_clustering"
        per_scen[s] = {p.name for p in d.iterdir() if p.is_dir()} if d.exists() else set()
        print(f"{s}: {len(per_scen[s])} datasets", flush=True)
    common = sorted(set.intersection(*per_scen.values())) if all(per_scen.values()) else []
    print("common datasets:", common, flush=True)
    reference = a.reference or (common[0] if common else None)

    counts, umap_rows = [], []
    for s in SCENARIOS:
        for ds in sorted(per_scen[s]):
            h5 = runs / s / "07_clustering" / ds / "clusters.h5ad"
            if not h5.exists():
                continue
            try:
                ad_obj = ad.read_h5ad(h5)
            except (OSError, KeyError, ValueError) as exc:
                print(f"  skip {s}/{ds}: {exc}", flush=True)
                continue
            col = next((c for c in ("leiden", "louvain", "cluster")
                        if c in ad_obj.obs.columns), None)
            if col is None:
                continue
            bc = _barcode(ad_obj)
            leiden = ad_obj.obs[col].astype(str).to_numpy()
            cmap = celltype_map(runs / s, ds)
            ct = np.array([cmap.get(b, UNLABELLED) for b in bc])
            n_lab = int((ct != UNLABELLED).sum())

            df = pd.DataFrame({"scenario": s, "dataset": ds,
                               "leiden": leiden, "celltype": ct})
            g = df.groupby(["scenario", "dataset", "leiden", "celltype"]).size()
            counts.append(g.reset_index(name="n_cells"))
            print(f"  {s}/{ds}: {len(df)} cells, {n_lab} labelled "
                  f"({100*n_lab/max(len(df),1):.1f}%), {len(set(leiden))} clusters",
                  flush=True)

            if ds == reference and "X_umap" in ad_obj.obsm:
                u = np.asarray(ad_obj.obsm["X_umap"])
                umap_rows.append(pd.DataFrame({
                    "scenario": s, "dataset": ds, "barcode": bc,
                    "umap_x": np.round(u[:, 0], 4), "umap_y": np.round(u[:, 1], 4),
                    "leiden": leiden, "celltype": ct,
                }))

    if counts:
        out = pd.concat(counts, ignore_index=True)
        out.to_csv(a.out / "fe_cluster_celltype_counts.tsv", sep="\t", index=False)
        print("wrote fe_cluster_celltype_counts.tsv", len(out), "rows", flush=True)
    if umap_rows:
        out = pd.concat(umap_rows, ignore_index=True)
        with gzip.open(a.out / "fe_umap_cells.tsv.gz", "wt") as fh:
            out.to_csv(fh, sep="\t", index=False)
        print("wrote fe_umap_cells.tsv.gz", len(out), "rows,",
              f"reference={reference}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
