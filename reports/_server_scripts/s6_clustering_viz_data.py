#!/usr/bin/env python
"""Export the real single-cell clustering data needed for the UMAP / Sankey /
cluster-size figures, for EVERY experiment in the corrected sweep.

Reads only real run artifacts on the analysis host:

  <run>/07_clustering/<dataset>/clusters.h5ad     -- obsm['X_umap'], obs['leiden']
  B1_cohort_full/B2_gex_celltyping/<ds>_pas_labeled.h5ad -- obs['celltype']

Cell-type labels exist only under the cohort run, but every experiment clusters
the same cells, so labels are joined onto the other runs **by barcode**.  The
join rate is recorded per (experiment, dataset) so a figure can never silently
present a poorly-matched panel as if it were complete.

Two output tiers, chosen so that "cover every experiment" does not mean reading
19 x 17 h5ads worth of coordinates:

1. ``cluster_celltype_counts.tsv`` -- the (experiment, dataset, leiden, celltype)
   contingency for ALL experiments and ALL datasets.  Small, and enough to drive
   every Sankey, every cluster-size bar and the composition panels.
2. ``umap_cells.tsv.gz`` -- per-cell UMAP coordinates, restricted to a
   representative dataset set, because coordinates are the only bulky part.
   Emitted for: every clustering variant on one shared reference dataset (so the
   variants are compared like for like), plus the cohort run on one dataset per
   stage.

Usage::

    python s6_clustering_viz_data.py [--jobs N]
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import anndata as ad
import numpy as np
import pandas as pd

SWEEP = Path("/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed")
RUNS = SWEEP / "runs"
B2 = RUNS / "B1_cohort_full" / "B2_gex_celltyping"
OUT = SWEEP / "analysis_extra"

#: Reference dataset for the like-for-like clustering-variant comparison. The
#: grid runs were executed on a different 6-dataset subset that does not include
#: the cohort reference, so they get their own; both are Normal-stage.
REFERENCE_DS = "GSM3516666-Normal"
GRID_REFERENCE_DS = "GSM3516675-Normal"
#: one dataset per stage, for the cohort UMAP panel
STAGE_DS = ["GSM3516666-Normal", "GSM3516662-StageIA",
            "GSM3516665-StageIVprimary", "GSM3516664-MetBone"]


def clean_celltype(raw: str) -> str:
    """`CELL_TYPES_WANSLEEBEN_HOGAN_2013:MACROPHAGE_M2` -> `Macrophage M2`."""
    if not isinstance(raw, str) or not raw:
        return "unlabelled"
    label = raw.split(":", 1)[1] if ":" in raw else raw
    label = re.sub(r"\s*\(.*?\)\s*", " ", label)          # drop parenthetical citations
    label = label.replace("_", " ").strip()
    label = re.sub(r"\bMARKERS?\b", "", label, flags=re.I).strip()
    label = re.sub(r"\s+", " ", label)
    fixed = {"AE2": "AE2", "DC": "DC", "NK": "NK", "M1": "M1", "M2": "M2",
             "CD8T": "CD8T", "TREG": "Treg", "PNECS": "PNECs", "KRT5 B": "KRT5-B"}
    parts = [fixed.get(w.upper(), w.capitalize()) for w in label.split()]
    return " ".join(parts) or "unlabelled"


def load_celltypes() -> dict[str, pd.Series]:
    """dataset -> Series(barcode -> clean cell type), from the real B2 outputs."""
    out: dict[str, pd.Series] = {}
    for p in sorted(B2.glob("*_pas_labeled.h5ad")):
        ds = p.name.replace("_pas_labeled.h5ad", "")
        try:
            a = ad.read_h5ad(p)
        except Exception as exc:
            print(f"  WARN {p.name}: {exc}", flush=True)
            continue
        if "celltype" not in a.obs:
            continue
        key = a.obs["barcode"].astype(str) if "barcode" in a.obs else pd.Series(a.obs_names, index=a.obs_names)
        s = pd.Series(a.obs["celltype"].astype(str).map(clean_celltype).values, index=key.values)
        out[ds] = s[~s.index.duplicated()]
    return out


def experiments() -> list[tuple[str, Path]]:
    exps = [("B1_cohort_full", RUNS / "B1_cohort_full")]
    exps += [(f"grid/{p.name}", p) for p in sorted((RUNS / "grid").iterdir()) if p.is_dir()]
    exps += [(f"reannotate/{p.name}", p) for p in sorted((RUNS / "reannotate").iterdir()) if p.is_dir()]
    return [(n, p) for n, p in exps if (p / "07_clustering").is_dir()]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    print("[celltypes] loading real B2 labels", flush=True)
    ct = load_celltypes()
    print(f"  {len(ct)} datasets labelled; "
          f"{sorted({v for s in ct.values() for v in s.unique()}).__len__()} distinct cell types", flush=True)

    counts_rows: list[dict] = []
    umap_frames: list[pd.DataFrame] = []

    for exp_name, exp_dir in experiments():
        cl = exp_dir / "07_clustering"
        datasets = sorted(d.name for d in cl.iterdir() if (d / "clusters.h5ad").exists())
        print(f"[{exp_name}] {len(datasets)} datasets", flush=True)
        for ds in datasets:
            try:
                a = ad.read_h5ad(cl / ds / "clusters.h5ad")
            except Exception as exc:
                print(f"  WARN {exp_name}/{ds}: {exc}", flush=True)
                continue
            if "leiden" not in a.obs:
                continue
            bc = a.obs["barcode"].astype(str).values if "barcode" in a.obs else np.asarray(a.obs_names, dtype=str)
            leiden = a.obs["leiden"].astype(str).values
            labels = ct.get(ds)
            celltype = (pd.Series(bc).map(labels).fillna("unlabelled").values
                        if labels is not None else np.full(len(bc), "unlabelled"))
            matched = float((celltype != "unlabelled").mean())

            df = pd.DataFrame({"leiden": leiden, "celltype": celltype})
            g = df.groupby(["leiden", "celltype"]).size().reset_index(name="n_cells")
            g.insert(0, "dataset", ds)
            g.insert(0, "experiment", exp_name)
            g["celltype_match_rate"] = round(matched, 4)
            g["n_cells_dataset"] = len(df)
            counts_rows.append(g)

            ref = GRID_REFERENCE_DS if exp_name.startswith("grid/") else REFERENCE_DS
            want_umap = (ds == ref) or (exp_name == "B1_cohort_full" and ds in STAGE_DS)
            if want_umap and "X_umap" in a.obsm:
                u = np.asarray(a.obsm["X_umap"])[:, :2]
                umap_frames.append(pd.DataFrame({
                    "experiment": exp_name, "dataset": ds, "barcode": bc,
                    "leiden": leiden, "celltype": celltype,
                    "umap1": u[:, 0].round(4), "umap2": u[:, 1].round(4),
                }))

    counts = pd.concat(counts_rows, ignore_index=True)
    counts.to_csv(args.out / "cluster_celltype_counts.tsv", sep="\t", index=False)
    print(f"wrote cluster_celltype_counts.tsv  {counts.shape} "
          f"({counts.experiment.nunique()} experiments)", flush=True)

    umap = pd.concat(umap_frames, ignore_index=True)
    with gzip.open(args.out / "umap_cells.tsv.gz", "wt") as fh:
        umap.to_csv(fh, sep="\t", index=False)
    print(f"wrote umap_cells.tsv.gz  {umap.shape} "
          f"({umap.experiment.nunique()} experiments, {umap.dataset.nunique()} datasets)", flush=True)
    print("DONE s6", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
