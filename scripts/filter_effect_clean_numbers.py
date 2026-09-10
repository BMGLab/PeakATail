#!/usr/bin/env python3
"""Clean (fisher/proportion-independent) filter-effect comparison numbers:
kept-vs-used PAS, clustering (n_clusters, ARI/AMI vs baseline + vs GEX),
%cells reassigned, and celltyping summary. Deliberately does NOT touch
B3_switch/{diff,length} (fisher/proportion, buggy pending fix).

Usage: python filter_effect_clean_numbers.py <FILTER_EFFECT_root>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

SCENARIOS = ["baseline", "atlas_filter", "annot_filter_3utr", "ip_filter"]


def load_obs(run_dir: Path, gsm: str) -> pd.DataFrame | None:
    h5 = run_dir / "07_clustering" / gsm / "clusters.h5ad"
    if not h5.exists():
        return None
    a = ad.read_h5ad(h5, backed="r")
    return a.obs[["leiden"]].copy()


def main(root: str) -> None:
    root = Path(root)
    runs = root / "runs"

    gsms = sorted(p.name for p in (runs / "baseline" / "07_clustering").iterdir())

    # 1. PAS kept-in-results vs used-for-clustering (from branch_manifest.json).
    pas_rows = []
    for sc in SCENARIOS:
        m = json.loads((runs / sc / "branch_manifest.json").read_text())
        pm = m["pas_labels_and_mask"]
        pas_rows.append({
            "scenario": sc,
            "n_kept_in_results": pm["n_kept_in_results"],
            "n_used_for_clustering": pm["n_used_for_clustering"],
            "n_excluded": pm["n_excluded_for_clustering"],
            "pct_excluded": round(100 * pm["n_excluded_for_clustering"] / pm["n_kept_in_results"], 2),
        })
    pas_df = pd.DataFrame(pas_rows)

    # 2. Per-GSM cluster counts + cross-scenario ARI/AMI vs baseline.
    obs_by_scenario = {sc: {} for sc in SCENARIOS}
    for sc in SCENARIOS:
        for gsm in gsms:
            obs = load_obs(runs / sc, gsm)
            if obs is not None:
                obs_by_scenario[sc][gsm] = obs

    cluster_rows = []
    for sc in SCENARIOS:
        n_clusters_list = []
        n_cells_list = []
        ari_list, ami_list, pct_reassigned_list = [], [], []
        for gsm in gsms:
            obs_sc = obs_by_scenario[sc].get(gsm)
            if obs_sc is None:
                continue
            n_clusters_list.append(obs_sc["leiden"].nunique())
            n_cells_list.append(len(obs_sc))
            obs_base = obs_by_scenario["baseline"].get(gsm)
            if obs_base is not None and sc != "baseline":
                common = obs_sc.index.intersection(obs_base.index)
                if len(common) > 2:
                    l_sc = obs_sc.loc[common, "leiden"].astype(str)
                    l_base = obs_base.loc[common, "leiden"].astype(str)
                    ari_list.append(adjusted_rand_score(l_base, l_sc))
                    ami_list.append(adjusted_mutual_info_score(l_base, l_sc))
                    # "% cells reassigned": cells NOT in the majority-overlap
                    # cluster pairing (best-effort Hungarian-free proxy: for
                    # each baseline cluster, the scenario cluster it overlaps
                    # with MOST is its "match"; any cell not in that matched
                    # pair for its baseline cluster counts as reassigned).
                    ct = pd.crosstab(l_base, l_sc)
                    match = ct.idxmax(axis=1)  # baseline_cluster -> best-matching scenario_cluster
                    matched_mask = l_sc.values == l_base.map(match).values
                    pct_reassigned_list.append(round(100 * (1 - matched_mask.mean()), 2))
        cluster_rows.append({
            "scenario": sc,
            "n_gsm": len(n_clusters_list),
            "n_clusters_mean": round(float(np.mean(n_clusters_list)), 2) if n_clusters_list else None,
            "n_cells_total": int(sum(n_cells_list)) if n_cells_list else None,
            "ARI_vs_baseline_mean": round(float(np.mean(ari_list)), 4) if ari_list else None,
            "AMI_vs_baseline_mean": round(float(np.mean(ami_list)), 4) if ami_list else None,
            "pct_cells_reassigned_vs_baseline_mean": (
                round(float(np.mean(pct_reassigned_list)), 2) if pct_reassigned_list else None
            ),
        })
    cluster_df = pd.DataFrame(cluster_rows)

    # 3. GEX concordance (PAS-leiden vs GEX celltype), already computed per scenario.
    gex_rows = []
    for sc in SCENARIOS:
        conc_path = runs / sc / "B2_gex_celltyping" / "concordance.csv"
        if not conc_path.exists():
            gex_rows.append({"scenario": sc, "n_gsm": 0})
            continue
        conc = pd.read_csv(conc_path)
        ari_cols = [c for c in conc.columns if c.startswith("ARI_")]
        ami_cols = [c for c in conc.columns if c.startswith("AMI_")]
        row = {"scenario": sc, "n_gsm": len(conc)}
        if ari_cols:
            row["ARI_vs_GEX_mean"] = round(float(conc[ari_cols].mean().mean()), 4)
        if ami_cols:
            row["AMI_vs_GEX_mean"] = round(float(conc[ami_cols].mean().mean()), 4)
        gex_rows.append(row)
    gex_df = pd.DataFrame(gex_rows)

    # 4. Celltyping summary: n distinct GEX celltypes labeled, per scenario/GSM.
    celltype_rows = []
    for sc in SCENARIOS:
        b2 = runs / sc / "B2_gex_celltyping"
        n_celltypes_per_gsm = []
        for gsm in gsms:
            h5 = b2 / f"{gsm}_pas_labeled.h5ad"
            if not h5.exists():
                continue
            a = ad.read_h5ad(h5, backed="r")
            col = "celltype" if "celltype" in a.obs.columns else None
            if col:
                n_celltypes_per_gsm.append(a.obs[col].nunique())
        celltype_rows.append({
            "scenario": sc,
            "n_gsm_typed": len(n_celltypes_per_gsm),
            "n_celltypes_mean": round(float(np.mean(n_celltypes_per_gsm)), 2) if n_celltypes_per_gsm else None,
        })
    celltype_df = pd.DataFrame(celltype_rows)

    out_dir = root / "logs" / "clean_numbers"
    out_dir.mkdir(parents=True, exist_ok=True)
    pas_df.to_csv(out_dir / "pas_kept_vs_used.csv", index=False)
    cluster_df.to_csv(out_dir / "clustering_ari_ami.csv", index=False)
    gex_df.to_csv(out_dir / "gex_concordance.csv", index=False)
    celltype_df.to_csv(out_dir / "celltyping_summary.csv", index=False)

    print("=== PAS kept-in-results vs used-for-clustering ===")
    print(pas_df.to_string(index=False))
    print("\n=== Clustering: n_clusters, ARI/AMI vs baseline, %cells reassigned ===")
    print(cluster_df.to_string(index=False))
    print("\n=== GEX concordance (PAS-leiden vs GEX celltype) ===")
    print(gex_df.to_string(index=False))
    print("\n=== Celltyping summary ===")
    print(celltype_df.to_string(index=False))
    print(f"\nWritten to: {out_dir}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "/mnt/ssd2/Laugney_Aligned/peakatail_experiments/FILTER_EFFECT_2026-08")
