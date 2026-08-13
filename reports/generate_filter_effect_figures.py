#!/usr/bin/env python
"""Figure 20 — filter-effect: what each PAS filter costs downstream.

Reads the clean per-scenario CSVs (computed without `fisher` or `proportion`,
so they are final) from `reports/cumulative_analysis/filter_effect/` and renders
the reserved slot 20.

    python reports/generate_filter_effect_figures.py \
        --in reports/cumulative_analysis/filter_effect \
        --out reports/figures/corrected

Slots 21-23 stay reserved for the differential, length and UMAP/Sankey panels,
which need the corrected diff/length engines.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

logger = logging.getLogger(__name__)

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight",
    "font.size": 9, "axes.grid": True, "grid.alpha": 0.25,
    "axes.spines.top": False, "axes.spines.right": False,
})

ORDER = ["baseline", "atlas_filter", "annot_filter_3utr", "ip_filter"]
LABEL = {
    "baseline": "baseline\n(keep all)",
    "atlas_filter": "atlas\nnon-match",
    "annot_filter_3utr": "3′UTR-only\nannotation",
    "ip_filter": "internal\npriming",
}
# The 3'UTR filter is the one that does damage, so it is the one coloured.
COLOR = {"baseline": "#9e9e9e", "atlas_filter": "#1f4e79",
         "annot_filter_3utr": "#c62828", "ip_filter": "#2e7d32"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="in_dir", required=True, type=Path)
    ap.add_argument("--out", dest="out_dir", required=True, type=Path)
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    pas = pd.read_csv(a.in_dir / "pas_kept_vs_used.csv")
    clu = pd.read_csv(a.in_dir / "clustering_ari_ami.csv")
    gex = pd.read_csv(a.in_dir / "gex_concordance.csv")
    d = pas.merge(clu, on="scenario").merge(gex, on="scenario")
    d["scenario"] = pd.Categorical(d["scenario"], ORDER, ordered=True)
    d = d.sort_values("scenario")
    base_cells = float(d.loc[d.scenario == "baseline", "n_cells_total"].iloc[0])
    d["pct_cells_lost"] = 100 * (base_cells - d["n_cells_total"]) / base_cells
    x = np.arange(len(d))
    cols = [COLOR[s] for s in d["scenario"]]
    labs = [LABEL[s] for s in d["scenario"]]

    fig, axes = plt.subplots(1, 4, figsize=(13.6, 3.5))

    ax = axes[0]
    ax.bar(x, d["pct_excluded"], color=cols)
    for xi, v in zip(x, d["pct_excluded"]):
        ax.text(xi, v + 1.5, f"{v:.1f}%", ha="center", fontsize=8)
    ax.set_ylim(0, 85)
    ax.set_ylabel("% of PAS excluded from clustering")
    ax.set_title("What the filter removes\n(all keep 400,253 in results)", fontsize=9.5)

    ax = axes[1]
    ax.bar(x - 0.2, d["ARI_vs_baseline_mean"], 0.4, color=cols, label="ARI")
    ax.bar(x + 0.2, d["AMI_vs_baseline_mean"], 0.4, color=cols, alpha=0.55, label="AMI")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("agreement with baseline clustering")
    ax.set_title("Does the clustering change?\n(solid ARI, pale AMI)", fontsize=9.5)

    ax = axes[2]
    ax.bar(x, d["ARI_vs_GEX_mean"], color=cols)
    base_gex = float(d.loc[d.scenario == "baseline", "ARI_vs_GEX_mean"].iloc[0])
    ax.axhline(base_gex, color="k", lw=0.9, ls="--")
    ax.annotate("baseline", xy=(len(d) - 0.6, base_gex), fontsize=7, va="bottom")
    ax.set_ylim(0, 0.55)
    ax.set_ylabel("ARI vs GEX cell type")
    ax.set_title("Does it cost biological signal?", fontsize=9.5)

    ax = axes[3]
    ax.bar(x - 0.2, d["pct_cells_reassigned_vs_baseline_mean"].fillna(0), 0.4,
           color=cols, label="reassigned")
    ax.bar(x + 0.2, d["pct_cells_lost"], 0.4, color=cols, alpha=0.5, hatch="//",
           label="dropped entirely")
    ax.set_ylabel("% of cells")
    ax.set_title("Cost to cells\n(solid reassigned, hatched dropped)", fontsize=9.5)
    ax.legend(fontsize=7)

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labs, fontsize=7.5)

    fig.suptitle(
        "Filter effect: restricting PAS to annotated 3′UTRs excludes 73% of calls and is the only "
        "filter that costs anything — atlas-nonmatch and internal-priming (~8% each) are nearly free.",
        fontsize=9.5, y=1.03,
    )
    fig.tight_layout()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(a.out_dir / f"20_filter_effect_clustering.{ext}")
        logger.info("wrote 20_filter_effect_clustering.%s", ext)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
