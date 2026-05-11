"""Bar chart of cluster sizes per dataset."""
from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_matplotlib


@register_viz_strategy
class ClusterSizesMatplotlib(VizStrategy):
    name = "cluster_sizes_matplotlib"
    plot_type = "cluster_sizes"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        adata, ds_id = data
        counts = Counter(adata.obs["leiden"]) if "leiden" in adata.obs else Counter()
        labels = sorted(counts, key=lambda x: int(x) if str(x).isdigit() else x)
        sizes = [counts[c] for c in labels]
        fig, ax = plt.subplots(figsize=(max(4, len(labels) * 0.5), 4))
        ax.bar(labels, sizes, color="#4c72b0")
        ax.set_xlabel("cluster")
        ax.set_ylabel("# cells")
        ax.set_title(f"Cluster sizes — {ds_id}")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for i, s in enumerate(sizes):
            ax.text(i, s, str(s), ha="center", va="bottom", fontsize=9)
        return save_matplotlib(fig, output_basepath)
