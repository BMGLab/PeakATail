from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import plotly.express as px

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_plotly
from ema.viz._meta import write_figure_meta


@register_viz_strategy
class ClusterSizesPlotly(VizStrategy):
    name = "cluster_sizes_plotly"
    plot_type = "cluster_sizes"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        adata, ds_id = data
        counts = Counter(adata.obs["leiden"]) if "leiden" in adata.obs else Counter()
        labels = sorted(counts, key=lambda x: int(x) if str(x).isdigit() else x)
        sizes = [counts[c] for c in labels]
        fig = px.bar(x=labels, y=sizes, labels={"x": "cluster", "y": "# cells"},
                     title=f"Cluster sizes — {ds_id}")
        fig.update_layout(template="plotly_white")
        paths = save_plotly(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "dataset_id": ds_id,
            "cluster_key": "leiden",
            "n_clusters": len(labels),
            "n_observations": int(adata.n_obs),
        })
        return paths
