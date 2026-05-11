"""Smoke test: gene_track_plotly renderer.

Picks a gene with >=3 PAS from the full_v8 run, builds a GenePanel,
renders through GeneTrackPlotly, and prints the meta.json content.

Run with:
    uv run python _smoke_gene_track_plotly.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

# Ensure ema is importable from working dir.
sys.path.insert(0, str(Path(__file__).parent))

from ema.viz._gene_track_helpers import build_gene_panel
# Importing triggers @register_viz_strategy
from ema.viz.gene_track_plotly import GeneTrackPlotly

RUN_DIR = Path("peakatail_runs/full_v8_2026-05-11_152746/per_dataset/only")
OUT_DIR = Path("peakatail_runs/full_v8_2026-05-11_152746/per_dataset/only/figs/gene_track")

def main() -> None:
    print("Loading clusters.h5ad ...")
    adata = ad.read_h5ad(RUN_DIR / "clusters.h5ad")

    print("Loading pasbed.bed ...")
    pasbed = pd.read_csv(
        RUN_DIR / "pasbed.bed",
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "pas_id", "score", "strand"],
    )

    # Pick a gene with >= 3 PAS.
    counts = adata.var["gene_id"].value_counts()
    gene_id = counts[counts >= 3].index[0]
    print(f"Selected gene: {gene_id} ({counts[gene_id]} PAS)")

    panel = build_gene_panel(gene_id, adata, pasbed, cluster_key="leiden")
    if panel is None:
        print("ERROR: build_gene_panel returned None")
        sys.exit(1)

    print(
        f"Panel: {len(panel.pas_ids)} PAS, "
        f"{len(panel.clusters)} clusters, "
        f"{panel.chrom}:{panel.start}-{panel.end} ({panel.strand})"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    basepath = OUT_DIR / f"gene_track_{gene_id}"

    strategy = GeneTrackPlotly()
    paths = strategy.render(panel, basepath)
    print(f"\nFiles written ({len(paths)}):")
    for p in paths:
        print(f"  {p}")

    meta_path = basepath.with_suffix(".meta.json")
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        print("\n--- meta.json ---")
        print(json.dumps(meta, indent=2))
    else:
        print("WARNING: meta.json not found")


if __name__ == "__main__":
    main()
