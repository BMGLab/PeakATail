"""Run this from the repo root to copy real figures into docs/assets/figures/."""
import shutil
from pathlib import Path

repo = Path(__file__).parent.parent.parent.parent
rundir = repo / "peakatail_runs/full_v8_2026-05-11_152746"
out = repo / "docs/assets/figures"
out.mkdir(parents=True, exist_ok=True)

copies = [
    (rundir / "figures/umap_default.png",           out / "umap_default.png"),
    (rundir / "figures/clusters_default.png",        out / "clusters_default.png"),
    (rundir / "figures/peak_qc_default.png",         out / "peak_qc_default.png"),
    (rundir / "figures/resource_timeline.png",       out / "resource_timeline.png"),
    (rundir / "switch_diff_2026-05-11_205015/figures/volcano_0_vs_4.png",
     out / "volcano_0_vs_4.png"),
    (rundir / "switch_length_2026-05-11_205045/figures/pdui_distribution.png",
     out / "pdui_distribution.png"),
    (rundir / "switch_length_2026-05-11_194149/figures/pdui_distribution.png",
     out / "proportion_distribution.png"),
    (rundir / "switch_geneview_2026-05-11_212341/figures/gene_ENSG00000103275.png",
     out / "gene_ENSG00000103275.png"),
    (rundir / "switch_diff_2026-05-11_205015/figures/figures_INDEX.md",
     out / "switch_diff_figures_INDEX.md"),
]

for src, dst in copies:
    if src.exists():
        shutil.copy2(src, dst)
        print(f"  OK  {dst.name}")
    else:
        print(f"  MISSING  {src}")

print(f"\nDone. {len(list(out.glob('*.png')))} PNG files in {out}")
