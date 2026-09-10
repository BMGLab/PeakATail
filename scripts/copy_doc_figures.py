"""Pre-build step: copy real figures from peakatail_runs/ into docs/assets/figures/.

Run from the repository root before `mkdocs build`:

    uv run python scripts/copy_doc_figures.py

This script is idempotent; it overwrites existing files.
"""
import shutil
import sys
from pathlib import Path

repo = Path(__file__).parent.parent
rundir = repo / "peakatail_runs/full_v8_2026-05-11_152746"
out = repo / "docs/assets/figures"
out.mkdir(parents=True, exist_ok=True)

copies = [
    (rundir / "figures/umap_default.png",
     out / "umap_default.png"),
    (rundir / "figures/clusters_default.png",
     out / "clusters_default.png"),
    (rundir / "figures/peak_qc_default.png",
     out / "peak_qc_default.png"),
    (rundir / "figures/resource_timeline.png",
     out / "resource_timeline.png"),
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

errors = 0
for src, dst in copies:
    if src.exists():
        shutil.copy2(src, dst)
        print(f"  OK      {dst.name}")
    else:
        print(f"  MISSING {src.relative_to(repo)}", file=sys.stderr)
        errors += 1

print(f"\n{len(copies) - errors}/{len(copies)} files copied to {out.relative_to(repo)}/")
if errors:
    sys.exit(1)
