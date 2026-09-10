"""Generate minimal placeholder PNG files for docs build.

This script is NOT the intended production path. Run scripts/copy_doc_figures.py
instead to copy the real figures from peakatail_runs/. This script exists as a
fallback for CI environments that don't have the run directory.

Usage (from repo root):
    uv run python docs/assets/figures/generate_placeholders.py
"""
import base64
import pathlib

# Minimal 1x1 white PNG (67 bytes), base64-encoded.
_MINIMAL_PNG_B64 = (
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8z8BQDwADhQGAWjR9awAAAABJRU5ErkJggg=="
)

out = pathlib.Path(__file__).parent

figures = [
    "umap_default.png",
    "clusters_default.png",
    "peak_qc_default.png",
    "resource_timeline.png",
    "volcano_0_vs_4.png",
    "pdui_distribution.png",
    "proportion_distribution.png",
    "gene_ENSG00000103275.png",
]

png_bytes = base64.b64decode(_MINIMAL_PNG_B64)
for name in figures:
    dst = out / name
    if not dst.exists():
        dst.write_bytes(png_bytes)
        print(f"  placeholder  {name}")
    else:
        print(f"  exists       {name}")

print("Done.")
