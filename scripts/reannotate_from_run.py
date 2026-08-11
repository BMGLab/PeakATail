#!/usr/bin/env python3
"""Branch a completed ``ema run`` into a new trim / clustering variant WITHOUT
re-running peak calling.

THIN SHIM: this script is now a CLI wrapper around
:func:`ema.reannotate.reannotate_run` — the same function backing the
``ema reannotate`` subcommand (``ema/cli/reannotate.py``).  It exists so
existing invocations of this script keep working unchanged; new callers
should prefer ``ema reannotate`` directly.

Peak calling (streaming the BAMs) is the expensive stage.  The "trim" —
``find_close(max_gene_distance, utr_multiplier, include_extended)`` — and
everything downstream of it (annotate → preprocess → cluster) is cheap and
depends only on artifacts a base run already wrote to disk:

    <base_run>/posbed.bed              (unified +strand PAS)
    <base_run>/negbed.bed              (unified -strand PAS)
    <base_run>/unified/concatenated.mtx        (PAS x cell counts)
    <base_run>/unified/concatenated_cbs.tsv    (namespaced barcodes)

``reannotate_run`` reuses the *exact* tested internals the pipeline uses
(``find_close`` + the per-dataset worker ``run_one_dataset_downstream``), so a
branch is behaviourally identical to having run ``ema run`` with those trim /
clustering parameters — it just skips peak calling.

This was the standalone form of what became ``ema reannotate``.  Validate
against one dataset before trusting a wide sweep: run a base ``ema run`` on
one GSM, then branch it here with the SAME trim params and confirm the
resulting clusters.h5ad matches the base run's.

Example
-------
    python scripts/reannotate_from_run.py \
        --base-run   /mnt/ssd2/.../peakatail_runs/A2_base_lambda_gradient \
        --gtf        /home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf \
        --out        /mnt/ssd2/.../peakatail_runs/A2_trim_d3000_ext \
        --max-gene-distance 3000 --include-extended \
        --resolution 1.0

Nothing is launched implicitly — this runs the downstream stages for the
datasets found in the base run and exits.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
log = logging.getLogger("reannotate")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-run", type=Path, required=True,
                    help="Completed `ema run` output dir to branch from.")
    ap.add_argument("--out", type=Path, required=True,
                    help="Fresh output dir for this branch (must not equal --base-run).")
    ap.add_argument("--gtf", type=Path, required=True, help="Same GTF as the base run.")

    # ── trim knobs (the whole point) ─────────────────────────────────────
    ap.add_argument("--max-gene-distance", type=int, default=5000,
                    help="TIER_3 distal cap (bp). Matters when --include-extended.")
    ap.add_argument("--utr-multiplier", type=float, default=2.0,
                    help="TIER_2 boundary = gene UTR length x this.")
    ap.add_argument("--include-extended", action="store_true",
                    help="Keep TIER_3 (distal/novel) PAS out to --max-gene-distance.")

    # ── clustering knobs (branch A3 from the same base peaks) ────────────
    ap.add_argument("--cluster-method", default="leiden_tfidf")
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--n-neighbors", type=int, default=None)
    ap.add_argument("--n-pcs", type=int, default=40)
    ap.add_argument("--n-svd-components", type=int, default=50)
    ap.add_argument("--n-top-hvg", type=int, default=2000,
                    help="HVG count (leiden_libsize only).")
    ap.add_argument("--random-seed", type=int, default=42)
    ap.add_argument("--tfidf-scale-factor", type=float, default=1e4,
                    help="Signac Method 1 TF-IDF scale factor (leiden_tfidf only).")
    ap.add_argument("--depth-corr-threshold", type=float, default=0.75,
                    help="Pearson |r| threshold for dropping depth-correlated LSI "
                         "components (leiden_tfidf only). 1.0 disables the filter.")
    ap.add_argument("--external-clusters", type=str, default=None,
                    help="Path to external cluster labels TSV (--cluster-method external).")

    # ── cell/PAS filters (match the base run's defaults unless overriding) ─
    # Defaults MUST match ema run's schema defaults (config_schema.py) so a
    # branch with unchanged params reproduces the base run's clustering:
    # min_read=1500, min_cells=3, min_pas_per_cell=50 (bridges to min_genes).
    ap.add_argument("--min-read", type=int, default=1500)
    ap.add_argument("--min-cells", type=int, default=3)
    ap.add_argument("--min-pas-per-cell", type=int, default=50)
    ap.add_argument("--threads", type=int, default=None,
                    help="Absolute worker ceiling for the downstream Pool.")
    args = ap.parse_args()

    from ema.reannotate import ReannotateError, reannotate_run

    try:
        manifest = reannotate_run(
            base_run=args.base_run,
            out=args.out,
            gtf=args.gtf,
            max_gene_distance=args.max_gene_distance,
            utr_multiplier=args.utr_multiplier,
            include_extended=args.include_extended,
            cluster_method=args.cluster_method,
            resolution=args.resolution,
            n_neighbors=args.n_neighbors,
            n_pcs=args.n_pcs,
            n_svd_components=args.n_svd_components,
            n_top_hvg=args.n_top_hvg,
            random_seed=args.random_seed,
            tfidf_scale_factor=args.tfidf_scale_factor,
            depth_corr_threshold=args.depth_corr_threshold,
            external_clusters=args.external_clusters,
            min_read=args.min_read,
            min_cells=args.min_cells,
            min_pas_per_cell=args.min_pas_per_cell,
            threads=args.threads,
        )
    except ReannotateError as e:
        sys.exit(f"ERROR: {e}")

    log.info("DONE — %d datasets clustered under %s", len(manifest["datasets"]), args.out)


if __name__ == "__main__":
    main()
