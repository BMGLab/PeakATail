"""`ema reannotate` — branch a completed ``ema run`` into a new trim / filter /
clustering variant WITHOUT re-running peak calling.

Peak calling (streaming the BAMs) is the expensive stage.  This subcommand
reuses the raw peak-call artifacts a base ``ema run`` already wrote
(``posbed.bed`` / ``negbed.bed`` / ``unified/concatenated*.mtx``) and re-runs
only the cheap downstream stages — trim (``find_close``) -> annotate ->
preprocess -> cluster — through the SAME tested internals ``ema run`` uses.
The result is a complete, chainable run dir: ``ema.data.Run.from_dir()`` loads
it directly and ``ema switch {diff,length,trend}`` can consume its
``07_clustering/<ds>/clusters.h5ad`` files.

All the actual work lives in :func:`ema.reannotate.reannotate_run` — this
module is a thin Click wrapper (mirrors ``ema/cli/merge.py`` /
``ema/cli/parse_gtf.py``: custom required path options, no ``--output``
timestamping, since ``--out`` here names the branch's own root directory).
"""
from __future__ import annotations

import logging
from pathlib import Path

import click

from ema.cli.common import common_options, parse_log_overrides

log = logging.getLogger(__name__)


@click.command(name="reannotate")
@common_options(include_output=False)
@click.option("--base-run", "base_run", required=True,
              type=click.Path(exists=True, file_okay=False, resolve_path=True),
              help="Completed `ema run` output dir to branch from.")
@click.option("--out", "out", required=True,
              type=click.Path(file_okay=False, resolve_path=True),
              help="Fresh output dir for this branch (must differ from --base-run).")
@click.option("--gtf", "gtf", required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True),
              help="Same GTF as the base run (or a new one to re-annotate against).")
# ── trim knobs (the whole point) ─────────────────────────────────────────
@click.option("--max-gene-distance", "max_gene_distance", type=int, default=5000,
              show_default=True,
              help="TIER_3 distal cap (bp). Matters when --include-extended.")
@click.option("--utr-multiplier", "utr_multiplier", type=float, default=2.0,
              show_default=True,
              help="TIER_2 boundary = gene UTR length x this.")
@click.option("--include-extended", "include_extended", is_flag=True, default=False,
              help="Keep TIER_3 (distal/novel) PAS out to --max-gene-distance.")
# ── clustering knobs (branch from the same base peaks) ───────────────────
@click.option("--cluster-method", "cluster_method", type=str, default="leiden_tfidf",
              show_default=True)
@click.option("--resolution", "resolution", type=float, default=1.0, show_default=True)
@click.option("--n-neighbors", "n_neighbors", type=int, default=None)
@click.option("--n-pcs", "n_pcs", type=int, default=40, show_default=True)
@click.option("--n-svd-components", "n_svd_components", type=int, default=50, show_default=True)
@click.option("--n-top-hvg", "n_top_hvg", type=int, default=2000, show_default=True,
              help="HVG count (leiden_libsize only).")
@click.option("--random-seed", "random_seed", type=int, default=42, show_default=True)
@click.option("--tfidf-scale-factor", "tfidf_scale_factor", type=float, default=1e4,
              show_default=True,
              help="Signac Method 1 TF-IDF scale factor (leiden_tfidf only).")
@click.option("--depth-corr-threshold", "depth_corr_threshold", type=float, default=0.75,
              show_default=True,
              help="Pearson |r| threshold for dropping depth-correlated LSI "
                   "components (leiden_tfidf only). 1.0 disables the filter.")
@click.option("--external-clusters", "external_clusters", type=click.Path(exists=True),
              default=None,
              help="Path to external cluster labels TSV (--cluster-method external).")
# ── cell/PAS filters (match the base run's defaults unless overriding) ───
# Defaults MUST match ema run's schema defaults (config_schema.py) so a
# branch with unchanged params reproduces the base run's clustering:
# min_read=1500, min_cells=3, min_pas_per_cell=50 (bridges to min_genes).
@click.option("--min-read", "min_read", type=int, default=1500, show_default=True)
@click.option("--min-cells", "min_cells", type=int, default=3, show_default=True)
@click.option("--min-pas-per-cell", "min_pas_per_cell", type=int, default=50, show_default=True)
def reannotate(**kwargs) -> None:
    """Branch a completed `ema run` into a new trim/cluster variant, skipping
    peak calling."""
    out_dir = Path(kwargs["out"])
    out_dir.mkdir(parents=True, exist_ok=True)

    from ema.logging_config import setup_logging, teardown_logging
    setup_logging(
        level=(kwargs["log_level"] or "INFO").split(",")[0],
        output_dir=out_dir,
        quiet=kwargs["quiet"],
        verbose_count=kwargs["verbose"],
        no_log_file=kwargs["no_log_file"],
        per_logger_overrides=parse_log_overrides(kwargs["log_level"]),
    )

    try:
        log.info("ema reannotate: %s -> %s", kwargs["base_run"], out_dir)
        if kwargs.get("threads") is not None:
            log.info("ResourceManager: --threads=%d (absolute ceiling)", kwargs["threads"])

        from ema.reannotate import ReannotateError, reannotate_run
        try:
            manifest = reannotate_run(
                base_run=kwargs["base_run"],
                out=out_dir,
                gtf=kwargs["gtf"],
                max_gene_distance=kwargs["max_gene_distance"],
                utr_multiplier=kwargs["utr_multiplier"],
                include_extended=kwargs["include_extended"],
                cluster_method=kwargs["cluster_method"],
                resolution=kwargs["resolution"],
                n_neighbors=kwargs["n_neighbors"],
                n_pcs=kwargs["n_pcs"],
                n_svd_components=kwargs["n_svd_components"],
                n_top_hvg=kwargs["n_top_hvg"],
                random_seed=kwargs["random_seed"],
                tfidf_scale_factor=kwargs["tfidf_scale_factor"],
                depth_corr_threshold=kwargs["depth_corr_threshold"],
                external_clusters=kwargs["external_clusters"],
                min_read=kwargs["min_read"],
                min_cells=kwargs["min_cells"],
                min_pas_per_cell=kwargs["min_pas_per_cell"],
                threads=kwargs["threads"],
            )
        except ReannotateError as e:
            raise click.ClickException(str(e))
        log.info(
            "ema reannotate: DONE — %d dataset(s) -> %s",
            len(manifest.get("datasets", [])), out_dir,
        )
    finally:
        teardown_logging()
