"""`ema merge` — merge BAM files (was `ema_merge`)."""
from __future__ import annotations

import logging

import click

from ema.cli.common import common_options, parse_log_overrides

log = logging.getLogger(__name__)


@click.command(name="merge")
@common_options(include_output=True, output_default="merged.bam")
@click.option("--bam-files", "-i", "bam_files", multiple=True, required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="BAM files to merge. Repeat -i for each.")
def merge(**kwargs) -> None:
    """Merge multiple BAM files into one sorted+indexed BAM."""
    from pathlib import Path
    from ema.logging_config import setup_logging

    out_path = Path(kwargs["output"])
    out_path.parent.mkdir(parents=True, exist_ok=True)

    setup_logging(
        level=(kwargs["log_level"] or "INFO").split(",")[0],
        output_dir=out_path.parent,
        quiet=kwargs["quiet"],
        verbose_count=kwargs["verbose"],
        no_log_file=kwargs["no_log_file"],
        per_logger_overrides=parse_log_overrides(kwargs["log_level"]),
    )

    log.info("ema merge: %d input BAMs -> %s", len(kwargs["bam_files"]), out_path)

    # Reuse existing merge_bam logic
    from ema.merge_bam.cli import run_merge  # will be promoted in Phase 9
    run_merge(
        bam_files=list(kwargs["bam_files"]),
        output=str(out_path),
        threads=kwargs["threads"] or 4,
    )
