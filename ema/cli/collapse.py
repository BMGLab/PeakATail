"""`ema collapse` — pool samtools-merge RG-suffix run tags in a unified matrix.

Corrects a ``before``-merge run whose barcodes were split by ``samtools merge``
RG suffixing, WITHOUT re-running peak calling. See :mod:`ema.datasets.collapse`.
"""
from __future__ import annotations

import logging

import click

from ema.cli.common import common_options, parse_log_overrides

log = logging.getLogger(__name__)


@click.command(name="collapse")
@common_options(include_output=False)
@click.option("--in-run", "in_run", required=True,
              type=click.Path(exists=True, file_okay=False),
              help="Completed before-merge run dir (has unified/concatenated.mtx).")
@click.option("--out-run", "out_run", required=True,
              type=click.Path(file_okay=False),
              help="Corrected run dir to create.")
def collapse(**kwargs) -> None:
    """Pool per-run cbs into per-library cells (RG-suffix collapse)."""
    from ema.logging_config import setup_logging, teardown_logging
    setup_logging(
        level=(kwargs["log_level"] or "INFO").split(",")[0],
        output_dir=None,
        quiet=kwargs["quiet"],
        verbose_count=kwargs["verbose"],
        no_log_file=True,
        per_logger_overrides=parse_log_overrides(kwargs["log_level"]),
    )
    try:
        from ema.datasets.collapse import collapse_run
        summary = collapse_run(kwargs["in_run"], kwargs["out_run"])
        log.info(
            "collapsed cells: %d -> %d (pooled %d per-run duplicates)",
            summary["n_cells_before"], summary["n_cells_after"], summary["n_pooled"],
        )
        log.info(
            "libraries: %d run-prefixes -> %d library-prefixes",
            summary["n_libraries_before"], summary["n_libraries_after"],
        )
        click.echo(f"corrected run written to {kwargs['out_run']}")
        if (summary["n_libraries_after"] > 40
                or summary["n_libraries_before"] == summary["n_libraries_after"]):
            log.warning(
                "library count did not collapse as expected — check the RG-suffix "
                "pattern against the actual cb prefixes."
            )
    finally:
        teardown_logging()
