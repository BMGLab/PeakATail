"""`ema switch trend` — ordered-stage APA trend from a PDUI long table.

Reads a per-cluster PDUI/proportion/entropy table (as produced by
``ema switch length``), treats the clusters as an ORDERED progression given by
``--stage-order``, and reports the trend (slope + Spearman monotonicity +
direction) overall and per gene. See :mod:`ema.switch_test.trend`.
"""
from __future__ import annotations

import json
import logging

import click

from ema.cli.common import common_options, parse_log_overrides

log = logging.getLogger(__name__)


@click.command(name="trend")
@common_options(include_output=True)
@click.option("--pdui", "pdui_tsv", required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="PDUI long TSV from `ema switch length` (has stage + value columns).")
@click.option("--stage-order", "stage_order", required=True,
              help="Comma-separated ordered stage labels, e.g. "
                   "'Normal,StageI,StageII,StageIII,StageIV,Met'.")
@click.option("--stage-col", "stage_col", default="cluster", show_default=True,
              help="Column holding the stage/cluster label.")
@click.option("--value-col", "value_col", default="pdui", show_default=True,
              help="Numeric value column (pdui/proportion/entropy).")
@click.option("--gene-col", "gene_col", default="gene_id", show_default=True,
              help="Gene column for the per-gene trend table (skipped if absent).")
def trend(**kwargs) -> None:
    """Aggregate an ordered-stage APA length trend."""
    from ema.logging_config import setup_logging, teardown_logging
    from pathlib import Path

    out_dir = Path(kwargs["output"] or ".")
    out_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(
        level=(kwargs["log_level"] or "INFO").split(",")[0],
        output_dir=str(out_dir),
        quiet=kwargs["quiet"],
        verbose_count=kwargs["verbose"],
        per_logger_overrides=parse_log_overrides(kwargs["log_level"]),
    )
    try:
        import pandas as pd
        from ema.switch_test.trend import pdui_trend, pdui_trend_by_gene

        stage_order = [s.strip() for s in kwargs["stage_order"].split(",") if s.strip()]
        if len(stage_order) < 2:
            raise click.UsageError("--stage-order needs at least 2 stages.")

        df = pd.read_csv(kwargs["pdui_tsv"], sep="\t")
        overall = pdui_trend(
            df, stage_order,
            stage_col=kwargs["stage_col"], value_col=kwargs["value_col"],
        )
        (out_dir / "length_trend.json").write_text(json.dumps(overall, indent=2, default=str))
        log.info(
            "overall %s trend: slope=%.4g spearman=%.3g direction=%s (over %d stages)",
            kwargs["value_col"], overall["slope"], overall["spearman"],
            overall["direction"], overall["n_stages"],
        )

        if kwargs["gene_col"] in df.columns:
            per_gene = pdui_trend_by_gene(
                df, stage_order,
                stage_col=kwargs["stage_col"], value_col=kwargs["value_col"],
                gene_col=kwargs["gene_col"],
            )
            per_gene.to_csv(out_dir / "length_trend_by_gene.tsv", sep="\t")
            log.info("wrote per-gene trend for %d genes", len(per_gene))
        else:
            log.warning(
                "no %r column — per-gene trend skipped (overall only).",
                kwargs["gene_col"],
            )
        click.echo(f"trend written to {out_dir}")
    finally:
        teardown_logging()
