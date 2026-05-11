"""`ema switch length` — 3'UTR shortening / lengthening quantification."""
from __future__ import annotations

import logging

import click

from ema.cli.common import (
    apply_yaml_to_kwargs,
    common_options,
    parse_log_overrides,
    resolve_subcommand_output_dir,
)
from ema.cli.defaults import DEFAULTS

log = logging.getLogger(__name__)


def _list_pdui_strategies() -> list[str]:
    try:
        from ema.quantification.strategies import list_pdui_strategies
        return list(list_pdui_strategies())
    except ImportError:
        return []


def _list_strategies_callback(ctx, param, value):
    if value:
        for s in _list_pdui_strategies():
            click.echo(s)
        ctx.exit(0)


@click.command(name="length")
@click.option("--list-strategies", is_flag=True, default=False,
              is_eager=True, expose_value=False,
              callback=_list_strategies_callback,
              help="Print available length strategies and exit.")
@common_options(output_default="switch_out")
@click.option("--h5ad", "-i", "h5ad", multiple=True, required=True,
              type=click.Path(exists=True, dir_okay=False))
@click.option("--gtf", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Required when --isoform-agg=per_isoform.")
@click.option("--cluster-pairs", "cluster_pairs", type=str, default=None)
@click.option("--cluster-key", "cluster_key", type=str, default="leiden")
@click.option("--strategy", "-s", "strategy", type=str, default="classic")
@click.option("--isoform-agg", "isoform_agg",
              type=click.Choice(["per_gene", "per_isoform"]),
              default=DEFAULTS["isoform-agg"], show_default=True,
              help="Aggregation level — must match strategy vocabulary "
                   "(per_gene collapses isoforms, per_isoform keeps them).")
@click.option("--isoform-collapse", "isoform_collapse",
              type=click.Choice(["none", "mean", "majority"]),
              default=DEFAULTS["isoform-collapse"], show_default=True,
              help="How to collapse multiple isoforms when --isoform-agg=per_gene "
                   "and the strategy tracks isoforms internally.")
@click.option("--pdui-pseudocount", "pdui_pseudocount", type=float,
              default=DEFAULTS["pdui-pseudocount"], show_default=True,
              help="Pseudocount added to counts before PDUI/entropy computation. "
                   "Default 0.0 (original behaviour). Use 1.0 to avoid NaN on "
                   "zero-count cells.")
@click.pass_context
def length(ctx: click.Context, **kwargs) -> None:
    """3'UTR shortening / lengthening quantification (PDUI variants)."""
    # Merge YAML config first; CLI flags (explicitly set) always win.
    apply_yaml_to_kwargs(ctx, kwargs)

    valid = _list_pdui_strategies()
    if valid and kwargs["strategy"] not in valid:
        raise click.BadParameter(
            f"Invalid --strategy {kwargs['strategy']!r}. Available: {', '.join(valid)}",
        )

    user_explicit_output = (
        ctx.get_parameter_source("output").name == "COMMANDLINE"
    )
    out_dir = resolve_subcommand_output_dir(
        kwargs["output"],
        user_explicit=user_explicit_output,
        source_paths=list(kwargs["h5ad"]),
        subdir="switch_length",
    )
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
        log.info("ema switch length: strategy=%s", kwargs["strategy"])
        from ema.switch_test.runner import run_length
        pdui_df, last_adata = run_length(
            h5ad_paths=list(kwargs["h5ad"]),
            gtf=kwargs["gtf"],
            output_dir=str(out_dir),
            cluster_pairs=kwargs["cluster_pairs"],
            cluster_key=kwargs["cluster_key"],
            strategy=kwargs["strategy"],
            isoform_agg=kwargs["isoform_agg"],
            isoform_collapse=kwargs["isoform_collapse"],
            threads=kwargs["threads"],
            pseudocount=kwargs["pdui_pseudocount"],
        )

        from ema.cli.common import parse_plot_engines
        from ema.viz.pipeline_hooks import render_switch_length_outputs
        render_switch_length_outputs(
            out_dir=out_dir,
            last_adata=last_adata,
            pdui_df=pdui_df,
            cluster_key=kwargs["cluster_key"],
            engines=parse_plot_engines(
                kwargs.get("plot_engine", "both"),
                kwargs.get("no_plots", False),
            ),
        )
    finally:
        teardown_logging()
