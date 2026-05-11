"""`ema switch diff` — differential APA test across cluster pairs."""
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


def _list_diff_strategies() -> list[str]:
    try:
        from ema.switch_test.strategies import list_diff_strategies
        return list(list_diff_strategies())
    except ImportError:
        return []


def _list_strategies_callback(ctx, param, value):
    if value:
        for s in _list_diff_strategies():
            click.echo(s)
        ctx.exit(0)


@click.command(name="diff")
@click.option("--list-strategies", is_flag=True, default=False,
              is_eager=True, expose_value=False,
              callback=_list_strategies_callback,
              help="Print available diff strategies and exit.")
@common_options(output_default="switch_out")
@click.option("--h5ad", "-i", "h5ad", multiple=True, required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="Per-dataset clusters.h5ad. Repeat for multi-dataset.")
@click.option("--pasbed", type=click.Path(exists=True, dir_okay=False),
              default=None, help="Optional PAS BED for context.")
@click.option("--gtf", type=click.Path(exists=True, dir_okay=False), default=None)
@click.option("--cluster-pairs", "cluster_pairs", type=str, default=None,
              help="`c1,c2;c3,c4` — limit to specific pairs.")
@click.option("--cluster-key", "cluster_key", type=str, default="leiden")
@click.option("--marker-top-n", "marker_top_n", type=int, default=DEFAULTS["marker-top-n"])
@click.option("--marker-method", "marker_method", type=str, default=DEFAULTS["marker-method"])
@click.option("--strategy", "-s", "strategy", type=str, default="fisher",
              show_default=True,
              help="Differential APA strategy (run --list-strategies to see).")
@click.option("--fdr", "fdr", type=float, default=DEFAULTS["fdr"])
@click.option("--per-worker-mb", "per_worker_mb", type=int, default=DEFAULTS["per-worker-mb"])
@click.option("--min-cells-per-group", "min_cells_per_group", type=int,
              default=DEFAULTS["min-cells-per-group"],
              help="Minimum cells per group for a PAS to enter differential testing.")
@click.option("--log2fc-thresh", "log2fc_thresh", type=float,
              default=DEFAULTS["log2fc-thresh"],
              help="log2 fold-change threshold drawn on the volcano plot. Default 1.0.")
@click.pass_context
def diff(ctx: click.Context, **kwargs) -> None:
    """Differential APA test (Fisher / NB regression) across cluster pairs."""
    # Merge YAML config first; CLI flags (explicitly set) always win.
    apply_yaml_to_kwargs(ctx, kwargs)

    valid = _list_diff_strategies()
    if valid and kwargs["strategy"] not in valid:
        raise click.BadParameter(
            f"Invalid --strategy {kwargs['strategy']!r}. Available: {', '.join(valid)}",
        )

    # When --output is left at its default and the input h5ads come from a
    # peakatail run dir, route output INSIDE that run dir so the run stays
    # self-contained: peakatail_runs/<run>/switch_diff_<ts>/
    user_explicit_output = (
        ctx.get_parameter_source("output").name == "COMMANDLINE"
    )
    out_dir = resolve_subcommand_output_dir(
        kwargs["output"],
        user_explicit=user_explicit_output,
        source_paths=list(kwargs["h5ad"]),
        subdir="switch_diff",
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
        log.info("ema switch diff: %d h5ad input(s); strategy=%s",
                 len(kwargs["h5ad"]), kwargs["strategy"])
        from ema.switch_test.runner import run_diff
        pair_results = run_diff(
            h5ad_paths=list(kwargs["h5ad"]),
            pasbed=kwargs["pasbed"],
            gtf=kwargs["gtf"],
            output_dir=str(out_dir),
            cluster_pairs=kwargs["cluster_pairs"],
            cluster_key=kwargs["cluster_key"],
            marker_top_n=kwargs["marker_top_n"],
            marker_method=kwargs["marker_method"],
            strategy=kwargs["strategy"],
            fdr=kwargs["fdr"],
            threads=kwargs["threads"],
            per_worker_mb=kwargs["per_worker_mb"],
            min_cells_per_group=kwargs["min_cells_per_group"],
        )

        # Visualisation lives in ema.viz.pipeline_hooks (one entry point per
        # CLI command).  Failures are warned, never raised.
        from ema.cli.common import parse_plot_engines
        from ema.viz.pipeline_hooks import render_switch_diff_outputs
        render_switch_diff_outputs(
            out_dir=out_dir,
            pair_results=pair_results,
            fdr=kwargs["fdr"],
            log2fc_thresh=kwargs["log2fc_thresh"],
            engines=parse_plot_engines(
                kwargs.get("plot_engine", "both"),
                kwargs.get("no_plots", False),
            ),
        )
    finally:
        teardown_logging()
