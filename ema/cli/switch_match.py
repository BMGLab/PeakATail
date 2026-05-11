"""`ema switch match` — cross-dataset cluster matching."""
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


def _list_match_strategies() -> list[str]:
    try:
        from ema.clustering.cross_dataset import list_match_strategies
        return list(list_match_strategies())
    except ImportError:
        return []


def _list_strategies_callback(ctx, param, value):
    if value:
        for s in _list_match_strategies():
            click.echo(s)
        ctx.exit(0)


@click.command(name="match")
@click.option("--list-strategies", is_flag=True, default=False,
              is_eager=True, expose_value=False,
              callback=_list_strategies_callback,
              help="Print available match strategies and exit.")
@common_options(output_default="switch_out")
@click.option("--h5ad", "-i", "h5ad", multiple=True, required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="Per-dataset clusters.h5ad files.")
@click.option("--strategy", "-s", "strategy", type=str, default=DEFAULTS["match-method"])
@click.option("--n-top-markers", "n_top_markers", type=int, default=DEFAULTS["n-top-markers"])
@click.option("--mnn-components", "mnn_components", type=int,
              default=DEFAULTS["mnn-components"],
              help="SVD components for MNN shared embedding (--strategy mnn). Default 30.")
@click.option("--mnn-k-neighbors", "mnn_k_neighbors", type=int,
              default=DEFAULTS["mnn-k-neighbors"],
              help="Nearest neighbours for MNN search (--strategy mnn). Default 10.")
@click.pass_context
def match(ctx: click.Context, **kwargs) -> None:
    """Cross-dataset cluster matching (marker_overlap / mnn / jaccard)."""
    # Merge YAML config first; CLI flags (explicitly set) always win.
    apply_yaml_to_kwargs(ctx, kwargs)

    valid = _list_match_strategies()
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
        subdir="switch_match",
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
        log.info("ema switch match: strategy=%s, %d datasets",
                 kwargs["strategy"], len(kwargs["h5ad"]))
        from ema.clustering.cross_dataset import get_match_strategy
        strat = get_match_strategy(
            kwargs["strategy"],
            mnn_components=kwargs["mnn_components"],
            mnn_k_neighbors=kwargs["mnn_k_neighbors"],
        )
        df = strat.match(
            list(kwargs["h5ad"]),
            [str(i) for i in range(len(kwargs["h5ad"]))],  # synthetic ds ids
            n_top_markers=kwargs["n_top_markers"],
            n_jobs=kwargs["threads"] or -1,
        )
        out_file = out_dir / "cluster_match.tsv"
        df.to_csv(out_file, sep="\t")
        log.info("Wrote %s (%d rows)", out_file, len(df))

        from ema.cli.common import parse_plot_engines
        from ema.viz.pipeline_hooks import render_switch_match_outputs
        render_switch_match_outputs(
            out_dir=out_dir,
            df=df,
            engines=parse_plot_engines(
                kwargs.get("plot_engine", "both"),
                kwargs.get("no_plots", False),
            ),
        )
    finally:
        teardown_logging()
