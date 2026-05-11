"""`ema switch match` — cross-dataset cluster matching."""
from __future__ import annotations

import logging

import click

from ema.cli.common import common_options, parse_log_overrides, resolve_output_dir
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
def match(**kwargs) -> None:
    """Cross-dataset cluster matching (marker_overlap / mnn / jaccard)."""
    valid = _list_match_strategies()
    if valid and kwargs["strategy"] not in valid:
        raise click.BadParameter(
            f"Invalid --strategy {kwargs['strategy']!r}. Available: {', '.join(valid)}",
        )

    out_dir = resolve_output_dir(kwargs["output"])
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
        strat = get_match_strategy(kwargs["strategy"])
        df = strat.match(
            list(kwargs["h5ad"]),
            [str(i) for i in range(len(kwargs["h5ad"]))],  # synthetic ds ids
            n_top_markers=kwargs["n_top_markers"],
            n_jobs=kwargs["threads"] or -1,
        )
        out_file = out_dir / "cluster_match.tsv"
        df.to_csv(out_file, sep="\t")
        log.info("Wrote %s (%d rows)", out_file, len(df))

        # Viz wire-in: render cluster_match_sankey + match_confidence figures
        try:
            from ema.cli.common import parse_plot_engines
            from ema.viz import render_all
            engines = parse_plot_engines(
                kwargs.get("plot_engine", "both"),
                kwargs.get("no_plots", False),
            )
            if engines and df is not None and not df.empty:
                figs_dir = out_dir / "figures"
                figs_dir.mkdir(parents=True, exist_ok=True)
                render_all(
                    plot_type="cluster_match_sankey",
                    data=df,
                    output_basepath=figs_dir / "cluster_match_sankey",
                    engines=engines,
                )
                render_all(
                    plot_type="match_confidence",
                    data=df,
                    output_basepath=figs_dir / "match_confidence",
                    engines=engines,
                )
                log.info("ema switch match: viz figures written to %s", figs_dir)
        except Exception as _viz_exc:
            log.warning("ema switch match: viz rendering failed (non-fatal): %s", _viz_exc)
    finally:
        teardown_logging()
