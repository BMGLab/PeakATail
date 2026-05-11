"""`ema switch diff` — differential APA test across cluster pairs."""
from __future__ import annotations

import logging

import click

from ema.cli.common import common_options, parse_log_overrides, resolve_output_dir
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
def diff(**kwargs) -> None:
    """Differential APA test (Fisher / NB regression) across cluster pairs."""
    valid = _list_diff_strategies()
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
        )

        # Viz wire-in: render volcano + diff_agreement figures
        try:
            from ema.cli.common import parse_plot_engines
            from ema.viz import render_all
            from pathlib import Path as _Path
            engines = parse_plot_engines(
                kwargs.get("plot_engine", "both"),
                kwargs.get("no_plots", False),
            )
            if engines and pair_results:
                figs_dir = out_dir / "figures"
                figs_dir.mkdir(parents=True, exist_ok=True)
                # Volcano: one plot per pair (use first pair for the main figure)
                for (c1, c2), df in pair_results.items():
                    pair_stem = f"volcano_{c1}_vs_{c2}" if c2 else "volcano_omnibus"
                    render_all(
                        plot_type="volcano",
                        data=df,
                        output_basepath=figs_dir / pair_stem,
                        engines=engines,
                    )
                # Diff agreement: build sig-PAS sets per pair and pass as dict
                # (only meaningful when multiple pairs exist)
                import pandas as _pd
                sig_sets: dict[str, set[str]] = {}
                for (c1, c2), df in pair_results.items():
                    key = f"{c1}_vs_{c2}" if c2 else "omnibus"
                    if "qvalue" in df.columns:
                        sig = set(
                            df.loc[df["qvalue"] < kwargs["fdr"], "pas_id"].astype(str)
                            if "pas_id" in df.columns
                            else df.index[df["qvalue"] < kwargs["fdr"]].astype(str)
                        )
                        sig_sets[key] = sig
                if len(sig_sets) >= 2:
                    render_all(
                        plot_type="diff_agreement",
                        data=sig_sets,
                        output_basepath=figs_dir / "diff_agreement",
                        engines=engines,
                    )
                log.info("ema switch diff: viz figures written to %s", figs_dir)
        except Exception as _viz_exc:
            log.warning("ema switch diff: viz rendering failed (non-fatal): %s", _viz_exc)
    finally:
        teardown_logging()
