"""`peakatail switch trend` — ordered-stage APA trend from a PDUI long table.

Reads a per-cluster PDUI/proportion/entropy table (as produced by
``peakatail switch length``), treats the clusters as an ORDERED progression given by
``--stage-order``, and reports the trend (slope + Spearman monotonicity +
direction) overall and per gene. See :mod:`ema.switch_test.trend`.

Two extra modes on top of the original single-table/``--value-col`` path:

* ``--metric {pdui,proportion,entropy}`` — a convenience selector that picks
  the right value column (and, with ``--length-dir``, the right source TSV)
  for one of ``peakatail switch length``'s three strategies. ``proportion`` is
  derived (per-(gene,cell) distal-PAS usage fraction) via
  :func:`ema.switch_test.trend.derive_distal_fraction`.
* ``--combine pdui,proportion,entropy`` — runs the per-gene trend for EACH
  listed metric and emits a consensus (shortening/lengthening/ambiguous)
  across them, written to ``trend_combined.tsv`` /
  ``trend_combined_summary.json`` instead of the single-metric outputs.

With neither ``--metric`` nor ``--combine`` passed, behavior is identical to
the original implementation (default metric = pdui, read from ``--pdui``
with ``--value-col`` verbatim).
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import click

from ema.cli.common import common_options, parse_log_overrides

log = logging.getLogger(__name__)

# strategy name -> (source filename under `peakatail switch length`, subdir to
# prefer when searching --length-dir, default value column on the raw table).
# "proportion" has no fixed value column: it is DERIVED per (gene, cell).
_METRIC_FILES = {
    "pdui": ("pdui_classic.tsv", "classic", "pdui"),
    "proportion": ("proportion.tsv", "proportion", None),
    "entropy": ("entropy_shannon.tsv", "shannon", "normalized_entropy"),
}


def _explicit(ctx: click.Context, name: str) -> bool:
    """True iff the user typed ``name`` on the command line (not a default)."""
    try:
        from click.core import ParameterSource
        return ctx.get_parameter_source(name) == ParameterSource.COMMANDLINE
    except Exception:
        return False


def _find_in_length_dir(length_dir: str, filename: str, subdir: str) -> str:
    """Locate ``filename`` under a `length/<celltype>/` directory.

    Prefers ``length_dir/<subdir>/**/filename`` (matches the
    ``length/<CT>/<strategy>/...`` layout ``peakatail switch length`` writes) and
    falls back to searching the whole ``length_dir`` if that subdir isn't
    present.
    """
    base = Path(length_dir)
    scoped = base / subdir
    search_root = scoped if scoped.is_dir() else base
    matches = sorted(search_root.rglob(filename))
    if not matches:
        raise click.UsageError(
            f"Could not find {filename!r} under {search_root} "
            f"(--length-dir {length_dir!r})."
        )
    if len(matches) > 1:
        log.warning(
            "multiple %s found under %s; using %s", filename, search_root, matches[0],
        )
    return str(matches[0])


def _resolve_metric_path(kwargs: dict, metric: str) -> str:
    """Resolve the source TSV for one metric from CLI kwargs.

    Precedence: the metric's own explicit ``--<metric>-tsv``/``--pdui`` >
    ``--pdui`` used as a generic "the input tsv" path (the pre-existing
    power-user convention, e.g. ``--pdui shannon.tsv --value-col
    normalized_entropy``) > ``--length-dir`` auto-discovery.
    """
    explicit_map = {
        "pdui": kwargs.get("pdui_tsv"),
        "proportion": kwargs.get("proportion_tsv"),
        "entropy": kwargs.get("entropy_tsv"),
    }
    path = explicit_map[metric]
    if path is None and metric != "pdui" and kwargs.get("pdui_tsv"):
        path = kwargs["pdui_tsv"]
    if path is None and kwargs.get("length_dir"):
        filename, subdir, _ = _METRIC_FILES[metric]
        path = _find_in_length_dir(kwargs["length_dir"], filename, subdir)
    if path is None:
        raise click.UsageError(
            f"metric={metric!r} needs a source TSV: pass --{metric}-tsv "
            "(or --pdui as a generic input path) or --length-dir <CT length dir>."
        )
    return path


def _load_metric_frame(kwargs: dict, metric: str, ctx: click.Context):
    """Return ``(df, stage_col, gene_col, value_col)`` ready for trend fitting."""
    import pandas as pd
    from ema.switch_test.trend import derive_distal_fraction

    path = _resolve_metric_path(kwargs, metric)
    raw = pd.read_csv(path, sep="\t")
    stage_col = kwargs["stage_col"]
    gene_col = kwargs["gene_col"]

    if metric == "proportion":
        derived = derive_distal_fraction(raw, gene_col=gene_col, stage_col=stage_col)
        return derived, "cluster", "gene_id", "value"

    _, _, default_value_col = _METRIC_FILES[metric]
    value_col = kwargs["value_col"] if _explicit(ctx, "value_col") else default_value_col
    return raw, stage_col, gene_col, value_col


@click.command(name="trend")
@click.pass_context
@common_options(include_output=True)
@click.option("--pdui", "pdui_tsv", required=False, default=None,
              type=click.Path(exists=True, dir_okay=False),
              help="Input long TSV from `peakatail switch length` (has stage + value "
                   "columns). Required for the default pdui metric unless "
                   "--length-dir is given; also usable as a generic 'the input "
                   "tsv' path for --metric entropy/proportion, and by every "
                   "metric listed in --combine that has no more specific "
                   "--<metric>-tsv / --length-dir.")
@click.option("--metric", "metric", default=None,
              type=click.Choice(["pdui", "proportion", "entropy"], case_sensitive=False),
              help="Convenience metric selector. Sets the value column (and, "
                   "with --length-dir, the source TSV) to sane defaults: "
                   "pdui -> 'pdui' on the classic table; entropy -> "
                   "'normalized_entropy' on the shannon table; proportion -> a "
                   "derived per-(gene,cell) distal-PAS-usage fraction from the "
                   "per-PAS proportion table. Unset (default) = pdui, matching "
                   "prior behavior exactly.")
@click.option("--length-dir", "length_dir", default=None,
              type=click.Path(exists=True, file_okay=False),
              help="A `length/<celltype>/` directory from `peakatail switch length` "
                   "(containing classic/proportion/shannon subdirs); used to "
                   "auto-locate a metric's source TSV when its explicit "
                   "--pdui/--proportion-tsv/--entropy-tsv is not given.")
@click.option("--proportion-tsv", "proportion_tsv", default=None,
              type=click.Path(exists=True, dir_okay=False),
              help="Per-PAS proportion TSV from `peakatail switch length -s "
                   "proportion` (gene_id/cell/rank/proportion/stage columns). "
                   "Used by --metric proportion and by --combine.")
@click.option("--entropy-tsv", "entropy_tsv", default=None,
              type=click.Path(exists=True, dir_okay=False),
              help="Shannon entropy TSV from `peakatail switch length -s shannon` "
                   "(gene_id/cell/normalized_entropy/stage columns). Used by "
                   "--metric entropy and by --combine.")
@click.option("--combine", "combine", default=None,
              help="Comma list of metrics to combine into a consensus trend, "
                   "e.g. 'pdui,proportion,entropy'. Each metric's source is "
                   "located the same way as --metric (its --<metric>-tsv or "
                   "--length-dir). Writes trend_combined.tsv + "
                   "trend_combined_summary.json instead of the single-metric "
                   "outputs.")
@click.option("--combine-min-agree", "combine_min_agree", type=int, default=None,
              help="Minimum number of VOTING metrics (see "
                   "--include-entropy-in-vote) that must agree on "
                   "increasing/decreasing for a gene's consensus to be called "
                   "'shortening'/'lengthening'; otherwise 'ambiguous'. Default: "
                   "all-but-one of the voting metrics.")
@click.option("--include-entropy-in-vote", "include_entropy_in_vote", is_flag=True,
              default=False,
              help="Count entropy's direction toward the shortening/"
                   "lengthening consensus vote. Off by default: rising entropy "
                   "means more even PAS usage, not cleanly 'lengthening' the "
                   "way a rising PDUI/distal-fraction is. Entropy is still "
                   "reported per gene either way.")
@click.option("--stage-order", "stage_order", required=True,
              help="Comma-separated ordered stage labels, e.g. "
                   "'Normal,StageI,StageII,StageIII,StageIV,Met'.")
@click.option("--stage-col", "stage_col", default="cluster", show_default=True,
              help="Column holding the stage/cluster label.")
@click.option("--value-col", "value_col", default="pdui", show_default=True,
              help="Numeric value column (pdui/proportion/entropy); ignored "
                   "for --metric/--combine proportion, which is derived.")
@click.option("--gene-col", "gene_col", default="gene_id", show_default=True,
              help="Gene column for the per-gene trend table (skipped if absent).")
def trend(ctx: click.Context, **kwargs) -> None:
    """Aggregate an ordered-stage APA length trend."""
    from ema.logging_config import setup_logging, teardown_logging

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
        from ema.switch_test.trend import (
            consensus_trend_by_gene,
            consensus_trend_summary,
            pdui_trend,
            pdui_trend_by_gene,
        )

        stage_order = [s.strip() for s in kwargs["stage_order"].split(",") if s.strip()]
        if len(stage_order) < 2:
            raise click.UsageError("--stage-order needs at least 2 stages.")

        if kwargs["combine"]:
            metrics = [m.strip().lower() for m in kwargs["combine"].split(",") if m.strip()]
            bad = [m for m in metrics if m not in _METRIC_FILES]
            if bad:
                raise click.UsageError(
                    f"--combine has unknown metric(s) {bad}; choose from "
                    f"{sorted(_METRIC_FILES)}."
                )
            if len(metrics) < 2:
                raise click.UsageError("--combine needs at least 2 metrics.")

            vote_metrics = [m for m in metrics if m != "entropy"
                             or kwargs["include_entropy_in_vote"]]
            min_agree = kwargs["combine_min_agree"]
            if min_agree is None:
                min_agree = max(1, len(vote_metrics) - 1)

            per_gene_by_metric: dict[str, pd.DataFrame] = {}
            for metric in metrics:
                df, stage_col, gene_col, value_col = _load_metric_frame(kwargs, metric, ctx)
                per_gene_by_metric[metric] = pdui_trend_by_gene(
                    df, stage_order,
                    stage_col=stage_col, value_col=value_col, gene_col=gene_col,
                )
                log.info("combine: %s trend computed for %d genes",
                          metric, len(per_gene_by_metric[metric]))

            combined = consensus_trend_by_gene(
                per_gene_by_metric, vote_metrics=vote_metrics, min_agree=min_agree,
            )
            combined.to_csv(out_dir / "trend_combined.tsv", sep="\t")
            summary = consensus_trend_summary(
                combined, metrics=metrics, vote_metrics=vote_metrics, min_agree=min_agree,
            )
            (out_dir / "trend_combined_summary.json").write_text(
                json.dumps(summary, indent=2, default=str)
            )
            log.info(
                "combined trend (%s; voting=%s; min_agree=%d): %d shortening, "
                "%d lengthening, %d ambiguous (of %d genes)",
                ",".join(metrics), ",".join(vote_metrics), min_agree,
                summary["n_shortening"], summary["n_lengthening"],
                summary["n_ambiguous"], summary["n_genes"],
            )
            click.echo(f"combined trend written to {out_dir}")
            return

        if kwargs["metric"]:
            metric = kwargs["metric"].lower()
            df, stage_col, gene_col, value_col = _load_metric_frame(kwargs, metric, ctx)
        else:
            # Exact original code path: default metric = pdui, --pdui required.
            if not kwargs["pdui_tsv"]:
                raise click.UsageError("--pdui is required (or pass --metric with --length-dir).")
            df = pd.read_csv(kwargs["pdui_tsv"], sep="\t")
            stage_col, gene_col, value_col = (
                kwargs["stage_col"], kwargs["gene_col"], kwargs["value_col"],
            )

        overall = pdui_trend(df, stage_order, stage_col=stage_col, value_col=value_col)
        (out_dir / "length_trend.json").write_text(json.dumps(overall, indent=2, default=str))
        log.info(
            "overall %s trend: slope=%.4g spearman=%.3g direction=%s (over %d stages)",
            value_col, overall["slope"], overall["spearman"],
            overall["direction"], overall["n_stages"],
        )

        if gene_col in df.columns:
            per_gene = pdui_trend_by_gene(
                df, stage_order,
                stage_col=stage_col, value_col=value_col, gene_col=gene_col,
            )
            per_gene.to_csv(out_dir / "length_trend_by_gene.tsv", sep="\t")
            log.info("wrote per-gene trend for %d genes", len(per_gene))
        else:
            log.warning("no %r column — per-gene trend skipped (overall only).", gene_col)
        click.echo(f"trend written to {out_dir}")
    finally:
        teardown_logging()
