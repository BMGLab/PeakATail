"""`ema run` — full pipeline subcommand.

This is a thin wrapper that:
1. Parses CLI flags (or loads YAML config + applies overrides)
2. Resolves the timestamped output dir
3. Sets up logging + ProgressManager
4. Calls into ema.main.run() with the resolved kwargs

The actual algorithm (peak calling → atlas snap → cluster → match) lives
in ema/main.py and is not changed by this CLI overhaul.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import click

from ema.cli.common import common_options, parse_log_overrides, resolve_output_dir
from ema.cli.defaults import DEFAULTS
from ema.cli.yaml_loader import load_run_yaml, RunYamlError

log = logging.getLogger(__name__)


def _list_strategies_and_exit() -> None:
    """Print every registered strategy across all 5 registries, then exit."""
    from ema.strategies import list_strategies as list_peak_strategies
    from ema.clustering.registry import list_clustering_strategies
    from ema.clustering.cross_dataset import list_match_strategies

    click.echo("Available strategies:")
    click.echo("  Peak calling   (--peak-strategy):   " + ", ".join(list_peak_strategies()))
    click.echo("  Clustering     (--cluster-method):  " + ", ".join(list_clustering_strategies()))
    click.echo("  Cluster match  (--match-method):    " + ", ".join(list_match_strategies()))
    sys.exit(0)


@click.command(name="run")
@common_options(output_default=DEFAULTS["output"])
# ─── inputs ─────────────────────────────────────────────────────────
@click.option("--bam-dir", "bam_dir", type=click.Path(exists=True),
              default=DEFAULTS["bam-dir"], help="Single-BAM convenience.")
@click.option("--bam-files", "bam_files", type=str, default=DEFAULTS["bam-files"],
              help="Comma-separated multi-BAM list.")
@click.option("--gtf", "gtf", type=click.Path(exists=True, dir_okay=False),
              default=DEFAULTS["gtf"], help="GTF file path.")
@click.option("--atlas", "atlas", type=click.Path(exists=True, dir_okay=False),
              default=DEFAULTS["atlas"], help="Reference PAS atlas BED.")
@click.option("--atlas-distance", "atlas_distance", type=int,
              default=DEFAULTS["atlas-distance"], help="Atlas snap distance (bp).")
# ─── read processing ────────────────────────────────────────────────
@click.option("--seq-len", "seq_len", type=int, default=DEFAULTS["seq-len"])
@click.option("--cb-len", "cb_len", type=int, default=DEFAULTS["cb-len"])
@click.option("--barcode-tag", "barcode_tag", type=str, default=DEFAULTS["barcode-tag"])
# ─── concurrency ────────────────────────────────────────────────────
@click.option("--bam-threads", "bam_threads", type=int, default=DEFAULTS["bam-threads"])
@click.option("--pipeline", "pipeline", is_flag=True, default=DEFAULTS["pipeline"])
@click.option("--batch-size", "batch_size", type=int, default=DEFAULTS["batch-size"])
@click.option("--tiles", "tiles", is_flag=True, default=DEFAULTS["tiles"])
@click.option("--tile-size", "tile_size", type=int, default=DEFAULTS["tile-size"])
@click.option("--tile-overlap", "tile_overlap", type=int, default=DEFAULTS["tile-overlap"])
# ─── peak calling ───────────────────────────────────────────────────
@click.option("--peak-strategy", "peak_strategy",
              type=str,  # validated lazily via _validate_strategy_choice
              default=DEFAULTS["peak-strategy"])
@click.option("--lambda-window", "lambda_window", type=int, default=DEFAULTS["lambda-window"])
@click.option("--lambda-method", "lambda_method", type=str, default=DEFAULTS["lambda-method"])
@click.option("--lambda-fold-change", "lambda_fold_change", type=float,
              default=DEFAULTS["lambda-fold-change"])
@click.option("--max-pas", "max_pas", type=int, default=DEFAULTS["max-pas"])
@click.option("--smoothing-window", "smoothing_window", type=int,
              default=DEFAULTS["smoothing-window"])
@click.option("--min-prominence", "min_prominence", type=float, default=DEFAULTS["min-prominence"])
@click.option("--dynamic-threshold", "dynamic_threshold", is_flag=True,
              default=DEFAULTS["dynamic-threshold"])
@click.option("--floor-threshold", "floor_threshold", type=int, default=DEFAULTS["floor-threshold"])
@click.option("--pas-gap", "pas_gap", type=int, default=DEFAULTS["pas-gap"])
# ─── filters ────────────────────────────────────────────────────────
@click.option("--ip-filter", "ip_filter", is_flag=True, default=DEFAULTS["ip-filter"])
@click.option("--genome-fasta", "genome_fasta", type=click.Path(exists=True),
              default=DEFAULTS["genome-fasta"])
@click.option("--annot-filter", "annot_filter", is_flag=True, default=DEFAULTS["annot-filter"])
@click.option("--ip-a-stretch", "ip_a_stretch", type=int, default=DEFAULTS["ip-a-stretch"])
@click.option("--min-pas-per-cell", "min_pas_per_cell", type=int,
              default=DEFAULTS["min-pas-per-cell"])
@click.option("--min-read", "min_read", type=int, default=DEFAULTS["min-read"])
@click.option("--min-cells", "min_cells", type=int, default=DEFAULTS["min-cells"])
# ─── annotation ─────────────────────────────────────────────────────
@click.option("--max-gene-distance", "max_gene_distance", type=int,
              default=DEFAULTS["max-gene-distance"])
@click.option("--utr-multiplier", "utr_multiplier", type=float,
              default=DEFAULTS["utr-multiplier"])
@click.option("--include-extended", "include_extended", is_flag=True,
              default=DEFAULTS["include-extended"])
# ─── clustering ─────────────────────────────────────────────────────
@click.option("--cluster-method", "cluster_method",
              type=str, default=DEFAULTS["cluster-method"])
@click.option("--resolution", "resolution", type=float, default=DEFAULTS["resolution"])
@click.option("--n-pcs", "n_pcs", type=int, default=DEFAULTS["n-pcs"])
@click.option("--external-clusters", "external_clusters", type=click.Path(exists=True),
              default=DEFAULTS["external-clusters"])
@click.option("--random-seed", "random_seed", type=int, default=DEFAULTS["random-seed"])
# ─── cross-dataset matching ─────────────────────────────────────────
@click.option("--match-method", "match_method",
              type=str, default=DEFAULTS["match-method"])
@click.option("--n-top-markers", "n_top_markers", type=int, default=DEFAULTS["n-top-markers"])
# ─── validation ─────────────────────────────────────────────────────
@click.option("--benchmark", "benchmark", is_flag=True, default=DEFAULTS["benchmark"])
@click.option("--validate-db", "validate_db", type=click.Path(exists=True),
              default=DEFAULTS["validate-db"])
# ─── helpers ────────────────────────────────────────────────────────
@click.option("--list-strategies", "list_strategies_flag", is_flag=True, default=False,
              help="Print available strategies and exit.")
def run(**kwargs) -> None:
    """Run the full PeakATail pipeline."""
    if kwargs.get("list_strategies_flag"):
        _list_strategies_and_exit()

    # Validate strategy choices NOW (after registries are populated by imports).
    _validate_strategy_choice(
        "peak-strategy", kwargs["peak_strategy"],
        _safe_list("ema.strategies", "list_strategies"),
    )
    _validate_strategy_choice(
        "cluster-method", kwargs["cluster_method"],
        _safe_list("ema.clustering.registry", "list_clustering_strategies"),
    )
    _validate_strategy_choice(
        "match-method", kwargs["match_method"],
        _safe_list("ema.clustering.cross_dataset", "list_match_strategies"),
    )

    # Load YAML if given; CLI flags override individual keys.
    cfg: dict = {}
    if kwargs["config"] is not None:
        try:
            cfg = load_run_yaml(kwargs["config"])
        except RunYamlError as e:
            raise click.BadParameter(str(e), param_hint="--config")

    # Lift --bam-dir into a synthetic single-dataset YAML if no datasets.
    if "datasets" not in cfg:
        if kwargs.get("bam_dir"):
            cfg["datasets"] = [{
                "id": "default", "merge_strategy": "none", "bams": [kwargs["bam_dir"]],
            }]
        elif kwargs.get("bam_files"):
            bams = [b.strip() for b in kwargs["bam_files"].split(",") if b.strip()]
            cfg["datasets"] = [{"id": "default", "merge_strategy": "none", "bams": bams}]
        else:
            raise click.UsageError(
                "No input. Provide --config <yaml>, --bam-dir <bam>, or --bam-files <a.bam,b.bam>."
            )

    # CLI flags override (only non-default values win).
    _apply_cli_overrides(cfg, kwargs)

    # Resolve output dir (with timestamp).
    out_dir = resolve_output_dir(kwargs["output"])
    out_dir.mkdir(parents=True, exist_ok=True)

    # Logging + progress.
    from ema.logging_config import setup_logging
    setup_logging(
        level=(kwargs["log_level"] or "INFO").split(",")[0],
        output_dir=out_dir,
        quiet=kwargs["quiet"],
        verbose_count=kwargs["verbose"],
        no_log_file=kwargs["no_log_file"],
        per_logger_overrides=parse_log_overrides(kwargs["log_level"]),
    )
    log.info("Output directory: %s", out_dir)

    # Warn loudly for flags that exist in the CLI surface but are not yet
    # wired into the pipeline body.  No silent fallbacks: if you set one of
    # these, we tell you it is a no-op so you know the result you are
    # looking at did NOT honour your request.
    _NOT_WIRED = {
        "ip_filter": "--ip-filter (internal-priming filter not integrated into `ema run`)",
        "genome_fasta": "--genome-fasta (only consumed by --ip-filter, currently no-op)",
        "annot_filter": "--annot-filter (annotation filter not integrated into `ema run`)",
        "ip_a_stretch": "--ip-a-stretch (only consumed by --ip-filter, currently no-op)",
        "benchmark": "--benchmark (run scripts/validate_strategies.py instead)",
        "validate_db": "--validate-db (run scripts/validate_strategies.py instead)",
    }
    for _k, _msg in _NOT_WIRED.items():
        _v = kwargs.get(_k)
        _default = DEFAULTS.get(_k.replace("_", "-"))
        if _v not in (None, False, _default):
            log.warning(
                "%s is currently a no-op; the value you supplied will NOT "
                "affect the output.", _msg,
            )

    from ema.progress import ProgressManager
    with ProgressManager(disable=kwargs["no_progress"]) as pm:
        # Hand off to the pipeline. All algorithm code is in ema.main.
        from ema.main import run as pipeline_run
        pipeline_run(cfg=cfg, out_dir=out_dir, progress=pm, **_pipeline_kwargs(kwargs))


def _safe_list(modname: str, fnname: str) -> list[str]:
    """Import the registry list-helper and call it. Returns [] if missing."""
    try:
        mod = __import__(modname, fromlist=[fnname])
        fn = getattr(mod, fnname)
        return list(fn())
    except (ImportError, AttributeError):
        return []


def _validate_strategy_choice(flag: str, value: str, valid: list[str]) -> None:
    if not valid:
        return  # registry not populated yet; skip
    if value not in valid:
        raise click.BadParameter(
            f"Invalid value for --{flag}: {value!r}. Available: {', '.join(valid)}",
            param_hint=f"--{flag}",
        )


def _apply_cli_overrides(cfg: dict, kwargs: dict) -> None:
    """Per spec: any non-None CLI flag overrides the corresponding YAML key."""
    mapping = {
        "gtf": "gtf", "atlas": "atlas", "atlas_distance": "atlas_distance",
        "seq_len": "seqlen", "cb_len": "cb_len", "barcode_tag": "barcode_tag",
        "min_read": "min_read", "min_cells": "min_cells",
        "min_pas_per_cell": "min_pas_per_cell", "pas_gap": "pas_gap",
    }
    for cli_key, yaml_key in mapping.items():
        v = kwargs.get(cli_key)
        if v is not None and v != DEFAULTS.get(cli_key.replace("_", "-")):
            cfg[yaml_key] = v


def _pipeline_kwargs(kwargs: dict) -> dict:
    """Translate Click kwargs into the keyword set ema.main.run() expects."""
    # Strip CLI-only keys
    drop = {
        "config", "output", "bam_dir", "bam_files",
        "verbose", "quiet", "log_level", "no_log_file", "no_progress",
        "list_strategies_flag",
        # Not-yet-wired flags (warned above) — strip so they don't end up
        # in args via _kwarg_to_args_map. Keeping the warn-loud-but-do-nothing
        # behavior visible in one place.
        "ip_filter", "genome_fasta", "annot_filter", "ip_a_stretch",
        "benchmark", "validate_db",
    }
    return {k: v for k, v in kwargs.items() if k not in drop}
