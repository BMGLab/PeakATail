"""`ema run` — full pipeline subcommand.

This module is intentionally thin. The pipeline body lives in
``ema.main`` and is not modified by the CLI overhaul.

Click options are generated from :class:`ema.cli.config_schema.RunConfig`.
The hand-rolled ``@click.option`` block was deleted in the centralisation
refactor; adding a new flag now means adding a single dataclass field.
The fields handled by ``common_options()`` (``--threads``, ``--config``,
``--output``, ``--verbose``, ``--quiet``, ``--log-level``,
``--no-log-file``, ``--no-progress``) are skipped from the auto-generator
to avoid double-registration.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import click

from ema.cli.common import common_options, parse_log_overrides, resolve_output_dir, parse_plot_engines
from ema.cli.config_schema import RunConfig, click_options_from_schema
from ema.cli.defaults import DEFAULTS
from ema.cli.yaml_loader import load_run_yaml, RunYamlError

log = logging.getLogger(__name__)


# Fields handled by common_options() or the bare-`ema` group; skip them
# from the schema-driven generator so the option is only registered once.
_COMMON_OPTION_FIELDS = (
    "config", "output", "threads",
    "verbose", "quiet", "log_level", "no_log_file", "no_progress",
)


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
@click_options_from_schema(RunConfig, skip=_COMMON_OPTION_FIELDS)
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

    # Lift --bam-dir / --bam-files into a synthetic single-dataset YAML
    # if no datasets list.
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

    # Identify which CLI params the user explicitly typed.  This drives the
    # schema-based override path: only user-set kwargs are forwarded to
    # pipeline_run(), so YAML values (or schema defaults) win for everything
    # the user did not touch on the command line.
    _ctx = click.get_current_context()
    _user_set = _user_supplied_params(_ctx)

    # Backstop defaults so the wizard / minimal CLI invocations do not
    # silently crash deep in pysam with "expected bytes, NoneType found"
    # when the user did not supply seqlen/cb_len/barcode_tag (these are
    # required by ema.countmatrix.read.read_check). No silent fallbacks:
    # a one-line WARNING tells the user the value being used.
    _BACKSTOP = {
        "seqlen": (150, "--seq-len / yaml `seqlen`"),
        "cb_len": (16, "--cb-len / yaml `cb_len`"),
        "barcode_tag": ("CB", "--barcode-tag / yaml `barcode_tag`"),
    }
    for _k, (_default, _flag) in _BACKSTOP.items():
        if cfg.get(_k) is None:
            cfg[_k] = _default
            log.warning(
                "%s not set; defaulting to %r. Set it explicitly to suppress this warning.",
                _flag, _default,
            )

    # Wire --threads into the ResourceManager singleton.  Without this,
    # the value the user supplied would silently be ignored: get_resource_manager()
    # only reads ema.config.args.threads (set by the legacy argparse shim, which
    # never sees Click flags).  Reset the singleton and re-create it with the
    # explicit ceiling so downstream get_n_jobs() calls honour it.
    # NOTE: the log.info() for this is emitted after setup_logging() below so
    # it actually reaches the per-run log file.
    if kwargs.get("threads") is not None:
        from ema.utils import reset_resource_manager
        from ema.utils.resource_manager import ResourceManager
        import ema.utils as _utils_mod
        reset_resource_manager()
        _utils_mod._RM_INSTANCE = ResourceManager(user_max_threads=kwargs["threads"])

    # Resolve output dir (with timestamp).
    # Precedence: CLI --output (when explicitly set) > YAML output_dir > default.
    # Click resolves --output to an absolute path, so compare basenames.
    _output_base = kwargs["output"]
    _output_base_name = Path(_output_base).name if _output_base else ""
    if _output_base_name == DEFAULTS["output"] and cfg.get("output_dir"):
        # CLI --output left at default but YAML provided one — honour YAML.
        _output_base = cfg["output_dir"]
    out_dir = resolve_output_dir(_output_base)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Logging + progress.
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
        log.info("Output directory: %s", out_dir)
        if kwargs.get("threads") is not None:
            log.info("ResourceManager: --threads=%d (absolute ceiling)", kwargs["threads"])

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

        # Resolve --plot-engine / --no-plots into an engine list for the pipeline.
        _plot_engines = parse_plot_engines(
            kwargs.get("plot_engine", "both"),
            no_plots=kwargs.get("no_plots", False),
        )

        from ema.progress import ProgressManager
        with ProgressManager(disable=kwargs["no_progress"]) as pm:
            # Hand off to the pipeline. All algorithm code is in ema.main.
            from ema.main import run as pipeline_run
            pipeline_run(
                cfg=cfg, out_dir=out_dir, progress=pm,
                plot_engines=_plot_engines,
                **_pipeline_kwargs(kwargs, user_set=_user_set),
            )
    finally:
        # Release the multiprocessing.Manager subprocess spawned by
        # setup_logging so the parent interpreter can shut down cleanly.
        teardown_logging()


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


def _user_supplied_params(ctx: click.Context | None) -> set[str]:
    """Return the set of Click parameter names the user typed on the CLI.

    Uses ``ctx.get_parameter_source(name) == ParameterSource.COMMANDLINE``.
    Without this, value-equals-default looks identical to user-typed-default
    and we silently drop overrides like ``--min-read 1500`` when 1500 is
    also the default.

    Returns an empty set when no Click context is active (e.g. the function
    is called from a unit test that didn't go through Click).  Callers
    should fall back to value comparison only when this set is empty.
    """
    if ctx is None:
        return set()
    try:
        from click.core import ParameterSource
    except ImportError:  # pragma: no cover -- Click >=8.0 ships ParameterSource
        return set()
    out: set[str] = set()
    for name in ctx.params:
        try:
            src = ctx.get_parameter_source(name)
        except Exception:
            continue
        if src == ParameterSource.COMMANDLINE:
            out.add(name)
    return out


def _pipeline_kwargs(kwargs: dict, user_set: set[str] | None = None) -> dict:
    """Translate Click kwargs into the keyword set ema.main.run() expects.

    Two responsibilities:
    1. Strip CLI-only keys (config, output, bam_dir, ...) and not-yet-wired
       flags (ip_filter, benchmark, ...) — they have no place in args.
    2. Drop any kwarg the user did NOT explicitly set on the CLI. Without
       this, ``--pas-gap``'s default of 100 would clobber a YAML
       ``pas_gap: 200`` in the bridge layer.

    The "user did not set it" check uses Click's ParameterSource (passed in
    via ``user_set``) so user-typed-default and Click-filled-default are
    distinguishable.  Without ``user_set`` (e.g. unit tests calling this
    directly), we fall back to a value-equals-default heuristic.
    """
    drop = {
        "config", "output", "bam_dir", "bam_files",
        "verbose", "quiet", "log_level", "no_log_file", "no_progress",
        "list_strategies_flag",
        # --threads is consumed by run() itself (wired into ResourceManager
        # before pipeline dispatch); no need to forward downstream.
        "threads",
        # Not-yet-wired flags (warned above) — strip so they don't end up
        # in args via _kwarg_to_args_map. Keeping the warn-loud-but-do-nothing
        # behavior visible in one place.
        "ip_filter", "genome_fasta", "annot_filter", "ip_a_stretch",
        "benchmark", "validate_db",
        # Plot flags are resolved into plot_engines before pipeline_run is called;
        # the raw strings should not be forwarded into the legacy arg bridge.
        "plot_engine", "plot_format", "no_plots",
    }
    out: dict = {}
    for k, v in kwargs.items():
        if k in drop:
            continue
        if user_set is not None:
            if k in user_set:
                out[k] = v
            continue
        # Fallback path: no Click context (e.g. unit test) -- best-effort.
        _default = DEFAULTS.get(k.replace("_", "-"))
        if v == _default:
            continue
        out[k] = v
    return out
