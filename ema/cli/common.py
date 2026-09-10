"""Decorator factory that applies the shared options to every subcommand.

Usage:
    @click.command()
    @common_options()
    def my_cmd(threads, output, verbose, quiet, log_level, no_log_file, no_progress, config, **kwargs):
        ...
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Callable

import click

from ema.cli.defaults import DEFAULTS

log = logging.getLogger(__name__)


def common_options(include_output: bool = True, output_default: str | None = None) -> Callable:
    """Apply --threads, -v/-vv, -q, --log-level, --no-log-file, --no-progress,
    --config, and (optionally) --output to a Click command.

    Args:
        include_output: If True, attach `--output / -o` (most subcommands have one).
        output_default: Override DEFAULTS['output'] for this command.
    """
    def decorator(fn: Callable) -> Callable:
        opts = [
            click.option(
                "--threads", "threads",
                type=int, default=DEFAULTS["threads"],
                help="Max parallel workers (auto-detected if not set). "
                     "Respected by ResourceManager as an absolute ceiling.",
            ),
            click.option(
                "-v", "--verbose", "verbose",
                count=True,
                help="Increase verbosity. -v = DEBUG for ema.*; -vv = DEBUG everywhere.",
            ),
            click.option(
                "-q", "--quiet", "quiet",
                is_flag=True, default=DEFAULTS["quiet"],
                help="WARNING and up only. Overrides --verbose.",
            ),
            click.option(
                "--log-level", "log_level",
                type=str, default=DEFAULTS["log-level"],
                help="Explicit logger level (DEBUG/INFO/WARNING/ERROR) "
                     "or `logger.name=LEVEL` (repeatable: comma-separated).",
            ),
            click.option(
                "--no-log-file", "no_log_file",
                is_flag=True, default=DEFAULTS["no-log-file"],
                help="Don't write peakatail_<ts>.log next to the outputs.",
            ),
            click.option(
                "--no-progress", "no_progress",
                is_flag=True, default=DEFAULTS["no-progress"],
                help="Suppress Rich progress bars.",
            ),
            click.option(
                "--config", "-c", "config",
                type=click.Path(exists=True, dir_okay=False, resolve_path=True),
                default=DEFAULTS["config"],
                help="YAML config; CLI flags override individual keys.",
            ),
        ]
        if include_output:
            opts.append(
                click.option(
                    "--output", "-o", "output",
                    # NOTE: resolve_path=False so plain names like ``-o emaout``
                    # stay relative — resolve_output_dir() then prefixes them
                    # with ``peakatail_runs/``.  With resolve_path=True the
                    # resolver sees an absolute path and skips the prefix,
                    # dumping outputs at the project root.
                    type=click.Path(file_okay=False, resolve_path=False),
                    default=output_default if output_default is not None else DEFAULTS["output"],
                    help="Output directory (timestamp suffix added automatically).",
                )
            )
        opts.extend([
            click.option(
                "--plot-engine", "plot_engine", type=str, default="matplotlib",
                help="Engines: 'matplotlib' (default), 'plotly', 'both', 'none', or comma list.",
            ),
            click.option(
                "--plot-format", "plot_format", type=str, default="all",
                help="Restrict output formats. Default 'all' = png+svg+html as appropriate. "
                     "Examples: 'svg' / 'png,svg' / 'html'.",
            ),
            click.option(
                "--no-plots", "no_plots", is_flag=True, default=False,
                help="Disable all plotting (alias for --plot-engine none).",
            ),
        ])
        for opt in reversed(opts):
            fn = opt(fn)
        return fn

    return decorator


def parse_log_overrides(spec: str | None) -> dict[str, str]:
    """Parse `--log-level` value into per-logger overrides.

    Examples:
        None              -> {}
        "DEBUG"           -> {"ema": "DEBUG"}
        "ema.cm=WARNING"  -> {"ema.cm": "WARNING"}
        "DEBUG,ema.cm=WARNING" -> {"ema": "DEBUG", "ema.cm": "WARNING"}
    """
    if not spec:
        return {}
    overrides: dict[str, str] = {}
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            name, lvl = part.split("=", 1)
            overrides[name.strip()] = lvl.strip().upper()
        else:
            overrides["ema"] = part.upper()
    return overrides


def apply_yaml_to_kwargs(
    ctx: click.Context,
    kwargs: dict[str, Any],
    *,
    warn_unknown: bool = True,
) -> dict[str, Any]:
    """Merge a YAML config file into *kwargs* for switch subcommands.

    This is the ``peakatail switch {diff,length,match}`` equivalent of the
    YAML-then-CLI merge that ``peakatail run`` performs.  It is intentionally
    separate from :func:`ema.cli.yaml_loader.load_run_yaml` because that
    function enforces a ``datasets:`` key which does not apply to the
    switch subcommands.

    Merge rules (mirrors the ``peakatail run`` contract):
    * If ``kwargs["config"]`` is ``None``, this is a no-op.
    * YAML keys are mapped to kwarg names via the :class:`RunConfig` schema
      (``yaml_key`` → field name).  Unknown YAML keys emit a WARNING.
    * Only kwargs the user did **not** supply on the CLI are overwritten.
      Detection uses ``ctx.get_parameter_source(name) == COMMANDLINE`` so
      an explicit ``--fdr 0.05`` is not silently eclipsed by a YAML
      ``fdr: 0.10``.

    Args:
        ctx: The active Click context (used for ParameterSource lookup).
        kwargs: Click-parsed kwargs dict (mutated in place and returned).
        warn_unknown: Emit a WARNING for YAML keys not in the schema.

    Returns:
        The same ``kwargs`` dict, mutated in place.
    """
    config_path = kwargs.get("config")
    if not config_path:
        return kwargs

    import yaml
    try:
        from click.core import ParameterSource
        _ps_commandline = ParameterSource.COMMANDLINE
    except ImportError:  # pragma: no cover
        _ps_commandline = None  # type: ignore[assignment]

    path = Path(str(config_path))
    try:
        with open(path) as _f:
            raw = yaml.safe_load(_f) or {}
    except OSError as exc:
        raise click.BadParameter(
            f"Cannot open YAML config {path}: {exc}", param_hint="--config"
        )
    if not isinstance(raw, dict):
        raise click.BadParameter(
            f"YAML config {path}: top-level must be a mapping", param_hint="--config"
        )

    # Build the yaml_key → field_name mapping from the schema.
    from ema.cli.config_schema import RunConfig, yaml_key_to_field_name
    yaml_to_field = yaml_key_to_field_name(RunConfig)

    # Identify which params the user explicitly typed on the CLI.
    explicitly_set: set[str] = set()
    if _ps_commandline is not None:
        for name in ctx.params:
            try:
                if ctx.get_parameter_source(name) == _ps_commandline:
                    explicitly_set.add(name)
            except Exception:
                pass

    for yaml_key, value in raw.items():
        field_name = yaml_to_field.get(yaml_key)
        if field_name is None:
            if warn_unknown:
                log.warning(
                    "YAML key %r is not recognised by the switch subcommand schema"
                    " — ignoring.",
                    yaml_key,
                )
            continue
        # Only apply if the kwarg exists in this subcommand's Click params
        # AND was not explicitly set by the user on the CLI.
        if field_name not in kwargs:
            continue
        if field_name in explicitly_set:
            continue
        kwargs[field_name] = value

    return kwargs


def parse_plot_engines(spec: str, no_plots: bool = False) -> list[str]:
    """Resolve --plot-engine string into a list of engine names. Empty list = disabled."""
    if no_plots or spec == "none":
        return []
    if spec == "both":
        return ["matplotlib", "plotly"]
    return [e.strip() for e in spec.split(",") if e.strip()]


_DEFAULT_PARENT_DIR = "peakatail_runs"


def _no_timestamp() -> bool:
    """True if PEAKATAIL_NO_TIMESTAMP disables the ``_<ts>`` output suffix.

    Lets a workflow engine (Nextflow/Snakemake) get deterministic output dirs
    so a downstream step can hand fixed paths to the next.
    """
    import os
    return os.environ.get("PEAKATAIL_NO_TIMESTAMP", "").strip().lower() in {"1", "true", "yes"}


def detect_run_dir(input_paths: list[str | Path] | tuple) -> Path | None:
    """Return the originating ``peakatail_runs/<run>/`` directory for inputs, if any.

    A subcommand like ``peakatail switch diff -i RUN/per_dataset/sampleA/clusters.h5ad``
    should write its outputs INSIDE ``RUN/`` so the run dir stays self-contained.
    This helper walks each input path upward looking for a directory whose
    parent is named ``peakatail_runs``.  All inputs must agree on the same run
    dir; otherwise ``None`` is returned and the caller falls back to the global
    ``peakatail_runs/<name>_<ts>/`` layout.

    Args:
        input_paths: Iterable of file paths the subcommand consumes (typically
            ``-i h5ad`` arguments).

    Returns:
        The shared run directory as a :class:`Path`, or ``None`` if the inputs
        don't all live in a recognisable run dir.
    """
    if not input_paths:
        return None
    run_dirs: set[Path] = set()
    for p in input_paths:
        path = Path(p).resolve()
        for parent in path.parents:
            if parent.parent.name == _DEFAULT_PARENT_DIR:
                run_dirs.add(parent)
                break
        else:
            return None  # this input is not under peakatail_runs/<run>/
    if len(run_dirs) != 1:
        return None  # mixed sources — can't pick a single home
    return next(iter(run_dirs))


def resolve_subcommand_output_dir(
    base: str | Path,
    *,
    user_explicit: bool,
    source_paths: list[str | Path] | tuple = (),
    subdir: str | None = None,
) -> Path:
    """Resolve a switch/downstream subcommand's output dir.

    When the user did NOT pass ``-o`` and all ``source_paths`` come from the
    same ``peakatail_runs/<run>/`` directory, the output is routed INSIDE that
    run dir as ``<run>/<subdir or base>_<ts>/``.  Otherwise it falls back to
    :func:`resolve_output_dir` (a fresh ``peakatail_runs/<base>_<ts>/``).

    Args:
        base: User-supplied output name (CLI ``--output`` or its default).
        user_explicit: ``True`` iff the user explicitly passed ``--output``.
            Detect via ``ctx.get_parameter_source("output").name == "COMMANDLINE"``.
        source_paths: Input files the subcommand consumes (``-i h5ad`` etc.).
        subdir: Optional override for the leaf directory name.  Defaults to
            ``base`` so e.g. ``base="switch_out"`` yields ``<run>/switch_out_<ts>/``.
    """
    if not user_explicit:
        run_dir = detect_run_dir(list(source_paths))
        if run_dir is not None:
            leaf_name = subdir if subdir else Path(base).name
            if _no_timestamp():
                return run_dir / leaf_name
            from datetime import datetime
            ts = datetime.now().strftime("%Y-%m-%d_%H%M%S")
            return run_dir / f"{leaf_name}_{ts}"
    return resolve_output_dir(base)


def resolve_output_dir(base: str | Path, parent: str | Path | None = None) -> Path:
    """Append a timestamp to the user-supplied output dir and nest under a parent.

    Default layout: ``peakatail_runs/<base>_<timestamp>/`` so all runs are
    grouped under one project-level dir instead of scattered at the repo root.

    If ``base`` is an absolute path or already lives inside an explicit parent
    (i.e. contains a path separator), the parent prefix is NOT applied — the
    user knows exactly where they want the run to go.

    Args:
        base: Output dir name from the user (CLI ``--output`` or YAML
            ``output_dir``). Plain name (e.g. ``"emaout"``) is nested under
            ``peakatail_runs/``; a path with separators is used as-is.
        parent: Override the default ``peakatail_runs`` parent dir. Use ``""``
            to disable nesting entirely.

    Timestamping can be disabled (for deterministic paths, e.g. under a
    Nextflow/Snakemake DAG that hands fixed paths between steps) by setting the
    environment variable ``PEAKATAIL_NO_TIMESTAMP=1``.  When disabled the leaf
    is exactly ``base`` with no ``_<ts>`` suffix, so the output dir equals the
    user-supplied path verbatim.

    Returns:
        A fresh ``Path``. Does NOT create the directory (caller decides when).
    """
    from datetime import datetime
    base = Path(base)
    if _no_timestamp():
        leaf = Path(base.name)
    else:
        ts = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        leaf = Path(f"{base.name}_{ts}")
    # Absolute paths or paths with explicit directories are passed through.
    if base.is_absolute() or len(base.parts) > 1:
        return base.parent / leaf
    parent_dir = _DEFAULT_PARENT_DIR if parent is None else str(parent)
    if not parent_dir:
        return leaf
    return Path(parent_dir) / leaf
