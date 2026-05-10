"""Decorator factory that applies the shared options to every subcommand.

Usage:
    @click.command()
    @common_options()
    def my_cmd(threads, output, verbose, quiet, log_level, no_log_file, no_progress, config, **kwargs):
        ...
"""
from __future__ import annotations

import functools
from pathlib import Path
from typing import Callable

import click

from ema.cli.defaults import DEFAULTS


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
                    type=click.Path(file_okay=False, resolve_path=True),
                    default=output_default if output_default is not None else DEFAULTS["output"],
                    help="Output directory (timestamp suffix added automatically).",
                )
            )
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


_DEFAULT_PARENT_DIR = "peakatail_runs"


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

    Returns:
        A fresh ``Path``. Does NOT create the directory (caller decides when).
    """
    from datetime import datetime
    ts = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    base = Path(base)
    leaf = Path(f"{base.name}_{ts}")
    # Absolute paths or paths with explicit directories are passed through.
    if base.is_absolute() or len(base.parts) > 1:
        return base.parent / leaf
    parent_dir = _DEFAULT_PARENT_DIR if parent is None else str(parent)
    if not parent_dir:
        return leaf
    return Path(parent_dir) / leaf
