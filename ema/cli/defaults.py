"""Single source of truth for every CLI flag's default value.

All defaults are now derived from :class:`ema.cli.config_schema.RunConfig`,
including switch-subcommand-only flags (``peakatail switch diff/length/match``).

Adding a flag:
    * Add a dataclass field on RunConfig in ``ema/cli/config_schema.py``.
      ``DEFAULTS`` picks it up automatically via :func:`cli_defaults`.
    * For switch-only fields, set ``skip_legacy_bridge=True`` and
      ``applies_to=frozenset({"switch_diff"})`` (or whichever subcommand).
"""
from __future__ import annotations

from typing import Any

from ema.cli.config_schema import RunConfig, cli_defaults


def _build_defaults() -> dict[str, Any]:
    """Build the unified defaults dict from the RunConfig schema."""
    return cli_defaults(RunConfig)


DEFAULTS: dict[str, object] = _build_defaults()


def get_default(key: str):
    """Return DEFAULTS[key], raising a clear KeyError on typos."""
    if key not in DEFAULTS:
        raise KeyError(
            f"No default registered for {key!r}. Add a dataclass field on "
            "RunConfig in ema/cli/config_schema.py and DEFAULTS picks it up."
        )
    return DEFAULTS[key]
