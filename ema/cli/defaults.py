"""Single source of truth for every CLI flag's default value.

The ``ema run`` defaults are derived from :class:`ema.cli.config_schema.RunConfig`
so the dataclass schema is the canonical source.  Sub-command-only flags
(``ema switch diff/length/match``, ``ema parse-gtf``) are listed inline
here because they are not part of the run-pipeline schema.

Adding a flag:
    * If it belongs to ``ema run``: add a dataclass field on RunConfig in
      ``ema/cli/config_schema.py`` and DEFAULTS picks it up automatically.
    * If it is sub-command-only: add it to ``_SUBCOMMAND_DEFAULTS`` below.
"""
from __future__ import annotations

from typing import Any

from ema.cli.config_schema import RunConfig, cli_defaults


# Sub-command-only flag defaults (NOT exposed on `ema run`).
# Keep this dict TINY -- the run pipeline schema is RunConfig.
_SUBCOMMAND_DEFAULTS: dict[str, Any] = {
    # `ema switch diff`
    "fdr": 0.05,
    "marker-method": "wilcoxon",
    "marker-top-n": 200,
    "per-worker-mb": 300,
    # `ema switch length` -- vocabulary MUST match strategy code
    # (per_gene/per_isoform), see commit "fix(switch length): align
    # --isoform-agg vocabulary with strategy code".
    "isoform-agg": "per_gene",
    "isoform-collapse": "none",
}


def _build_defaults() -> dict[str, Any]:
    """Combine schema-derived run defaults with subcommand-only defaults.

    Run defaults (``cli_defaults(RunConfig)``) take precedence over
    subcommand defaults if they overlap -- they are the authoritative
    source.
    """
    out: dict[str, Any] = dict(_SUBCOMMAND_DEFAULTS)
    out.update(cli_defaults(RunConfig))
    return out


DEFAULTS: dict[str, object] = _build_defaults()


def get_default(key: str):
    """Return DEFAULTS[key], raising a clear KeyError on typos."""
    if key not in DEFAULTS:
        raise KeyError(
            f"No default registered for {key!r}. Add it to "
            "ema/cli/config_schema.py::RunConfig (for ema run) or "
            "ema/cli/defaults.py::_SUBCOMMAND_DEFAULTS (for ema switch *)."
        )
    return DEFAULTS[key]
