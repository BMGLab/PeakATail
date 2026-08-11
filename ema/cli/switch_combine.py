"""`ema switch combine` — stitch stage-labelled h5ads for grouped switch analysis.

Combines several clustered ``clusters.h5ad`` files (one per dataset/stage) into
per-group ``<slug>.h5ad`` files whose ``obs[group_key]`` carries the stage
label, ready for ``ema switch {diff,length} --cluster-key <group_key>``.
Optionally splits by an existing cell-type obs column. See
:mod:`ema.switch_test.combine`.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import click

from ema.cli.common import common_options, parse_log_overrides

log = logging.getLogger(__name__)


def _parse_input(spec: str) -> tuple[str, str]:
    """Parse a ``LABEL=PATH`` input spec."""
    if "=" not in spec:
        raise click.UsageError(
            f"--input must be LABEL=PATH (e.g. Normal=run1/clusters.h5ad), got {spec!r}"
        )
    label, path = spec.split("=", 1)
    label, path = label.strip(), path.strip()
    if not label or not path:
        raise click.UsageError(f"--input LABEL and PATH must both be non-empty, got {spec!r}")
    return label, path


@click.command(name="combine")
@common_options(output_default="switch_combined")
@click.option("--input", "-i", "inputs", multiple=True, required=True,
              help="LABEL=PATH for each stage-labelled clusters.h5ad "
                   "(repeatable, >=2 required).")
@click.option("--group-key", "group_key", default="stage", show_default=True,
              help="obs column to stamp with the input LABEL.")
@click.option("--split-key", "split_key", default=None,
              help="Optional obs column (e.g. 'celltype') to split by before combining.")
@click.option("--min-cells", "min_cells", type=int, default=1, show_default=True,
              help="Drop a group value with fewer than this many cells.")
def combine(**kwargs) -> None:
    """Combine stage-labelled h5ads into per-group h5ads for switch analysis."""
    from ema.logging_config import setup_logging, teardown_logging

    out_dir = Path(kwargs["output"] or "switch_combined")
    out_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(
        level=(kwargs["log_level"] or "INFO").split(",")[0],
        output_dir=str(out_dir),
        quiet=kwargs["quiet"],
        verbose_count=kwargs["verbose"],
        per_logger_overrides=parse_log_overrides(kwargs["log_level"]),
    )
    try:
        from ema.switch_test.combine import combine_to_dir

        parsed = [_parse_input(s) for s in kwargs["inputs"]]
        if len(parsed) < 2:
            raise click.UsageError("need >= 2 --input LABEL=PATH entries to contrast.")

        manifest = combine_to_dir(
            parsed, out_dir,
            group_key=kwargs["group_key"],
            split_key=kwargs["split_key"],
            min_cells=kwargs["min_cells"],
        )
        (out_dir / "combine_manifest.json").write_text(json.dumps(manifest, indent=2))
        if not manifest:
            log.warning(
                "no combined groups written — every split spanned < 2 groups with "
                ">= %d cells. Check --split-key / --min-cells.", kwargs["min_cells"],
            )
        for split_val, info in manifest.items():
            log.info("combined %s: %d cells over groups %s",
                     split_val, info["n_cells"], sorted(info["groups"]))
        click.echo(f"combined h5ads written to {out_dir} ({len(manifest)} groups)")
    finally:
        teardown_logging()
