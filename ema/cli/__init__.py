"""Click root group for PeakATail.

The single `ema` entry point. Subcommands:
    ema run                  — full pipeline
    ema switch {diff,length,match} — per-cluster analyses
    ema merge                — merge BAMs
    ema parse-gtf            — pre-warm GTF cache
    ema wizard               — interactive setup
    (bare `ema`)             — same as `ema wizard`
"""
from __future__ import annotations

import sys
from importlib import metadata as md

import click


def _get_version() -> str:
    try:
        return md.version("peakatail")
    except md.PackageNotFoundError:
        return "0.0.0-dev"


@click.group(
    name="ema",
    invoke_without_command=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)
@click.version_option(version=_get_version(), package_name="peakatail")
@click.pass_context
def main(ctx: click.Context) -> None:
    """PeakATail — single-cell poly(A) site detection and APA analysis."""
    if ctx.invoked_subcommand is None:
        # Bare `ema` → wizard
        from ema.cli import wizard  # local import: lets test monkeypatch
        rc = wizard.run()
        sys.exit(rc)


# Subcommand registration happens after import. We import the leaves at
# the bottom so they can register themselves on `main`.
from ema.cli import run as _run  # noqa: E402,F401
from ema.cli import switch as _switch  # noqa: E402,F401
from ema.cli import merge as _merge  # noqa: E402,F401
from ema.cli import parse_gtf as _parse_gtf  # noqa: E402,F401
from ema.cli import wizard as _wizard_mod  # noqa: E402,F401

main.add_command(_run.run)
main.add_command(_switch.switch)
main.add_command(_merge.merge)
main.add_command(_parse_gtf.parse_gtf)
main.add_command(_wizard_mod.wizard)


# ---------------------------------------------------------------------------
# Backward-compatibility shim for ema/config.py
# ---------------------------------------------------------------------------
# ema/config.py does ``args = cli()`` at module load time.  The 80-line
# argparse parser this used to invoke was the original ``ema`` CLI; it has
# been superseded by the Click subcommands above.  The remaining
# responsibility of ``cli()`` is therefore minimal:
#
#   1. Return an object with the same attributes the legacy pipeline body
#      reads (e.g. ``args.strategy``, ``args.bam_threads``, ...) so module
#      load doesn't AttributeError.
#   2. Honour ``--config <yaml>`` if provided on the actual command line so
#      a power user can still bootstrap the legacy globals.
#
# Defaults are derived from RunConfig -- the schema is the single source.
# Subsequent ``RunConfig.apply_to_legacy_globals()`` calls (from
# ema.main.run) will override these starter values with the actual
# user-supplied configuration.
# ---------------------------------------------------------------------------
def cli():
    """Build the legacy ``args`` namespace from the RunConfig schema.

    Returns:
        argparse.Namespace -- carries every attribute the legacy globals
        and pipeline body read.

    Notes:
        * Unknown CLI arguments (e.g. pytest's ``-x``) are ignored via
          ``parse_known_args``; only ``--config`` is honoured.
        * The Click CLI does not call this -- it is only loaded by the
          import-time chain ``ema.config -> ema.cli.cli()``.  The clean
          fix (lazy config) is tracked separately; this shim keeps the
          legacy contract while we transition.
    """
    import argparse

    from ema.cli.config_schema import RunConfig, field_specs

    ns = argparse.Namespace()

    # Seed every attribute the schema knows about.
    specs = field_specs(RunConfig)
    from dataclasses import fields as _fields
    for f in _fields(RunConfig):
        spec = specs[f.name]
        # Map field -> legacy args attribute name.
        attr = f.name
        if spec.legacy_args_attr:
            attr = spec.legacy_args_attr
        setattr(ns, attr, f.default)

    # Add the legacy attribute names that aren't 1:1 schema fields.
    ns.datasets = []
    ns.gtf_dir = None
    ns.bam_dir = None
    ns.atlas = None
    ns.atlas_distance = 50
    ns.seqlen = None
    ns.cb_len = None
    ns.barcode_tag = None
    ns.min_pas_per_cell = 50

    # Honour the small handful of legacy argparse flags some module-import
    # smoke tests still rely on (test_region_fetch.py patches sys.argv with
    # the original flag names so peak_calling can seed variable_config at
    # import time).  Everything else flows through the Click subcommands.
    parser = argparse.ArgumentParser(prog="ema", add_help=False)
    parser.add_argument("--config", dest="config", type=str, default=None)
    parser.add_argument("--bamDir", dest="_bam_dir", type=str, default=None)
    parser.add_argument("--sequenceLen", dest="_seqlen", type=int, default=None)
    parser.add_argument("--CellBarcodeLen", dest="_cb_len", type=int, default=None)
    parser.add_argument("--BarcodeTag", dest="_barcode_tag", type=str, default=None)
    parsed, _ = parser.parse_known_args()
    if parsed._bam_dir is not None:
        ns.bam_dir = parsed._bam_dir
        ns.datasets = [{"id": "default", "merge_strategy": "none",
                        "bams": [parsed._bam_dir]}]
    if parsed._seqlen is not None:
        ns.seqlen = parsed._seqlen
    if parsed._cb_len is not None:
        ns.cb_len = parsed._cb_len
    if parsed._barcode_tag is not None:
        ns.barcode_tag = parsed._barcode_tag
    if parsed.config:
        import yaml
        with open(parsed.config) as f:
            cfg = yaml.safe_load(f) or {}
        # Bridge the YAML through RunConfig so legacy aliases / defaults
        # are handled in one place.
        rc = RunConfig.from_yaml_dict(cfg)
        for f in _fields(RunConfig):
            spec = specs[f.name]
            attr = spec.legacy_args_attr or f.name
            value = getattr(rc, f.name)
            setattr(ns, attr, value)
        # YAML-only attributes.
        if "datasets" in cfg:
            ns.datasets = cfg["datasets"]
        if "gtf" in cfg:
            ns.gtf_dir = cfg["gtf"]
        if "seqlen" in cfg:
            ns.seqlen = cfg["seqlen"]
        if "cb_len" in cfg:
            ns.cb_len = cfg["cb_len"]
        if "barcode_tag" in cfg:
            ns.barcode_tag = cfg["barcode_tag"]
        if "atlas" in cfg:
            ns.atlas = cfg["atlas"]
        if "atlas_distance" in cfg:
            ns.atlas_distance = cfg["atlas_distance"]
        if "min_pas_per_cell" in cfg:
            ns.min_pas_per_cell = cfg["min_pas_per_cell"]

    return ns
