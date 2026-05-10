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
# Backward-compatibility shim for ema/config.py (Phase 9 will remove this)
# ema/config.py does `from ema.cli import cli; args = cli()` at module level.
# Until Phase 9 rewrites config.py, we provide a cli() that returns defaults.
# ---------------------------------------------------------------------------
def cli():  # noqa: D103
    """Legacy argparse shim — returns a Namespace with all default values.

    Phase 9 (legacy YAML loader extraction) will remove this shim once
    ema/config.py no longer calls argparse at import time.
    """
    import argparse
    import sys

    parser = argparse.ArgumentParser(prog="ema", add_help=False)
    parser.add_argument("--config", dest="config", type=str, default=None)
    parser.add_argument("--bamDir", dest="bam_dir", type=str, default=None)
    parser.add_argument("--sequenceLen", dest="seqlen", type=int, default=None)
    parser.add_argument("--CellBarcodeLen", dest="cb_len", type=int, default=None)
    parser.add_argument("--BarcodeTag", type=str, dest="barcode_tag", default=None)
    parser.add_argument("--gtfDir", dest="gtf_dir", type=str, default=None)
    parser.add_argument("--cell_combinations", type=str, default=None)
    parser.add_argument("--bamFiles", dest="bam_files", type=str, default=None)
    parser.add_argument("--threads", dest="threads", type=int, default=None)
    parser.add_argument("--bam-threads", dest="bam_threads", type=int, default=4)
    parser.add_argument("--pipeline", action="store_true", default=False)
    parser.add_argument("--batch-size", dest="batch_size", type=int, default=10000)
    parser.add_argument("--tiles", action="store_true", default=False)
    parser.add_argument("--tile-size", dest="tile_size", type=int, default=25_000_000)
    parser.add_argument("--tile-overlap", dest="tile_overlap", type=int, default=10_000)
    parser.add_argument("--strategy", type=str, default="original")
    parser.add_argument("--lambda-window", dest="lambda_window", type=int, default=5000)
    parser.add_argument("--lambda-method", dest="lambda_method", type=str, default="median")
    parser.add_argument("--max-pas", dest="max_pas", type=int, default=5)
    parser.add_argument("--smoothing-window", dest="smoothing_window", type=int, default=50)
    parser.add_argument("--min-prominence", dest="min_prominence", type=float, default=5.0)
    parser.add_argument("--benchmark", action="store_true", default=False)
    parser.add_argument("--validate-db", dest="validate_db", type=str, default=None)
    parser.add_argument("--internal-priming-filter", dest="ip_filter", action="store_true", default=False)
    parser.add_argument("--genome-fasta", dest="genome_fasta", type=str, default=None)
    parser.add_argument("--annotation-filter", dest="annot_filter", action="store_true", default=False)
    parser.add_argument("--ip-a-stretch", dest="ip_a_stretch", type=int, default=6)
    parser.add_argument("--min-pas-per-cell", dest="min_pas_per_cell", type=int, default=50)
    parser.add_argument("--max-gene-distance", dest="max_gene_distance", type=int, default=5000)
    parser.add_argument("--utr-multiplier", dest="utr_multiplier", type=float, default=2.0)
    parser.add_argument("--include-extended", dest="include_extended", action="store_true", default=False)
    parser.add_argument("--dynamic-threshold", dest="dynamic_threshold", action="store_true", default=False)
    parser.add_argument("--floor-threshold", dest="floor_threshold", type=int, default=3)
    parser.add_argument("--lambda-fold-change", dest="lambda_fold_change", type=float, default=2.0)
    parser.add_argument("--atlas", type=str, default=None)
    parser.add_argument("--atlas-distance", dest="atlas_distance", type=int, default=50)
    parser.add_argument("--clustering-method", dest="clustering_method", type=str, default="leiden_tfidf")
    parser.add_argument("--resolution", type=float, default=1.0)
    parser.add_argument("--n-pcs", dest="n_pcs", type=int, default=40)
    parser.add_argument("--external-clusters", dest="external_clusters", type=str, default=None)
    parser.add_argument("--random-seed", dest="random_seed", type=int, default=42)
    parser.add_argument("--pdui-method", dest="pdui_method", type=str, default="classic")
    parser.add_argument("--pdui-isoform-agg", dest="pdui_isoform_agg", type=str, default="per_gene")
    parser.add_argument("--pdui-isoform-collapse", dest="pdui_isoform_collapse", type=str, default="none")
    parser.add_argument("--diff-method", dest="diff_method", type=str, default="fisher")
    parser.add_argument("--cluster-match-method", dest="cluster_match_method", type=str, default="marker_overlap")
    parser.add_argument("--n-top-markers", dest="n_top_markers", type=int, default=50)

    # parse_known_args ignores unknown args (pytest args) without exiting
    args, _ = parser.parse_known_args()

    # Load YAML config if provided
    if args.config:
        import yaml
        with open(args.config) as f:
            cfg = yaml.safe_load(f)
        args.datasets = cfg.get("datasets", [])
        for key in [
            "seqlen", "cb_len", "barcode_tag", "min_read", "min_cells",
            "min_pas_per_cell", "pas_gap", "gtf_dir", "atlas", "atlas_distance",
            "pdui_method", "pdui_isoform_agg", "pdui_isoform_collapse",
            "diff_method", "cluster_match_method", "n_top_markers",
        ]:
            yaml_key = "gtf" if key == "gtf_dir" else key
            if yaml_key in cfg:
                setattr(args, key, cfg[yaml_key])
    elif args.bam_dir:
        args.datasets = [{"id": "default", "merge_strategy": "none", "bams": [args.bam_dir]}]
    else:
        args.datasets = []

    return args
