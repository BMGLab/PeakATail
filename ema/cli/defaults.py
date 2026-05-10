"""Single source of truth for every CLI flag's default value.

Both Click commands and the interactive wizard read from this dict so
the two surfaces never drift. To add a new flag, add it here and reference
DEFAULTS[key] in the @click.option(... default=...) and the wizard prompt.
"""
from __future__ import annotations

# Each entry: key = canonical flag name (kebab-case without leading --),
# value = default value used by both Click and wizard.
DEFAULTS: dict[str, object] = {
    # ─── input / output ───
    "config": None,
    "output": "emaout",
    "bam-dir": None,
    "bam-files": None,
    "gtf": None,
    "atlas": None,
    "atlas-distance": 50,

    # ─── read processing ───
    "seq-len": None,
    "cb-len": None,
    "barcode-tag": None,

    # ─── concurrency ───
    "threads": None,
    "bam-threads": 4,
    "pipeline": False,
    "batch-size": 10000,
    "tiles": False,
    "tile-size": None,  # auto via ResourceManager.get_tile_size
    "tile-overlap": 10000,

    # ─── peak calling ───
    "peak-strategy": "original",
    "lambda-window": 5000,
    "lambda-method": "median",
    "lambda-fold-change": 2.0,
    "max-pas": 5,
    "smoothing-window": 50,
    "min-prominence": 5.0,
    "dynamic-threshold": False,
    "floor-threshold": 3,
    "pas-gap": 100,

    # ─── filters ───
    "ip-filter": False,
    "genome-fasta": None,
    "annot-filter": False,
    "ip-a-stretch": 6,
    "min-pas-per-cell": 50,
    "min-read": 1500,
    "min-cells": 3,

    # ─── annotation ───
    "max-gene-distance": 5000,
    "utr-multiplier": 2.0,
    "include-extended": False,

    # ─── clustering ───
    "cluster-method": "leiden_tfidf",
    "resolution": 1.0,
    "n-pcs": 40,
    "external-clusters": None,
    "random-seed": 42,

    # ─── cross-dataset matching ───
    "match-method": "marker_overlap",
    "n-top-markers": 50,

    # ─── switch diff ───
    "fdr": 0.05,
    "marker-method": "wilcoxon",
    "marker-top-n": 200,
    "per-worker-mb": 300,

    # ─── switch length ───
    # Vocabulary must match ema/quantification/strategies/base.py:
    #   AggregationMode = Literal["per_isoform", "per_gene"]
    #   IsoformCollapseMode = Literal["none", "mean", "majority"]
    # Earlier values "gene" / "weighted" silently exercised the WRONG branch
    # because strategy code only matched the per_* tokens.
    "isoform-agg": "per_gene",
    "isoform-collapse": "none",

    # ─── validation ───
    "benchmark": False,
    "validate-db": None,

    # ─── logging / UX ───
    "log-level": None,
    "no-log-file": False,
    "no-progress": False,
    "verbose": 0,    # -v count
    "quiet": False,  # -q
}


def get_default(key: str):
    """Return DEFAULTS[key], raising a clear KeyError on typos."""
    if key not in DEFAULTS:
        raise KeyError(
            f"No default registered for {key!r}. Add it to ema/cli/defaults.py "
            "so Click flags and the wizard stay in sync."
        )
    return DEFAULTS[key]
