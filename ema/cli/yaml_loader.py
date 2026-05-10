"""YAML config loading + schema validation for `ema run`.

Schema is a strict superset of today's main.py YAML reader:
- `datasets` (required, list of {id, merge_strategy, bams})
- `gtf`, `output_dir`, `seqlen`, `cb_len`, `barcode_tag`
- `min_read`, `min_cells`, `min_pas_per_cell`, `pas_gap`
- `atlas`, `atlas_distance`
- `cluster_match_method`, `n_top_markers`

Dead keys (no longer used by ema run) emit a one-shot warning telling
the user where to find the equivalent: `ema switch length` / `ema switch diff`.
Unknown keys also warn (don't fail).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import yaml

log = logging.getLogger(__name__)


class RunYamlError(ValueError):
    """Raised when a YAML config is malformed beyond warning."""


_VALID_MERGE_STRATEGIES = {"before", "after", "none"}

# Live keys: every YAML key the wizard or `ema run` may emit.  Any key
# missing here triggers a false "unknown YAML key" warning even though
# the value is honoured downstream.  Part B will derive this set from
# the centralised schema (RunConfig); pre-Part-B we hand-extend it to
# cover every accepted key:
#   - every key the wizard writes (threads, output_dir, cluster_method,
#     resolution, n_pcs, random_seed, external_clusters, peak_strategy
#     and its hyperparameters, ...)
#   - every key main.py reads (`peak_strategy`, `bam_threads`, `tiles`,
#     ...)
_LIVE_KEYS = {
    # Inputs / outputs
    "datasets", "gtf", "output_dir", "atlas", "atlas_distance",
    # Read processing
    "seqlen", "cb_len", "barcode_tag",
    # Concurrency / runtime
    "threads", "bam_threads", "pipeline", "batch_size",
    "tiles", "tile_size", "tile_overlap",
    # Peak calling
    "peak_strategy", "lambda_window", "lambda_method",
    "lambda_fold_change", "max_pas", "smoothing_window",
    "min_prominence", "dynamic_threshold", "floor_threshold",
    "pas_gap",
    # Filters
    "ip_filter", "genome_fasta", "annot_filter", "ip_a_stretch",
    "min_pas_per_cell", "min_read", "min_cells", "min_genes",
    # Annotation
    "max_gene_distance", "utr_multiplier", "include_extended",
    # Clustering
    "cluster_method", "resolution", "n_pcs",
    "external_clusters", "random_seed",
    # Cross-dataset matching
    "match_method", "cluster_match_method", "n_top_markers",
}

# Dead keys: were used by old `ema` but moved to `ema switch ...`
_DEAD_KEYS = {
    "pdui_method": "ema switch length",
    "pdui_isoform_agg": "ema switch length",
    "pdui_isoform_collapse": "ema switch length",
    "diff_method": "ema switch diff",
}


def load_run_yaml(path: str | Path) -> dict[str, Any]:
    """Load + validate a run-mode YAML config. Returns the loaded dict.

    Raises RunYamlError on missing required keys or invalid values.
    Logs WARNINGS for dead keys (preserves them in returned dict so any
    backward-compat caller can still read them) and unknown keys.
    """
    path = Path(path)
    with open(path) as f:
        cfg = yaml.safe_load(f) or {}
    if not isinstance(cfg, dict):
        raise RunYamlError(f"{path}: top-level YAML must be a mapping")

    # 1. required + valid datasets
    if "datasets" not in cfg or not cfg["datasets"]:
        raise RunYamlError(f"{path}: missing required `datasets:` list")
    for i, ds in enumerate(cfg["datasets"]):
        if not isinstance(ds, dict):
            raise RunYamlError(f"{path}: datasets[{i}] must be a mapping")
        if "id" not in ds:
            raise RunYamlError(f"{path}: datasets[{i}] missing `id`")
        ms = ds.get("merge_strategy", "none")
        if ms not in _VALID_MERGE_STRATEGIES:
            raise RunYamlError(
                f"{path}: datasets[{i}].merge_strategy={ms!r} — "
                f"must be one of {sorted(_VALID_MERGE_STRATEGIES)}"
            )
        if "bams" not in ds or not ds["bams"]:
            raise RunYamlError(f"{path}: datasets[{i}] missing `bams:` list")

    # 2. dead-key warnings
    for k in cfg:
        if k in _DEAD_KEYS:
            log.warning(
                "YAML key %r is no longer used by `ema run` — moved to %r. Ignored.",
                k, _DEAD_KEYS[k],
            )
        elif k not in _LIVE_KEYS:
            log.warning("YAML key %r is unknown — ignoring.", k)

    return cfg
