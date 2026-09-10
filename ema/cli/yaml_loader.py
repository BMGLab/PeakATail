"""YAML config loading + schema validation for `peakatail run`.

Schema is a strict superset of today's main.py YAML reader:
- `datasets` (required, list of {id, merge_strategy, bams})
- `gtf`, `output_dir`, `seqlen`, `cb_len`, `barcode_tag`
- `min_read`, `min_cells`, `min_pas_per_cell`, `pas_gap`
- `atlas`, `atlas_distance`
- `cluster_match_method`, `n_top_markers`

Dead keys (no longer used by peakatail run) emit a one-shot warning telling
the user where to find the equivalent: `peakatail switch length` / `peakatail switch diff`.
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


def _build_live_keys() -> set[str]:
    """Derive the live YAML key set from the RunConfig schema.

    Adds two keys the schema does not model:
      * ``datasets`` -- dict-shaped, validated structurally below.
      * any future schema gaps caught by tests/test_config_schema.py.
    """
    from ema.cli.config_schema import RunConfig, yaml_keys_from_schema
    return {"datasets"} | yaml_keys_from_schema(RunConfig)


# Single source of truth -- never hand-edit this set; modify the schema.
_LIVE_KEYS = _build_live_keys()

# Dead keys: were used by old `ema` but moved to `peakatail switch ...`
_DEAD_KEYS = {
    "pdui_method": "peakatail switch length",
    "pdui_isoform_agg": "peakatail switch length",
    "pdui_isoform_collapse": "peakatail switch length",
    "diff_method": "peakatail switch diff",
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
    _seen_ids: set[str] = set()
    for i, ds in enumerate(cfg["datasets"]):
        if not isinstance(ds, dict):
            raise RunYamlError(f"{path}: datasets[{i}] must be a mapping")
        if "id" not in ds:
            raise RunYamlError(f"{path}: datasets[{i}] missing `id`")
        # Barcode-namespace safety: cells are keyed as "<id>_<barcode>" and the
        # indexer splits on the FIRST '_'. An id containing '_' would swallow
        # part of the barcode → every cell collapses to one; a duplicate id
        # would merge two libraries' cells (same barcode = same cell) — a
        # cross-library collision. Reject both structurally.
        _id = ds["id"]
        if "_" in str(_id):
            raise RunYamlError(
                f"{path}: datasets[{i}].id={_id!r} contains '_'. Dataset ids "
                f"namespace cell barcodes and must not contain underscores "
                f"(use '-', e.g. 'GSM123-StageI')."
            )
        if _id in _seen_ids:
            raise RunYamlError(
                f"{path}: duplicate dataset id {_id!r}. Each id is one library; "
                f"reusing an id would merge distinct libraries' cells (same "
                f"barcode would collide across libraries)."
            )
        _seen_ids.add(_id)
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
                "YAML key %r is no longer used by `peakatail run` — moved to %r. Ignored.",
                k, _DEAD_KEYS[k],
            )
        elif k not in _LIVE_KEYS:
            log.warning("YAML key %r is unknown — ignoring.", k)

    return cfg
