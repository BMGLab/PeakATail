"""Round-trip the canonical cluster map back into each dataset's h5ad obs (B6).

Bug B6: ``canonical_cluster_map.tsv`` was written but never read back into
``clusters.h5ad``. Cross-sample comparisons then silently compared per-sample
``leiden`` labels — which are NOT comparable across datasets (cluster 0 in
sample A is unrelated to cluster 0 in sample B). Writing the shared
``canonical_cluster`` into ``obs`` makes cross-sample joins correct.
"""
from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)


def build_canonical_lookup(match_df) -> dict[tuple[str, str], int]:
    """Build ``(dataset_id, original_cluster) -> canonical_cluster`` from match_df.

    ``original_cluster`` keys are coerced to ``str`` so they join against a
    stringified ``obs['leiden']`` regardless of dtype (categorical/int/str).
    """
    required = {"dataset_id", "original_cluster", "canonical_cluster"}
    missing = required - set(match_df.columns)
    if missing:
        raise KeyError(f"match_df missing columns for canonical round-trip: {missing}")
    lookup: dict[tuple[str, str], int] = {}
    for row in match_df.itertuples(index=False):
        lookup[(str(row.dataset_id), str(row.original_cluster))] = int(row.canonical_cluster)
    return lookup


def annotate_h5ad_with_canonical(
    h5ad_path: str | Path,
    dataset_id: str,
    lookup: dict[tuple[str, str], int],
    *,
    leiden_key: str = "leiden",
    out_key: str = "canonical_cluster",
) -> int:
    """Add ``obs[out_key]`` to one h5ad by mapping its per-sample leiden labels.

    Returns the number of cells that received a (non-NA) canonical id. Cells
    whose leiden label has no entry in ``lookup`` get ``-1`` (an explicit
    "unmatched" sentinel — never silently NA). Writes the h5ad in place.
    """
    import anndata as ad  # heavy import kept local

    adata = ad.read_h5ad(h5ad_path)
    if leiden_key not in adata.obs.columns:
        log.warning(
            "canonical round-trip: %s has no obs[%r]; skipping.",
            h5ad_path, leiden_key,
        )
        return 0
    leiden = adata.obs[leiden_key].astype(str)
    canon = leiden.map(lambda lab: lookup.get((str(dataset_id), lab), -1)).astype(int)
    adata.obs[out_key] = canon.values
    n_assigned = int((canon != -1).sum())
    adata.write(h5ad_path)
    log.info(
        "canonical round-trip: %s -> obs[%r] (%d/%d cells matched)",
        h5ad_path, out_key, n_assigned, adata.n_obs,
    )
    return n_assigned


def write_canonical_clusters(
    match_df,
    h5ad_by_dataset: list[tuple[str | Path, str]],
    *,
    leiden_key: str = "leiden",
    out_key: str = "canonical_cluster",
) -> dict[str, int]:
    """Write ``canonical_cluster`` into every dataset's obs from ``match_df``.

    ``h5ad_by_dataset`` is a list of ``(h5ad_path, dataset_id)``. Returns a
    ``{dataset_id: n_cells_matched}`` summary. Never raises on a single missing
    h5ad — logs and continues.
    """
    lookup = build_canonical_lookup(match_df)
    summary: dict[str, int] = {}
    for h5ad_path, dataset_id in h5ad_by_dataset:
        p = Path(h5ad_path)
        if not p.exists():
            log.warning("canonical round-trip: %s missing; skipping %r.", p, dataset_id)
            summary[dataset_id] = 0
            continue
        summary[dataset_id] = annotate_h5ad_with_canonical(
            p, dataset_id, lookup, leiden_key=leiden_key, out_key=out_key,
        )
    return summary
