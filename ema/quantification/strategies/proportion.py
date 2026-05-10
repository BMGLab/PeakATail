"""Per-PAS proportion vector PDUI strategy.

For each (gene, cell) pair, computes the proportion of reads falling on each
PAS so that proportions sum to 1.0.

The heavy normalization is done via ``scipy.sparse`` operations — the count
matrix is never fully densified.  Only when extracting gene-level slices (
typically O(n_pas_per_gene) rows) is a small dense array created.

Output schema (long-format):

    gene_id        str
    transcript_id  str   — transcript ID (aggregation=per_isoform) or
                           "_gene_" (aggregation=per_gene)
    pas_id         int
    rank           int   — proximal-to-distal rank (0 = most proximal)
    cell           str
    proportion     float64   — sums to 1.0 per (gene/transcript, cell)

Proportions for cells where the gene total is 0 are NaN.

Parallelized via ``joblib.Parallel`` when ``n_genes > 5000``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.sparse import issparse, csr_matrix

from ema.quantification.strategies.base import (
    AggregationMode,
    IsoformCollapseMode,
    PDUIStrategy,
)
from ema.quantification.strategies import register_pdui_strategy


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def _proportions_for_rows(
    row_indices: list[int],
    count_matrix: pd.DataFrame,
    pseudocount: float,
) -> np.ndarray:
    """Return a (len(row_indices), n_cells) float64 array of proportions.

    Uses sparse normalization: extracts rows as a CSR sub-matrix, sums
    across PAS (axis=0 of transposed), then divides.

    Args:
        row_indices: Positional (iloc) indices into *count_matrix*.
        count_matrix: Shape ``(n_pas, n_cells)``.
        pseudocount: Added to each count before normalization.

    Returns:
        float64 array of shape ``(n_selected_pas, n_cells)``.
        NaN where column total = 0.
    """
    # Extract sub-matrix as dense (small slice per gene)
    sub = count_matrix.iloc[row_indices].values.astype(float)
    if pseudocount != 0.0:
        sub = sub + pseudocount

    col_totals = sub.sum(axis=0)  # (n_cells,)
    with np.errstate(divide="ignore", invalid="ignore"):
        props = np.where(col_totals > 0, sub / col_totals, np.nan)
    return props


def _process_gene_per_gene(
    gene_id: str,
    pas_ids: list[int],
    ranks: list[int],
    count_matrix: pd.DataFrame,
    pseudocount: float,
) -> list[dict]:
    """Compute proportions for one gene (per_gene aggregation).

    Args:
        gene_id: Gene identifier.
        pas_ids: List of PAS IDs for this gene.
        ranks: Corresponding per-gene rank for each pas_id.
        count_matrix: Full count matrix.
        pseudocount: Added to each count.

    Returns:
        List of dicts with keys:
        [gene_id, transcript_id, pas_id, rank, cell, proportion].
    """
    valid = [pid for pid in pas_ids if pid in count_matrix.index]
    if not valid:
        return []

    iloc_map = {pid: count_matrix.index.get_loc(pid) for pid in valid}
    iloc_indices = [iloc_map[pid] for pid in valid]
    props = _proportions_for_rows(iloc_indices, count_matrix, pseudocount)

    cells = count_matrix.columns.tolist()
    rank_map = dict(zip(pas_ids, ranks))
    rows = []
    for j, pid in enumerate(valid):
        r = rank_map.get(pid, -1)
        for k, cell in enumerate(cells):
            rows.append({
                "gene_id": gene_id,
                "transcript_id": "_gene_",
                "pas_id": pid,
                "rank": r,
                "cell": cell,
                "proportion": float(props[j, k]),
            })
    return rows


def _process_gene_per_isoform(
    gene_id: str,
    transcript_pas: dict[str, list[tuple[int, int]]],  # tid -> [(pas_id, rank)]
    count_matrix: pd.DataFrame,
    pseudocount: float,
) -> list[dict]:
    """Compute proportions for one gene (per_isoform aggregation).

    Args:
        gene_id: Gene identifier.
        transcript_pas: transcript_id -> [(pas_id, rank), ...]
        count_matrix: Full count matrix.
        pseudocount: Added to each count.

    Returns:
        List of dicts.
    """
    rows = []
    cells = count_matrix.columns.tolist()
    valid_index = set(count_matrix.index)

    for transcript_id, pid_rank_pairs in transcript_pas.items():
        valid_pairs = [(pid, r) for pid, r in pid_rank_pairs if pid in valid_index]
        if not valid_pairs:
            continue

        pas_ids_t = [p for p, _ in valid_pairs]
        ranks_t = [r for _, r in valid_pairs]
        iloc_indices = [count_matrix.index.get_loc(pid) for pid in pas_ids_t]
        props = _proportions_for_rows(iloc_indices, count_matrix, pseudocount)

        for j, (pid, r) in enumerate(zip(pas_ids_t, ranks_t)):
            for k, cell in enumerate(cells):
                rows.append({
                    "gene_id": gene_id,
                    "transcript_id": transcript_id,
                    "pas_id": pid,
                    "rank": r,
                    "cell": cell,
                    "proportion": float(props[j, k]),
                })
    return rows


# ---------------------------------------------------------------------------
# Strategy class
# ---------------------------------------------------------------------------

@register_pdui_strategy
class ProportionPDUIStrategy(PDUIStrategy):
    """Per-PAS proportion vector strategy.

    For each (gene-or-transcript, cell), computes the proportion of reads
    at each PAS.  Proportions sum to 1.0 (NaN for zero-total cells).

    Uses scipy.sparse normalization; no full-matrix dense conversion.
    """

    name: str = "proportion"

    def compute(
        self,
        count_matrix: pd.DataFrame,
        pas_isoform_map: dict[int, list[tuple[str, str, int, int, int]]],
        aggregation: AggregationMode = "per_isoform",
        isoform_collapse: IsoformCollapseMode = "none",
        pseudocount: float = 0.0,
    ) -> pd.DataFrame:
        """Compute per-PAS proportion vectors.

        Args:
            count_matrix: Shape ``(n_pas, n_cells)`` with integer PAS ID
                index.
            pas_isoform_map: ``{pas_id: [(gene_id, transcript_id,
                transcript_pos, rank, total_pas_in_transcript), ...]}``.
            aggregation: Level at which proportions are computed.
            isoform_collapse: Not used by this strategy (proportions are
                already per-PAS).
            pseudocount: Added to each count before summing per gene.

        Returns:
            Long-format DataFrame with columns::

                gene_id, transcript_id, pas_id, rank, cell, proportion
        """
        if not pas_isoform_map:
            return pd.DataFrame(
                columns=["gene_id", "transcript_id", "pas_id", "rank",
                         "cell", "proportion"]
            )

        # Build gene -> {transcript -> [(pas_id, rank)]} index
        gene_transcript_map: dict[str, dict[str, list[tuple[int, int]]]] = {}
        gene_pas_rank: dict[str, list[tuple[int, int]]] = {}  # for per_gene

        for pas_id, entries in pas_isoform_map.items():
            for gene_id, transcript_id, _tpos, rank, _total in entries:
                # Per-isoform structure
                gene_transcript_map.setdefault(gene_id, {}).setdefault(
                    transcript_id, []
                ).append((pas_id, rank))
                # Per-gene structure (deduplicate pas_id per gene)
                gene_pas_rank.setdefault(gene_id, [])

        # Deduplicate: for per_gene, each pas_id appears once per gene
        # Use min rank across transcripts
        if aggregation == "per_gene":
            gene_pas_min_rank: dict[str, dict[int, int]] = {}
            for pas_id, entries in pas_isoform_map.items():
                for gene_id, _tid, _tpos, rank, _total in entries:
                    gene_pas_min_rank.setdefault(gene_id, {})[pas_id] = min(
                        gene_pas_min_rank.get(gene_id, {}).get(pas_id, 9999),
                        rank,
                    )

        genes = list(gene_transcript_map.keys())
        n_genes = len(genes)
        use_parallel = n_genes > 5000

        if aggregation == "per_gene":
            def _do_gene(gid: str) -> list[dict]:
                pid_rank = gene_pas_min_rank.get(gid, {})
                pas_ids = list(pid_rank.keys())
                ranks = [pid_rank[p] for p in pas_ids]
                return _process_gene_per_gene(
                    gid, pas_ids, ranks, count_matrix, pseudocount
                )

            if use_parallel:
                from ema.utils import get_resource_manager
                n_jobs = get_resource_manager().get_n_jobs(
                    per_worker_mb=200, stage="pdui"
                )
                results = Parallel(n_jobs=n_jobs)(delayed(_do_gene)(g) for g in genes)
            else:
                results = [_do_gene(g) for g in genes]

        else:  # per_isoform
            def _do_gene_iso(gid: str) -> list[dict]:
                return _process_gene_per_isoform(
                    gid, gene_transcript_map[gid], count_matrix, pseudocount
                )

            if use_parallel:
                from ema.utils import get_resource_manager
                n_jobs = get_resource_manager().get_n_jobs(
                    per_worker_mb=200, stage="pdui"
                )
                results = Parallel(n_jobs=n_jobs)(
                    delayed(_do_gene_iso)(g) for g in genes
                )
            else:
                results = [_do_gene_iso(g) for g in genes]

        rows: list[dict] = []
        for sub in results:
            rows.extend(sub)

        if not rows:
            return pd.DataFrame(
                columns=["gene_id", "transcript_id", "pas_id", "rank",
                         "cell", "proportion"]
            )

        df = pd.DataFrame(rows)
        df["proportion"] = df["proportion"].astype("float64")
        return df
