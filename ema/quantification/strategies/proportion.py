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
    cells: np.ndarray | None = None,
) -> dict | None:
    """Compute proportions for one gene (per_gene aggregation).

    Args:
        gene_id: Gene identifier.
        pas_ids: List of PAS IDs for this gene.
        ranks: Corresponding per-gene rank for each pas_id.
        count_matrix: Full count matrix.
        pseudocount: Added to each count.
        cells: Pre-extracted cell index as object array. Passed in by the
            caller to avoid rebuilding it for every gene.

    Returns:
        Dict of per-column numpy arrays (each length n_pas * n_cells) or
        ``None`` when no valid PAS for this gene. Returning column arrays
        instead of a list of per-row dicts is the difference between
        ~250 bytes/row and ~50 bytes/row at this scale — the dict-of-rows
        materialisation peaked at 10 GB RSS on real datasets and tripped
        the OOM-killer.
    """
    valid = [pid for pid in pas_ids if pid in count_matrix.index]
    if not valid:
        return None

    iloc_indices = [count_matrix.index.get_loc(pid) for pid in valid]
    props = _proportions_for_rows(iloc_indices, count_matrix, pseudocount)
    # props: (n_pas, n_cells), float64

    if cells is None:
        cells = np.asarray(count_matrix.columns.tolist(), dtype=object)
    n_pas, n_cells = props.shape
    n_total = n_pas * n_cells

    rank_map = dict(zip(pas_ids, ranks))
    ranks_arr = np.fromiter(
        (rank_map.get(p, -1) for p in valid), dtype=np.int64, count=n_pas,
    )
    pas_arr = np.asarray(valid, dtype=np.int64)

    return {
        "gene_id": np.full(n_total, gene_id, dtype=object),
        "transcript_id": np.full(n_total, "_gene_", dtype=object),
        "pas_id": np.repeat(pas_arr, n_cells),
        "rank": np.repeat(ranks_arr, n_cells),
        "cell": np.tile(cells, n_pas),
        "proportion": props.reshape(n_total).astype(np.float64, copy=False),
    }


def _process_gene_per_isoform(
    gene_id: str,
    transcript_pas: dict[str, list[tuple[int, int]]],  # tid -> [(pas_id, rank)]
    count_matrix: pd.DataFrame,
    pseudocount: float,
    cells: np.ndarray | None = None,
) -> dict | None:
    """Compute proportions for one gene (per_isoform aggregation).

    Returns column-oriented numpy arrays (see :func:`_process_gene_per_gene`
    for the rationale).  Returns ``None`` when no transcript yields a valid
    PAS for this gene.
    """
    if cells is None:
        cells = np.asarray(count_matrix.columns.tolist(), dtype=object)
    n_cells = len(cells)
    valid_index = set(count_matrix.index)

    # Per-transcript chunks; we concatenate them into one gene-level chunk.
    chunks: list[dict] = []
    for transcript_id, pid_rank_pairs in transcript_pas.items():
        valid_pairs = [(pid, r) for pid, r in pid_rank_pairs if pid in valid_index]
        if not valid_pairs:
            continue

        pas_ids_t = [p for p, _ in valid_pairs]
        ranks_t = [r for _, r in valid_pairs]
        iloc_indices = [count_matrix.index.get_loc(pid) for pid in pas_ids_t]
        props = _proportions_for_rows(iloc_indices, count_matrix, pseudocount)
        n_pas = props.shape[0]
        n_total = n_pas * n_cells

        chunks.append({
            "gene_id": np.full(n_total, gene_id, dtype=object),
            "transcript_id": np.full(n_total, transcript_id, dtype=object),
            "pas_id": np.repeat(np.asarray(pas_ids_t, dtype=np.int64), n_cells),
            "rank": np.repeat(np.asarray(ranks_t, dtype=np.int64), n_cells),
            "cell": np.tile(cells, n_pas),
            "proportion": props.reshape(n_total).astype(np.float64, copy=False),
        })

    if not chunks:
        return None
    if len(chunks) == 1:
        return chunks[0]
    return {
        col: np.concatenate([c[col] for c in chunks])
        for col in ("gene_id", "transcript_id", "pas_id", "rank",
                    "cell", "proportion")
    }


# ---------------------------------------------------------------------------
# Strategy class
# ---------------------------------------------------------------------------

@register_pdui_strategy
class ProportionPDUIStrategy(PDUIStrategy):
    """Per-PAS proportion vector strategy.

    For each (gene-or-transcript, cell), computes the proportion of reads
    at each PAS.  Proportions sum to 1.0 (NaN for zero-total cells).

    Uses scipy.sparse normalization; no full-matrix dense conversion.

    Tunable hyperparameters:
        pseudocount (default 0.0): Added to each count before per-gene
            normalization.  The default 0.0 preserves original behaviour;
            set to e.g. 1.0 to avoid NaN for zero-total cells.
            CLI: ``--pdui-pseudocount`` / YAML: ``pdui_pseudocount``.
        aggregation (default "per_isoform"): Whether proportions are computed
            per gene or per isoform.  CLI: ``--isoform-agg`` /
            YAML: ``isoform_agg``.
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

        # Pre-extract cell index ONCE; passed to every gene worker so we
        # don't pay 22 k list rebuilds.
        cells_arr = np.asarray(count_matrix.columns.tolist(), dtype=object)

        if aggregation == "per_gene":
            def _do_gene(gid: str) -> dict | None:
                pid_rank = gene_pas_min_rank.get(gid, {})
                pas_ids = list(pid_rank.keys())
                ranks = [pid_rank[p] for p in pas_ids]
                return _process_gene_per_gene(
                    gid, pas_ids, ranks, count_matrix, pseudocount,
                    cells=cells_arr,
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
            def _do_gene_iso(gid: str) -> dict | None:
                return _process_gene_per_isoform(
                    gid, gene_transcript_map[gid], count_matrix, pseudocount,
                    cells=cells_arr,
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

        # Drop empty-gene results and assemble column-wise.  Going through
        # a long-format dict-per-row list peaks at ~10 GB on real datasets;
        # column-wise concat is bounded by the final DataFrame footprint.
        chunks = [r for r in results if r is not None]
        if not chunks:
            return pd.DataFrame(
                columns=["gene_id", "transcript_id", "pas_id", "rank",
                         "cell", "proportion"]
            )
        df = pd.DataFrame({
            col: np.concatenate([c[col] for c in chunks])
            for col in ("gene_id", "transcript_id", "pas_id", "rank",
                        "cell", "proportion")
        })
        df["proportion"] = df["proportion"].astype("float64", copy=False)
        return df
