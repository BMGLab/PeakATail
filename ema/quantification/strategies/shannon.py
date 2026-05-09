"""Shannon entropy PDUI strategy.

For each (gene-or-transcript, cell), computes:

    H      = -Σ p_i · log2(p_i)          (where p_i = count_i / total_counts)
    H_norm = H / log2(N)                  (N = number of PAS with p_i > 0)

Boundary cases:
    - All reads on one PAS (delta distribution): H = 0.0, H_norm = 0.0
    - Uniform distribution (p_i = 1/N for all i): H = log2(N), H_norm = 1.0
    - Single PAS gene (N=1): H = 0.0, H_norm = NaN (log2(1) = 0, undefined)
    - Zero-total cell: H = NaN, H_norm = NaN

Output schema (long-format):

    gene_id             str
    transcript_id       str   — transcript ID or "_gene_" sentinel
    cell                str
    entropy             float64  — H in bits
    normalized_entropy  float64  — H / log2(N), in [0, 1]
    n_pas               int      — number of PAS in the gene/transcript

Parallelized via ``joblib.Parallel`` when ``n_genes > 5000``.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from ema.quantification.strategies.base import (
    AggregationMode,
    IsoformCollapseMode,
    PDUIStrategy,
)
from ema.quantification.strategies import register_pdui_strategy


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def _shannon_from_counts(
    counts_matrix: np.ndarray,
    pseudocount: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Compute per-cell Shannon entropy from a (n_pas, n_cells) count array.

    Args:
        counts_matrix: float64 array of shape ``(n_pas, n_cells)``.
        pseudocount: Added to each count before computing proportions.

    Returns:
        Tuple of:
            - entropy: float64 array of shape ``(n_cells,)`` — H in bits.
            - normalized_entropy: float64 array of shape ``(n_cells,)`` —
              H / log2(N) where N = n_pas.  NaN when N <= 1 or total = 0.
            - n_pas: int — number of PAS (rows in counts_matrix).
    """
    n_pas, n_cells = counts_matrix.shape
    counts = counts_matrix.astype(float)
    if pseudocount != 0.0:
        counts = counts + pseudocount

    col_totals = counts.sum(axis=0)  # (n_cells,)

    with np.errstate(divide="ignore", invalid="ignore"):
        props = np.where(col_totals > 0, counts / col_totals, 0.0)
        # H = -Σ p * log2(p), treat 0*log2(0)=0
        log_props = np.where(props > 0, np.log2(props), 0.0)
        entropy = -(props * log_props).sum(axis=0)  # (n_cells,)

    # Set entropy to NaN for zero-total cells
    entropy = np.where(col_totals > 0, entropy, np.nan)

    # Normalized entropy
    if n_pas > 1:
        log2_n = math.log2(n_pas)
        normalized_entropy = np.where(col_totals > 0, entropy / log2_n, np.nan)
    else:
        normalized_entropy = np.full(n_cells, np.nan)

    return entropy, normalized_entropy, n_pas


def _process_gene_per_gene(
    gene_id: str,
    pas_ids: list[int],
    count_matrix: pd.DataFrame,
    pseudocount: float,
) -> list[dict]:
    """Compute entropy for one gene (per_gene aggregation).

    Args:
        gene_id: Gene identifier.
        pas_ids: List of PAS IDs belonging to this gene.
        count_matrix: Full count matrix.
        pseudocount: Added per cell count.

    Returns:
        List of dicts with keys [gene_id, transcript_id, cell, entropy,
        normalized_entropy, n_pas].
    """
    valid = [pid for pid in pas_ids if pid in count_matrix.index]
    if not valid:
        return []

    sub = count_matrix.loc[valid].values.astype(float)  # (n_valid, n_cells)
    entropy, norm_entropy, n_pas = _shannon_from_counts(sub, pseudocount)
    cells = count_matrix.columns.tolist()

    return [
        {
            "gene_id": gene_id,
            "transcript_id": "_gene_",
            "cell": cell,
            "entropy": float(entropy[i]),
            "normalized_entropy": float(norm_entropy[i]),
            "n_pas": n_pas,
        }
        for i, cell in enumerate(cells)
    ]


def _process_gene_per_isoform(
    gene_id: str,
    transcript_pas: dict[str, list[int]],  # tid -> [pas_id, ...]
    count_matrix: pd.DataFrame,
    pseudocount: float,
) -> list[dict]:
    """Compute entropy for one gene (per_isoform aggregation).

    Args:
        gene_id: Gene identifier.
        transcript_pas: transcript_id -> list of pas_ids.
        count_matrix: Full count matrix.
        pseudocount: Added per cell count.

    Returns:
        List of dicts.
    """
    rows = []
    cells = count_matrix.columns.tolist()
    valid_index = set(count_matrix.index)

    for transcript_id, pas_ids in transcript_pas.items():
        valid = [pid for pid in pas_ids if pid in valid_index]
        if not valid:
            continue

        sub = count_matrix.loc[valid].values.astype(float)
        entropy, norm_entropy, n_pas = _shannon_from_counts(sub, pseudocount)

        for i, cell in enumerate(cells):
            rows.append({
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "cell": cell,
                "entropy": float(entropy[i]),
                "normalized_entropy": float(norm_entropy[i]),
                "n_pas": n_pas,
            })

    return rows


# ---------------------------------------------------------------------------
# Strategy class
# ---------------------------------------------------------------------------

@register_pdui_strategy
class ShannonPDUIStrategy(PDUIStrategy):
    """Shannon entropy of PAS usage per gene per cell.

    H = 0 when all reads go to one PAS (most specific).
    H = log2(N) for uniform usage across N PAS (least specific).
    H_norm = H / log2(N) normalizes to [0, 1].
    """

    name: str = "shannon"

    def compute(
        self,
        count_matrix: pd.DataFrame,
        pas_isoform_map: dict[int, list[tuple[str, str, int, int, int]]],
        aggregation: AggregationMode = "per_isoform",
        isoform_collapse: IsoformCollapseMode = "none",
        pseudocount: float = 0.0,
    ) -> pd.DataFrame:
        """Compute Shannon entropy of PAS usage.

        Args:
            count_matrix: Shape ``(n_pas, n_cells)`` with integer PAS ID
                index.
            pas_isoform_map: ``{pas_id: [(gene_id, transcript_id,
                transcript_pos, rank, total_pas_in_transcript), ...]}``.
            aggregation: ``"per_gene"`` or ``"per_isoform"``.
            isoform_collapse: Not used — entropy is computed per-isoform or
                per-gene natively.
            pseudocount: Added to each count before computing proportions.

        Returns:
            Long-format DataFrame with columns::

                gene_id, transcript_id, cell, entropy,
                normalized_entropy, n_pas
        """
        if not pas_isoform_map:
            return pd.DataFrame(
                columns=["gene_id", "transcript_id", "cell",
                         "entropy", "normalized_entropy", "n_pas"]
            )

        # Build gene -> {transcript -> [pas_ids]} index
        gene_transcript_pas: dict[str, dict[str, list[int]]] = {}
        gene_pas_set: dict[str, set[int]] = {}

        for pas_id, entries in pas_isoform_map.items():
            for gene_id, transcript_id, _tpos, _rank, _total in entries:
                gene_transcript_pas.setdefault(gene_id, {}).setdefault(
                    transcript_id, []
                ).append(pas_id)
                gene_pas_set.setdefault(gene_id, set()).add(pas_id)

        genes = list(gene_transcript_pas.keys())
        n_genes = len(genes)
        use_parallel = n_genes > 5000

        if aggregation == "per_gene":
            def _do_gene(gid: str) -> list[dict]:
                return _process_gene_per_gene(
                    gid, list(gene_pas_set[gid]), count_matrix, pseudocount
                )

            if use_parallel:
                results = Parallel(n_jobs=-1)(delayed(_do_gene)(g) for g in genes)
            else:
                results = [_do_gene(g) for g in genes]

        else:  # per_isoform
            def _do_gene_iso(gid: str) -> list[dict]:
                return _process_gene_per_isoform(
                    gid, gene_transcript_pas[gid], count_matrix, pseudocount
                )

            if use_parallel:
                results = Parallel(n_jobs=-1)(
                    delayed(_do_gene_iso)(g) for g in genes
                )
            else:
                results = [_do_gene_iso(g) for g in genes]

        rows: list[dict] = []
        for sub in results:
            rows.extend(sub)

        if not rows:
            return pd.DataFrame(
                columns=["gene_id", "transcript_id", "cell",
                         "entropy", "normalized_entropy", "n_pas"]
            )

        df = pd.DataFrame(rows)
        df["entropy"] = df["entropy"].astype("float64")
        df["normalized_entropy"] = df["normalized_entropy"].astype("float64")
        df["n_pas"] = df["n_pas"].astype(int)
        return df
