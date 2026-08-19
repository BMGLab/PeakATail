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
    rank           int   — proximal-to-distal rank of the PAS within the
                           gene/transcript, as supplied by *pas_isoform_map*
                           (see :mod:`ema.quantification.pas_to_isoform`).
    cell           str
    proportion     float64   — sums to 1.0 per (gene/transcript, cell)
    reads_at_pas   float64   — raw read count behind ``proportion`` (never
                               includes ``pseudocount``)
    total_reads_gene / total_reads_transcript
                   float64   — the REAL summed read count for the gene (or
                               transcript) in that cell — i.e. the sum of
                               ``reads_at_pas`` across the gene's PAS, never
                               inflated by ``pseudocount``.

(gene/transcript, cell) pairs with ZERO real coverage — the gene has no
reads at all in that cell — are OMITTED entirely rather than emitting a
synthetic ``proportion = 1/n_PAS`` placeholder. This matters in particular
when ``pseudocount > 0``: naively adding the same pseudocount to every PAS
before normalizing would make an all-zero cell's proportions collapse to a
uniform ``1/n_PAS`` prior (and its "total" would equal ``n_PAS *
pseudocount`` — a PAS *count*, not a read count) — this is real biological
noise, not signal, and previously flooded ~98% of rows on production data,
making every downstream trend a padding artifact. ``pseudocount`` is still
applied to the proportion *ratio* for cells that DO have real coverage
(the original opt-in smoothing knob is unchanged there); it is only
disallowed from manufacturing coverage that does not exist.

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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return per-PAS proportions plus the raw read counts behind them.

    Args:
        row_indices: Positional (iloc) indices into *count_matrix*.
        count_matrix: Shape ``(n_pas, n_cells)``.
        pseudocount: Added to each count before normalization (applied only
            to the proportion *ratio* — never to ``raw_totals``, see below).

    Returns:
        ``(props, reads, raw_totals, covered)``:
          * ``props``      — (n_selected_pas, n_cells) proportion. NaN only
            if ``pseudocount == 0`` AND the column has zero raw coverage
            (callers must still drop those columns — see ``covered``).
          * ``reads``      — (n_selected_pas, n_cells) raw counts, BEFORE
            pseudocount.
          * ``raw_totals`` — (n_cells,) sum of RAW reads across selected PAS
            (the real per-cell, per-gene read depth — never includes
            ``pseudocount``, so it can't be mistaken for a PAS count).
          * ``covered``    — (n_cells,) boolean, True where ``raw_totals >
            0``, i.e. the gene/transcript genuinely has reads in that cell.
            ``pseudocount`` must never manufacture coverage: a
            zero-coverage column always has ``covered=False`` regardless of
            ``pseudocount``, and callers MUST drop those columns rather
            than emit a padded ``1/n_PAS`` proportion.

        Returning reads + raw_totals so the row-wise TSV can carry the raw
        evidence behind each proportion (``reads_at_pas`` and
        ``total_reads_gene``).  Researchers can then audit a 0.8
        proportion as e.g. "8 of 10 reads", not just trust the ratio.
    """
    reads = count_matrix.iloc[row_indices].values.astype(float)
    raw_totals = reads.sum(axis=0)  # (n_cells,) — TRUE coverage, no pseudocount
    covered = raw_totals > 0
    sub = reads + pseudocount if pseudocount != 0.0 else reads
    col_totals = sub.sum(axis=0)  # (n_cells,) — denominator for the ratio only
    with np.errstate(divide="ignore", invalid="ignore"):
        props = np.where(col_totals > 0, sub / col_totals, np.nan)
    return props, reads, raw_totals, covered


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
        Dict of per-column numpy arrays (each length n_pas * n_covered_cells)
        or ``None`` when there's no valid PAS for this gene, or no cell has
        real coverage (gene total reads > 0). Zero-coverage (gene, cell)
        pairs are OMITTED — never padded with a synthetic ``1/n_PAS``
        proportion (see module docstring). Returning column arrays instead
        of a list of per-row dicts is the difference between ~250 bytes/row
        and ~50 bytes/row at this scale — the dict-of-rows materialisation
        peaked at 10 GB RSS on real datasets and tripped the OOM-killer.
    """
    valid = [pid for pid in pas_ids if pid in count_matrix.index]
    if not valid:
        return None

    iloc_indices = [count_matrix.index.get_loc(pid) for pid in valid]
    props, reads, raw_totals, covered = _proportions_for_rows(
        iloc_indices, count_matrix, pseudocount,
    )
    if not covered.any():
        return None

    if cells is None:
        cells = np.asarray(count_matrix.columns.tolist(), dtype=object)
    n_pas = props.shape[0]

    # Restrict to cells with real (gene, cell) coverage — drop the rest
    # entirely rather than emitting a padded 1/n_PAS placeholder.
    props = props[:, covered]
    reads = reads[:, covered]
    raw_totals = raw_totals[covered]
    cells = cells[covered]
    n_cells = int(covered.sum())
    n_total = n_pas * n_cells

    rank_map = dict(zip(pas_ids, ranks))
    ranks_arr = np.fromiter(
        (rank_map.get(p, -1) for p in valid), dtype=np.int64, count=n_pas,
    )
    pas_arr = np.asarray(valid, dtype=np.int64)
    totals_tiled = np.tile(raw_totals.astype(np.float64), n_pas)

    return {
        "gene_id": np.full(n_total, gene_id, dtype=object),
        "transcript_id": np.full(n_total, "_gene_", dtype=object),
        "pas_id": np.repeat(pas_arr, n_cells),
        "rank": np.repeat(ranks_arr, n_cells),
        "cell": np.tile(cells, n_pas),
        "proportion": props.reshape(n_total).astype(np.float64, copy=False),
        "reads_at_pas": reads.reshape(n_total).astype(np.float64, copy=False),
        "total_reads_gene": totals_tiled,
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
    for the rationale).  Returns ``None`` when no transcript yields any
    (transcript, cell) pair with real coverage. Zero-coverage
    (transcript, cell) pairs are OMITTED per-transcript (each transcript's
    coverage is judged independently, since different transcripts of the
    same gene may not share the same PAS set).
    """
    if cells is None:
        cells = np.asarray(count_matrix.columns.tolist(), dtype=object)
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
        props, reads, raw_totals, covered = _proportions_for_rows(
            iloc_indices, count_matrix, pseudocount,
        )
        if not covered.any():
            continue
        n_pas = props.shape[0]

        props_c = props[:, covered]
        reads_c = reads[:, covered]
        raw_totals_c = raw_totals[covered]
        cells_c = cells[covered]
        n_cells_c = int(covered.sum())
        n_total = n_pas * n_cells_c
        totals_tiled = np.tile(raw_totals_c.astype(np.float64), n_pas)

        chunks.append({
            "gene_id": np.full(n_total, gene_id, dtype=object),
            "transcript_id": np.full(n_total, transcript_id, dtype=object),
            "pas_id": np.repeat(np.asarray(pas_ids_t, dtype=np.int64), n_cells_c),
            "rank": np.repeat(np.asarray(ranks_t, dtype=np.int64), n_cells_c),
            "cell": np.tile(cells_c, n_pas),
            "proportion": props_c.reshape(n_total).astype(np.float64, copy=False),
            "reads_at_pas": reads_c.reshape(n_total).astype(np.float64, copy=False),
            "total_reads_transcript": totals_tiled,
        })

    if not chunks:
        return None
    if len(chunks) == 1:
        return chunks[0]
    # NOTE: previously this concatenation silently DROPPED "reads_at_pas"
    # and "total_reads_transcript" for genes with >1 transcript chunk (the
    # column tuple below used to stop at "proportion"), which would raise a
    # KeyError one level up in compute() the first time such a gene was
    # encountered.  All 8 columns must survive the multi-chunk concat.
    return {
        col: np.concatenate([c[col] for c in chunks])
        for col in ("gene_id", "transcript_id", "pas_id", "rank", "cell",
                    "proportion", "reads_at_pas", "total_reads_transcript")
    }


# ---------------------------------------------------------------------------
# Strategy class
# ---------------------------------------------------------------------------

@register_pdui_strategy
class ProportionPDUIStrategy(PDUIStrategy):
    """Per-PAS proportion vector strategy.

    For each (gene-or-transcript, cell), computes the proportion of reads
    at each PAS.  Proportions sum to 1.0 over the emitted rows.
    Zero-coverage (gene/transcript, cell) pairs are OMITTED — not padded
    with a synthetic ``1/n_PAS`` value, regardless of ``pseudocount``
    (see module docstring).

    Uses scipy.sparse normalization; no full-matrix dense conversion.

    Tunable hyperparameters:
        pseudocount (default 0.0): Added to the proportion ratio's
            numerator/denominator for cells that already have real
            coverage.  The default 0.0 preserves original behaviour for
            those cells; set to e.g. 1.0 to smooth per-PAS proportions
            when some individual PAS (but not the whole gene) have zero
            reads. It can no longer manufacture coverage for a genuinely
            zero-read (gene, cell) pair — those rows are omitted outright.
            CLI: ``--pdui-pseudocount`` / YAML: ``pdui_pseudocount``.
        aggregation (default "per_isoform"): Whether proportions are computed
            per gene or per isoform.  CLI: ``--isoform-agg`` /
            YAML: ``isoform_agg``.
    """

    name: str = "proportion"
    output_filename: str = "proportion.tsv"
    score_column: str = "proportion"

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
        # Column order: identifiers, score, raw read evidence.
        # The denominator column differs by aggregation:
        #   per_gene     -> total_reads_gene
        #   per_isoform  -> total_reads_transcript
        if aggregation == "per_gene":
            col_order = ("gene_id", "transcript_id", "pas_id", "rank", "cell",
                         "proportion", "reads_at_pas", "total_reads_gene")
        else:
            col_order = ("gene_id", "transcript_id", "pas_id", "rank", "cell",
                         "proportion", "reads_at_pas", "total_reads_transcript")
        if not chunks:
            return pd.DataFrame(columns=list(col_order))
        df = pd.DataFrame({
            col: np.concatenate([c[col] for c in chunks])
            for col in col_order
        })
        df["proportion"] = df["proportion"].astype("float64", copy=False)
        return df
