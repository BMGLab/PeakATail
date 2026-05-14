"""Shared helpers for the gene-track visualisation.

The gene-track viz renders, for one gene:

  * gene structure (isoforms, when a GTF is supplied)
  * per-cluster PAS coverage (reads / cell normalised)
  * per-cluster within-gene PAS proportions

This module collects the data-prep work shared between the matplotlib
and plotly back-ends so neither has to duplicate the per-cluster
aggregation or GTF parsing.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class GenePanel:
    """All data needed to render one gene's track view.

    Attributes:
        gene_id: Ensembl / public gene identifier.
        chrom, start, end: Genomic span covering the gene's PAS (BED 0-based,
            half-open).
        strand: '+' or '-'.  Used to orient the track left->right by direction
            of transcription.
        pas_ids: Integer PAS identifiers (one per peak).
        pas_positions: Genomic position of each PAS (same length as pas_ids).
        clusters: Cluster labels in display order (sorted natural).
        n_cells_per_cluster: Number of cells per cluster (same order).
        reads: (n_clusters, n_pas) sum of reads per PAS per cluster.
        reads_per_cell: (n_clusters, n_pas) reads / n_cells (depth-normalised).
        proportions: (n_clusters, n_pas) within-gene proportion (rows ~ sum to 1
            for clusters where the gene was expressed).
        isoforms: Optional list of (transcript_id, [(exon_start, exon_end), ...])
            tuples for the gene structure track.
    """
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str
    pas_ids: list[int]
    pas_positions: list[int]
    clusters: list[str]
    n_cells_per_cluster: list[int]
    reads: np.ndarray
    reads_per_cell: np.ndarray
    proportions: np.ndarray
    isoforms: list[tuple[str, list[tuple[int, int]]]] = field(default_factory=list)
    # Actual genomic span of each PAS region (BED half-open).  Same length
    # and order as ``pas_ids`` / ``pas_positions``.  When the merger widens
    # a PAS, ``pas_ends - pas_starts`` reflects that.  The rendering layer
    # uses these to draw bars at their real width rather than a fixed
    # 1%-of-gene-span placeholder.
    pas_starts: list[int] = field(default_factory=list)
    pas_ends: list[int] = field(default_factory=list)
    # Human-readable gene symbol from the GTF (e.g. "CLIC2", "DPYD").
    # When available it's used in figure titles in front of the Ensembl ID
    # so the reader sees ``CLIC2 (ENSG00000155962) — chrX:...`` instead of
    # the opaque Ensembl identifier alone.  Empty when no GTF supplied or
    # the gene wasn't found.
    gene_name: str = ""


# ---------------------------------------------------------------------------
# Top-N gene ranking from a switch_diff result set
# ---------------------------------------------------------------------------

def rank_top_genes(
    pair_results: dict[tuple[str, str | None], pd.DataFrame],
    n: int = 5,
    qvalue_col: str = "qvalue",
    log2fc_col: str = "log2fc",
    gene_id_col: str = "gene_id",
) -> list[str]:
    """Pick the top-N most-shifted genes across all cluster pairs.

    Score per (pas_id, pair) = ``-log10(q) * |log2fc|`` (volcano score).
    Score per gene = max over its PAS across all pairs.

    Args:
        pair_results: dict from (cluster1, cluster2) to result DataFrame.
        n: Number of genes to return.
        qvalue_col, log2fc_col, gene_id_col: column names in the result frames.

    Returns:
        list of gene IDs (length <= n).  Empty if no result has gene_id.
    """
    score_by_gene: dict[str, float] = {}
    for _, df in pair_results.items():
        if df is None or df.empty or gene_id_col not in df.columns:
            continue
        if qvalue_col not in df.columns or log2fc_col not in df.columns:
            continue
        # Replace inf/-inf in log2fc so they don't dominate the ranking.
        lf = pd.to_numeric(df[log2fc_col], errors="coerce").replace(
            [np.inf, -np.inf], np.nan
        )
        q = pd.to_numeric(df[qvalue_col], errors="coerce")
        q = q.where(q > 0, 1e-300)
        score = (-np.log10(q.fillna(1.0))) * lf.abs().fillna(0.0)
        for gene, s in zip(df[gene_id_col].astype(str), score):
            if not gene or gene == "nan":
                continue
            prev = score_by_gene.get(gene, -np.inf)
            if s > prev:
                score_by_gene[gene] = float(s)
    return sorted(score_by_gene, key=lambda g: score_by_gene[g], reverse=True)[:n]


# ---------------------------------------------------------------------------
# Gene panel assembly
# ---------------------------------------------------------------------------

def build_gene_panel(
    gene_id: str,
    adata: Any,
    pasbed: pd.DataFrame,
    cluster_key: str = "leiden",
    gene_id_col: str = "gene_id",
    isoforms: list[tuple[str, list[tuple[int, int]]]] | None = None,
    gene_name: str = "",
) -> GenePanel | None:
    """Build a :class:`GenePanel` for one gene.

    Args:
        gene_id: Gene to extract.
        adata: AnnData with ``var[gene_id_col]`` and ``obs[cluster_key]``.
            ``X`` is the cells x PAS count matrix (sparse or dense).
        pasbed: BED-format DataFrame with at least columns
            ``[chrom, start, end, pas_id, strand]``.  Indexed by pas_id (str)
            or carrying a ``pas_id`` column.
        cluster_key: Column in ``obs`` for cluster labels.
        gene_id_col: Column in ``var`` carrying gene identifiers.
        isoforms: Optional pre-loaded isoform structure.

    Returns:
        :class:`GenePanel` or ``None`` if the gene has no PAS in the AnnData.
    """
    if gene_id_col not in adata.var.columns or cluster_key not in adata.obs.columns:
        return None
    var = adata.var
    mask = var[gene_id_col].astype(str) == str(gene_id)
    if not mask.any():
        return None
    pas_ids_var = var.index[mask].astype(str).tolist()
    pas_ids = [int(p) for p in pas_ids_var if str(p).isdigit()]
    if not pas_ids:
        return None

    # Look up coords (pasbed indexed by pas_id, str-keyed).
    if "pas_id" in pasbed.columns:
        pb = pasbed.set_index(pasbed["pas_id"].astype(str))
    else:
        pb = pasbed.copy()
        pb.index = pb.index.astype(str)
    coords = pb.reindex([str(p) for p in pas_ids]).dropna(subset=["chrom"])
    if coords.empty:
        return None
    # Re-derive pas_ids in coord-order so positions align with coverage cols.
    pas_ids = [int(p) for p in coords.index.tolist()]
    pas_starts = [int(row["start"]) for _, row in coords.iterrows()]
    pas_ends = [int(row["end"]) for _, row in coords.iterrows()]
    pas_positions = [(s + e) // 2 for s, e in zip(pas_starts, pas_ends)]
    chrom = str(coords["chrom"].iloc[0])
    strand = str(coords["strand"].iloc[0]) if "strand" in coords.columns else "+"
    g_start = int(coords["start"].min())
    g_end = int(coords["end"].max())

    # Extract per-PAS counts (cells x len(pas_ids)).
    import scipy.sparse as _sp
    var_pos_by_id = {str(v): i for i, v in enumerate(var.index.astype(str))}
    col_idx = [var_pos_by_id[str(p)] for p in pas_ids]
    X = adata.X[:, col_idx]
    if _sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)  # (n_cells, n_pas)

    # Per-cluster aggregation.
    clusters = sorted(
        adata.obs[cluster_key].astype(str).unique(),
        key=lambda x: int(x) if x.isdigit() else x,
    )
    cluster_labels = adata.obs[cluster_key].astype(str).values
    n_clusters = len(clusters)
    n_pas = len(pas_ids)
    reads = np.zeros((n_clusters, n_pas), dtype=np.float64)
    n_cells = np.zeros(n_clusters, dtype=np.int64)
    for i, c in enumerate(clusters):
        rows = cluster_labels == c
        n_cells[i] = int(rows.sum())
        if n_cells[i] > 0:
            reads[i, :] = X[rows, :].sum(axis=0)
    rpc = np.where(n_cells[:, None] > 0, reads / np.maximum(n_cells[:, None], 1), 0.0)
    row_totals = reads.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        proportions = np.where(row_totals > 0, reads / row_totals, np.nan)

    return GenePanel(
        gene_id=str(gene_id),
        gene_name=str(gene_name),
        chrom=chrom,
        start=g_start,
        end=g_end,
        strand=strand,
        pas_ids=pas_ids,
        pas_positions=pas_positions,
        pas_starts=pas_starts,
        pas_ends=pas_ends,
        clusters=clusters,
        n_cells_per_cluster=n_cells.tolist(),
        reads=reads,
        reads_per_cell=rpc,
        proportions=proportions,
        isoforms=isoforms or [],
    )


# ---------------------------------------------------------------------------
# Isoform structure helpers
# ---------------------------------------------------------------------------

def load_gene_name_from_gtf(gtf_path: Path, gene_id: str) -> str:
    """Return the ``gene_name`` attribute for *gene_id* from the GTF (or '').

    Streams the GTF, stopping at the first matching ``gene`` feature.  Used
    to populate :attr:`GenePanel.gene_name` so figure titles can show a
    human-readable symbol next to the Ensembl ID.

    Args:
        gtf_path: Path to a GTF annotation file.
        gene_id: Ensembl gene ID to look up.

    Returns:
        The gene_name attribute, or ``""`` if not found / file unreadable.
    """
    needle = f'gene_id "{gene_id}"'
    try:
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                if needle not in line:
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                # Match in any line that references this gene_id; gene_name
                # is identical across all rows for one gene so any hit
                # suffices.
                name = _gtf_attr(parts[8], "gene_name")
                if name:
                    return name
    except OSError as exc:
        log.warning("load_gene_name_from_gtf: could not read %s: %s", gtf_path, exc)
    return ""


def load_isoforms_for_gene(
    gtf_path: Path,
    gene_id: str,
    feature_type: str = "exon",
) -> list[tuple[str, list[tuple[int, int]]]]:
    """Return ``[(transcript_id, [(exon_start, exon_end), ...]), ...]``.

    Streams the GTF once per call.  For one-off invocation that's fine;
    callers rendering many genes should batch via a parsed cache.
    """
    isoforms: dict[str, list[tuple[int, int]]] = {}
    gene_id_str = str(gene_id)
    try:
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != feature_type:
                    continue
                attrs = parts[8]
                if gene_id_str not in attrs:
                    continue
                # crude attribute parse
                gid = _gtf_attr(attrs, "gene_id")
                if gid != gene_id_str:
                    continue
                tid = _gtf_attr(attrs, "transcript_id") or "unknown"
                isoforms.setdefault(tid, []).append((int(parts[3]) - 1, int(parts[4])))
    except OSError as exc:
        log.warning("load_isoforms_for_gene: could not read %s: %s", gtf_path, exc)
        return []
    return [(t, sorted(exons)) for t, exons in isoforms.items()]


def _gtf_attr(attr_field: str, key: str) -> str:
    """Extract a single attribute value from a GTF attribute field."""
    needle = f'{key} "'
    i = attr_field.find(needle)
    if i < 0:
        return ""
    j = attr_field.find('"', i + len(needle))
    if j < 0:
        return ""
    return attr_field[i + len(needle):j]
