"""Proximal-Distal Usage Index (PDUI) calculation.

PDUI quantifies 3'UTR usage:
    PDUI = distal_count / (proximal_count + distal_count)

- PDUI = 0.0: all reads at proximal PAS (3'UTR shortening)
- PDUI = 1.0: all reads at distal PAS (3'UTR lengthening)
- delta-PDUI between conditions: APA switching

For genes with >2 PAS, proximal = first PAS (closest to stop codon)
and distal = last PAS (furthest from stop codon). Strand is respected:
    + strand: proximal = min coordinate, distal = max coordinate
    - strand: proximal = max coordinate, distal = min coordinate
"""

import numpy as np
import pandas as pd
from scipy.io import mmread
from typing import Optional, Tuple


def load_pas_data(matrix_path: str,
                  cellbarcodes_path: str,
                  pas_gene_path: str,
                  pasbed_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load PAS count matrix and coordinate data.

    Args:
        matrix_path: Path to filtered count matrix (MTX format).
        cellbarcodes_path: Path to cell barcodes TSV.
        pas_gene_path: Path to PAS-gene mapping TSV.
        pasbed_path: Path to PAS BED file with coordinates.

    Returns:
        Tuple of (count_matrix, pas_info):
            - count_matrix: DataFrame (PAS x cells) with counts.
            - pas_info: DataFrame with columns [chrom, start, end,
              pas_id, score, strand, gene_id].
    """
    # Load sparse matrix and metadata
    sparse_matrix = mmread(matrix_path).tocsr()
    cellbarcodes = pd.read_csv(cellbarcodes_path, sep='\t', header=None)[0].values
    pas_gene = pd.read_csv(pas_gene_path, sep='\t', header=None, names=['pas_id', 'gene_id'])

    # Load BED file with coordinates
    pasbed = pd.read_csv(
        pasbed_path, sep='\t', header=None,
        names=['chrom', 'start', 'end', 'pas_id', 'score', 'strand'],
    )

    # Merge PAS coordinates with gene assignments
    pas_info = pasbed.merge(pas_gene, on='pas_id', how='inner')

    # Only keep PAS IDs within matrix range (sparse PAS IDs)
    n_rows = sparse_matrix.shape[0]
    pas_info = pas_info[(pas_info['pas_id'] >= 1) & (pas_info['pas_id'] <= n_rows)].copy()

    # Extract only annotated rows from sparse matrix (avoids full dense conversion)
    row_indices = pas_info['pas_id'].values - 1  # 0-based
    sub_matrix = sparse_matrix[row_indices, :].toarray()

    count_matrix = pd.DataFrame(
        sub_matrix,
        index=pas_info['pas_id'].values,
        columns=cellbarcodes[:sparse_matrix.shape[1]],
    )

    return count_matrix, pas_info


def calculate_pdui(count_matrix: pd.DataFrame,
                   pas_info: pd.DataFrame,
                   pseudocount: float = 0.0) -> pd.DataFrame:
    """Calculate PDUI for each gene across all cells/clusters.

    For each gene with >=2 PAS:
    1. Sort PAS by genomic coordinate (respecting strand)
    2. Identify proximal (closest to stop codon) and distal (furthest)
    3. PDUI = distal / (proximal + distal)

    Args:
        count_matrix: DataFrame (PAS x cells/clusters) with counts.
            Index must be PAS identifiers matching pas_info.
        pas_info: DataFrame with columns [chrom, start, end, pas_id,
            score, strand, gene_id].
        pseudocount: Small value added to avoid division by zero
            (default: 0.0, NaN for zero-count genes).

    Returns:
        DataFrame (genes x cells/clusters) with PDUI values.
        Genes with <2 PAS are excluded.
    """
    # Group PAS by gene
    gene_groups = pas_info.groupby('gene_id')

    pdui_results = {}

    for gene_id, gene_pas in gene_groups:
        if len(gene_pas) < 2:
            continue

        # Sort PAS by coordinate
        gene_pas_sorted = gene_pas.sort_values('start')
        strand = gene_pas_sorted['strand'].iloc[0]

        # Identify proximal and distal PAS based on strand
        if strand == '+':
            # + strand: proximal = min coord (closest to stop codon)
            proximal_id = gene_pas_sorted.iloc[0]['pas_id']
            distal_id = gene_pas_sorted.iloc[-1]['pas_id']
        else:
            # - strand: proximal = max coord (closest to stop codon)
            proximal_id = gene_pas_sorted.iloc[-1]['pas_id']
            distal_id = gene_pas_sorted.iloc[0]['pas_id']

        # Get counts for proximal and distal PAS
        if proximal_id not in count_matrix.index or distal_id not in count_matrix.index:
            continue

        proximal_counts = count_matrix.loc[proximal_id].values.astype(float)
        distal_counts = count_matrix.loc[distal_id].values.astype(float)

        # Calculate PDUI
        total = proximal_counts + distal_counts + pseudocount
        with np.errstate(divide='ignore', invalid='ignore'):
            pdui = np.where(total > 0, distal_counts / total, np.nan)

        pdui_results[gene_id] = pdui

    if len(pdui_results) == 0:
        return pd.DataFrame(columns=count_matrix.columns)

    pdui_df = pd.DataFrame(
        pdui_results,
        index=count_matrix.columns,
    ).T
    pdui_df.index.name = 'gene_id'

    return pdui_df


def calculate_delta_pdui(pdui_matrix: pd.DataFrame,
                         cluster1: str,
                         cluster2: str) -> pd.Series:
    """Calculate delta-PDUI between two clusters.

    delta-PDUI = PDUI(cluster2) - PDUI(cluster1)
    Positive = lengthening in cluster2 relative to cluster1.
    Negative = shortening in cluster2 relative to cluster1.

    Args:
        pdui_matrix: DataFrame (genes x clusters) with PDUI values.
        cluster1: Name of first cluster column.
        cluster2: Name of second cluster column.

    Returns:
        Series of delta-PDUI values per gene.
    """
    delta = pdui_matrix[cluster2] - pdui_matrix[cluster1]
    delta.name = f'delta_pdui_{cluster1}_vs_{cluster2}'
    return delta


def calculate_pdui_per_cluster(count_matrix: pd.DataFrame,
                               pas_info: pd.DataFrame,
                               cluster_labels: pd.Series,
                               pseudocount: float = 0.0) -> pd.DataFrame:
    """Calculate PDUI on pseudo-bulk (cluster-aggregated) counts.

    Aggregates per-cell counts by cluster, then computes PDUI
    for each gene in each cluster.

    Args:
        count_matrix: DataFrame (PAS x cells) with counts.
        pas_info: DataFrame with PAS coordinates and gene mapping.
        cluster_labels: Series mapping cell barcodes to cluster IDs.
        pseudocount: Pseudocount for PDUI calculation (default: 0.0).

    Returns:
        DataFrame (genes x clusters) with PDUI values.
    """
    # Aggregate counts by cluster
    common_cells = count_matrix.columns.intersection(cluster_labels.index)
    cm = count_matrix[common_cells]
    labels = cluster_labels[common_cells]

    cluster_counts = cm.T.groupby(labels).sum().T

    return calculate_pdui(cluster_counts, pas_info, pseudocount=pseudocount)
