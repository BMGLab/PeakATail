"""Fisher's exact test for differential PAS usage between clusters.

Tests whether the proportion of reads at each PAS (relative to all PAS
in the same gene) differs significantly between two clusters.

Includes FDR correction (Benjamini-Hochberg) across all tests.
"""

from scipy.stats import fisher_exact, false_discovery_control
import pandas as pd
import numpy as np


def fishertest(selected_cells: pd.DataFrame,
               result_dir: str,
               level: str = 'Ensemble_ID',
               fdr_threshold: float = 0.05) -> pd.DataFrame:
    """Perform Fisher's exact test on selected cells with FDR correction.

    For each gene with multiple PAS, tests whether the proportion of
    reads at each PAS differs between two clusters.

    Args:
        selected_cells: DataFrame with multi-index (gene, PAS) and
            two columns (one per cluster) of aggregated counts.
        result_dir: Path to save the results TSV file.
        level: Index level to group by (default: 'Ensemble_ID').
        fdr_threshold: FDR q-value threshold for significance
            (default: 0.05). Used only for reporting, all results
            are saved.

    Returns:
        DataFrame with all test results including q-values.
    """
    cluster1, cluster2 = selected_cells.columns

    # Collect all results in memory first (for FDR correction)
    all_results = []
    for gene, group in selected_cells.groupby(level=level):
        gene_results = _process_gene(gene, group, cluster1, cluster2)
        all_results.extend(gene_results)

    if len(all_results) == 0:
        # No testable genes
        results_df = pd.DataFrame(
            columns=['gene', 'pas', 'cluster1', 'cluster2',
                     'p_value', 'q_value', 'odds_ratio', 'delta_proportion']
        )
        results_df.to_csv(result_dir, sep='\t', index=False)
        return results_df

    # Build results DataFrame
    results_df = pd.DataFrame(all_results, columns=[
        'gene', 'pas', 'cluster1', 'cluster2',
        'p_value', 'odds_ratio', 'delta_proportion',
    ])

    # Apply FDR correction (Benjamini-Hochberg)
    pvalues = results_df['p_value'].values
    qvalues = false_discovery_control(pvalues, method='bh')
    results_df['q_value'] = qvalues

    # Reorder columns
    results_df = results_df[[
        'gene', 'pas', 'cluster1', 'cluster2',
        'p_value', 'q_value', 'odds_ratio', 'delta_proportion',
    ]]

    # Sort by q-value
    results_df = results_df.sort_values('q_value')

    # Save to file
    results_df.to_csv(result_dir, sep='\t', index=False)

    # Report summary
    n_significant = (results_df['q_value'] < fdr_threshold).sum()
    n_total = len(results_df)
    n_genes = results_df['gene'].nunique()
    print(
        f"Fisher test: {n_significant}/{n_total} significant PAS "
        f"(q < {fdr_threshold}) across {n_genes} genes "
        f"[{cluster1} vs {cluster2}]"
    )

    return results_df


def _process_gene(gene: str, group: pd.DataFrame,
                  cluster1: str, cluster2: str) -> list:
    """Process a single gene: test each PAS against the rest.

    For each PAS in the gene, constructs a 2x2 contingency table:
        [[PAS_count_c1, PAS_count_c2],
         [other_PAS_count_c1, other_PAS_count_c2]]

    Args:
        gene: Gene identifier.
        group: DataFrame with counts for all PAS in this gene.
        cluster1: Name of first cluster.
        cluster2: Name of second cluster.

    Returns:
        List of tuples: (gene, pas, cluster1, cluster2, pvalue,
        odds_ratio, delta_proportion).
    """
    results = []

    for pas in group.index.unique():
        pas_counts = group.loc[pas]
        other_counts = group.drop(index=pas).sum()

        # Ensure we have scalar values (handle single-row case)
        if isinstance(pas_counts, pd.DataFrame):
            pas_counts = pas_counts.iloc[0]

        table = np.array([
            [pas_counts[cluster1], pas_counts[cluster2]],
            [other_counts[cluster1], other_counts[cluster2]],
        ])

        # Skip if table has zero margins (no counts in a cluster)
        if table.sum() == 0:
            continue
        if table[:, 0].sum() == 0 or table[:, 1].sum() == 0:
            continue

        odds_ratio, pvalue = fisher_exact(table, alternative='two-sided')

        # Compute delta proportion (proportion of PAS in c1 - proportion in c2)
        total_c1 = table[:, 0].sum()
        total_c2 = table[:, 1].sum()
        prop_c1 = table[0, 0] / total_c1 if total_c1 > 0 else 0.0
        prop_c2 = table[0, 1] / total_c2 if total_c2 > 0 else 0.0
        delta_prop = prop_c1 - prop_c2

        # Extract the PAS identifier from the multi-index
        pas_id = pas[1] if isinstance(pas, tuple) else pas

        results.append((
            gene, pas_id, cluster1, cluster2,
            pvalue, odds_ratio, delta_prop,
        ))

    return results
