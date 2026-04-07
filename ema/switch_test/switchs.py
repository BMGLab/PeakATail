"""APA switch test orchestrator.

Coordinates the differential APA analysis pipeline:
1. Load filtered count matrix and group by cluster labels
2. For each cluster pair, run Fisher's exact test with FDR correction
3. Optionally compute PDUI for each cluster pair
"""

import os
import pandas as pd
from .pop_apa import pop_apa
from .groupby import groupby_cb
from .fishertest import fishertest
from ema.quantification.pdui import (
    load_pas_data, calculate_pdui_per_cluster, calculate_delta_pdui,
)


def apa_switch(cell_combine: list,
               fdr_threshold: float = 0.05,
               compute_pdui: bool = True,
               output_dir: str = None) -> dict:
    """Perform APA switch analysis on the given cell combinations.

    Args:
        cell_combine: List of cluster pairs to compare. Each element
            is a list of two cluster IDs, e.g. [['0', '1'], ['0', '2']].
        fdr_threshold: FDR threshold for significance (default: 0.05).
        compute_pdui: Whether to compute PDUI (default: True).
        output_dir: Output directory for results. If None, uses
            directory_config.switch_dir.

    Returns:
        Dictionary mapping cluster pair tuples to result DataFrames.
    """
    from ema.config import directory_config

    if output_dir is None:
        output_dir = directory_config.switch_dir
    os.makedirs(output_dir, exist_ok=True)

    # Load and aggregate count matrix
    matrix = pop_apa()
    grouped_matrix = groupby_cb(matrix)

    # Optionally compute PDUI
    pdui_matrix = None
    if compute_pdui:
        try:
            count_matrix, pas_info = load_pas_data(
                matrix_path=directory_config.filterd_matrix,
                cellbarcodes_path=directory_config.filtered_cb,
                pas_gene_path=directory_config.pas_geneid,
                pasbed_path=directory_config.pasbed,
            )
            cluster_labels = pd.read_csv(
                directory_config.cluster_labels, header=None,
            )
            cluster_series = pd.Series(
                cluster_labels.iloc[:, 1].values,
                index=cluster_labels.iloc[:, 0].values,
            ).astype(str)

            pdui_matrix = calculate_pdui_per_cluster(
                count_matrix, pas_info, cluster_series,
            )

            # Save PDUI matrix
            pdui_path = os.path.join(output_dir, 'pdui_matrix.tsv')
            pdui_matrix.to_csv(pdui_path, sep='\t')
            print(f"PDUI matrix saved to: {pdui_path}")
            print(
                f"PDUI computed for {pdui_matrix.shape[0]} genes "
                f"across {pdui_matrix.shape[1]} clusters"
            )
        except Exception as e:
            print(f"Warning: PDUI calculation failed: {e}")
            print("Continuing with Fisher test only.")
            pdui_matrix = None

    # Run Fisher test for each cluster pair
    all_results = {}
    # Ensure column types match (groupby produces int columns, pairs may be strings)
    available_cols = set(grouped_matrix.columns)
    for clusters in cell_combine:
        cluster1, cluster2 = clusters[0], clusters[1]
        # Try both string and int versions of cluster IDs
        cols = []
        for c in [cluster1, cluster2]:
            if c in available_cols:
                cols.append(c)
            elif int(c) in available_cols:
                cols.append(int(c))
            elif str(c) in available_cols:
                cols.append(str(c))
            else:
                print(f"Warning: cluster '{c}' not found in grouped matrix columns: {sorted(available_cols)}")
                cols.append(c)
        selected_cells = grouped_matrix[cols]

        result_path = os.path.join(
            output_dir,
            f"fisherresults_{cluster1}_{cluster2}.tsv",
        )

        result_df = fishertest(
            selected_cells=selected_cells,
            result_dir=result_path,
            fdr_threshold=fdr_threshold,
        )

        all_results[(cluster1, cluster2)] = result_df

        # Compute delta-PDUI for this pair if available
        if pdui_matrix is not None:
            c1_str, c2_str = str(cluster1), str(cluster2)
            if c1_str in pdui_matrix.columns and c2_str in pdui_matrix.columns:
                delta_pdui = calculate_delta_pdui(pdui_matrix, c1_str, c2_str)
                delta_path = os.path.join(
                    output_dir,
                    f"delta_pdui_{cluster1}_{cluster2}.tsv",
                )
                delta_pdui.to_csv(delta_path, sep='\t', header=True)
                print(
                    f"Delta-PDUI saved: {delta_path} "
                    f"({(delta_pdui.abs() > 0.1).sum()} genes with "
                    f"|delta-PDUI| > 0.1)"
                )

    return all_results


if __name__ == "__main__":
    apa_switch()
