from .pop_apa import pop_apa
from .groupby import groupby_cb
from .fishertest import fishertest
from ema.config import directory_config

def apa_switch(cell_combine: list) -> None:
    """
    Perform APA switch analysis on the given cell combinations.

    Args:
        cell_combine (list): List of cell combinations. The list must be of size 2*N, where N is the number of pairs.
                             Each pair should be a list of two elements representing the cell indices.

    Returns:
        None
    """
    matrix = pop_apa()
    matrix = groupby_cb(matrix)

    for clusters in cell_combine:
        selected_cells = matrix[clusters]
        print(clusters)
        result_dir = directory_config.get_fisher_dir(directory_config,cluster1=clusters[0], cluster2=clusters[1])
        fishertest(selected_cells=selected_cells, result_dir=result_dir)

if __name__ == "__main__":
    apa_switch()