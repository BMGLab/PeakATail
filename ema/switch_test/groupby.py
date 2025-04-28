import pandas as pd
from ema.config import directory_config
import pandas as pd

def groupby_cb(count_matrix: pd.DataFrame, 
               cluster_labels_dir: str = directory_config.cluster_labels
               ) -> pd.DataFrame:
    """
    Group the count matrix by cluster labels.

    Args:
        count_matrix (pd.DataFrame): The count matrix to be grouped.
        cluster_labels_dir (str, optional): The directory of the cluster labels file. 
                                            Defaults to directory_config.cluster_labels.

    Returns:
        pd.DataFrame: The grouped count matrix.
    """
    cluster_labels = pd.read_csv(cluster_labels_dir, sep=",", header=None)
    count_matrix = count_matrix.transpose()
    cell_intersect = count_matrix.index.intersection(cluster_labels.iloc[:, 0])
    cluster_labels = cluster_labels[cluster_labels.iloc[:, 0].isin(cell_intersect)]

    multi_index = pd.MultiIndex.from_frame(cluster_labels)

    order = [multi_index.levels[1].name, multi_index.levels[0].name]

    count_matrix.index = multi_index.reorder_levels(order=order)
    count_matrix = count_matrix.transpose()

    count_matrix.columns = pd.MultiIndex.from_frame(count_matrix.columns.to_frame().fillna(''))
    count_matrix = count_matrix.groupby(level=1, axis=1).sum()
    return count_matrix