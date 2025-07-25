from multiprocessing import Pool
from scipy.stats import fisher_exact, false_discovery_control
import pandas as pd
import numpy as np


def benjamini_hochberg_test(result_file: str, out_file: str):
    
    df = pd.read_csv(result_file, sep="\t", header=None)
    orig_p = df.iloc[:, 4].to_numpy()

    adj_p = false_discovery_control(orig_p, method="bh")

    df.loc[:, 5] = adj_p

    df.to_csv(out_file, sep="\t", header=False, index=False)

     
def fishertest(selected_cells:pd.DataFrame, result_dir:str, level:str='Ensemble_ID')->None:
    """
    Perform Fisher's exact test on selected cells.

    Parameters:
    selected_cells (pd.DataFrame): DataFrame containing the selected cells.
    result_dir (str): Directory to save the results.
    level (str, optional): Level to group the selected cells by. Defaults to 'Ensemble_ID'.

    Returns:
    None
    """
    cluster1, cluster2 = selected_cells.columns
    with Pool() as p:
        all_results = p.starmap(process_gene, 
                        [(gene, group, cluster1, cluster2) for gene, group in selected_cells.groupby(level= level)])
        
    with open(result_dir, 'w') as f:
        f.write(''.join(sum(all_results, [])))


def process_gene(gene:str, group:pd.DataFrame, cluster1:str, cluster2:str)->list:
    """
    Process a gene and calculate Fisher's exact test p-value for each unique index in the group DataFrame.

    Parameters:
    gene (str): The name of the gene being processed.
    group (pd.DataFrame): The DataFrame containing the data for the gene.
    cluster1 (str): The name of the first cluster.
    cluster2 (str): The name of the second cluster.

    Returns:
    list: A list of strings representing the results of the Fisher's exact test for each unique index in the group DataFrame.
    """
    results = []

    for pas in group.index.unique():
        table = np.array([group.loc[pas], group.drop(index=pas).sum()])
        res = fisher_exact(table, alternative='two-sided')
        results.append(f"{gene}\t{pas[1]}\t{cluster1}\t{cluster2}\t{res[1]}\n")
    return results

if __name__ == "__main__":
    fishertest()

    #benjamini_hochberg_test("../../../fisherresults_3_4.tsv","../../../fisherresults_adjusted.tsv")
    #benjamini_hochberg_test("./test_results.txt","./test_results_adjusted.txt")
