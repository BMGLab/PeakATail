import anndata as ad
import scanpy as sc
import pandas as pd
from .config import directory_config

def clustering(adata:ad.AnnData(), outputpath=directory_config.cluster_labels,
                result_file=directory_config.cluster_output):
    """TODO
    """
    sc.pp.normalize_total(adata=adata, target_sum= 1e4)
    sc.pp.log1p(adata)
    
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable]#TODO not sure

    sc.pp.scale(adata, max_value=10)

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    #sc.pl.paga(adata, plot=False)
    sc.tl.umap(adata) #init_pos='paga'
    #sc.tl.umap(adata)
    sc.tl.louvain(adata)
    #sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
    #adata.write(result_file)
    cluster_labels = adata.obs['louvain']
    cluster_labels.sort_index()
    cluster_labels.to_csv(outputpath, header=True, index=True)


