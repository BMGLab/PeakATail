import pandas as df
import anndata as ad
import scipy.io as sci
from ema.matrixfilter import filter_cb, make_dataframe, preprocessing, filtered_cb_list
from ema.countmatrix.indexing import cb_total
from ema.config import directory_config, variable_config, filter_config
from ema.clustering import clustering



#peak_calling(True, bedfilepath=directory_config.negbed, matrixpath=directory_config.negmatrixpath)
#peak_calling(False, bedfilepath=directory_config.posbed, matrixpath=directory_config.posmatrixpath)






def intersolo(solomatrixpath="/home/user/D/BAMdata/proje/ProjectEMA/test/testdata/Solomatrix/matrix.mtx"
              ,solo_cbpath="/home/user/D/BAMdata/proje/ProjectEMA/test/testdata/Solomatrix/barcodes.tsv"):
    solo_list = open(solo_cbpath, "r").readlines()
    solo_df = make_dataframe(matrixpath=solomatrixpath, collist=solo_list)
    solo_data = preprocessing(df=solo_df, min_cells=3, min_genes=200)
    filter_cb(negativematrixpath=  "/home/user/elyar/emaout/negmatrix.mtx",
              positivematrixpath="/home/user/elyar/emaout/posmatrix.mtx",
              min_read=2000)
    
    ema_df = make_dataframe(collist=filtered_cb_list)
    ema_data = preprocessing(df=ema_df, min_cells=3, min_genes=200)

    common_col = solo_data.columns.intersection(ema_data.columns)
    intersectema = ema_data.loc[:,common_col].transpose()
    intersectemaadat = ad.AnnData(X=intersectema)
    intersectsolo = solo_data.loc[:,common_col].transpose()
    intersectsoloadata = ad.AnnData(X=intersectsolo)
    clustering(adata=intersectemaadat, 
               outputpath="/home/user/D/BAMdata/proje/ProjectEMA/test/testdata/emalabels.csv")
    
    clustering(adata=intersectsoloadata, 
               outputpath= "/home/user/D/BAMdata/proje/ProjectEMA/test/testdata/sololabels.csv")
    
if __name__ == "__main__":
    intersolo()



    



