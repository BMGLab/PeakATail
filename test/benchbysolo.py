import pandas as df
import anndata as ad
import scipy.io as sci
from ema.matrixfilter import filter_cb, make_dataframe, preprocessing, filtered_cb_list
from ema.countmatrix.indexing import cb_total
from ema.config import directory_config, variable_config, filter_config
from ema.clustering.clustering import clustering
from ema.countmatrix.peackcalling import peak_calling

""" cb =open("/home/user/D/BAMdata/proje/ProjectEMA/ema/cbtotal.py", "w")
peak_calling(True, bedfilepath=directory_config.negbed, matrixpath=directory_config.negmatrixpath, bamfile_dir="/home/user/D/BAMdata/proje/ProjectEMA/test/testdata/srr8325981/Aligned.sortedByCoord.out.bam")
peak_calling(False, bedfilepath=directory_config.posbed, matrixpath=directory_config.posmatrixpath, bamfile_dir="/home/user/D/BAMdata/proje/ProjectEMA/test/testdata/srr8325981/Aligned.sortedByCoord.out.bam")


cb.write(f"cb_dic={cb_total}") """


def intersolo(solomatrixpath="/home/user/D/BAMdata/proje/ProjectEMA/test/testdata/srr8325981/filtered/matrix.mtx"
              ,solo_cbpath="/home/user/D/BAMdata/proje/ProjectEMA/test/testdata/srr8325981/filtered/barcodes.tsv"):
    solo_list = open(solo_cbpath, "r").readlines()
    solo_list = [line.strip() for line in solo_list]
    solo_df = make_dataframe(matrixpath=solomatrixpath, collist=solo_list)
    solo_data = preprocessing(df=solo_df, min_cells=3, min_genes=200)
    filter_cb(negativematrixpath=  "/home/user/D/BAMdata/proje/ProjectEMA/emaout1/negmatrix.mtx",
              positivematrixpath="/home/user/D/BAMdata/proje/ProjectEMA/emaout1/posmatrix.mtx",
              min_read=2000)
    
    ema_df = make_dataframe(collist=filtered_cb_list)
    ema_data = preprocessing(df=ema_df, min_cells=3, min_genes=200)


    common_col = solo_data.columns.intersection(ema_data.columns)
    print("thisssssssssss" + common_col)
    intersectema = ema_data.loc[:,common_col].transpose()
    print(intersectema)
    intersectemaadat = ad.AnnData(X=intersectema)
    intersectsolo = solo_data.loc[:,common_col].transpose()
    intersectsoloadata = ad.AnnData(X=intersectsolo)
    #print(intersectema)
    #print(intersectsolo)
    clustering(adata=intersectemaadat, 
               outputpath="/home/user/D/BAMdata/proje/ProjectEMA/test/testdata/srr8325981/emalabels.csv")
    
    clustering(adata=intersectsoloadata, 
               outputpath= "/home/user/D/BAMdata/proje/ProjectEMA/test/testdata/srr8325981/sololabels.csv")
    
if __name__ == "__main__":
    intersolo()


    



