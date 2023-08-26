from ema.matrixfilter import filter_cb, make_dataframe, preprocessing, filtered_cb_list
from ema.clustering import clustering
filter_cb()

print("cb ok")
matrixfram = make_dataframe()
adata = preprocessing(df=matrixfram)
clustering(adata=adata)
