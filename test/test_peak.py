from ema.countmatrix.peackcalling import peak_calling, directory_config
from ema.countmatrix.indexing import cb_total
peak_calling(True, bedfilepath=directory_config.negbed, matrixpath=directory_config.negmatrixpath)
peak_calling(False, bedfilepath=directory_config.posbed, matrixpath=directory_config.posmatrixpath)
sve = open("cbtotal.py", "w")
sve.write(f"cb_dic = {cb_total}")