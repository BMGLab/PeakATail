from ema.countmatrix.peackcalling import peak_calling
from ema.countmatrix.indexing import cb_total
peak_calling(True, bedfilepath="negbed.bed", matrixpath="negmatrix.mtx")
peak_calling(False, bedfilepath="posbed.bed", matrixpath="posmatrix.mtx")
sve = open("cbtotal.py", "w")
sve.write(f"cb_dic = {cb_total}")