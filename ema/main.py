from .countmatrix.peackcalling import peak_calling
from .matrixfilter import default_filtiration
import multiprocessing

def main():
    peak_calling(True)
    peak_calling(False)
    default_filtiration()





