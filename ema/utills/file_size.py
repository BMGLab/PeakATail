import pysam as ps
import os   

def get_bam_size(bam_file:str):

    return os.path.getsize(bam_file)