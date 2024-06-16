import pysam as ps
import os   

def get_bam_size(bam_file:str):

    with ps.AlignmentFile(bam_file, 'rb') as bam:
        bam.seek(0, os.SEEK_END)
        return bam.tell()