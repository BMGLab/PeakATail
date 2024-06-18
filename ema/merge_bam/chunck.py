from ema.utills.file_size import get_bam_size
import pysam as ps
import os
import multiprocessing as mp


def process_chuncks(file_dir:str, chunck_size:int):
    sample_id = os.path.splitext(os.path.basename(file_dir))[0]
    chuncks = calculate_chunck(file_dir, chunck_size, sample_id)

    with mp.Pool(mp.cpu_count()) as pool:
        for reads in pool.starmap(tag_chunck, chuncks):
            yield from reads


def calculate_chunck(file_dir:str, chunck_size:int, sample_id:str) -> list:
    file_size = get_bam_size(file_dir)
    chuncks = [(file_dir, sample_id,start, min(start + chunck_size, file_size)) for start in range(0, file_size, chunck_size)]
    return chuncks

def tag_chunck(file_dir:str,  sample_id:str, start:int, end:int, tag_symbol:str="RG") -> list:

    reads = []
    with ps.AlignmentFile(file_dir, "rb") as bam:
        bam.seek(start)

        for read in bam:
            if bam.tell() > end:
                break
            
            read.set_tag(tag_symbol, sample_id)
            reads.append(read)
        
    return reads