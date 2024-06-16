from ema.config import directory_config
from ema.merge_bam.chunck import process_chuncks
from ema.utills.csvpars import csv_to_list
import pysam as ps
import os

def merge_bam(bam_files:list, temp_merge:str="temp.bam", output_dir:str="emamerged.Aligned.sortedByCoord.out.bam", chunck_size:int=10**6) -> None:
    "TODO ADD args form config"

    if not bam_files:
        print("No BAM files found in the file.")
        return

    with ps.AlignmentFile(temp_merge, "wb", header=ps.AlignmentFile(bam_files[0]).header) as temp_bam:
        for bam_file in bam_files:
            for read in process_chuncks(bam_file, chunck_size):
                temp_bam.write(read)

    ps.sort("-o", output_dir, temp_merge)

    print(f"merge and sort {len(bam_files)} BAM files in {csv_to_list} to {output_dir}")


if __name__ == "__main__":
    merge_bam()