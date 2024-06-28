import pysam as ps 

def merge(bam_files:list,threads:int=4, output:str="Aligned.sortedByCoord.merged.out.bam") -> None:

    ps.merge("-r", "-p", "--threads", f"{threads}", "-o", output, *bam_files)


if __name__ == "__main__":
    merge()
