import pysam as ps 

def merge(bam_files:list,threads:int=4, output:str="Aligned.sortedByCoord.merged.out.bam") -> None:
    temp = "temp.bam"
    ps.merge("-r", "-p", "--threads", f"{threads}", "-o", temp, *bam_files)
    ps.sort("-o", temp, output)


if __name__ == "__main__":
    merge()
