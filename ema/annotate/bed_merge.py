import pybedtools
from ema.config import directory_config as dc
def bed_merge(bed_file= dc.endbed, 
             chro_index=0,
             start_index=1,
             end_index=2,
             strand_index=5,
             transcript_id = 3,
             symbol_index=4,
             ) -> None:
   
        
    bed = pybedtools.BedTool(bed)
    
    sorted_bed = bed.sort()
    
    # Merge overlapping intervals and retain IDs and all columns after the ID column
    retained_index_tuple = tuple( map( lambda n:n+1 , [transcript_id , symbol_index , strand_index]) )
    merged_bed = sorted_bed.merge(c=retained_index_tuple, o="collapse")
    
    
    merged_bed.saveas(bed_file)


    

if __name__ == "__main__":
    bed_merge()