
from ema.config import directory_config as dc

def bed_merge(bed_file= dc.endbed, 
              output_bed_file = dc.utr_bed,
             chro_index=0,
             start_index=1,
             end_index=2,
             strand_index=5,
             transcript_id = 3,
             symbol_index=4,
             ) -> None:
    bed = open(bed_file)
    output_bed = open(output_bed_file , "w")
    prev_chro = "0"
    string_start = 0
    string_end =0
    id_list = []
    symbol_list = []
    strand_list = []
    def listToString(myList):
        resultStr = ""
        for x in myList:
            resultStr += x + ","
        return resultStr[:-1]

    for record in bed:
        split_record = record.split()
        chro , start , end , id , symbol , strand = split_record[chro_index] ,int(split_record[start_index]) , int(split_record[end_index]) , split_record[transcript_id] , split_record[symbol_index] , split_record[strand_index]
        if string_start == 0:
            string_start = start 
        if string_end == 0:
            string_end = end 
        if chro == prev_chro:
            id_list.append(id)
            symbol_list.append(symbol)
            strand_list.append(strand)
            if start < string_end:
                if end > string_end:
                    string_end = end
                
            else:

                output_bed.write(f"{chro}\t{string_start}\t{string_end}\t{listToString(id_list)}\t{listToString(symbol_list)}\t{listToString(strand_list)}\n")
                string_start = start
                string_end = end
                id_list = []
                symbol_list = []
                strand_list = []
            
        else:
            output_bed.write(f"{chro}\t{string_start}\t{string_end}\t{listToString(id_list)}\t{listToString(symbol_list)}\t{listToString(strand_list)}\n")
            string_start = 0
            string_end =0
            id_list = []
            symbol_list = []
            strand_list = []
        prev_chro = chro

    
                
