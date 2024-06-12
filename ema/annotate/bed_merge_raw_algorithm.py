
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
    string_start , string_end = 0 , 0
    id_list , symbol_list , strand_list= [],[],[]

    def list_to_string(myList):
        resultStr = ""
        for x in myList:
            resultStr += x + ","
        return resultStr[:-1]
    
    def write_to_output(bed,record,string_start ,string_end  , id_list , symbol_list , strand_list , start = 0 , end = 0):
        bed.write(record)
        string_start = start
        string_end = end
        id_list , symbol_list , strand_list= [],[],[]

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
                output_record = f"{chro}\t{string_start}\t{string_end}\t{list_to_string(id_list)}\t{list_to_string(symbol_list)}\t{list_to_string(strand_list)}\n"
                write_to_output(output_bed,output_record,string_start ,string_end  , id_list , symbol_list , strand_list , start, end)
                
            
        else:
            output_record = f"{chro}\t{string_start}\t{string_end}\t{list_to_string(id_list)}\t{list_to_string(symbol_list)}\t{list_to_string(strand_list)}\n"
            write_to_output(output_bed,output_record,string_start ,string_end  , id_list , symbol_list , strand_list)
                
        prev_chro = chro

    
                
