import time 
bed = open(r"testBed.bed")
output_bed = open(r"output_bed.bed" , "w")
ttime = time.time()
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

for r in bed:
    rr = r.split()
    chro , start , end , id , symbol , strand = rr[0] ,int(rr[1]) , int(rr[2]) , rr[3] , rr[4] , rr[5]
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

print(time.time() - ttime)

            
