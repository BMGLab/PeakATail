
import pandas as pd

import random
def ratio_score(bed_file):#TODO change the parameter from bed to dataframe
    df = pd.read_csv(bed_file, sep='\t', header=None)

    grouped = df.groupby(5)#TODO add .filter
    for key, group in grouped:
        if len(group) > 1:
            continue
        print(f"Group '{key}':")
        print(group)
        group_c = (((group[1] + group[2])/2) - group[8])
        print(group_c)
        
        break




#To erase:start
def add_random_coverage_to_example_data(bed_file_input , bed_file_output):
    with open(bed_file_input, 'r') as input_file , open(bed_file_output, "w") as output_file:
        for row in input_file:
            written_row = row.rstrip("\n")
            output_file.write(f"{written_row}\t{random.randint(20,100)}\n")
        
def add_random_gene_range_to_example_data(bed_file_input , bed_file_output):
    with open(bed_file_output , "w") as file:
        df = pd.read_csv(bed_file_input, sep='\t', header=None)

        grouped = df.groupby(5)#TODO add .filter
        for key, group in grouped:
            
            start = min(group[1])-100
            end = max(group[2]) + 100
            
            for x in range(len(group)):
                row_as_string = '\t'.join(group.iloc[x].astype(str))
                file.write(f"{row_as_string}\t{start}\t{end}\n")
#to erase:end
      
        

ratio_score(r"example_bed_2.bed")
#add_random_coverage_to_example_data1(r"example_bed.bed" , r"example_bed_2.bed")




