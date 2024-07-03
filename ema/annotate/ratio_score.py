
import pandas as pd
from multiprocessing import Pool

def ratio_score(bed_file , result_dir):#TODO change the parameter from bed to dataframe
    df = pd.read_csv(bed_file, sep='\t', header=None)#TODO name the columns based on correct file example
    grouped = df.groupby(5)#TODO add .filter()
    with Pool() as p:
        all_results = p.starmap(process_pas_group, 
                        [(group,) for key,group in grouped])
        
    with open(result_dir, 'w') as f:
        f.write(''.join(sum(all_results, [])))

def process_pas_group(group):
    if len(group) < 1:
        return
    result = []
    group_distance = (((group[1] + group[2])/2) - group[8])
    group_c = group_distance.sum()
    group_ratio_scores = group_distance / group_c
    group['group_ratio_scores'] = group_ratio_scores
    for x in range(len(group)):
        row_as_string = '\t'.join(group.iloc[x].astype(str))
        result.append(row_as_string+"\n") 
    return result

if __name__ == "__main__":
    ratio_score(r"example_bed_2.bed" , r"result.bed")




