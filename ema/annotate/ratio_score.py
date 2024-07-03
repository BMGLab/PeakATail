
import pandas as pd
from multiprocessing import Pool

def ratio_score(df , result_dir):

    df.columns= ["chromosomes" , "PAS_start" , "PAS_end" , " " , "ID" , "name" , "strand" , "coverage" , "gene_start", "gene_end"]#TODO name the columns based on correct file example
    grouped = df.groupby("ID")#TODO add .filter() 
    with Pool() as p:
        all_results = p.starmap(process_pas_group, 
                        [(group,) for key,group in grouped])
        
    with open(result_dir, 'w') as f:
        f.write(''.join(sum(all_results, [])))

def process_pas_group(group):
    if len(group) < 1:
        return
    result = []
    group_distance = (((group["PAS_start"] + group["PAS_end"])/2) - group["gene_start"])
    group_c = group_distance.sum()
    group_ratio_scores = group_distance / group_c
    group['group_ratio_scores'] = group_ratio_scores
    for x in range(len(group)):
        row_as_string = '\t'.join(group.iloc[x].astype(str))
        result.append(row_as_string+"\n") 
    return result

def bed_to_df(bed_file):
    df = pd.read_csv(bed_file, sep='\t', header=None)
    return df

if __name__ == "__main__":
    ratio_score()




