from ema.config import directory_config
import pybedtools
import pandas as pd

def find_close(posbed_dir=directory_config.posbed,
             negbe_dir=directory_config.negbed,
             genomebed_dir=directory_config.endbed,
             annotatedbed_dir=directory_config.annotatedbed,
             mergebed=directory_config.pasbed,
             pas_col=3,
             gene_col=9
             ) -> pd.DataFrame():
    
    #TODO this 5 line will be other module 
    posbed = pybedtools.BedTool(posbed_dir)
    negbed = pybedtools.BedTool(negbe_dir)
    genomebed = pybedtools.BedTool(genomebed_dir).sort()

    pasbed = negbed.cat(posbed, postmerge=False).sort()
    pasbed.saveas(mergebed)
    #find closest endpoint of  agene for a pas and annotated it for that
    annotaded = pasbed.intersect(genomebed, D="ref", t="first", s=True, loj=True)
    annotaded.saveas(annotatedbed_dir)
    
    annotated_frame = pd.read_csv(annotatedbed_dir, delimiter="\t", header=None)
    # getthat distance are more than 5000, 13 = distance column
    #condition = abs(annotated_frame.iloc[:,12]) <= 5000
    # set gene_id column to NA
    #annotated_frame.iloc[condition, 9] = "NA"
    annotated_frame.to_csv(annotatedbed_dir, sep="\t", header=False, index=False)
    # get columns that contain pas_id and gene_id
    genes_frame = annotated_frame.iloc[:,[pas_col,gene_col]]
    genes_frame = genes_frame.sort_values(by=pas_col, kind="quicksort")
    # remove that pas could not annotate to any gene and bedtools returned "."
    genes_frame = genes_frame.set_index(genes_frame.columns[0])
    mask = genes_frame != "."
    genes_frame = genes_frame[mask].dropna()

    return genes_frame

if __name__ == "__main__":
    find_close()