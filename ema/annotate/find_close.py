import pybedtools
import pandas as pd

def find_close(posbed_dir:str,
             negbe_dir:str,
             genomebed_dir:str,
             annotatedbed_dir:str,
             mergebed:str,
             features:str,
             pas_col=3,
             gene_col=9
             ) -> pd.DataFrame():
    
    #TODO this 5 line will be other module 
    posbed = pybedtools.BedTool(posbed_dir)
    negbed = pybedtools.BedTool(negbe_dir)
    genomebed = pybedtools.BedTool(genomebed_dir)

    pasbed = negbed.cat(posbed, postmerge=False).sort()
    pasbed.saveas(mergebed)

    annotaded = pasbed.closet(genomebed, D="ref", t="first", s=True)
    annotaded.saveas(annotatedbed_dir)

    annotated_frame = pd.read_csv(annotatedbed_dir, delimiter="\t", header=None)
    genes_frame = annotated_frame[:,[pas_col,gene_col]]
    genes_frame.sort_values(by=0, kind="quicksort")
    genes_frame.to_csv(features, sep="\t", index=False)
    return genes_frame

if __name__ == "__main__":
    find_close()