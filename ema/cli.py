import argparse

def cli():

    parser = argparse.ArgumentParser(prog="ema")
    
    parser.add_argument("--bamFile", dest="bam_dir", type=str, required=True,
                        help="a sorted bam file (an output of CellRanger or Star Solo)")
    
    print("bam")

    parser.add_argument("--sequenceLen", dest="seqlen", type=int, required=True,
                        help="the sequencing read length (in *_2.fastq)")
    
    print("seq")

    parser.add_argument("--CellBarcodeLen", dest="cb_len", type=int, required=True)

    parser.add_argument("--gap", dest="pas_gap", type=int)
    parser.add_argument("--min_read", default = 2000, dest="min_read", type=int)
    parser.add_argument("--min_cells", default =3, dest="min_cells", type=int)
    parser.add_argument("--min_genes", default = 200, dest="min_genes", type=int)

    parser.add_argument("--BarcodeTag", type=str, dest="barcode_tag", required=True)

    parser.add_argument("--gtfFile", dest="gtf_dir", type=str, required=True)
    parser.add_argument("--cell_combinations", type=str, required=False)
    parser.add_argument("--threads", dest="threads", type=int, required=False)

    return parser.parse_args()
    
    
if __name__ == "__main__":
    cli()