import argparse

def cli():

    parser = argparse.ArgumentParser(prog="ema")

    parser.add_argument("--bamDir", dest="bam_dir", type=str, 
                        help="directory of bamfile examplebam.bam")
    print("bam")

    parser.add_argument("--sequenceLen", dest="seqlen", type=int)
    print("seq")

    parser.add_argument("--CellBarcodeLen", dest="cb_len", type=int)

    parser.add_argument("--gap", dest="pas_gap", type=int)
    parser.add_argument("--min_read", default = 2000, dest="min_read", type=int)
    parser.add_argument("--min_cells", default =3, dest="min_cells", type=int)
    parser.add_argument("--min_genes", default = 200, dest="min_genes", type=int)

    parser.add_argument("--BarcodeTag", type=str, dest="barcode_tag")

    parser.add_argument("--gtfDir", dest="gtf_dir", type=str)
    parser.add_argument("--cell_combinations", type=str, required=False)
    parser.add_argument("--bamFiles", dest="bam_files", type=str, required=False)
    parser.add_argument("--threads", dest="threads", type=int, required=False)

    return parser.parse_args()
    
    
if __name__ == "__main__":
    cli()