import argparse

def cli():

    parser = argparse.ArgumentParser(prog="ema")

    parser.add_argument("--bamDir", dest="bam_dir", type=str, required=True,
                        help="directory of bamfile examplebam.bam")
    print("bam")

    parser.add_argument("--sequenceLen", dest="seqlen", type=int, required=True)
    print("seq")

    parser.add_argument("--CellBarcodeLen", dest="cb_len", type=int, required=True)

    parser.add_argument("--BarcodeTag", type=str, dest="barcode_tag", required=True)
    return parser.parse_args()
    
    
if __name__ == "__main__":
    cli()