import argparse
from ema.config import directory_config, variable_config, filter_config
from ema.main import main, peak_calling, filter_cb

def cli():

    parser = argparse.ArgumentParser(prog="ema", 
                                     discription="Peak calling with base on pas points")

    parser.add_argument("--bamDir", dest="bam_dir", type=str, required=True,
                        help="directory of bamfile examplebam.bam")

    parser.add_argument("--sequenceLen", dest="seqlen", type=int, required=True)

    parser.add_argument("--CellBarcodeLen", dest="cb_len", type=int, required=True)

    parser.add_argument("--BarcodeTag", type=str, dest="barcode_tag", required=True)
    args = parser.parse_args()
    directory_config.bam_dir = args.bam_dir
    variable_config.seqlen = args.seqlen
    variable_config.cb_len = args.cb_len
    variable_config.barcode_tag = args.barcode_tag
    main()