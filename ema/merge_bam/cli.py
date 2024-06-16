import argparse
from ema.merge_bam.merge import merge_bam

def cli():
    parser = argparse.ArgumentParser(prog="ema-merge")

    parser.add_argument("--BAMcsv", dest="csv_directory", type=str, required=True,
                        help="directory of csv file containing paths to BAM files")
    
    args = parser.parse_args()
    merge_bam(bam_files=args.csv_directory)
