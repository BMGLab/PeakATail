import argparse
from ema.utills.csvpars import csv_to_list
from ema.merge_bam.merge import merge_bam

def cli():
    parser = argparse.ArgumentParser(prog="ema_merge")

    parser.add_argument("--BAMcsv", dest="csv_directory", type=str, required=True,
                        help="directory of csv file containing paths to BAM files")
    
    args = parser.parse_args()
    bam_list = csv_to_list(args.csv_directory)
    print(bam_list[0])
    merge_bam(bam_files=bam_list[0])
