import argparse
import yaml
from ema.merge_bam.merge import merge

def cli():

    parser = argparse.ArgumentParser(prog="ema_merge")

    parser.add_argument("--bamFiles", dest="bam_files", type=str)
    parser.add_argument("--threads", dest="threads", type=int, required=False)
    
    args = parser.parse_args()

    with open(args.bam_files, "r") as yml:
        bam_list = yaml.safe_load(yml)

    if not args.threads:
        merge(bam_files=bam_list)
    
    else:
        merge(bam_files=bam_list, threads=args.threads)
