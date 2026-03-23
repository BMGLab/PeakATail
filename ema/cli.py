import argparse

def cli():

    parser = argparse.ArgumentParser(prog="ema")

    parser.add_argument("--bamDir", dest="bam_dir", type=str, 
                        help="directory of bamfile examplebam.bam")
    print("bam")

    parser.add_argument("--sequenceLen", dest="seqlen", type=int)
    print("seq")

    parser.add_argument("--CellBarcodeLen", dest="cb_len", type=int)

    parser.add_argument("--BarcodeTag", type=str, dest="barcode_tag")

    parser.add_argument("--gtfDir", dest="gtf_dir", type=str)
    parser.add_argument("--cell_combinations", type=str, required=False)
    parser.add_argument("--bamFiles", dest="bam_files", type=str, required=False)
    parser.add_argument("--threads", dest="threads", type=int, required=False)

    # Strategy selection and parameters
    parser.add_argument("--strategy", type=str, default="original",
                        choices=["original", "lambda_poisson", "lambda_gradient", "sierra_iterative"],
                        help="Peak finding strategy (default: original)")
    parser.add_argument("--lambda-window", dest="lambda_window", type=int, default=5000,
                        help="Flanking window size in bp for lambda estimation (default: 5000)")
    parser.add_argument("--lambda-method", dest="lambda_method", type=str, default="median",
                        choices=["median", "percentile20", "global"],
                        help="Method for computing background lambda (default: median)")
    parser.add_argument("--max-pas", dest="max_pas", type=int, default=5,
                        help="Maximum PAS per peak region (default: 5)")
    parser.add_argument("--smoothing-window", dest="smoothing_window", type=int, default=50,
                        help="Smoothing window for gradient strategy in bp (default: 50)")
    parser.add_argument("--min-prominence", dest="min_prominence", type=float, default=5.0,
                        help="Minimum prominence for gradient peaks (default: 5.0)")
    parser.add_argument("--benchmark", action="store_true", default=False,
                        help="Run benchmark mode: compare all strategies")
    parser.add_argument("--validate-db", dest="validate_db", type=str, default=None,
                        help="Path to PolyASite/PolyA_DB BED for validation")

    # Post-processing filters (all optional, enable/disable)
    parser.add_argument("--internal-priming-filter", dest="ip_filter", action="store_true", default=False,
                        help="Enable internal priming filter (requires --genome-fasta)")
    parser.add_argument("--genome-fasta", dest="genome_fasta", type=str, default=None,
                        help="Path to genome FASTA for internal priming filter")
    parser.add_argument("--annotation-filter", dest="annot_filter", action="store_true", default=False,
                        help="Enable annotation region filter (keep only peaks in annotated regions)")
    parser.add_argument("--ip-a-stretch", dest="ip_a_stretch", type=int, default=6,
                        help="Minimum consecutive A's for internal priming (default: 6)")

    # Matrix filtering parameters
    parser.add_argument("--min-pas-per-cell", dest="min_pas_per_cell", type=int, default=50,
                        help="Minimum PAS per cell for filtering (default: 50)")

    # Dynamic lambda-based threshold parameters
    parser.add_argument("--dynamic-threshold", dest="dynamic_threshold", action="store_true", default=False,
                        help="Use dynamic lambda-based threshold instead of fixed (default: off)")
    parser.add_argument("--floor-threshold", dest="floor_threshold", type=int, default=3,
                        help="Minimum peak height threshold when --dynamic-threshold is active (default: 3)")
    parser.add_argument("--lambda-fold-change", dest="lambda_fold_change", type=float, default=2.0,
                        help="Fold change above background lambda for dynamic threshold (default: 2.0)")

    return parser.parse_args()
    
    
if __name__ == "__main__":
    cli()