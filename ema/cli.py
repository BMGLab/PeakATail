import argparse

def cli():

    parser = argparse.ArgumentParser(prog="ema")

    parser.add_argument("--config", dest="config", type=str, default=None,
                        help="Path to YAML configuration file")

    parser.add_argument("--bamDir", dest="bam_dir", type=str,
                        help="directory of bamfile examplebam.bam")

    parser.add_argument("--sequenceLen", dest="seqlen", type=int)

    parser.add_argument("--CellBarcodeLen", dest="cb_len", type=int)

    parser.add_argument("--BarcodeTag", type=str, dest="barcode_tag")

    parser.add_argument("--gtfDir", dest="gtf_dir", type=str)
    parser.add_argument("--cell_combinations", type=str, required=False)
    parser.add_argument("--bamFiles", dest="bam_files", type=str, required=False)
    parser.add_argument("--threads", dest="threads", type=int, required=False)
    parser.add_argument("--bam-threads", dest="bam_threads", type=int, default=4,
                        help="Threads for pysam BGZF block decompression (default: 4)")

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

    # Gene annotation / find_close parameters
    parser.add_argument("--max-gene-distance", dest="max_gene_distance", type=int, default=5000,
                        help="Maximum distance (bp) for PAS-to-gene assignment (default: 5000)")
    parser.add_argument("--utr-multiplier", dest="utr_multiplier", type=float, default=2.0,
                        help="Multiplier for 3'UTR length to define TIER_2 boundary (default: 2.0)")
    parser.add_argument("--include-extended", dest="include_extended", action="store_true", default=False,
                        help="Include TIER_3 PAS (within max-gene-distance but beyond UTR x multiplier)")

    # Dynamic lambda-based threshold parameters
    parser.add_argument("--dynamic-threshold", dest="dynamic_threshold", action="store_true", default=False,
                        help="Use dynamic lambda-based threshold instead of fixed (default: off)")
    parser.add_argument("--floor-threshold", dest="floor_threshold", type=int, default=3,
                        help="Minimum peak height threshold when --dynamic-threshold is active (default: 3)")
    parser.add_argument("--lambda-fold-change", dest="lambda_fold_change", type=float, default=2.0,
                        help="Fold change above background lambda for dynamic threshold (default: 2.0)")

    # Atlas merge strategy parameters
    parser.add_argument("--atlas", type=str, default=None,
                        help="Path to reference PAS atlas BED for atlas merge strategy")
    parser.add_argument("--atlas-distance", dest="atlas_distance", type=int, default=50,
                        help="Max snap distance in bp for atlas strategy (default: 50)")

    # Clustering parameters
    parser.add_argument("--clustering-method", dest="clustering_method", type=str, default="leiden_tfidf",
                        choices=["leiden_tfidf", "leiden_libsize", "external"],
                        help="Clustering strategy (default: leiden_tfidf)")
    parser.add_argument("--resolution", type=float, default=1.0,
                        help="Leiden clustering resolution (default: 1.0)")
    parser.add_argument("--n-pcs", dest="n_pcs", type=int, default=40,
                        help="Number of dimensions for neighbor computation (default: 40)")
    parser.add_argument("--external-clusters", dest="external_clusters", type=str, default=None,
                        help="Path to pre-computed cluster labels CSV (for --clustering-method external)")
    parser.add_argument("--random-seed", dest="random_seed", type=int, default=42,
                        help="Random seed for reproducibility (default: 42)")

    # PDUI strategy parameters (multi-PAS isoform-aware)
    parser.add_argument("--pdui-method", dest="pdui_method", type=str, default="classic",
                        help="Comma-separated PDUI methods: classic, proportion, shannon (default: classic)")
    parser.add_argument("--pdui-isoform-agg", dest="pdui_isoform_agg", type=str,
                        default="per_gene", choices=["per_gene", "per_isoform"],
                        help="PDUI aggregation level (default: per_gene)")
    parser.add_argument("--pdui-isoform-collapse", dest="pdui_isoform_collapse", type=str,
                        default="none", choices=["none", "mean", "majority"],
                        help="When per_isoform: how to collapse isoforms for gene-level summary (default: none)")

    # Differential APA strategy
    parser.add_argument("--diff-method", dest="diff_method", type=str, default="fisher",
                        choices=["fisher", "nb_pairwise", "nb_multi"],
                        help="Differential APA test method (default: fisher)")

    # Cross-dataset cluster matching strategy
    parser.add_argument("--cluster-match-method", dest="cluster_match_method", type=str,
                        default="marker_overlap",
                        choices=["marker_overlap", "mnn", "jaccard"],
                        help="Cross-dataset cluster matching strategy (default: marker_overlap)")
    parser.add_argument("--n-top-markers", dest="n_top_markers", type=int, default=50,
                        help="Top-N marker PAS for cluster matching (default: 50)")

    args = parser.parse_args()

    # Load YAML config if provided
    if args.config:
        import yaml
        with open(args.config) as f:
            cfg = yaml.safe_load(f)
        args.datasets = cfg.get('datasets', [])
        # Override CLI defaults with YAML values if present
        for key in ['seqlen', 'cb_len', 'barcode_tag', 'min_read', 'min_cells', 'min_pas_per_cell',
                    'pas_gap', 'gtf_dir', 'atlas', 'atlas_distance',
                    'pdui_method', 'pdui_isoform_agg', 'pdui_isoform_collapse',
                    'diff_method', 'cluster_match_method', 'n_top_markers']:
            yaml_key = 'gtf' if key == 'gtf_dir' else key
            if yaml_key in cfg:
                setattr(args, key, cfg[yaml_key])
    elif args.bam_dir:
        args.datasets = [{"id": "default", "merge_strategy": "none", "bams": [args.bam_dir]}]
    else:
        args.datasets = []

    return args
    
    
if __name__ == "__main__":
    cli()