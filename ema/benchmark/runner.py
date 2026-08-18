import logging
import os
import time
import json
import pandas as pd
from datetime import datetime
from typing import List, Optional, Dict

from ema.strategies import get_strategy, list_strategies

log = logging.getLogger(__name__)


def run_benchmark(bam_path: str, gtf_path: str,
                  strategies: Optional[List[str]] = None,
                  reference_db: Optional[str] = None,
                  output_dir: str = "benchmark_results",
                  seqlen: int = 98, cb_len: int = 16,
                  barcode_tag: str = "CB",
                  point_strand: bool = False,
                  restricted_reference_bed: Optional[str] = None,
                  null_shuffle: bool = False,
                  null_n_seeds: int = 3,
                  null_gene_body_bed: Optional[str] = None) -> pd.DataFrame:
    """Run multiple strategies on the same data and collect metrics.

    Args:
        bam_path: Path to BAM file
        gtf_path: Path to GTF file
        strategies: List of strategy names (default: all registered)
        reference_db: Path to reference PAS BED (PolyASite/PolyA_DB) for validation
        output_dir: Directory for benchmark outputs
        seqlen: Sequence length parameter
        cb_len: Cell barcode length
        barcode_tag: BAM barcode tag
        point_strand: Also score point-mode/strand-matched (adds
            ``point_precision_*`` / ``reference_coverage_*`` columns).
            Default False keeps the legacy output unchanged.
        restricted_reference_bed: Optional cohort-restricted reference BED for
            ``recall_restricted_*`` columns (implies point_strand).
        null_shuffle: Also score a shuffled-null baseline (adds
            ``null_precision_*`` columns with the across-seed mean).
        null_n_seeds: Number of shuffle seeds for the null baseline.
        null_gene_body_bed: Optional gene-body BED constraining null placement.

    Returns:
        DataFrame with metrics per strategy
    """
    if strategies is None:
        strategies = list_strategies()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(output_dir, f"benchmark_{timestamp}")
    os.makedirs(run_dir, exist_ok=True)

    results = []
    bed_files = {}

    for strategy_name in strategies:
        log.info("Running strategy: %s", strategy_name)

        strategy = get_strategy(strategy_name)
        strategy_dir = os.path.join(run_dir, strategy_name)
        os.makedirs(strategy_dir, exist_ok=True)

        pos_bed = os.path.join(strategy_dir, "posbed.bed")
        neg_bed = os.path.join(strategy_dir, "negbed.bed")
        pos_mtx = os.path.join(strategy_dir, "posmatrix.mtx")
        neg_mtx = os.path.join(strategy_dir, "negmatrix.mtx")
        combined_bed = os.path.join(strategy_dir, "combined.bed")

        # Import here to avoid circular imports with config
        from ema.countmatrix.peackcalling import peak_calling
        from ema.countmatrix.peak import Peak
        from ema.countmatrix.indexing import reset_index

        # Reset state between runs
        Peak.pasnumber = 0
        reset_index()

        start_time = time.time()

        peak_calling(True, bedfilepath=neg_bed, matrixpath=neg_mtx,
                     bamfile_dir=bam_path, strategy=strategy)

        # Reset for second strand
        Peak.pasnumber = 0
        reset_index()

        peak_calling(False, bedfilepath=pos_bed, matrixpath=pos_mtx,
                     bamfile_dir=bam_path, strategy=strategy)

        runtime = time.time() - start_time

        # Combine BED files
        with open(combined_bed, 'w') as out:
            for bed_path in [pos_bed, neg_bed]:
                if os.path.exists(bed_path):
                    with open(bed_path) as f:
                        out.write(f.read())

        bed_files[strategy_name] = combined_bed

        # Count results
        n_peaks = 0
        chroms = set()
        strand_plus = 0
        strand_minus = 0

        if os.path.exists(combined_bed):
            with open(combined_bed) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        n_peaks += 1
                        chroms.add(parts[0])
                        if parts[5] == '+':
                            strand_plus += 1
                        else:
                            strand_minus += 1

        result = {
            "strategy": strategy_name,
            "n_peaks": n_peaks,
            "n_chromosomes": len(chroms),
            "strand_plus": strand_plus,
            "strand_minus": strand_minus,
            "runtime_sec": round(runtime, 2),
            "params": json.dumps(strategy.get_params()),
            "bed_path": combined_bed
        }

        # Database validation if reference provided
        if reference_db and os.path.exists(reference_db):
            from ema.benchmark.metrics import compute_metrics
            metrics = compute_metrics(
                combined_bed, reference_db,
                point_strand=point_strand,
                restricted_reference_bed=restricted_reference_bed,
                null_shuffle=null_shuffle,
                null_n_seeds=null_n_seeds,
                null_gene_body_bed=null_gene_body_bed)

            # Add metrics at each cutoff (legacy interval/strand-agnostic)
            for cutoff, m in metrics["cutoffs"].items():
                result[f"precision_{cutoff}bp"] = m["precision"]
                result[f"recall_{cutoff}bp"] = m["recall"]
                result[f"f1_{cutoff}bp"] = m["f1"]

            # Point-mode strand-matched metrics (only present when requested)
            ps = metrics.get("point_strand")
            if ps:
                for cutoff, m in ps["cutoffs"].items():
                    result[f"point_precision_{cutoff}bp"] = m["precision"]
                    result[f"reference_coverage_{cutoff}bp"] = m["reference_coverage"]
                    if "recall_restricted" in m:
                        result[f"recall_restricted_{cutoff}bp"] = m["recall_restricted"]

            # Shuffled-null baseline (across-seed mean; only when requested)
            null = metrics.get("null")
            if null:
                for cutoff, m in null["mean"]["cutoffs"].items():
                    result[f"null_precision_{cutoff}bp"] = m["precision"]
                    result[f"null_reference_coverage_{cutoff}bp"] = m["reference_coverage"]

        results.append(result)
        log.info("  Peaks: %d | Runtime: %.1fs | Chromosomes: %d", n_peaks, runtime, len(chroms))

    # Save results
    df = pd.DataFrame(results)
    csv_path = os.path.join(run_dir, "summary.csv")
    df.to_csv(csv_path, index=False)
    log.info("Results saved to: %s", csv_path)

    # Save metadata
    meta = {
        "timestamp": timestamp,
        "bam_path": bam_path,
        "gtf_path": gtf_path,
        "strategies": strategies,
        "reference_db": reference_db
    }
    with open(os.path.join(run_dir, "metadata.json"), 'w') as f:
        json.dump(meta, f, indent=2)

    return df


def compare_results(results_df: pd.DataFrame, output_dir: str = ".") -> str:
    """Generate a text comparison report from benchmark results.

    Args:
        results_df: DataFrame from run_benchmark()
        output_dir: Where to save the report

    Returns:
        Path to generated report
    """
    report_lines = []
    report_lines.append("=" * 70)
    report_lines.append("PeakATail Strategy Benchmark Report")
    report_lines.append("=" * 70)
    report_lines.append("")

    # Summary table
    summary_cols = ["strategy", "n_peaks", "runtime_sec"]
    f1_cols = [c for c in results_df.columns if c.startswith("f1_")]

    report_lines.append("Strategy Comparison:")
    report_lines.append("-" * 50)

    for _, row in results_df.iterrows():
        report_lines.append(f"\n  {row['strategy']}:")
        report_lines.append(f"    Peaks found: {row['n_peaks']}")
        report_lines.append(f"    Runtime: {row['runtime_sec']}s")
        report_lines.append(f"    Strands: {row['strand_plus']}+ / {row['strand_minus']}-")

        for col in f1_cols:
            cutoff = col.replace("f1_", "").replace("bp", "")
            p_col = f"precision_{cutoff}bp"
            r_col = f"recall_{cutoff}bp"
            if p_col in row and r_col in row:
                report_lines.append(
                    f"    @{cutoff}bp: P={row[p_col]:.3f} R={row[r_col]:.3f} F1={row[col]:.3f}"
                )

    report_lines.append("")
    report_lines.append("=" * 70)

    report_text = "\n".join(report_lines)
    report_path = os.path.join(output_dir, "benchmark_report.txt")
    with open(report_path, 'w') as f:
        f.write(report_text)

    log.info("%s", report_text)
    return report_path
