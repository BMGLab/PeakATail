import os
import json
import pandas as pd
from typing import Optional


def generate_report(results_df: pd.DataFrame, output_dir: str,
                    reference_db: Optional[str] = None) -> str:
    """Generate a complete benchmark report directory.

    Args:
        results_df: DataFrame from run_benchmark()
        output_dir: Directory for report files
        reference_db: Path to reference database used

    Returns:
        Path to report directory
    """
    os.makedirs(output_dir, exist_ok=True)

    # CSV summary
    results_df.to_csv(os.path.join(output_dir, "summary.csv"), index=False)

    # Per-strategy details
    for _, row in results_df.iterrows():
        strategy_name = row["strategy"]
        detail = {
            "strategy": strategy_name,
            "n_peaks": int(row["n_peaks"]),
            "runtime_sec": float(row["runtime_sec"]),
            "strand_plus": int(row["strand_plus"]),
            "strand_minus": int(row["strand_minus"]),
            "n_chromosomes": int(row["n_chromosomes"]),
            "params": json.loads(row["params"]) if isinstance(row["params"], str) else row["params"]
        }

        # Add precision/recall if available
        metrics = {}
        for col in row.index:
            if col.startswith("precision_") or col.startswith("recall_") or col.startswith("f1_"):
                try:
                    metrics[col] = float(row[col])
                except (ValueError, TypeError):
                    pass
        if metrics:
            detail["metrics"] = metrics

        detail_path = os.path.join(output_dir, f"{strategy_name}_detail.json")
        with open(detail_path, 'w') as f:
            json.dump(detail, f, indent=2)

    # Comparison table (markdown)
    md_path = os.path.join(output_dir, "comparison.md")
    _write_markdown_table(results_df, md_path)

    return output_dir


def _write_markdown_table(df: pd.DataFrame, path: str):
    """Write a markdown comparison table."""
    cols = ["strategy", "n_peaks", "runtime_sec"]
    metric_cols = sorted([c for c in df.columns if c.startswith("f1_")])

    with open(path, 'w') as f:
        f.write("# PeakATail Strategy Benchmark\n\n")

        # Basic stats
        header = "| Strategy | Peaks | Runtime (s) |"
        sep = "|----------|-------|-------------|"
        for mc in metric_cols:
            cutoff = mc.replace("f1_", "").replace("bp", "")
            header += f" F1@{cutoff}bp |"
            sep += "---------|"

        f.write(header + "\n")
        f.write(sep + "\n")

        for _, row in df.iterrows():
            line = f"| {row['strategy']} | {row['n_peaks']} | {row['runtime_sec']} |"
            for mc in metric_cols:
                val = row.get(mc, "N/A")
                if isinstance(val, float):
                    line += f" {val:.3f} |"
                else:
                    line += f" {val} |"
            f.write(line + "\n")
