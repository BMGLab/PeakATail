import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt
from typing import Dict, List, Optional
import pandas as pd


def plot_precision_by_distance(results_df: pd.DataFrame,
                                output_path: str = "precision_by_distance.png"):
    """Plot precision vs distance cutoff for each strategy.

    Primary diagnostic plot for APA analysis: "How does precision improve
    as we widen the matching window?"
    """
    fig, ax = plt.subplots(figsize=(10, 7))

    strategies = results_df["strategy"].unique()
    colors = ['#2ca02c', '#1f77b4', '#ff7f0e', '#7f7f7f']

    cutoff_cols = sorted(
        [c for c in results_df.columns if c.startswith("precision_")],
        key=lambda c: int(c.replace("precision_", "").replace("bp", ""))
    )
    cutoffs = [int(c.replace("precision_", "").replace("bp", "")) for c in cutoff_cols]

    for idx, strategy in enumerate(strategies):
        row = results_df[results_df["strategy"] == strategy].iloc[0]
        color = colors[idx % len(colors)]

        precisions = []
        valid_cutoffs = []
        for col, cutoff in zip(cutoff_cols, cutoffs):
            val = row.get(col)
            if pd.notna(val):
                precisions.append(val)
                valid_cutoffs.append(cutoff)

        if not precisions:
            continue

        ax.plot(valid_cutoffs, precisions, 'o-', color=color, label=strategy,
                markersize=9, linewidth=2.5)

        for cutoff, p in zip(valid_cutoffs, precisions):
            ax.annotate(f'{p:.1%}', (cutoff, p), textcoords="offset points",
                        xytext=(6, 8), fontsize=9, fontweight='bold', color=color)

    # Target lines
    precision_targets = {50: 0.70, 200: 0.80, 1000: 0.90, 5000: 0.95}
    for cutoff, target in precision_targets.items():
        ax.axhline(y=target, color='red', linestyle='--', linewidth=0.8, alpha=0.4)
        ax.text(max(cutoffs) * 1.05, target, f'target {target:.0%}',
                fontsize=8, color='red', alpha=0.6, va='center')

    ax.set_xscale('log')
    ax.set_xlabel("Distance Cutoff (bp, log scale)", fontsize=13)
    ax.set_ylabel("Precision", fontsize=13)
    ax.set_title("Precision by Distance Cutoff — Are Our PAS Real?", fontsize=15, fontweight='bold')
    ax.legend(loc="upper left", fontsize=11, framealpha=0.9)
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.2)
    ax.set_xticks(cutoffs)
    ax.set_xticklabels([f'{c}bp' for c in cutoffs], fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Precision-by-distance plot saved to: {output_path}")


def plot_pr_curves(results_df: pd.DataFrame,
                   output_path: str = "pr_curves.png"):
    """Plot traditional precision-recall curves for each strategy."""
    fig, ax = plt.subplots(figsize=(10, 7))

    strategies = results_df["strategy"].unique()
    colors = ['#2ca02c', '#1f77b4', '#ff7f0e', '#7f7f7f']

    for idx, strategy in enumerate(strategies):
        row = results_df[results_df["strategy"] == strategy].iloc[0]
        color = colors[idx % len(colors)]

        precisions = []
        recalls = []
        cutoffs = []

        for col in sorted(results_df.columns):
            if col.startswith("precision_"):
                cutoff = int(col.replace("precision_", "").replace("bp", ""))
                p_val = row.get(col)
                r_col = f"recall_{cutoff}bp"
                r_val = row.get(r_col)

                if pd.notna(p_val) and pd.notna(r_val):
                    precisions.append(p_val)
                    recalls.append(r_val)
                    cutoffs.append(cutoff)

        if precisions:
            ax.plot(recalls, precisions, 'o-', color=color, label=strategy,
                    markersize=9, linewidth=2.5)

            for p, r, c in zip(precisions, recalls, cutoffs):
                ax.annotate(f'{c}bp', (r, p), textcoords="offset points",
                           xytext=(6, 6), fontsize=9, color=color)

    ax.set_xlabel("Recall", fontsize=13)
    ax.set_ylabel("Precision", fontsize=13)
    ax.set_title("Precision-Recall Curves by Strategy", fontsize=15, fontweight='bold')
    ax.legend(loc="upper right", fontsize=11)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"PR curves saved to: {output_path}")


def plot_peak_count_comparison(results_df: pd.DataFrame,
                                output_path: str = "peak_counts.png"):
    """Bar chart comparing peak counts across strategies."""
    # Only show runtime panel if we have real runtime data
    has_runtime = results_df["runtime_sec"].sum() > 0

    if has_runtime:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    else:
        fig, ax_peaks = plt.subplots(figsize=(8, 6))
        axes = [ax_peaks]

    strategies = results_df["strategy"].tolist()
    n_peaks = results_df["n_peaks"].tolist()
    colors = ['#2ca02c', '#1f77b4', '#ff7f0e', '#7f7f7f'][:len(strategies)]

    # Peak counts
    bars = axes[0].bar(strategies, n_peaks, color=colors, edgecolor='white', linewidth=0.5)
    axes[0].set_ylabel("Number of Peaks", fontsize=12)
    axes[0].set_title("Peaks Found per Strategy", fontsize=14, fontweight='bold')
    axes[0].tick_params(axis='x', rotation=20)
    for bar, count in zip(bars, n_peaks):
        axes[0].text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                    f'{count:,}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    axes[0].grid(True, alpha=0.2, axis='y')

    # Runtime (only if data exists)
    if has_runtime:
        runtimes = results_df["runtime_sec"].tolist()
        bars = axes[1].bar(strategies, runtimes, color=colors, edgecolor='white', linewidth=0.5)
        axes[1].set_ylabel("Runtime (seconds)", fontsize=12)
        axes[1].set_title("Runtime per Strategy", fontsize=14, fontweight='bold')
        axes[1].tick_params(axis='x', rotation=20)
        for bar, rt in zip(bars, runtimes):
            axes[1].text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                        f'{rt:.1f}s', ha='center', va='bottom', fontsize=11)
        axes[1].grid(True, alpha=0.2, axis='y')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Peak count comparison saved to: {output_path}")


def plot_distance_histogram(distances: Dict[str, List[int]],
                            output_path: str = "distance_hist.png",
                            max_distance: int = 5000):
    """Histogram of distances to nearest reference PAS per strategy.

    Uses two panels: close-range (0-200bp) and full-range (0-max_distance)
    with log Y scale to show both the spike and the tail.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    colors = ['#2ca02c', '#1f77b4', '#ff7f0e', '#7f7f7f']

    # Panel 1: Close-range (0-200bp) — linear Y, fine bins
    bins_close = np.arange(0, 210, 5)
    for idx, (name, dists) in enumerate(distances.items()):
        color = colors[idx % len(colors)]
        filtered = [d for d in dists if d <= 200]
        if filtered:
            axes[0].hist(filtered, bins=bins_close, alpha=0.6, label=name,
                        color=color, edgecolor='white', linewidth=0.3)

    axes[0].set_xlabel("Distance to Nearest Known PAS (bp)", fontsize=12)
    axes[0].set_ylabel("Count", fontsize=12)
    axes[0].set_title("Close Range (0-200bp)", fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.2, axis='y')

    # Panel 2: Full range (0-max_distance) — log Y, wider bins
    bins_full = np.arange(0, max_distance + 100, 50)
    for idx, (name, dists) in enumerate(distances.items()):
        color = colors[idx % len(colors)]
        filtered = [d for d in dists if d <= max_distance]
        if filtered:
            axes[1].hist(filtered, bins=bins_full, alpha=0.6, label=name,
                        color=color, edgecolor='white', linewidth=0.3)

    axes[1].set_xlabel("Distance to Nearest Known PAS (bp)", fontsize=12)
    axes[1].set_ylabel("Count (log scale)", fontsize=12)
    axes[1].set_title(f"Full Range (0-{max_distance}bp)", fontsize=13, fontweight='bold')
    axes[1].set_yscale('log')
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.2, axis='y')

    plt.suptitle("Distance Distribution to Reference PAS (PolyASite 2.0)",
                 fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Distance histogram saved to: {output_path}")


def plot_metrics_table(results_df: pd.DataFrame,
                       output_path: str = "metrics_table.png"):
    """Render key metrics as a readable table figure."""
    # Select only the most important columns
    base_cols = ["strategy", "n_peaks"]
    precision_cols = sorted(
        [c for c in results_df.columns if c.startswith("precision_")],
        key=lambda c: int(c.replace("precision_", "").replace("bp", ""))
    )

    display_cols = base_cols + precision_cols
    display_cols = [c for c in display_cols if c in results_df.columns]
    display_df = results_df[display_cols].copy()

    # Rename columns for readability
    rename = {"strategy": "Strategy", "n_peaks": "Peaks"}
    for col in precision_cols:
        cutoff = col.replace("precision_", "").replace("bp", "")
        rename[col] = f"P@{cutoff}bp"
    display_df = display_df.rename(columns=rename)

    # Format precision as percentage
    for col in display_df.columns:
        if col.startswith("P@"):
            display_df[col] = display_df[col].apply(
                lambda x: f"{x:.1%}" if pd.notna(x) else "N/A"
            )
    display_df["Peaks"] = display_df["Peaks"].apply(lambda x: f"{x:,}")

    n_rows = len(display_df)
    n_cols = len(display_df.columns)

    fig, ax = plt.subplots(figsize=(n_cols * 1.4, n_rows * 0.7 + 1.5))
    ax.axis('off')

    table = ax.table(
        cellText=display_df.values,
        colLabels=display_df.columns,
        cellLoc='center',
        loc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.0, 2.0)

    # Style header
    for j in range(n_cols):
        table[0, j].set_facecolor('#2E75B6')
        table[0, j].set_text_props(color='white', fontweight='bold', fontsize=11)

    # Alternate row colors
    for i in range(1, n_rows + 1):
        bg = '#F2F2F2' if i % 2 == 0 else '#FFFFFF'
        for j in range(n_cols):
            table[i, j].set_facecolor(bg)

    # Highlight best precision per column
    for j, col in enumerate(display_df.columns):
        if col.startswith("P@"):
            vals = []
            for i in range(n_rows):
                try:
                    vals.append(float(display_df.iloc[i][col].replace('%', '')) / 100)
                except (ValueError, AttributeError):
                    vals.append(0)
            best_idx = np.argmax(vals)
            table[best_idx + 1, j].set_text_props(fontweight='bold', color='#006400')

    ax.set_title("PeakATail Strategy Benchmark — Precision Summary",
                 fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Metrics table saved to: {output_path}")


def plot_venn_overlap(sets: Dict[str, set], output_path: str = "venn_overlap.png"):
    """Venn diagram showing PAS overlap between 2-3 strategies."""
    n = len(sets)
    names = list(sets.keys())

    if n < 2 or n > 3:
        print(f"Venn diagram requires 2-3 strategies, got {n}. Skipping.")
        return

    fig, ax = plt.subplots(figsize=(10, 8))

    try:
        from matplotlib_venn import venn2, venn3

        if n == 2:
            v = venn2([sets[names[0]], sets[names[1]]],
                      set_labels=names, ax=ax)
        else:
            v = venn3([sets[names[0]], sets[names[1]], sets[names[2]]],
                      set_labels=names, ax=ax)

        # Style the labels
        for text in ax.texts:
            text.set_fontsize(11)

        ax.set_title("PAS Detection Overlap Between Strategies",
                     fontsize=14, fontweight='bold')

    except ImportError:
        ax.axis('off')
        text_lines = ["PAS Overlap (matplotlib-venn not installed)\n"]
        for i, n1 in enumerate(names):
            for n2 in names[i+1:]:
                overlap = len(sets[n1] & sets[n2])
                only_1 = len(sets[n1] - sets[n2])
                only_2 = len(sets[n2] - sets[n1])
                text_lines.append(f"{n1} \u2229 {n2}: {overlap:,} shared")
                text_lines.append(f"  {n1} only: {only_1:,}, {n2} only: {only_2:,}")
        ax.text(0.5, 0.5, "\n".join(text_lines), transform=ax.transAxes,
               fontsize=13, va='center', ha='center', family='monospace')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Venn diagram saved to: {output_path}")
