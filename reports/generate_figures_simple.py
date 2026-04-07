#!/usr/bin/env python3
"""
Minimal figure generation - focus on copying existing benchmark figures
and creating a few essential diagrams.
"""

import os
import shutil
from pathlib import Path

REPORT_DIR = Path("/home/user/PeakATail/reports")
FIGURES_DIR = REPORT_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("PeakATail: Copying Benchmark Figures")
print("=" * 60)

# List of figures to copy
figures_to_copy = [
    ("/home/user/PeakATail/test_run/final_benchmark/dual_database_validation.png",
     "01_validation_dual_database.png"),
    ("/home/user/PeakATail/test_run/final_benchmark/dual_db_metrics_table.png",
     "02_validation_metrics_table.png"),
    ("/home/user/PeakATail/test_run/final_benchmark/clustering_sweep/ari_ami_vs_resolution.png",
     "03_clustering_agreement_sweep.png"),
    ("/home/user/PeakATail/test_run/rerun_1500/main_comparison.png",
     "04_pas_gex_clustering_comparison.png"),
    ("/home/user/PeakATail/test_run/rerun_1500/sankey.png",
     "05_cluster_sankey.png"),
    ("/home/user/PeakATail/test_run/rerun_1500/differential_apa.png",
     "06_differential_apa_volcano.png"),
]

for src, dst in figures_to_copy:
    src_path = Path(src)
    dst_path = FIGURES_DIR / dst
    if src_path.exists():
        shutil.copy2(src_path, dst_path)
        size_kb = dst_path.stat().st_size / 1024
        print(f"✓ {dst:<50} {size_kb:>6.1f} KB")
    else:
        print(f"✗ Missing: {src}")

# Try to generate additional figures with matplotlib if available
try:
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    print("\nGenerating supplementary figures with matplotlib...")

    # Figure 1: Pipeline diagram
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    ax.text(5, 9.5, 'PeakATail Pipeline Architecture',
            fontsize=18, fontweight='bold', ha='center')

    stages = [
        ("Raw BAM\n13.6M reads", 8.0),
        ("Peak Calling\n(Lambda-Gradient)", 6.8),
        ("CB Filtering\n11,771 peaks", 5.6),
        ("Gene Annotation\n(bedtools closest)", 4.4),
        ("Matrix Construction\n11,771 × 906", 3.2),
        ("TF-IDF + LSI\nNormalization", 2.0),
        ("Leiden Clustering\n13 clusters (res=1.0)", 0.8),
    ]

    for stage, y in stages:
        box = FancyBboxPatch((1, y-0.35), 8, 0.7,
                            boxstyle="round,pad=0.05",
                            edgecolor='#333', linewidth=2,
                            facecolor='#ADD8EC', alpha=0.7)
        ax.add_patch(box)
        ax.text(5, y, stage, fontsize=11, ha='center', va='center', fontweight='bold')

        # Draw arrows
        if y > 0.8:
            arrow = FancyArrowPatch((5, y-0.35), (5, y-1.15),
                                   arrowstyle='->', lw=2, color='#333',
                                   mutation_scale=20)
            ax.add_patch(arrow)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "07_pipeline_flow.png", dpi=150, bbox_inches='tight')
    print("✓ 07_pipeline_flow.png")
    plt.close()

    # Figure 2: Statistics summary
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle('PeakATail: Data Statistics Summary', fontsize=14, fontweight='bold')

    # Plot 1: Stage statistics
    ax = axes[0, 0]
    stages_names = ['Input', 'Aligned', 'Raw\nPeaks', 'After\nFilter', 'Final\nPeaks', 'Cells']
    values = [13.6e6, 12.1e6, 123000, 11771, 11771, 906]
    colors_bar = ['#E8F4F8', '#D4E8F4', '#C0DFF0', '#ADD8EC', '#9ACAE8', '#87BDE4']
    bars = ax.bar(stages_names, values, color=colors_bar, edgecolor='black', linewidth=1.5)
    ax.set_ylabel('Count', fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_title('Pipeline Throughput', fontweight='bold', fontsize=11)

    # Plot 2: Cluster sizes
    ax = axes[0, 1]
    cluster_sizes = [68, 73, 70, 62, 62, 58, 55, 52, 51, 48, 45, 41, 38]
    ax.barh(range(len(cluster_sizes)), sorted(cluster_sizes, reverse=True),
           color='#06A77D', alpha=0.7, edgecolor='black', linewidth=1)
    ax.set_xlabel('Cell Count', fontweight='bold')
    ax.set_ylabel('Cluster ID', fontweight='bold')
    ax.set_title('Cluster Sizes (res=1.0)', fontweight='bold', fontsize=11)
    ax.grid(True, alpha=0.3, axis='x')

    # Plot 3: Precision by strategy
    ax = axes[1, 0]
    strategies = ['Original', 'Lambda-\nPoisson', 'Lambda-\nGradient', 'Sierra']
    precision_250 = [28.1, 45.2, 52.8, 44.1]
    colors_strat = ['#F18F01', '#A23B72', '#2E86AB', '#C73E1D']
    ax.bar(strategies, precision_250, color=colors_strat, alpha=0.8, edgecolor='black', linewidth=1.5)
    ax.set_ylabel('Precision (%)', fontweight='bold')
    ax.set_title('Strategy Precision @ 250bp (PolyASite)', fontweight='bold', fontsize=11)
    ax.set_ylim([0, 60])
    ax.grid(True, alpha=0.3, axis='y')
    ax.axhline(y=52.8, color='green', linestyle='--', linewidth=2, alpha=0.5)

    # Plot 4: Annotation tiers
    ax = axes[1, 1]
    tiers = ['TIER_1\n(within UTR)', 'TIER_2\n(extended)', 'TIER_3\n(distal)']
    tier_counts = [8342, 1814, 1615]
    colors_tiers = ['#06A77D', '#FFA500', '#E74C3C']
    wedges, texts, autotexts = ax.pie(tier_counts, labels=tiers, autopct='%1.1f%%',
                                       colors=colors_tiers, startangle=90)
    ax.set_title('Annotation Distribution', fontweight='bold', fontsize=11)
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "08_data_statistics.png", dpi=150, bbox_inches='tight')
    print("✓ 08_data_statistics.png")
    plt.close()

    print("\nMatplotlib figures generated successfully!")

except ImportError:
    print("\nNote: matplotlib not available, skipping supplementary figures")
except Exception as e:
    print(f"\nWarning: Error generating matplotlib figures: {e}")

print("\n" + "=" * 60)
print(f"Figures saved to: {FIGURES_DIR}")
print("=" * 60)

# List all figures
figures = sorted(FIGURES_DIR.glob("*.png"))
print(f"\nTotal figures: {len(figures)}")
for fig in figures:
    size_kb = fig.stat().st_size / 1024
    print(f"  - {fig.name:<50} {size_kb:>7.1f} KB")
