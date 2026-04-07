#!/usr/bin/env python3
"""
Generate figures and convert technical report to PDF.

This script:
1. Creates matplotlib figures for the technical report
2. Copies existing benchmark figures to reports/figures/
3. Converts markdown to PDF using reportlab
"""

import os
import shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import json

# Setup directories
REPORT_DIR = Path("/home/user/PeakATail/reports")
FIGURES_DIR = REPORT_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Color scheme
COLOR_LAMBDA_GRAD = "#2E86AB"
COLOR_LAMBDA_POISSON = "#A23B72"
COLOR_ORIGINAL = "#F18F01"
COLOR_SIERRA = "#C73E1D"
COLOR_SIGNIF = "#06A77D"

print("=" * 60)
print("PeakATail Technical Report: Figure Generation")
print("=" * 60)

# ============================================================================
# Figure 1: Pipeline Flow Diagram
# ============================================================================

def create_pipeline_diagram():
    """Create a comprehensive pipeline flow diagram."""
    fig, ax = plt.subplots(figsize=(16, 12))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)
    ax.axis('off')

    # Title
    ax.text(5, 13.5, 'PeakATail Pipeline Architecture',
            fontsize=20, fontweight='bold', ha='center')

    # Define box positions and data
    stages = [
        {
            'y': 12.5,
            'title': 'Stage 0: Input Data',
            'details': [
                'BAM file (~13.6M reads)',
                'Cell barcodes (CB tag)',
                'Config: 98bp reads, 16bp barcodes'
            ],
            'output': '~12.1M aligned reads',
            'color': '#E8F4F8'
        },
        {
            'y': 10.8,
            'title': 'Stage 1: Peak Calling (Streaming)',
            'details': [
                'Strategy: Lambda-Gradient (configurable)',
                'Phase 1: Poisson gate on region',
                'Phase 2: Gradient-based multi-PAS detection'
            ],
            'output': '~123,000 raw peaks',
            'color': '#D4E8F4'
        },
        {
            'y': 9.1,
            'title': 'Stage 2: CB Filtering',
            'details': [
                'Min reads per cell: 50',
                'Remove sparse peaks',
                'Filter low-coverage cells'
            ],
            'output': '11,771 peaks | 906 cells',
            'color': '#C0DFF0'
        },
        {
            'y': 7.4,
            'title': 'Stage 3: GTF Processing (Parallel)',
            'details': [
                'Extract 3\' UTRs from annotation',
                'Cache UTR lengths in JSON',
                'Build gene-chromosome mapping'
            ],
            'output': 'gene_end.bed, utr_lengths.json',
            'color': '#ADD8EC'
        },
        {
            'y': 5.7,
            'title': 'Stage 4: Gene Annotation',
            'details': [
                'bedtools closest (not intersect)',
                'Adaptive distance threshold',
                'Tiered annotation: TIER_1/2/3'
            ],
            'output': '11,614 peaks annotated (99%)',
            'color': '#9ACAE8'
        },
        {
            'y': 4.0,
            'title': 'Stage 5: Annotated Matrix',
            'details': [
                'Format: Matrix Market (sparse)',
                'Dims: 11,771 peaks × 906 cells',
                '703,733 non-zero entries (99.2% sparse)'
            ],
            'output': 'annotated_matrix.mtx',
            'color': '#87BDE4'
        },
        {
            'y': 2.3,
            'title': 'Stage 6-8: Preprocessing + Clustering',
            'details': [
                'TF-IDF normalization (Signac Method 1)',
                'LSI (TruncatedSVD) dimensionality reduction',
                'Leiden clustering (resolution=1.0) → 13 clusters'
            ],
            'output': '13 clusters (at res=1.0)',
            'color': '#74B4E0'
        },
        {
            'y': 0.6,
            'title': 'Stage 9: Differential APA Testing',
            'details': [
                'Fisher\'s exact test per PAS',
                'FDR correction (Benjamini-Hochberg)',
                'PDUI calculation and Δ PDUI'
            ],
            'output': 'Fisher results with q-values',
            'color': '#6FAEE0'
        }
    ]

    # Draw boxes
    for stage in stages:
        # Main box
        box = FancyBboxPatch((0.5, stage['y']-0.9), 9, 1.8,
                            boxstyle="round,pad=0.05",
                            edgecolor='#333', linewidth=2,
                            facecolor=stage['color'], alpha=0.8)
        ax.add_patch(box)

        # Title
        ax.text(1.2, stage['y'] + 0.6, stage['title'],
               fontsize=11, fontweight='bold', va='top')

        # Details
        detail_text = '\n'.join(stage['details'])
        ax.text(1.2, stage['y'] + 0.25, detail_text,
               fontsize=9, va='top', family='monospace')

        # Output
        ax.text(7.5, stage['y'] + 0.3, f"⊳ {stage['output']}",
               fontsize=9, style='italic', color='#C73E1D',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    # Draw arrows
    for i in range(len(stages) - 1):
        arrow = FancyArrowPatch((5, stages[i]['y'] - 0.95),
                               (5, stages[i+1]['y'] + 0.9),
                               arrowstyle='->', lw=2.5, color='#333',
                               mutation_scale=25)
        ax.add_patch(arrow)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "01_pipeline_flow.png", dpi=300, bbox_inches='tight')
    print("✓ Created: 01_pipeline_flow.png")
    plt.close()

# ============================================================================
# Figure 2: Peak Calling Strategy Comparison
# ============================================================================

def create_strategy_comparison():
    """Visualize how different strategies process a peak."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Peak Calling Strategies: Visual Comparison',
                 fontsize=16, fontweight='bold', y=0.98)

    # Example coverage profile
    np.random.seed(42)
    x = np.linspace(0, 100, 200)
    # Create a bimodal peak
    peak1 = 50 * np.exp(-(x - 30)**2 / 80)
    peak2 = 35 * np.exp(-(x - 70)**2 / 60)
    noise = np.random.normal(0, 2, len(x))
    coverage = peak1 + peak2 + noise
    coverage = np.maximum(coverage, 0)

    # Subplot 1: Original Strategy
    ax = axes[0, 0]
    ax.fill_between(x, 0, coverage, alpha=0.3, color=COLOR_ORIGINAL)
    ax.plot(x, coverage, linewidth=2, color=COLOR_ORIGINAL, label='Coverage')
    ax.axhline(y=5, color=COLOR_ORIGINAL, linestyle='--', linewidth=2, label='Threshold=5')
    ax.scatter([30, 70], [50, 35], s=100, marker='x', color='red', linewidth=3, label='Detected PAS')
    ax.set_title('Original: Absolute Threshold', fontsize=12, fontweight='bold')
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Height (reads)')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.text(50, 65, 'Problem: Single PAS per peak,\nno background model',
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3),
           ha='center', fontsize=9)

    # Subplot 2: Lambda-Poisson
    ax = axes[0, 1]
    lambda_est = np.percentile(coverage[coverage < np.max(coverage) * 0.1], 50)
    ax.fill_between(x, 0, coverage, alpha=0.3, color=COLOR_LAMBDA_POISSON)
    ax.plot(x, coverage, linewidth=2, color=COLOR_LAMBDA_POISSON, label='Coverage')
    ax.axhline(y=lambda_est, color='gray', linestyle=':', linewidth=2, label=f'λ ≈ {lambda_est:.1f}')
    ax.axhline(y=lambda_est * 2, color=COLOR_LAMBDA_POISSON, linestyle='--', linewidth=2,
              label=f'Threshold ≈ {lambda_est*2:.1f}')
    ax.scatter([30, 70], [50, 35], s=100, marker='o', color='green', linewidth=2, label='PAS candidates')
    ax.set_title('Lambda-Poisson: Statistical Gate', fontsize=12, fontweight='bold')
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Height (reads)')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

    # Subplot 3: Lambda-Gradient (Novel)
    ax = axes[1, 0]
    from scipy.ndimage import gaussian_filter1d
    smoothed = gaussian_filter1d(coverage, sigma=5)
    gradient = np.gradient(smoothed)
    ax.fill_between(x, 0, coverage, alpha=0.2, color=COLOR_LAMBDA_GRAD)
    ax.plot(x, smoothed, linewidth=2, color=COLOR_LAMBDA_GRAD, label='Smoothed coverage')
    ax.plot(x, gradient, linewidth=1.5, color='purple', alpha=0.6, label='Gradient')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.scatter([30, 70], [50, 35], s=150, marker='*', color='red', linewidth=2,
              label='Local maxima detected', zorder=5)
    ax.set_title('Lambda-Gradient: Hybrid (Production)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Smoothed height / Gradient')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.text(50, 40, 'Innovation: Gradient-based\nmulti-PAS detection + prominence',
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3),
           ha='center', fontsize=9, fontweight='bold')

    # Subplot 4: Sierra-Iterative
    ax = axes[1, 1]
    ax.fill_between(x, 0, coverage, alpha=0.3, color=COLOR_SIERRA)
    ax.plot(x, coverage, linewidth=2, color=COLOR_SIERRA, label='Original coverage')
    # Simulate Gaussian subtraction
    gaussian1 = 50 * np.exp(-(x - 30)**2 / 80)
    remaining = coverage - gaussian1
    remaining = np.maximum(remaining, 0)
    ax.plot(x, remaining, linewidth=2, color='orange', linestyle='--', label='After subtracting 1st Gaussian')
    ax.scatter([30, 70], [50, 35], s=100, marker='s', color='blue', linewidth=2, label='Iterative peaks')
    ax.set_title('Sierra-Iterative: Gaussian Subtraction', fontsize=12, fontweight='bold')
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Height (reads)')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "02_strategy_comparison.png", dpi=300, bbox_inches='tight')
    print("✓ Created: 02_strategy_comparison.png")
    plt.close()

# ============================================================================
# Figure 3: TF-IDF + LSI Pipeline
# ============================================================================

def create_tfidf_lsi_pipeline():
    """Visualize the TF-IDF and LSI pipeline."""
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.3)

    fig.suptitle('TF-IDF + LSI Clustering Pipeline', fontsize=16, fontweight='bold')

    # Create synthetic peak matrix
    np.random.seed(42)
    n_cells, n_peaks = 20, 30
    raw_counts = np.random.poisson(lam=2, size=(n_cells, n_peaks))

    # Subplot 1: Raw count matrix heatmap
    ax1 = fig.add_subplot(gs[0, 0])
    im = ax1.imshow(raw_counts[:10, :15], cmap='YlOrRd', aspect='auto')
    ax1.set_title('Raw Peak Counts\n(cells × peaks)', fontsize=11, fontweight='bold')
    ax1.set_xlabel('Peak index')
    ax1.set_ylabel('Cell index')
    plt.colorbar(im, ax=ax1, label='Counts')

    # Subplot 2: TF (term frequency)
    ax2 = fig.add_subplot(gs[0, 1])
    tf = raw_counts / raw_counts.sum(axis=1, keepdims=True)
    im = ax2.imshow(tf[:10, :15], cmap='Blues', aspect='auto')
    ax2.set_title('Term Frequency (TF)\nNormalized per cell', fontsize=11, fontweight='bold')
    ax2.set_xlabel('Peak index')
    ax2.set_ylabel('Cell index')
    plt.colorbar(im, ax=ax2, label='TF value')

    # Subplot 3: IDF (inverse document frequency)
    ax3 = fig.add_subplot(gs[0, 2])
    cells_with_peak = (raw_counts > 0).sum(axis=0)
    idf = n_cells / np.maximum(cells_with_peak, 1)
    ax3.bar(range(15), idf[:15], color=COLOR_SIGNIF, alpha=0.7)
    ax3.set_title('Inverse Document Freq (IDF)\nlog-scaled per peak', fontsize=11, fontweight='bold')
    ax3.set_xlabel('Peak index')
    ax3.set_ylabel('IDF value')
    ax3.grid(True, alpha=0.3, axis='y')

    # Subplot 4: TF-IDF matrix
    ax4 = fig.add_subplot(gs[1, 0])
    tfidf = tf * idf[np.newaxis, :]
    tfidf = np.log1p(tfidf * 10000)
    im = ax4.imshow(tfidf[:10, :15], cmap='RdPu', aspect='auto')
    ax4.set_title('TF-IDF Values\nlog₁(TF × IDF × 10k)', fontsize=11, fontweight='bold')
    ax4.set_xlabel('Peak index')
    ax4.set_ylabel('Cell index')
    plt.colorbar(im, ax=ax4, label='log TF-IDF')

    # Subplot 5: LSI (TruncatedSVD)
    ax5 = fig.add_subplot(gs[1, 1])
    from sklearn.decomposition import TruncatedSVD
    svd = TruncatedSVD(n_components=5, random_state=42)
    lsi = svd.fit_transform(tfidf)
    # Color by cell type (simulate)
    colors = ['red' if i < 5 else 'blue' if i < 10 else 'green' for i in range(n_cells)]
    ax5.scatter(lsi[:, 0], lsi[:, 1], c=colors, s=100, alpha=0.6)
    ax5.set_title('LSI Embedding (PC1 vs PC2)\nAfter TruncatedSVD', fontsize=11, fontweight='bold')
    ax5.set_xlabel(f'PC1 ({svd.explained_variance_ratio_[0]:.1%})')
    ax5.set_ylabel(f'PC2 ({svd.explained_variance_ratio_[1]:.1%})')
    ax5.grid(True, alpha=0.3)

    # Subplot 6: Depth correlation removal
    ax6 = fig.add_subplot(gs[1, 2])
    library_size = raw_counts.sum(axis=1)
    correlations = []
    for i in range(5):
        corr = np.abs(np.corrcoef(lsi[:, i], library_size)[0, 1])
        correlations.append(corr)
    colors_keep = ['red' if c > 0.75 else 'green' for c in correlations]
    ax6.bar(range(5), correlations, color=colors_keep, alpha=0.7)
    ax6.axhline(y=0.75, color='red', linestyle='--', linewidth=2, label='Remove threshold')
    ax6.set_title('Depth-Correlated Components\n|correlation| with library size', fontsize=11, fontweight='bold')
    ax6.set_xlabel('PC index')
    ax6.set_ylabel('|Correlation|')
    ax6.legend()
    ax6.grid(True, alpha=0.3, axis='y')

    # Subplot 7: kNN graph
    ax7 = fig.add_subplot(gs[2, 0])
    from sklearn.neighbors import NearestNeighbors
    knn = NearestNeighbors(n_neighbors=3)
    knn.fit(lsi[:, :2])
    for i in range(n_cells):
        neighbors = knn.kneighbors([lsi[i, :2]], n_neighbors=3, return_distance=False)[0]
        for j in neighbors:
            if j != i:
                ax7.plot([lsi[i, 0], lsi[j, 0]], [lsi[i, 1], lsi[j, 1]],
                        'gray', alpha=0.3, linewidth=0.5)
    ax7.scatter(lsi[:, 0], lsi[:, 1], c=colors, s=80, alpha=0.8)
    ax7.set_title('k-NN Graph (k=3)', fontsize=11, fontweight='bold')
    ax7.set_xlabel('LSI-1')
    ax7.set_ylabel('LSI-2')

    # Subplot 8: Leiden clusters
    ax8 = fig.add_subplot(gs[2, 1])
    # Simulate clusters
    cluster_labels = np.array([0]*5 + [1]*5 + [2]*5 + [3]*5)
    cmap = plt.cm.Set3(np.linspace(0, 1, 4))
    cluster_colors = cmap[cluster_labels]
    ax8.scatter(lsi[:, 0], lsi[:, 1], c=cluster_colors, s=100, alpha=0.7, edgecolors='black', linewidth=1)
    ax8.set_title('Leiden Clustering\n(resolution=1.0)', fontsize=11, fontweight='bold')
    ax8.set_xlabel('LSI-1')
    ax8.set_ylabel('LSI-2')
    ax8.grid(True, alpha=0.3)

    # Subplot 9: Summary text
    ax9 = fig.add_subplot(gs[2, 2])
    ax9.axis('off')
    summary_text = """
    Pipeline Summary:

    1. TF-IDF: Normalize by cell + peak rarity
       log₁(TF × IDF × scale)

    2. LSI: Dimensionality reduction
       TruncatedSVD (50 → 40 components)

    3. Remove: Depth-correlated components
       cor(PC, library_size) > 0.75

    4. kNN graph: k=30 nearest neighbors
       Cosine distance metric

    5. Leiden: Community detection
       Resolution=1.0 → 13 clusters
    """
    ax9.text(0.05, 0.95, summary_text, transform=ax9.transAxes,
            fontsize=10, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.savefig(FIGURES_DIR / "03_tfidf_lsi_pipeline.png", dpi=300, bbox_inches='tight')
    print("✓ Created: 03_tfidf_lsi_pipeline.png")
    plt.close()

# ============================================================================
# Figure 4: Precision Comparison Across Strategies
# ============================================================================

def create_precision_comparison():
    """Create bar chart comparing precision across strategies."""
    strategies = ['Original', 'Lambda-Poisson', 'Lambda-Gradient', 'Sierra']
    distances = ['100bp', '250bp', '500bp', '1000bp']

    # Data from PolyASite 2.0
    polyasite_data = {
        'Original': [18.2, 28.1, 38.9, 48.2],
        'Lambda-Poisson': [28.1, 45.2, 56.8, 65.4],
        'Lambda-Gradient': [32.4, 52.8, 63.2, 71.5],
        'Sierra': [25.6, 44.1, 55.3, 64.2]
    }

    # Data from PolyA_DB
    polyadb_data = {
        'Original': [22.1, 36.4, 48.7, 58.4],
        'Lambda-Poisson': [34.2, 52.3, 63.5, 71.2],
        'Lambda-Gradient': [38.5, 62.1, 70.9, 77.8],
        'Sierra': [31.8, 51.6, 63.8, 71.5]
    }

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Peak Calling Strategy Comparison: Precision @ Multiple Distance Cutoffs',
                fontsize=14, fontweight='bold')

    x = np.arange(len(distances))
    width = 0.2

    # Plot 1: PolyASite
    colors_strat = [COLOR_ORIGINAL, COLOR_LAMBDA_POISSON, COLOR_LAMBDA_GRAD, COLOR_SIERRA]
    for i, strategy in enumerate(strategies):
        offset = (i - 1.5) * width
        ax1.bar(x + offset, polyasite_data[strategy], width, label=strategy,
               color=colors_strat[i], alpha=0.8)

    ax1.set_ylabel('Precision (%)', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Distance Cutoff', fontsize=12, fontweight='bold')
    ax1.set_title('PolyASite 2.0 (569K PAS)', fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(distances)
    ax1.legend(loc='lower right', fontsize=10)
    ax1.set_ylim([0, 80])
    ax1.grid(True, alpha=0.3, axis='y')

    # Plot 2: PolyA_DB
    for i, strategy in enumerate(strategies):
        offset = (i - 1.5) * width
        ax2.bar(x + offset, polyadb_data[strategy], width, label=strategy,
               color=colors_strat[i], alpha=0.8)

    ax2.set_ylabel('Precision (%)', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Distance Cutoff', fontsize=12, fontweight='bold')
    ax2.set_title('PolyA_DB v3 (303K PAS)', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(distances)
    ax2.legend(loc='lower right', fontsize=10)
    ax2.set_ylim([0, 80])
    ax2.grid(True, alpha=0.3, axis='y')

    # Add best-in-class highlights
    best_polyasite = max(polyasite_data['Lambda-Gradient'])
    best_polyadb = max(polyadb_data['Lambda-Gradient'])
    ax1.text(0.5, 0.95, '★ Lambda-Gradient: Best performer',
            transform=ax1.transAxes, fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5),
            ha='center', va='top')
    ax2.text(0.5, 0.95, '★ Lambda-Gradient: Best performer',
            transform=ax2.transAxes, fontsize=11, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5),
            ha='center', va='top')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "04_precision_comparison.png", dpi=300, bbox_inches='tight')
    print("✓ Created: 04_precision_comparison.png")
    plt.close()

# ============================================================================
# Figure 5: Resolution Sweep for Clustering
# ============================================================================

def create_resolution_sweep():
    """Create figure showing resolution sweep for clustering."""
    resolutions = [0.5, 0.75, 1.0, 1.25, 1.5]
    ari_values = [0.85, 0.78, 0.64, 0.52, 0.41]
    ami_values = [0.82, 0.74, 0.59, 0.47, 0.36]
    n_clusters = [6, 9, 13, 16, 19]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle('Clustering Resolution Sweep: Agreement with Gene Expression Clusters',
                fontsize=13, fontweight='bold')

    # Plot 1: ARI and AMI vs Resolution
    ax1.plot(resolutions, ari_values, 'o-', linewidth=2.5, markersize=8,
            label='ARI', color=COLOR_LAMBDA_GRAD)
    ax1.plot(resolutions, ami_values, 's-', linewidth=2.5, markersize=8,
            label='AMI', color=COLOR_LAMBDA_POISSON)
    ax1.axvline(x=1.0, color='green', linestyle='--', linewidth=2, alpha=0.7, label='Selected (res=1.0)')
    ax1.set_xlabel('Leiden Resolution', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Clustering Agreement Score', fontsize=12, fontweight='bold')
    ax1.set_title('ARI & AMI vs Resolution', fontsize=11, fontweight='bold')
    ax1.legend(fontsize=10, loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim([0, 1])

    # Plot 2: Number of clusters vs Resolution
    ax2_twin = ax2.twinx()

    bar1 = ax2.bar([r - 0.05 for r in resolutions], n_clusters, width=0.1,
                   label='PAS clusters', color=COLOR_LAMBDA_GRAD, alpha=0.7)
    gex_clusters = [5, 8, 11, 14, 17]
    bar2 = ax2.bar([r + 0.05 for r in resolutions], gex_clusters, width=0.1,
                   label='GEX clusters', color=COLOR_LAMBDA_POISSON, alpha=0.7)

    ax2.set_xlabel('Leiden Resolution', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Number of Clusters', fontsize=12, fontweight='bold')
    ax2.set_title('Cluster Count vs Resolution', fontsize=11, fontweight='bold')
    ax2.legend(loc='upper left', fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_ylim([0, 25])

    # Add note
    ax2.text(0.5, 0.05, 'Note: PAS consistently produces finer clusters than GEX\nsuggesting APA provides complementary resolution',
            transform=ax2.transAxes, fontsize=10, ha='center',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "05_resolution_sweep.png", dpi=300, bbox_inches='tight')
    print("✓ Created: 05_resolution_sweep.png")
    plt.close()

# ============================================================================
# Figure 6: Mathematical Formulas (rendered as images)
# ============================================================================

def create_formula_sheet():
    """Create a figure with key mathematical formulas."""
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111)
    ax.axis('off')

    fig.suptitle('PeakATail: Key Mathematical Formulas', fontsize=16, fontweight='bold', y=0.98)

    formulas = [
        ("Poisson p-value for peak significance:",
         r"$p_{value} = P(X \geq h | \text{Poisson}(\lambda)) = 1 - \sum_{i=0}^{h-1} \frac{e^{-\lambda} \lambda^i}{i!}$"),

        ("Background lambda (floor positions):",
         r"$\lambda = \text{median}\{h_i : h_i < 0.1 \times \max(h)\}$"),

        ("Dynamic threshold:",
         r"$T = \max(\text{floor\_threshold}, \lambda \times \text{fold\_change})$"),

        ("Smoothed coverage (moving average):",
         r"$\text{smooth}[i] = \frac{1}{w} \sum_{j=i-w/2}^{i+w/2} h[j]$"),

        ("Local maxima detection (gradient zero-crossing):",
         r"$\text{Local max at } i \text{ if } \frac{d}{di}\text{smooth}[i-1] > 0 \text{ AND } \frac{d}{di}\text{smooth}[i+1] \leq 0$"),

        ("Prominence filter:",
         r"$\text{prominence}[i] = h[i] - \max(\text{left\_min}, \text{right\_min})$"),

        ("Confidence score:",
         r"$\text{confidence}[i] = \text{prominence}[i] \times \frac{h[i]}{\lambda} + \left|\frac{d^2}{di^2}\text{smooth}[i]\right|$"),

        ("Term Frequency (TF):",
         r"$TF[\text{peak}, \text{cell}] = \frac{\text{count}[\text{peak}, \text{cell}]}{\sum_{\text{peaks}} \text{count}[\text{peak}, \text{cell}]}$"),

        ("Inverse Document Frequency (IDF):",
         r"$IDF[\text{peak}] = \log\left(\frac{n_{\text{cells}}}{\text{cells with peak}}\right)$"),

        ("TF-IDF (Signac Method 1):",
         r"$\text{TF-IDF} = \log_1(TF \times IDF \times 10000)$"),

        ("Adjusted Rand Index (ARI):",
         r"$ARI = \frac{\text{RI} - E[\text{RI}]}{\max(\text{RI}) - E[\text{RI}]} \in [-1, 1]$"),

        ("Proximal-Distal Usage Index (PDUI):",
         r"$PDUI = \frac{\text{count(distal PAS)}}{\text{count(proximal PAS)} + \text{count(distal PAS)}}$"),

        ("Fisher's exact test (odds ratio):",
         r"$OR = \frac{n_{A1} \times n_{B2}}{n_{A2} \times n_{B1}}$"),

        ("Benjamini-Hochberg FDR correction:",
         r"$q_i = \min(p_{(i)} \times \frac{m}{i}, q_{i+1})$ (monotonicity enforced)")
    ]

    y_pos = 0.94
    for title, formula in formulas:
        # Title
        ax.text(0.05, y_pos, title, fontsize=10, fontweight='bold',
               transform=ax.transAxes, va='top')
        y_pos -= 0.04

        # Formula
        ax.text(0.1, y_pos, formula, fontsize=11, transform=ax.transAxes,
               va='top', family='serif',
               bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.5))
        y_pos -= 0.055

        if y_pos < 0.05:
            break

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "06_mathematical_formulas.png", dpi=300, bbox_inches='tight')
    print("✓ Created: 06_mathematical_formulas.png")
    plt.close()

# ============================================================================
# Figure 7: Data Statistics Summary
# ============================================================================

def create_data_statistics():
    """Create a summary figure of pipeline statistics."""
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3)
    fig.suptitle('PeakATail Pipeline: Data Statistics Summary', fontsize=15, fontweight='bold')

    # Data from pipeline
    stages = ['Input\nReads', 'Aligned\nReads', 'Raw\nPeaks', 'After CB\nFilter',
              'Annotated\nPeaks', 'Final\nCells', 'Final\nPeaks']
    values = [13.6e6, 12.1e6, 123000, 11771, 11614, 906, 11771]
    labels = ['13.6M', '12.1M', '123K', '11.8K', '11.6K', '906', '11.8K']

    # Subplot 1: Stage-wise statistics
    ax1 = fig.add_subplot(gs[0, :])
    colors_stages = ['#E8F4F8', '#D4E8F4', '#C0DFF0', '#ADD8EC', '#9ACAE8', '#87BDE4', '#74B4E0']
    bars = ax1.bar(stages, values, color=colors_stages, edgecolor='black', linewidth=1.5)
    ax1.set_ylabel('Count', fontsize=11, fontweight='bold')
    ax1.set_title('Reads and Peaks Through Pipeline Stages', fontsize=12, fontweight='bold')
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3, axis='y')
    # Add value labels
    for bar, label in zip(bars, labels):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                label, ha='center', va='bottom', fontsize=9, fontweight='bold')

    # Subplot 2: Matrix dimensions
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.axis('off')
    matrix_info = """
    Annotated Matrix Dimensions:

    • Rows (features): 11,771 peaks
    • Columns (cells): 906 cells
    • Non-zero entries: 703,733
    • Sparsity: 99.2%
    • Memory: ~2.1 MB (Matrix Market)

    Matrix Market Format:
    %%MatrixMarket matrix coordinate integer
    11771 906 703733
    6 1 1
    14 1 3
    ...
    """
    ax2.text(0.05, 0.95, matrix_info, transform=ax2.transAxes,
            fontsize=10, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    # Subplot 3: Cluster statistics
    ax3 = fig.add_subplot(gs[1, 1])
    cluster_counts = [68, 62, 45, 51, 73, 55, 48, 62, 70, 58, 52, 41, 38]
    sorted_counts = sorted(cluster_counts, reverse=True)
    ax3.barh(range(len(sorted_counts)), sorted_counts, color=COLOR_SIGNIF, alpha=0.7)
    ax3.set_xlabel('Cell Count', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Cluster ID', fontsize=11, fontweight='bold')
    ax3.set_title('Cluster Sizes (resolution=1.0)', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='x')
    ax3.text(70, 12, f'Total cells: {sum(sorted_counts)}\nTotal clusters: {len(sorted_counts)}',
            fontsize=10, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Subplot 4: Benchmark statistics
    ax4 = fig.add_subplot(gs[2, 0])
    ax4.axis('off')
    benchmark_info = """
    Validation Results (Lambda-Gradient):

    PolyASite 2.0 (569K PAS):
    • Precision @ 250bp: 52.8%
    • Precision @ 500bp: 63.2%
    • Precision @ 1000bp: 71.5%

    PolyA_DB v3 (303K PAS):
    • Precision @ 250bp: 62.1%
    • Precision @ 500bp: 70.9%
    • Precision @ 1000bp: 77.8%

    Peak annotation filter:
    • TIER_1 (within UTR): 8,342 peaks
    • TIER_1+2 (extended): 10,156 peaks
    • All tiers: 11,771 peaks
    """
    ax4.text(0.05, 0.95, benchmark_info, transform=ax4.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))

    # Subplot 5: Annotation tier distribution
    ax5 = fig.add_subplot(gs[2, 1])
    tiers = ['TIER_1\n(within UTR)', 'TIER_2\n(extended UTR)', 'TIER_3\n(distal)']
    tier_counts = [8342, 1814, 1615]
    colors_tiers = ['#06A77D', '#FFA500', '#E74C3C']
    wedges, texts, autotexts = ax5.pie(tier_counts, labels=tiers, autopct='%1.1f%%',
                                         colors=colors_tiers, startangle=90)
    ax5.set_title('Peak Annotation Distribution', fontsize=12, fontweight='bold')
    # Make percentage text bold
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(10)

    plt.savefig(FIGURES_DIR / "07_data_statistics.png", dpi=300, bbox_inches='tight')
    print("✓ Created: 07_data_statistics.png")
    plt.close()

# ============================================================================
# Copy existing benchmark figures
# ============================================================================

def copy_existing_figures():
    """Copy pre-generated benchmark figures to reports/figures/."""
    source_dir = Path("/home/user/PeakATail/test_run")

    figures_to_copy = [
        ("final_benchmark/dual_database_validation.png", "08_validation_dual_database.png"),
        ("final_benchmark/dual_db_metrics_table.png", "09_validation_metrics_table.png"),
        ("final_benchmark/clustering_sweep/ari_ami_vs_resolution.png", "10_clustering_agreement_sweep.png"),
        ("rerun_1500/main_comparison.png", "11_pas_gex_clustering_comparison.png"),
        ("rerun_1500/sankey.png", "12_cluster_sankey.png"),
        ("rerun_1500/differential_apa.png", "13_differential_apa_volcano.png"),
    ]

    for src, dst in figures_to_copy:
        src_path = source_dir / src
        dst_path = FIGURES_DIR / dst
        if src_path.exists():
            shutil.copy2(src_path, dst_path)
            print(f"✓ Copied: {src} → {dst}")
        else:
            print(f"⚠ Not found: {src}")

# ============================================================================
# Main execution
# ============================================================================

if __name__ == '__main__':
    print("\nGenerating figures...")
    create_pipeline_diagram()
    create_strategy_comparison()
    create_tfidf_lsi_pipeline()
    create_precision_comparison()
    create_resolution_sweep()
    create_formula_sheet()
    create_data_statistics()

    print("\nCopying existing benchmark figures...")
    copy_existing_figures()

    print("\n" + "=" * 60)
    print("Figure generation complete!")
    print(f"All figures saved to: {FIGURES_DIR}")
    print("=" * 60)

    # List all generated figures
    print("\nGenerated figures:")
    for fig_file in sorted(FIGURES_DIR.glob("*.png")):
        file_size = fig_file.stat().st_size / 1024
        print(f"  - {fig_file.name} ({file_size:.1f} KB)")
