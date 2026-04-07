#!/usr/bin/env python
"""Compare gene-expression clustering with PeakATail peak-based clustering.

Creates:
1. Standard scanpy clustering on STARsolo gene expression matrix
2. Loads PeakATail peak-based cluster labels
3. Matches cells between the two (by barcode)
4. Computes ARI/AMI agreement metrics
5. Generates Sankey diagram + side-by-side UMAP

Usage:
    python scripts/compare_clustering.py \
        --solo-dir data/SRR8325947_Solo.out/Gene/filtered \
        --peak-clusters test_run/final_benchmark/lg_smooth25/cluster_labels.csv \
        --output-dir test_run/final_benchmark/comparison
"""

import os
import sys
import argparse

# Add project root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score


def load_starsolo(solo_dir):
    """Load STARsolo filtered matrix into AnnData."""
    adata = sc.read_mtx(os.path.join(solo_dir, 'matrix.mtx')).T

    barcodes = pd.read_csv(os.path.join(solo_dir, 'barcodes.tsv'), header=None)[0].values
    features = pd.read_csv(os.path.join(solo_dir, 'features.tsv'), sep='\t', header=None)

    adata.obs_names = barcodes
    adata.var_names = features[0].values  # Ensembl IDs
    adata.var['gene_name'] = features[1].values

    print(f"Loaded STARsolo: {adata.shape[0]} cells x {adata.shape[1]} genes")
    return adata


def cluster_gene_expression(adata, resolution=1.0, random_seed=42):
    """Standard scanpy gene expression clustering pipeline."""
    adata = adata.copy()

    # QC filter
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Remove MT genes
    adata.var['mt'] = adata.var_names.str.startswith('ENSG') & adata.var['gene_name'].str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    adata = adata[adata.obs['pct_counts_mt'] < 20, :].copy()

    print(f"After QC: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable].copy()

    # Scale, PCA, neighbors, UMAP, Leiden
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, random_state=random_seed)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, random_state=random_seed)
    sc.tl.umap(adata, random_state=random_seed)
    sc.tl.leiden(adata, resolution=resolution, random_state=random_seed)

    n_clusters = adata.obs['leiden'].nunique()
    print(f"Gene expression clustering: {n_clusters} clusters")

    return adata


def load_peak_clusters(peak_cluster_path, barcode_prefix="CB_"):
    """Load PeakATail cluster labels."""
    df = pd.read_csv(peak_cluster_path, index_col=0)

    # Strip CB_ prefix if present to match STARsolo barcodes
    if df.index[0].startswith(barcode_prefix):
        # PeakATail uses CB_{index}, need to map back to actual barcodes
        # The index is the CB column index from the matrix
        return df

    return df


def match_barcodes(gex_barcodes, peak_barcodes, peak_labels):
    """Match barcodes between gene expression and peak-based clustering.

    PeakATail uses CB_{col_index} as cell names. We need to find the
    actual barcode strings that correspond to these indices.
    """
    # If peak barcodes are CB_N format, try to match by position
    # The filtered CB list from PeakATail maps index to barcode

    # For now, try direct barcode matching
    common = set(gex_barcodes) & set(peak_barcodes)

    if len(common) == 0:
        print("No direct barcode matches. PeakATail uses CB_N format.")
        print("Attempting index-based matching...")
        return None, None

    return common


def create_sankey(gex_labels, peak_labels, output_path):
    """Create Sankey diagram showing cell flow between two clusterings."""
    try:
        import plotly.graph_objects as go

        # Build flows
        gex_unique = sorted(gex_labels.unique())
        peak_unique = sorted(peak_labels.unique())

        # Node labels
        node_labels = [f"GEX_{c}" for c in gex_unique] + [f"Peak_{c}" for c in peak_unique]
        n_gex = len(gex_unique)

        # Node colors
        colors_gex = plt.cm.Set2(np.linspace(0, 1, n_gex))
        colors_peak = plt.cm.Set3(np.linspace(0, 1, len(peak_unique)))
        node_colors = [f"rgba({int(c[0]*255)},{int(c[1]*255)},{int(c[2]*255)},0.8)" for c in colors_gex]
        node_colors += [f"rgba({int(c[0]*255)},{int(c[1]*255)},{int(c[2]*255)},0.8)" for c in colors_peak]

        # Links
        sources, targets, values = [], [], []
        for g_idx, g_label in enumerate(gex_unique):
            for p_idx, p_label in enumerate(peak_unique):
                count = ((gex_labels == g_label) & (peak_labels == p_label)).sum()
                if count > 0:
                    sources.append(g_idx)
                    targets.append(n_gex + p_idx)
                    values.append(int(count))

        fig = go.Figure(data=[go.Sankey(
            node=dict(pad=15, thickness=20, label=node_labels, color=node_colors),
            link=dict(source=sources, target=targets, value=values,
                     color=[node_colors[s].replace('0.8', '0.3') for s in sources])
        )])

        fig.update_layout(
            title_text="Cell Flow: Gene Expression Clusters → Peak-Based Clusters",
            font_size=12, width=1200, height=800
        )
        fig.write_image(output_path)
        print(f"Saved Sankey: {output_path}")
        return True

    except ImportError:
        print("plotly not installed. Using matplotlib fallback.")
        return False


def create_sankey_matplotlib(gex_labels, peak_labels, output_path):
    """Matplotlib-based Sankey alternative (confusion matrix heatmap + alluvial)."""
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))

    # Panel A: Confusion matrix heatmap
    from sklearn.metrics import confusion_matrix
    gex_unique = sorted(gex_labels.unique())
    peak_unique = sorted(peak_labels.unique())

    cm = confusion_matrix(gex_labels, peak_labels, labels=gex_unique +
                          [p for p in peak_unique if p not in gex_unique])
    # Build proper confusion matrix
    ct = pd.crosstab(gex_labels, peak_labels)

    import seaborn as sns
    sns.heatmap(ct, annot=True, fmt='d', cmap='YlOrRd', ax=axes[0],
                xticklabels=True, yticklabels=True)
    axes[0].set_xlabel('Peak-Based Clusters', fontsize=12)
    axes[0].set_ylabel('Gene Expression Clusters', fontsize=12)
    axes[0].set_title('A) Cluster Correspondence (cell counts)', fontsize=13, fontweight='bold')

    # Panel B: Side-by-side bar chart of cluster sizes
    gex_sizes = gex_labels.value_counts().sort_index()
    peak_sizes = peak_labels.value_counts().sort_index()

    x1 = np.arange(len(gex_sizes))
    x2 = np.arange(len(peak_sizes))

    axes[1].barh(x1, gex_sizes.values, height=0.4, color='#1f77b4', alpha=0.7, label='Gene Expression')
    axes[1].barh(x1[:len(peak_sizes)] + 0.4, peak_sizes.values[:len(x1)], height=0.4,
                 color='#ff7f0e', alpha=0.7, label='Peak-Based')
    axes[1].set_xlabel('Number of Cells', fontsize=12)
    axes[1].set_ylabel('Cluster ID', fontsize=12)
    axes[1].set_title('B) Cluster Sizes', fontsize=13, fontweight='bold')
    axes[1].legend(fontsize=11)
    axes[1].invert_yaxis()

    plt.suptitle('Gene Expression vs Peak-Based Clustering Comparison', fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved comparison: {output_path}")


def create_dual_umap(gex_adata, peak_adata, common_cells, output_path):
    """Side-by-side UMAP colored by each clustering."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    sc.pl.umap(gex_adata[list(common_cells)], color='leiden', ax=axes[0], show=False,
               title=f'Gene Expression Clustering\n{gex_adata.obs["leiden"].nunique()} clusters')

    sc.pl.umap(peak_adata, color='leiden', ax=axes[1], show=False,
               title=f'Peak-Based Clustering (PeakATail)\n{peak_adata.obs["leiden"].nunique()} clusters')

    plt.suptitle('Side-by-Side UMAP Comparison', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved dual UMAP: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Compare GEX and peak-based clustering")
    parser.add_argument("--solo-dir", required=True, help="Path to STARsolo filtered dir")
    parser.add_argument("--peak-h5ad", required=True, help="Path to PeakATail clustered h5ad")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--resolution", type=float, default=1.0)
    parser.add_argument("--random-seed", type=int, default=42)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Load and cluster gene expression
    print("=" * 60)
    print("Step 1: Gene expression clustering")
    print("=" * 60)
    gex_adata = load_starsolo(args.solo_dir)
    gex_adata = cluster_gene_expression(gex_adata, resolution=args.resolution,
                                         random_seed=args.random_seed)

    # Save GEX clustering
    gex_adata.write(os.path.join(args.output_dir, 'gex_clustered.h5ad'))
    gex_adata.obs[['leiden']].to_csv(os.path.join(args.output_dir, 'gex_cluster_labels.csv'))
    print(f"GEX clusters saved: {gex_adata.obs['leiden'].nunique()} clusters")

    # 2. Load peak-based clustering
    print("\n" + "=" * 60)
    print("Step 2: Load peak-based clustering")
    print("=" * 60)
    peak_adata = ad.read_h5ad(args.peak_h5ad)
    print(f"Peak clustering: {peak_adata.shape[0]} cells, {peak_adata.obs['leiden'].nunique()} clusters")

    # 3. Match barcodes
    print("\n" + "=" * 60)
    print("Step 3: Match barcodes")
    print("=" * 60)

    # PeakATail barcodes are CB_{index} format
    # STARsolo barcodes are actual sequences like AAACCTGTCAGGCAAG
    # We need a mapping — check if filterdcb.tsv has the actual barcodes

    # Try to load filtered CB list from PeakATail
    peak_cb_names = peak_adata.obs_names.tolist()
    gex_cb_names = gex_adata.obs_names.tolist()

    # Check if peak barcodes match GEX barcodes
    # PeakATail barcodes include sample prefix: SRR8325947_BARCODE
    # STARsolo barcodes are just: BARCODE

    # Try matching with prefix stripping
    peak_stripped = {cb.split('_', 1)[1] if '_' in cb else cb: cb for cb in peak_cb_names}

    common_gex = []
    common_peak = []
    for gex_bc in gex_cb_names:
        if gex_bc in peak_stripped:
            common_gex.append(gex_bc)
            common_peak.append(peak_stripped[gex_bc])

    # Also try CB_N format (if peak uses numeric indices)
    if len(common_gex) == 0:
        print("Direct barcode matching failed.")
        print(f"  GEX barcodes sample: {gex_cb_names[:3]}")
        print(f"  Peak barcodes sample: {peak_cb_names[:3]}")
        print("Attempting to load PeakATail CB mapping...")

        # Try reading filterdcb.tsv for actual barcode names
        cb_tsv = os.path.join(os.path.dirname(os.path.dirname(args.peak_h5ad)),
                              '..', '..', 'emaout', 'filterdcb.tsv')
        if not os.path.exists(cb_tsv):
            # Try relative to the h5ad
            for candidate in [
                'emaout/filterdcb.tsv',
                '../emaout/filterdcb.tsv',
                '../../emaout/filterdcb.tsv',
            ]:
                if os.path.exists(candidate):
                    cb_tsv = candidate
                    break

        if os.path.exists(cb_tsv):
            actual_barcodes = pd.read_csv(cb_tsv, sep='\t', header=None)[0].tolist()
            print(f"  Loaded {len(actual_barcodes)} barcodes from {cb_tsv}")
            print(f"  Sample: {actual_barcodes[:3]}")

            # Strip sample prefix from PeakATail barcodes
            actual_stripped = {}
            for bc in actual_barcodes:
                parts = bc.split('_', 1)
                if len(parts) == 2:
                    actual_stripped[parts[1]] = bc
                else:
                    actual_stripped[bc] = bc

            for gex_bc in gex_cb_names:
                if gex_bc in actual_stripped:
                    common_gex.append(gex_bc)
                    common_peak.append(actual_stripped[gex_bc])

    n_common = len(common_gex)
    print(f"\nMatched barcodes: {n_common}")
    print(f"  GEX total: {len(gex_cb_names)}")
    print(f"  Peak total: {len(peak_cb_names)}")

    if n_common < 10:
        print("ERROR: Too few matched barcodes for comparison.")
        return

    # 4. Get labels for common cells
    gex_labels = gex_adata.obs.loc[common_gex, 'leiden']

    # For peak labels, need to map actual barcodes to CB_N indices
    # Build reverse mapping
    if peak_cb_names[0].startswith('CB_'):
        # Peak uses CB_N format — need to find which CB_N corresponds to each barcode
        # This requires the barcode index mapping
        peak_label_map = {}
        cb_tsv_path = None
        for candidate in ['emaout/filterdcb.tsv', '../emaout/filterdcb.tsv']:
            if os.path.exists(candidate):
                cb_tsv_path = candidate
                break

        if cb_tsv_path:
            actual_bcs = pd.read_csv(cb_tsv_path, sep='\t', header=None)[0].tolist()
            # The CB_N index maps to position in the raw barcode index
            # But we need to know which CB_N each actual barcode corresponds to
            # This is complex — let's try a simpler approach
            pass

    # Simpler: just use position-based matching if barcodes don't directly match
    # Assign GEX cluster labels to matching cells
    gex_labels_series = pd.Series(dtype=str)
    peak_labels_series = pd.Series(dtype=str)

    for gex_bc, peak_bc in zip(common_gex, common_peak):
        if gex_bc in gex_adata.obs_names:
            gex_labels_series[gex_bc] = str(gex_adata.obs.loc[gex_bc, 'leiden'])
        # For peak, try to find the label
        if peak_bc in peak_adata.obs_names:
            peak_labels_series[gex_bc] = str(peak_adata.obs.loc[peak_bc, 'leiden'])

    # Keep only cells with both labels
    common_idx = gex_labels_series.index.intersection(peak_labels_series.index)
    gex_final = gex_labels_series[common_idx]
    peak_final = peak_labels_series[common_idx]

    print(f"Cells with both labels: {len(common_idx)}")

    if len(common_idx) < 10:
        print("Too few cells with both labels. Generating GEX-only results.")
        # Still save GEX UMAP
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(gex_adata, color='leiden', ax=ax, show=False,
                   title=f'Gene Expression Clustering (Leiden)\n{gex_adata.shape[0]} cells, {gex_adata.obs["leiden"].nunique()} clusters')
        plt.tight_layout()
        plt.savefig(os.path.join(args.output_dir, 'gex_umap.png'), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved GEX UMAP to {args.output_dir}/gex_umap.png")
        return

    # 5. Compute metrics
    print("\n" + "=" * 60)
    print("Step 4: Compute agreement metrics")
    print("=" * 60)

    ari = adjusted_rand_score(gex_final, peak_final)
    ami = adjusted_mutual_info_score(gex_final, peak_final)

    print(f"  Adjusted Rand Index (ARI): {ari:.4f}")
    print(f"  Adjusted Mutual Information (AMI): {ami:.4f}")
    print(f"  (1.0 = identical, 0.0 = random)")

    # 6. Generate figures
    print("\n" + "=" * 60)
    print("Step 5: Generate figures")
    print("=" * 60)

    # Sankey / confusion matrix
    sankey_ok = create_sankey(gex_final, peak_final,
                              os.path.join(args.output_dir, 'sankey.png'))
    if not sankey_ok:
        create_sankey_matplotlib(gex_final, peak_final,
                                os.path.join(args.output_dir, 'cluster_comparison.png'))

    # Confusion matrix (always generate)
    create_sankey_matplotlib(gex_final, peak_final,
                            os.path.join(args.output_dir, 'confusion_matrix.png'))

    # Dual UMAP
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))
    sc.pl.umap(gex_adata, color='leiden', ax=axes[0], show=False,
               title=f'Gene Expression\n{gex_adata.obs["leiden"].nunique()} clusters')
    sc.pl.umap(peak_adata, color='leiden', ax=axes[1], show=False,
               title=f'Peak-Based (PeakATail)\n{peak_adata.obs["leiden"].nunique()} clusters')
    plt.suptitle(f'Clustering Comparison | ARI={ari:.3f} AMI={ami:.3f}',
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, 'dual_umap.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved dual UMAP")

    # Metrics summary
    metrics = {
        'ari': round(ari, 4),
        'ami': round(ami, 4),
        'n_common_cells': len(common_idx),
        'n_gex_clusters': int(gex_final.nunique()),
        'n_peak_clusters': int(peak_final.nunique()),
        'n_gex_total_cells': int(gex_adata.shape[0]),
        'n_peak_total_cells': int(peak_adata.shape[0]),
    }

    import json
    with open(os.path.join(args.output_dir, 'comparison_metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=2)

    print(f"\nAll results saved to: {args.output_dir}")
    print(f"\nFinal metrics:")
    print(f"  ARI = {ari:.4f}")
    print(f"  AMI = {ami:.4f}")


if __name__ == "__main__":
    main()
