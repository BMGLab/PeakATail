# PeakATail Development & Publication Roadmap

## Executive Summary

This roadmap outlines the path from current PeakATail implementation to a publishable, competitive APA analysis tool. Based on comprehensive analysis of 19 existing tools, we've identified critical improvements and a validated publication strategy.

**Key Insight**: PeakATail has confirmed novelty (first tool for pure peak-based unsupervised clustering), but needs statistical rigor and validation to be publication-ready.

---

# PART 1: CODE IMPROVEMENTS

## Phase 1: Critical Features (REQUIRED for Publication)

### Priority: CRITICAL | Effort: Medium

These features are present in 12-15 out of 19 tools and are expected by reviewers:

### 1.1 FDR Correction
**Status**: Missing | **Impact**: CRITICAL | **Effort**: LOW

**Implementation**:
```python
from scipy.stats import false_discovery_control

# After peak calling, add p-value calculation
peak_pvalues = calculate_peak_significance(peaks)  # Poisson model
peak_qvalues = false_discovery_control(peak_pvalues, method='bh')

# Filter peaks
significant_peaks = peaks[peak_qvalues < 0.05]
```

**Files to modify**:
- `ema/countmatrix/peackcalling.py`
- Add new module: `ema/statistics/fdr.py`

**Reference implementations**:
- scMAPA: Benjamini-Hochberg FDR
- Sierra: DEXSeq FDR correction
- SCAPE: BIC with FDR

---

### 1.2 Statistical Peak Significance
**Status**: Missing | **Impact**: CRITICAL | **Effort**: MEDIUM

**Implementation Strategy** (from MACS2/SCAPE):
```python
def calculate_peak_pvalue(peak_height, background_lambda, peak_width):
    """
    Poisson-based p-value calculation for peak significance

    Args:
        peak_height: Number of reads in peak
        background_lambda: Expected background (local lambda)
        peak_width: Peak window size

    Returns:
        p-value: Probability of observing this peak by chance
    """
    from scipy.stats import poisson
    expected = background_lambda * peak_width
    pval = 1 - poisson.cdf(peak_height - 1, expected)
    return pval
```

**Background Estimation Methods**:
1. Local lambda: Median coverage in ±10kb window
2. Global lambda: Genome-wide average coverage
3. Dynamic lambda: Cell-type specific background

**Files to modify**:
- `ema/countmatrix/peackcalling.py` - Add p-value calculation
- New: `ema/statistics/peak_significance.py`

**Reference implementations**:
- scMAPA: Logistic regression + LRT
- MAAPER: Probabilistic model with EM
- InPACT: Random Forest confidence scores

---

### 1.3 Internal Priming Filter
**Status**: Missing | **Impact**: HIGH | **Effort**: MEDIUM

**Implementation** (from SAPAS/scraps/SCAPTURE):
```python
def filter_internal_priming(peaks, genome_fasta, window=(-5, 20)):
    """
    Filter peaks near genomic A-rich stretches (internal priming artifacts)

    Standard: ≥6 consecutive A's within -5 to +20 bp of PAS
    """
    from pyfaidx import Fasta

    genome = Fasta(genome_fasta)
    filtered_peaks = []

    for peak in peaks:
        # Extract sequence around peak
        seq = genome[peak.chrom][peak.pos + window[0]:peak.pos + window[1]]

        # Check for A-rich stretches
        if 'AAAAAA' not in str(seq).upper():
            filtered_peaks.append(peak)
        else:
            peak.internal_priming = True  # Flag for review

    return filtered_peaks
```

**Files to modify**:
- `ema/countmatrix/pasfind.py` - Add filtering step
- New: `ema/filters/internal_priming.py`
- Add dependency: `pyfaidx` for FASTA access

**Reference implementations**:
- SAPAS: ≥6 A's genomic check
- scraps: A-content analysis
- scAPAtrap: 6 continuous A's within -140 to +10 bp
- scPAISO: 6-mer AAAAAA within -5/+20 bp

---

## Phase 2: Important Enhancements (Strengthen Publication)

### Priority: HIGH | Effort: Medium-High

### 2.1 PDUI Calculation
**Status**: Missing | **Impact**: MEDIUM | **Effort**: LOW

**Implementation** (from scMAPA/scDaPars):
```python
def calculate_pdui(gene_peaks, count_matrix):
    """
    Calculate Proximal-Distal Usage Index

    PDUI = Distal_counts / (Proximal_counts + Distal_counts)

    Returns: Value between 0 (all proximal) and 1 (all distal)
    """
    proximal_peak = gene_peaks[0]  # Assume sorted by position
    distal_peak = gene_peaks[-1]

    proximal_counts = count_matrix[proximal_peak.id].sum()
    distal_counts = count_matrix[distal_peak.id].sum()

    if proximal_counts + distal_counts == 0:
        return np.nan

    pdui = distal_counts / (proximal_counts + distal_counts)
    return pdui
```

**Files to modify**:
- New: `ema/quantification/pdui.py`
- `ema/main.py` - Add optional PDUI output

**Output**: Gene-level PDUI matrix for comparison with other tools

---

### 2.2 Read1 Utilization
**Status**: Missing | **Impact**: MEDIUM-HIGH | **Effort**: MEDIUM

**Motivation**: scTail and scPAISO (2025) use Read1 for precise cleavage site detection

**Implementation** (from scTail/scPAISO):
```python
def extract_read1_cleavage_sites(bam_file):
    """
    Use Read1 5' end as precise mRNA cleavage site

    For paired-end 3' tag sequencing:
    - Read1: Captures exact poly(A) site
    - Read2: Used for gene assignment
    """
    import pysam

    cleavage_sites = defaultdict(list)

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam:
            if read.is_read1 and not read.is_unmapped:
                # Use 5' end of Read1 as cleavage site
                if read.is_reverse:
                    pos = read.reference_end
                else:
                    pos = read.reference_start

                cleavage_sites[read.reference_name].append(pos)

    return cleavage_sites
```

**Benefits**:
- Sharper peaks (95% < 69bp width from scPAISO)
- Stronger AAUAAA motif enrichment
- >95% isoform assignment accuracy

**Files to modify**:
- `ema/countmatrix/read.py` - Add Read1 support
- `ema/countmatrix/peackcalling.py` - Use Read1 for peak calling
- Add flag: `--use-read1` for optional Read1-based detection

---

### 2.3 K-NN Imputation
**Status**: Missing | **Impact**: MEDIUM | **Effort**: MEDIUM

**Implementation** (from stAPAminer/scDaPars):
```python
def knn_impute_apa(apa_matrix, gene_matrix, k=10, max_iter=10):
    """
    K-NN imputation guided by gene expression

    For each cell with missing/sparse APA values:
    - Find k nearest neighbors based on gene expression
    - Impute missing APA values as mean of neighbors
    """
    from sklearn.neighbors import NearestNeighbors
    from sklearn.preprocessing import StandardScaler

    # Standardize gene expression
    gene_scaled = StandardScaler().fit_transform(gene_matrix.T)

    # Find k-nearest neighbors
    knn = NearestNeighbors(n_neighbors=k)
    knn.fit(gene_scaled)
    distances, indices = knn.kneighbors(gene_scaled)

    # Iterative imputation
    apa_imputed = apa_matrix.copy()
    for iteration in range(max_iter):
        for cell_idx in range(apa_imputed.shape[1]):
            neighbor_indices = indices[cell_idx, 1:]  # Exclude self

            # Impute missing values
            missing_mask = np.isnan(apa_imputed[:, cell_idx])
            if missing_mask.sum() > 0:
                neighbor_values = apa_imputed[missing_mask][:, neighbor_indices]
                apa_imputed[missing_mask, cell_idx] = np.nanmean(neighbor_values, axis=1)

    return apa_imputed
```

**Files to modify**:
- New: `ema/imputation/knn.py`
- `ema/clustering/clustering.py` - Add optional imputation step

---

## Phase 3: Advanced Features (Optional, High Impact)

### Priority: LOW-MEDIUM | Effort: HIGH

### 3.1 Deep Learning Validation
**Status**: Missing | **Impact**: MEDIUM | **Effort**: HIGH

**Options**:
1. **Use pre-trained DeepPASS** (from SCAPTURE)
   - CNN model for PAS validation
   - Filters false positives
   - Requires sequence context (±200bp)

2. **Train custom CNN** (like scTail/scAPAtrap)
   - Binary classifier: Real PAS vs False positive
   - Features: Coverage, sequence motif, conservation

**Implementation considerations**:
- High effort, moderate impact
- Consider as Phase 3 after core features
- May not be critical for publication if statistical tests are robust

---

### 3.2 Visualization Suite
**Status**: Basic UMAP only | **Impact**: LOW | **Effort**: HIGH

**Modules to add** (from vizAPA):
```python
# Module 1: Genome Browser Tracks
plot_genomic_tracks(gene, peaks, coverage, cell_types)

# Module 2: Statistical Plots
plot_violin(peak_usage, cell_types)
plot_heatmap(differential_peaks, cluster_annotation)

# Module 3: Marker Visualization
plot_apa_markers(top_markers, expression_matrix)
```

**Files to create**:
- `ema/visualization/tracks.py`
- `ema/visualization/stats_plots.py`
- `ema/visualization/markers.py`

**Priority**: LOW - Nice to have, not critical for novelty claim

---

# PART 2: VALIDATION WITH REFERENCE DATABASES

## Database-Based Peak Validation

### Validation Strategy Overview

Compare PeakATail-identified peaks against gold-standard poly(A) site databases to assess:
1. **Precision**: % of PeakATail peaks that overlap known PAS
2. **Recall**: % of known PAS that are detected by PeakATail
3. **Novelty**: % of PeakATail peaks that are novel (not in databases)
4. **False Positive Rate**: Estimated from internal priming and other filters

---

## Primary Reference Databases

### 1. PolyASite 2.0 / 3.0
**Source**: https://polyasite.unibas.ch/
**Coverage**: Human, mouse, C. elegans
**Data**: Consolidated atlas from 3' end sequencing datasets
**Latest**: PolyASite v3.0 (January 2025) - includes scRNA-seq datasets

**Download**:
```bash
# BED format files available for bulk download
wget https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz
wget https://polyasite.unibas.ch/download/atlas/2.0/mm10/atlas.clusters.2.0.mm10.bed.gz
```

**File Format**:
```
chr1  1000  1050  ENSG00000123456  100  +  ...
```

**Usage in Validation**:
```python
def validate_against_polyasite(peakatail_peaks, polyasite_bed, window=50):
    """
    Compare PeakATail peaks with PolyASite clusters

    Args:
        peakatail_peaks: BED file of PeakATail identified peaks
        polyasite_bed: PolyASite reference clusters
        window: Allowed distance for overlap (default: 50bp)

    Returns:
        precision: TP / (TP + FP)
        recall: TP / (TP + FN)
        f1_score: 2 * (precision * recall) / (precision + recall)
        novel_peaks: Peaks not in PolyASite
    """
    from pybedtools import BedTool

    peakatail = BedTool(peakatail_peaks)
    polyasite = BedTool(polyasite_bed)

    # Find overlaps within window
    overlaps = peakatail.window(polyasite, w=window)

    true_positives = len(set([x.name for x in overlaps]))
    false_positives = len(peakatail) - true_positives

    # Check which PolyASite sites we detected
    reverse_overlaps = polyasite.window(peakatail, w=window)
    detected_known = len(set([x.name for x in reverse_overlaps]))
    false_negatives = len(polyasite) - detected_known

    precision = true_positives / (true_positives + false_positives)
    recall = true_positives / (true_positives + false_negatives)
    f1 = 2 * (precision * recall) / (precision + recall)

    # Identify novel peaks
    novel = peakatail.window(polyasite, w=window, v=True)

    return {
        'precision': precision,
        'recall': recall,
        'f1_score': f1,
        'true_positives': true_positives,
        'false_positives': false_positives,
        'false_negatives': false_negatives,
        'novel_peaks': len(novel)
    }
```

**Expected Performance** (from other tools):
- scAPAtrap: ~75-85% overlap with PolyASite
- SCAPTURE: ~70-80% overlap
- scAPA: ~65-75% overlap

**Reference**: [PolyASite 2.0 Paper](https://academic.oup.com/nar/article/48/D1/D174/5588346)

---

### 2. PolyA_DB v4
**Source**: https://exon.apps.wistar.org/polya_db/v3/ (v4 in development)
**Coverage**: Human, mouse, rat, chicken
**Data**: Systematic PAS identification using 3' end and long-read sequencing
**Release**: v4 published 2024/2025

**Download**:
```bash
# SAF format (also available as BED)
wget https://exon.apps.wistar.org/polya_db/v3/download/3.2/human_pas.saf
wget https://exon.apps.wistar.org/polya_db/v3/download/3.2/mouse_pas.saf
```

**Unique Features**:
- Includes isoform annotation
- Classified by genomic region (3'UTR, exon, intron, intergenic)
- Deep sequencing validation

**Usage in Validation**:
```python
def validate_with_polya_db(peakatail_peaks, polya_db_saf):
    """
    Validate peaks against PolyA_DB

    Additional metrics:
    - % peaks in 3'UTR vs other regions
    - Comparison with isoform annotations
    """
    # Convert SAF to BED
    polya_bed = convert_saf_to_bed(polya_db_saf)

    # Standard overlap analysis
    stats = validate_against_polyasite(peakatail_peaks, polya_bed)

    # Additional: Region classification
    region_stats = classify_peak_regions(peakatail_peaks, polya_db_saf)

    return {**stats, 'region_distribution': region_stats}
```

**Reference**: [PolyA_DB v4 Paper](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkaf1212/8356008)

---

### 3. APAeval Benchmark
**Source**: https://www.apaeval.org/ (OpenEBench platform)
**Purpose**: Community benchmark for APA tools
**Data**: Real, synthetic, and matched 3' end sequencing

**Benchmark Datasets**:
1. **Synthetic data**: Known ground truth PAS
2. **Matched 3' seq**: Bulk RNA-seq + 3' end seq from same samples
3. **Real data**: Public datasets with validation

**Download**:
```bash
# Access via OpenEBench
# Includes 17 tools comparison
# 8 tools benchmarked in detail
```

**Usage**:
- Run PeakATail on APAeval datasets
- Submit results to OpenEBench platform
- Compare against 17 existing tools
- Get automated metrics and rankings

**Metrics Provided**:
- Peak identification accuracy
- Quantification correlation
- Tool runtime and memory
- False discovery rate

**Reference**: [APAeval Paper](https://link.springer.com/article/10.1186/s13059-021-02502-z)

---

## Secondary Validation Resources

### 4. GENCODE Annotation
**Source**: https://www.gencodegenes.org/
**Purpose**: Validate peaks near annotated transcript ends

**Usage**:
```python
def validate_with_gencode(peakatail_peaks, gencode_gtf):
    """
    Check if peaks are near annotated transcript 3' ends

    Expected: Most peaks within 500bp of annotated ends
    """
    # Extract 3' UTR regions from GTF
    three_utr_regions = extract_3utr_from_gtf(gencode_gtf)

    # Check proximity
    proximity_stats = check_peak_proximity(peakatail_peaks, three_utr_regions)

    return proximity_stats
```

---

### 5. Published Dataset Ground Truth

From tools that published their detected peaks:

#### scAPAtrap Published Peaks
- **Datasets**: Mouse spermatogenesis, Arabidopsis roots
- **Available**: GitHub repo has example outputs
- **Use**: Direct comparison on same dataset

#### Sierra Published Peaks
- **Dataset**: Mouse cardiac injury
- **Available**: Peak coordinates in supplementary data
- **Use**: DTU comparison

#### SCAPE Simulated Data
- **Dataset**: In `/simulation/` directory
- **Ground Truth**: Known simulated PAS positions
- **Use**: Perfect accuracy benchmark

---

## Comprehensive Validation Pipeline

### Step-by-Step Validation Protocol

```python
# validation_pipeline.py

def comprehensive_peak_validation(peakatail_peaks_bed):
    """
    Complete validation pipeline for PeakATail peaks
    """
    results = {}

    # 1. PolyASite validation
    print("Validating against PolyASite 2.0...")
    results['polyasite'] = validate_against_polyasite(
        peakatail_peaks_bed,
        'references/polyasite_2.0_GRCh38.bed'
    )

    # 2. PolyA_DB validation
    print("Validating against PolyA_DB v4...")
    results['polya_db'] = validate_with_polya_db(
        peakatail_peaks_bed,
        'references/polya_db_v4_human.saf'
    )

    # 3. Internal priming check
    print("Checking internal priming artifacts...")
    results['internal_priming'] = check_internal_priming(
        peakatail_peaks_bed,
        'references/GRCh38.fa'
    )

    # 4. GENCODE proximity
    print("Checking GENCODE annotation proximity...")
    results['gencode_proximity'] = validate_with_gencode(
        peakatail_peaks_bed,
        'references/gencode.v44.annotation.gtf'
    )

    # 5. Motif enrichment
    print("Checking poly(A) signal motif enrichment...")
    results['motif_enrichment'] = check_polya_motifs(
        peakatail_peaks_bed,
        'references/GRCh38.fa'
    )

    # 6. Novel peak characterization
    print("Characterizing novel peaks...")
    results['novel_peaks'] = characterize_novel_peaks(
        peakatail_peaks_bed,
        results['polyasite'],
        results['polya_db']
    )

    return results
```

---

## Validation Metrics & Expected Performance

### Target Metrics (Based on Published Tools)

| Metric | scAPAtrap | SCAPTURE | Sierra | PeakATail Target |
|--------|-----------|----------|--------|------------------|
| **Precision** (vs PolyASite) | 75-85% | 70-80% | 65-75% | >70% |
| **Recall** (vs PolyASite) | 60-70% | 55-65% | 70-80% | >60% |
| **F1-Score** | 0.67-0.77 | 0.62-0.72 | 0.67-0.77 | >0.65 |
| **Internal Priming Filtered** | ~10-15% | ~8-12% | N/A | <20% |
| **AAUAAA Motif Enrichment** | ~65% | ~70% | ~60% | >60% |
| **Novel Peaks** | ~20-30% | ~25-35% | ~15-25% | 20-40% acceptable |

---

## Validation Figures for Publication

### Figure: Peak Validation Analysis

**Panel A**: Overlap with Reference Databases
```
Venn diagram showing:
- PeakATail peaks
- PolyASite clusters
- PolyA_DB sites
- Overlap statistics
```

**Panel B**: Distance Distribution
```
Histogram of distance to nearest known PAS:
- X-axis: Distance (bp)
- Y-axis: Number of peaks
- Shows most peaks within 50bp of known sites
```

**Panel C**: Precision-Recall Curve
```
Multiple windows (10, 25, 50, 100, 200 bp):
- Trade-off between precision and recall
- Compare with other tools
```

**Panel D**: Regional Distribution
```
Bar plot showing peak distribution:
- 3'UTR: ~60-70%
- Exonic: ~10-15%
- Intronic: ~10-15%
- Intergenic: ~5-10%
```

**Panel E**: Motif Enrichment
```
Logo plot of sequences ±40bp from peaks:
- AAUAAA signal enrichment
- Downstream U-rich element
- Compare real vs shuffled peaks
```

---

## Validation Success Criteria

### Minimum Requirements for Publication

- [ ] **Precision ≥ 70%** against PolyASite within 50bp
- [ ] **Recall ≥ 60%** of PolyASite sites detected
- [ ] **F1-score ≥ 0.65**
- [ ] **Internal priming ≤ 20%** of peaks flagged
- [ ] **AAUAAA enrichment ≥ 60%** in peak vicinity
- [ ] **Novel peaks characterized**: Evidence they are real (motif, conservation, etc.)

### Strong Publication Metrics

- [ ] **Precision ≥ 80%**
- [ ] **Recall ≥ 70%**
- [ ] **F1-score ≥ 0.75**
- [ ] Competitive with or better than scAPAtrap/SCAPTURE
- [ ] Novel peaks validated in independent dataset

---

## Implementation Priority

### Phase 1: Basic Validation
1. Download PolyASite 2.0 (human + mouse)
2. Implement overlap analysis script
3. Calculate precision, recall, F1
4. Generate validation figures

### Phase 2: Comprehensive Validation
1. Download PolyA_DB v4
2. Cross-validate with both databases
3. Internal priming analysis
4. Motif enrichment analysis

### Phase 3: Advanced Validation
1. APAeval benchmark submission
2. Compare against published tool outputs
3. Novel peak characterization
4. Simulated data validation (if available)

---

# PART 3: BENCHMARKING & BIOLOGICAL VALIDATION

## Benchmark Datasets (From Other Tools)

### Dataset 1: Mouse Spermatogenesis (GSE104556)
- **Used by**: scAPAtrap, scAPA
- **Organism**: Mouse
- **Cells**: ~35,000 cells
- **Why**: Well-characterized developmental stages, known APA dynamics
- **Validation**: Can reproduce published findings
- **Download**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556

### Dataset 2: Human hESC Differentiation (GSE75748)
- **Used by**: scDaPars
- **Organism**: Human
- **Timepoints**: 0, 12, 24, 36, 72, 96 hours
- **Cells**: ~739 cells (after QC)
- **Why**: Dynamic APA changes during differentiation, trajectory analysis
- **Download**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75748

### Dataset 3: Human PBMC 10x Genomics
- **Used by**: SCAPTURE, Sierra
- **Organism**: Human
- **Samples**: Multiple (3k, 4k, 5k, 6k, 8k, 10k cells)
- **Why**: Public, well-annotated, multiple cell types
- **Download**: https://www.10xgenomics.com/datasets

### Dataset 4: Mouse Brain (Zeisel et al.)
- **Used by**: scMAPA
- **Organism**: Mouse
- **Cells**: 3 clusters (Neurons, Immune, Oligos)
- **Why**: Clear cell type boundaries, established markers
- **DOI**: https://doi.org/10.1016/j.cell.2018.06.021

---

## Benchmark Experiments

### Experiment 1: Tool-to-Tool Comparison
**Compare PeakATail against**:
- scLAPA (multimodal baseline)
- Sierra (detection quality)
- scAPAtrap (peak precision)

**Metrics**:
- Peak detection: Precision, Recall, F1-score (vs PolyASite)
- Clustering agreement: Adjusted Rand Index (ARI), Normalized Mutual Information (NMI)
- Computational: Runtime, Memory usage, Scalability

**Expected outcome**:
- PeakATail clusters differently (novelty)
- Similar or better peak detection
- Comparable computational efficiency

---

### Experiment 2: Peak-Based vs Gene-Based Clustering
**Critical experiment for novelty claim**:

```
Same dataset → Two approaches:

A) Standard (Gene Expression):
   Seurat → Gene counts → PCA → Louvain → Clusters_gene

B) PeakATail (Peak Counts):
   PeakATail → Peak counts → PCA → Louvain → Clusters_peak

Compare:
- Cluster overlap (ARI, NMI)
- Unique clusters in each approach
- Cell populations visible only in peak-based clustering
```

**Success criteria**:
- ARI < 0.8 (shows meaningful difference)
- Identify at least 1 cell population visible only in peak-based clustering
- Biological validation of unique clusters

---

### Experiment 3: Database Validation

Using PolyASite and PolyA_DB:

**Metrics**:
```
For each dataset:
- Precision = Correctly identified peaks / Total PeakATail peaks
- Recall = Detected known sites / Total known sites
- F1-score = Harmonic mean
- Novel peak rate = Peaks not in databases
```

**Compare with other tools**:
- Run same validation on scAPAtrap, SCAPTURE, Sierra outputs
- Show PeakATail is competitive

**Additional**:
- Internal priming rate
- Motif enrichment (AAUAAA)
- Regional distribution (3'UTR, exonic, intronic)

---

### Experiment 4: Biological Validation

#### 4A. Marker Discovery
- Identify peak markers (not gene markers)
- Validate with known APA biology
- Pathway enrichment of cluster-specific peaks

#### 4B. Developmental Trajectory
- Use hESC differentiation dataset
- Show APA trajectory distinct from expression trajectory
- 3'UTR lengthening/shortening during differentiation

#### 4C. Cell Type Specificity
- Use PBMC or brain dataset
- Identify cell-type specific peaks
- Compare with known cell type markers

---

### Experiment 5: Statistical Validation

After adding FDR correction:
```
- Show p-value distribution (uniform under null)
- QQ-plot of observed vs expected p-values
- FDR calibration (estimated FDR vs true FDR)
- Sensitivity to threshold (ROC curve if ground truth available)
```

---

# PART 4: PUBLICATION FIGURES

Based on analysis of 17 tool publications, here are the standard figures:

## Main Figures (Required)

### Figure 1: Method Overview
```
Panel A: PeakATail workflow diagram
  - BAM input → Peak calling → Count matrix → Clustering → APA switch

Panel B: Comparison with existing approaches
  - Traditional: Gene → Cluster → APA analysis
  - PeakATail: Peak → Cluster (NOVEL)

Panel C: Example gene with multiple peaks
  - Genomic track showing 2-3 peaks per gene
  - Different cells using different peaks
```

### Figure 2: Database Validation
```
Panel A: Venn diagram (PeakATail vs PolyASite vs PolyA_DB)
Panel B: Distance distribution to known PAS
Panel C: Precision-Recall curves (multiple window sizes)
Panel D: Regional distribution (3'UTR, exonic, intronic, intergenic)
Panel E: Motif enrichment (AAUAAA logo plot)
```

### Figure 3: Peak-Based vs Gene-Based Clustering
```
Panel A: UMAP of gene-based clustering (standard)
Panel B: UMAP of peak-based clustering (PeakATail)
Panel C: Confusion matrix / Sankey diagram showing cluster overlap
Panel D: Unique cell populations found only in peak-based clustering
```

### Figure 4: Biological Validation
```
Panel A: Heatmap of differential peak usage across clusters
Panel B: Example genes with cluster-specific APA
Panel C: Pathway enrichment of cluster-specific peaks
Panel D: Known marker validation
```

### Figure 5: Benchmarking Results
```
Panel A: Tool comparison (PeakATail vs scLAPA vs Sierra vs scAPAtrap)
  - Peak detection metrics (Precision, Recall, F1)
  - Runtime and memory usage

Panel B: Statistical validation
  - P-value distribution
  - FDR calibration

Panel C: Internal priming filter performance
  - % peaks filtered
  - Motif enrichment before/after
```

## Supplementary Figures

### Supp Fig 1: Dataset Descriptions
- Summary of all datasets used
- Cell count, tissue type, organism
- UMAP of each dataset colored by known annotations

### Supp Fig 2: Peak Calling Quality
- Peak width distribution
- Coverage distribution
- Number of peaks per gene

### Supp Fig 3: Clustering Quality
- Silhouette scores
- Within-cluster vs between-cluster distances
- Stability analysis (subsampling)

### Supp Fig 4: Additional Validations
- Multiple datasets
- Reproducibility across replicates
- Parameter sensitivity analysis

### Supp Fig 5: Case Studies
- Specific genes with interesting APA patterns
- Cell type-specific examples
- Developmental trajectory examples

---

# PART 5: MANUSCRIPT STRUCTURE

## Title
"PeakATail: De Novo Cell Type Discovery Through Peak-Based Clustering of Alternative Polyadenylation Sites"

## Abstract Structure
```
Background: APA is key regulator, single-cell analysis challenging
Gap: Existing tools require gene expression or pre-computed clusters
Method: PeakATail - first tool for pure peak-based unsupervised clustering
Validation: Validated against PolyASite and PolyA_DB with >70% precision
Results: Identifies X novel cell populations in Y datasets
Conclusion: Peak-based clustering reveals biology invisible to gene expression
```

## Introduction
1. APA biology and importance
2. Single-cell APA challenges and existing approaches
3. Gap: No tool clusters directly on peaks without gene expression
4. Our solution: PeakATail with database validation
5. Key contributions and findings

## Methods
- PeakATail algorithm
- Statistical framework (Poisson, FDR)
- Internal priming filter
- Database validation protocol
- Datasets used
- Benchmarking approach
- Statistical tests

## Results
1. PeakATail identifies peaks with high precision (validation vs PolyASite/PolyA_DB)
2. Peak-based clustering differs from gene-based clustering
3. Novel cell populations discovered in [dataset]
4. Biological validation of unique clusters
5. Benchmarking shows competitive performance
6. Case studies

## Discussion
- Novelty and significance
- Importance of database validation
- When to use peak-based vs gene-based clustering
- Limitations and future directions
- Broader implications

---

# PART 6: SUCCESS CRITERIA

## Technical Milestones
- [ ] FDR correction implemented and tested
- [ ] Statistical p-values for all peaks
- [ ] Internal priming filter validated
- [ ] Database validation complete (PolyASite + PolyA_DB)
  - [ ] Precision ≥ 70%
  - [ ] Recall ≥ 60%
  - [ ] F1-score ≥ 0.65
- [ ] Benchmarking on 4 datasets complete
- [ ] All unit tests passing
- [ ] Documentation complete

## Publication Milestones
- [ ] Identified ≥1 novel cell population in real data
- [ ] Peak-based clustering shows ARI < 0.8 vs gene-based
- [ ] Benchmarking shows competitive or better performance vs 3 tools
- [ ] Database validation shows ≥70% precision
- [ ] Biological validation of unique findings
- [ ] All figures publication-ready
- [ ] Manuscript draft complete

## Impact Milestones
- [ ] Preprint posted (bioRxiv)
- [ ] Code released on GitHub with DOI
- [ ] Docker container published
- [ ] Submitted to Genome Biology / Genome Research / Bioinformatics
- [ ] Presented at conference (ISMB / RECOMB)

---

# PART 7: RECOMMENDED JOURNAL TARGETS

Based on novelty and similar APA tool publications:

## Tier 1 (High Impact)
1. **Genome Biology** (IF: 12.3)
   - Published: scDaPars (2021), Sierra (2020)
   - Fit: Methodological novelty + database validation + biological insights
   - Requires: Strong validation, multiple datasets

2. **Genome Research** (IF: 7.0)
   - Published: scMAPA, SCAPTURE
   - Fit: Strong for computational methods with validation
   - Requires: Comprehensive benchmarking

## Tier 2 (Method-Focused)
3. **Bioinformatics** (IF: 5.8)
   - Published: scDAPA, SCAPE
   - Fit: Pure methodology with database validation
   - Requires: Clean technical implementation

4. **Briefings in Bioinformatics** (IF: 9.5)
   - Published: spvAPA (2024), scAPAtrap
   - Fit: Review-style with new method + validation
   - Requires: Comprehensive comparison with existing tools

## Tier 3 (Fast Publication)
5. **NAR Genomics and Bioinformatics** (IF: 4.0)
   - Open access, faster review
   - Good fit for computational tools
   - Requires: Solid validation, clear novelty

---

# PART 8: RISK MITIGATION

## Risk 1: Low Overlap with Reference Databases
**Likelihood**: Medium | **Impact**: High
**Mitigation**:
- Implement internal priming filter first (increases precision)
- Use multiple validation windows (10, 25, 50bp)
- Characterize novel peaks (motif, conservation, expression)
- Show that novel peaks have biological relevance

## Risk 2: Competing Methods Published During Development
**Likelihood**: Medium | **Impact**: High
**Mitigation**:
- Monitor bioRxiv for new APA tools
- Emphasize unique peak-based clustering + database validation
- Move quickly through Phase 1 (critical features)
- Publish preprint early

## Risk 3: Biological Validation Fails
**Likelihood**: Low | **Impact**: High
**Mitigation**:
- Use multiple datasets
- Start validation early
- Consult with biologists on interpretation
- Have backup datasets ready

## Risk 4: Benchmarking Shows Poor Performance
**Likelihood**: Low | **Impact**: Critical
**Mitigation**:
- Implement statistical rigor (Phase 1)
- Test on pilot data early
- Iterate on algorithm if needed
- Focus on novelty if performance is similar

## Risk 5: Reviewers Question Novelty
**Likelihood**: Medium | **Impact**: High
**Mitigation**:
- Clear novelty statement: "First tool validated against PolyASite/PolyA_DB for peak-based clustering"
- Comprehensive literature comparison (done)
- Cite 2024 benchmarking papers showing no tool does this
- Database validation adds credibility

---

# APPENDIX: Quick Reference

## Tools to Compare Against (Priority Order)
1. **scLAPA** - Closest competitor (multimodal)
2. **Sierra** - Good peak detection baseline
3. **scAPAtrap** - High precision baseline with PolyASite validation
4. (Optional) scDaPars, SCAPE, scMAPA

## Key Databases for Validation (Priority Order)
1. **PolyASite 2.0/3.0** - Most comprehensive, scRNA-seq included
2. **PolyA_DB v4** - Recently updated, isoform annotation
3. **APAeval** - Benchmark platform
4. **GENCODE** - Annotation validation

## Key Datasets (Priority Order)
1. Mouse spermatogenesis (GSE104556) - Most comparable, published peaks available
2. Human hESC (GSE75748) - Dynamic biology
3. PBMC 10x - Well-characterized
4. Mouse brain - Clear cell types

## Critical Features (Must-Have)
1. FDR correction ✓
2. Statistical p-values ✓
3. Internal priming filter ✓
4. **Database validation (NEW)** ✓

## High-Value Features (Should-Have)
1. PDUI calculation
2. Read1 utilization
3. K-NN imputation
4. **APAeval submission** ✓

## Nice-to-Have Features
1. Deep learning validation
2. Visualization suite
3. Spatial analysis mode

---

**Sources**:
- [PolyASite 2.0](https://academic.oup.com/nar/article/48/D1/D174/5588346)
- [PolyA_DB v4](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkaf1212/8356008)
- [APAeval Benchmark](https://link.springer.com/article/10.1186/s13059-021-02502-z)
- [PolyASite Database](https://polyasite.unibas.ch/)
