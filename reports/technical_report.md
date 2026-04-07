# PeakATail Technical Report
## Methods, Strategies, and Preliminary Results

**March 2026**

> **Note**: All results use test data (SRR8325947, ~13.6M reads, human GRCh38). Validation on deeper sequencing data is planned. The aim is to present methods and approaches. Clustering and differential APA sections are proof-of-concept under active development.

---

## 1. Pipeline Architecture

![Pipeline Architecture](figures/fig_pipeline_flow.png)

The pipeline runs in 9 stages. GTF processing runs in a parallel thread during peak calling (cached — only parses once per GTF file). The streaming bisect architecture processes reads one-by-one without loading the full BAM into memory.

### Data dimensions at each stage (test dataset)

| Stage | Input | Output |
|-------|-------|--------|
| BAM input | 13.6M reads | 11.7M uniquely mapped |
| Peak calling (lambda_gradient) | 11.7M reads | 17,108 raw peaks (neg+pos) |
| CB filtering (min_read=1500) | 69,879 unique barcodes | 906 cells passing filter |
| GTF processing | 60,676 gene features | 19,555 genes with 3'UTR lengths |
| Gene annotation | 17,108 peaks | 11,771 annotated peaks (TIER 1+2+3) |
| Matrix | 17,108 x 906 sparse | 648,674 non-zero entries (99.6% sparse) |
| Clustering (Leiden res=1.0) | 906 cells x 16,448 PAS | 12 clusters |
| Differential APA | 2,223 multi-PAS genes | 1,165 significant PAS (q<0.05) |

---

## 2. Peak Calling Strategies

Four strategies implemented as a pluggable registry. Each receives a `Peak` object with coverage profile (`peak_list`) and per-position cell barcode tracking (`cb_positions`), returns a list of PAS positions with boundaries.

### 2.1 Original (Baseline)

Fixed height threshold with `pasfind()` heuristic:

1. Streaming: reads processed one-by-one, `bisect.insort()` maintains sorted coverage array
2. When `height >= 5` (fixed), signal starts
3. `pasfind()` finds first position where coverage exceeds 5% of maximum
4. Returns single PAS per peak region

**Limitation**: No statistical test, no multi-PAS, threshold not data-adaptive.

### 2.2 Lambda-Poisson

MACS2-inspired with local background estimation:

1. **Background lambda** from floor positions (coverage < 10% of peak max):
   - `lambda = median({h_i : h_i < 0.1 * max(h)})`

2. **Poisson p-value** tests if summit is significantly above background:
   - `p = P(X >= height | Poisson(lambda)) = 1 - CDF(height-1, lambda)`

3. **Boundary walk** from summit until coverage drops to lambda level
4. Returns single PAS with proper (start, end) boundaries

### 2.3 Lambda-Gradient (Hybrid — Novel Approach)

Two-phase: Poisson validates the region, gradient finds multiple PAS.

**Phase 1 — Region significance gate:**
- Compute local lambda from floor positions
- Poisson test on regional maximum. If `p >= 0.05`, reject region

**Phase 2 — Multi-PAS localization:**

1. **Smooth** coverage with moving average (window=50bp):
   - `smooth[j] = (1/w) * sum(h[j-w/2 : j+w/2])`

2. **Gradient** (first derivative): `gradient = diff(smooth)`

3. **Zero-crossings** where gradient goes positive to negative = local maxima = PAS candidates

4. **Per-PAS Poisson test**: Each candidate individually tested against lambda

5. **Prominence filter**:
   - `prominence[i] = h[i] - max(left_valley, right_valley)`
   - Reject if prominence < min_prominence

6. **Confidence score**:
   - `confidence = prominence * (h / lambda) + |second_derivative|`

7. **Partitioning**: Peak region divided by midpoints between summits. Each PAS gets own cell barcode counts from `cb_positions` — no reads lost.

### 2.4 Sierra-Iterative

Based on Sierra (Patrick et al., Genome Biology 2020):

1. Find maximum coverage position -> PAS_1
2. Zero out coverage within +/-300bp
3. Find new maximum -> PAS_2
4. Repeat until remaining max < min_height or max_pas reached
5. Partition by midpoints (same as lambda-gradient)

### Strategy visual comparison (Figure 1)

![Strategy explained](figures/fig_strategy_explained.png)

*Figure 1: Visual comparison of all 4 strategies on a simulated bimodal peak. (A) Original finds only 1 PAS via 5% threshold. (B) Lambda-Poisson adds statistical significance and boundary walk but still finds 1 PAS. (C) Lambda-Gradient (our hybrid) finds BOTH PAS via gradient zero-crossings — the smoothed coverage derivative changes sign at each local maximum. (D) Sierra-Iterative finds both via sequential subtraction.*

### Mathematical methods (Figure 2)

![Math methods](figures/fig_math_methods.png)

*Figure 2: (A) Lambda estimation: floor positions (green, < 10% of max) give the background level; peak positions (red) are excluded. Lambda = median of floor. (B) Poisson test: observed height 18 vs Poisson(lambda=3) gives p = 3.5e-09 — highly significant. (C) Gradient zero-crossing: where the first derivative crosses zero from positive to negative = local maximum = PAS candidate.*

### Peak counts per strategy (Figure 3)

![Peak counts](figures/fig_peak_counts.png)

*Figure 3: Total peaks found per strategy. Lambda-gradient and lambda-poisson find ~17K (statistically filtered). Sierra ~39K (aggressive multi-PAS). Original ~33K (no filter).*

---

## 3. Background Lambda Estimation

### Why Poisson for peak calling?

At peak calling, reads from all cells are pooled. Pooling reduces overdispersion, making Poisson (variance = mean) reasonable. Lambda is estimated from inter-peak floor:

`lambda = median({h_i : h_i < 0.1 * max_height})`

Robust to peaks because peaks are sparse outliers — median is insensitive when peaks < 50% of window.

### Dynamic streaming threshold

A `collections.deque` tracks recent read positions in sliding window (default 10kb):

```
local_lambda = len(deque) / window_size
threshold = max(floor_threshold, int(local_lambda * fold_change))
```

Replaces fixed threshold=5 with position-aware threshold adapting to local coverage.

### Why NOT Poisson for differential APA?

Per-cell scRNA-seq counts are overdispersed (variance >> mean). Fisher's exact test on pseudo-bulk is anti-conservative. Negative Binomial is the correct model (planned).

---

## 4. Gene Annotation

### Why `bedtools closest` instead of `intersect`

`intersect` only finds OVERLAPPING genes. A PAS 100bp downstream of a 3'UTR boundary would be silently dropped. `closest` finds the NEAREST gene for every PAS and reports distance.

### Adaptive UTR threshold

Fixed distance thresholds are wrong for half the genes. Using per-gene 3'UTR lengths from GTF:

| Tier | Condition | Count |
|------|-----------|-------|
| TIER_1 | distance <= UTR_length | 8,342 peaks (70.9%) |
| TIER_2 | distance <= UTR_length x 2.0 | 1,814 peaks (15.4%) |
| TIER_3 | distance <= 5000bp max | 1,615 peaks (13.7%) |
| INTERGENIC | distance > 5000bp | Removed |

### Annotation filter impact (Figure 4)

![Annotation impact](figures/fig_annotation_impact.png)

*Figure 4: (A) Precision @50bp before/after annotation filter. Annotation improves precision +8-13 points across all strategies. (B) ~22% of peaks removed as intergenic noise.*

---

## 5. Validation Against PAS Databases

Validated against two independent references:
- **PolyASite 2.0**: 569,005 PAS (3' end sequencing, Ensembl GRCh38)
- **PolyA_DB v4**: 303,312 PAS (3' end + long-read sequencing)

Primary metric: **precision** ("of PAS we found, what fraction are near a known PAS?"). Tested at 50, 100, 200, 500, 1000, 5000bp cutoffs.

### Dual database validation (Figure 5)

![Dual DB validation](figures/fig_dual_db_validation.png)

*Figure 5: Precision curves for all strategies against both databases. Lambda-gradient (orange) highest precision. Results concordant across databases.*

### Precision summary (Figure 6)

![Dual DB table](figures/fig_dual_db_table.png)

*Figure 6: Full precision table. Lambda-gradient: 40.6% @50bp (PolyASite), 33.0% @50bp (PolyA_DB). Full pipeline with matrix filtering reaches 50.7% @50bp and 94.4% @5kbp.*

### Precision by distance (Figure 7)

![Precision v3](figures/fig_precision_v3.png)

*Figure 7: All strategies with target lines. Lambda-gradient separates from baseline at tight cutoffs, indicating better PAS positioning.*

### Distance distribution (Figure 8)

![Distance hist](figures/fig_distance_hist.png)

*Figure 8: Left (0-200bp linear): spike at 0-10bp = precise PAS detection. Right (0-5000bp log): exponential decay confirms real signal near genes.*

### Strategy overlap (Figure 9)

![Venn](figures/fig_venn_overlap.png)

*Figure 9: Lambda-gradient shares 14,925 PAS with sierra but only 231 with original. Statistical strategies find fundamentally different peaks.*

### Parameter sensitivity (Figure 10)

![Metrics table](figures/fig_metrics_table.png)

*Figure 16: Lambda-gradient precision varies 31.5%-32.9% across all tested parameter combinations. The algorithm is robust and not parameter-sensitive.*

---

## 6. Clustering (Under Development)

### Method: TF-IDF + LSI + Leiden

Peak data is sparse and binary-ish (like scATAC-seq), so we use the scATAC-seq approach instead of standard gene expression normalization.

![TF-IDF + LSI Pipeline](figures/fig_tfidf_pipeline.png)

**TF-IDF** (Signac Method 1):
- TF = count / total_counts_per_cell
- IDF = n_cells / cells_with_peak
- Result = log1p(TF * IDF * 10000)

**LSI**: TruncatedSVD on sparse TF-IDF matrix. First component removed (depth-correlated).

**Cosine distance** for k-NN: appropriate for sparse high-dimensional data.

**Leiden**: Replaces deprecated Louvain. Guarantees connected communities.

### Resolution sweep (Figure 11)

![ARI sweep](figures/fig_ari_ami_sweep.png)

*Figure 15: ARI/AMI vs gene expression clusters across resolutions. TF-IDF (orange) beats library-size (blue) at every resolution. Best: res=1.0, 12 clusters, ARI=0.869.*

### Clustering comparison (Figure 12)

![Main comparison](figures/fig_main_comparison.png)

*Figure 16: (A) Gene expression UMAP. (B) Peak-based UMAP. (C) Confusion matrix (ARI=0.869). (D) GEX cluster composition. (E) Resolution sweep. (F) GEX clusters split by APA.*

### Sankey diagram (Figure 13)

![Sankey](figures/fig_sankey.png)

*Figure 15: Cell flow from gene expression (left) to peak-based clusters (right). 687 common cells, mostly 1:1 mapping.*

### Cluster sizes (Figure 14)

![Cluster sizes](figures/fig_cluster_sizes.png)

*Figure 16: Cluster size distribution for 5 strategies. Power-law shape consistent with real biology.*

### APA subclusters (Figure 15)

![Subclusters](figures/fig_subcluster.png)

*Figure 15: GEX clusters split by APA. GEX_6 splits into 3 APA subclusters — cells with identical gene expression but different polyadenylation.*

---

## 7. Differential APA — Switch Test (Under Development)

### Method

For each gene with >=2 PAS, Fisher's exact test on 2x2 table per PAS. BH FDR correction across all tests.

**PDUI** (Proximal-Distal Usage Index):
- + strand: proximal = min coord, distal = max coord
- - strand: proximal = max coord, distal = min coord
- `PDUI = distal / (proximal + distal)`
- PDUI=0: 3'UTR shortening, PDUI=1: lengthening

### Results

| Pair | Tests | Significant (q<0.05) | Strong (q<0.05, |dP|>0.2) |
|------|-------|---------------------|--------------------------|
| 0 vs 1 | 5,603 | 207 | 195 |
| 0 vs 2 | 6,309 | 49 | 44 |
| 1 vs 2 | 5,933 | 230 | 191 |
| 0 vs 3 | 6,548 | 107 | 84 |
| 1 vs 3 | 6,294 | 572 | 481 |
| **Total** | **30,687** | **1,165** | **995** |

174 genes significant in multiple comparisons (robust). 65% of paired PAS have complementary delta-proportions (biological consistency).

### Differential APA figure (Figure 16)

![Differential APA](figures/fig_differential_apa.png)

*Figure 16: (A) PDUI heatmap — 2,223 genes x 12 clusters. Red=distal, blue=proximal. (B) Delta-PDUI distribution per pair. (C) Volcano: 1,165 significant PAS.*

**Caveat**: P-values inflated ~2x (pseudo-bulk Fisher). NB regression planned.

---

## 8. Summary

### What works

| Component | Key metric |
|-----------|------------|
| Peak calling (lambda_gradient) | 40.6% precision @50bp |
| Gene annotation (adaptive UTR) | 99% assignment, tiered confidence |
| Dual-database validation | Concordant on PolyASite + PolyA_DB |
| TF-IDF + LSI clustering | ARI=0.869 vs gene expression |
| Fisher + FDR + PDUI | 1,165 significant PAS, 2,223 PDUI genes |

### What needs improvement

| Issue | Plan |
|-------|------|
| No internal priming filter | Needs genome FASTA |
| Fisher p-values inflated ~2x | Implement NB regression |
| Single dataset | Test on >=3 datasets |
| No experimental validation | qPCR/RNA-FISH on top hits |

---

## Appendix: Figure Index

| # | Description | File |
|---|-------------|------|
| 1 | Peak counts per strategy | `fig_peak_counts.png` |
| 16 | Annotation filter impact | `fig_annotation_impact.png` |
| 15 | Dual database validation | `fig_dual_db_validation.png` |
| 16 | Precision metrics table | `fig_dual_db_table.png` |
| 15 | Precision by distance | `fig_precision_v3.png` |
| 16 | Distance histogram | `fig_distance_hist.png` |
| 15 | Strategy overlap Venn | `fig_venn_overlap.png` |
| 16 | Parameter sensitivity | `fig_metrics_table.png` |
| 15 | ARI/AMI resolution sweep | `fig_ari_ami_sweep.png` |
| 16 | 6-panel clustering comparison | `fig_main_comparison.png` |
| 15 | Sankey flow diagram | `fig_sankey.png` |
| 16 | Cluster sizes | `fig_cluster_sizes.png` |
| 15 | APA subclusters | `fig_subcluster.png` |
| 16 | Differential APA (PDUI + volcano) | `fig_differential_apa.png` |
