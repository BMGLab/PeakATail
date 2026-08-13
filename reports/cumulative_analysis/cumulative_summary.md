# PeakATail Cumulative Sweep Analysis

Sweep root: `/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed`

## 1. Peak-calling strategy comparison
| run | strategy | n_PAS | F1@50bp | F1@100bp | F1@500bp | F1@1000bp | atlas match |
|---|---|---|---|---|---|---|---|
| si_annotate | sierra_iterative | 43,035 | 0.086 | 0.095 | 0.146 | 0.196 | 77.9% |
| lg_annotate | lambda_gradient | 22,633 | 0.058 | 0.063 | 0.093 | 0.125 | 82.0% |
| lg_ip_off | lambda_gradient | 22,633 | 0.058 | 0.063 | 0.093 | 0.125 | 82.0% |
| lg_ip_filter | lambda_gradient | 20,609 | 0.053 | 0.058 | 0.086 | 0.116 | 82.0% |
| lp_annotate | lambda_poisson | 21,474 | 0.034 | 0.039 | 0.072 | 0.106 | 84.3% |

## 2. Internal-priming / fasta-filter axis
| ip_mode | n_PAS | n_PAS dropped vs off | % flagged internal-priming |
|---|---|---|---|
| annotate | 22,633 | 0 | 0.0% |
| filter | 20,609 | 2,024 | 8.94% |
| off | 22,633 | 0 | 0.0% |

## 3. Trim sensitivity
| branch | max_gene_distance | utr_multiplier | include_extended | n_PAS | n_cells_total | mean n_clusters |
|---|---|---|---|---|---|---|
| A2_trim_d1000_ext | 1000 | 2.0 | True | 17,758 | 55,647 | 14.88 |
| A2_trim_d2000_ext | 2000 | 2.0 | True | 17,829 | 55,662 | 14.82 |
| A2_trim_d3000_ext | 3000 | 2.0 | True | 17,915 | 55,676 | 14.82 |
| A2_trim_mult1.5 | 5000 | 1.5 | False | 16,454 | 55,411 | 14.94 |
| A2_trim_d5000_ext | 5000 | 2.0 | True | 18,097 | 55,697 | 15.06 |
| A2_trim_default | 5000 | 2.0 | False | 16,500 | 55,422 | 14.82 |
| A2_trim_mult3.0 | 5000 | 3.0 | False | 16,575 | 55,440 | 14.82 |
| A2_trim_d10000_ext | 10000 | 2.0 | True | 18,512 | 55,734 | 14.88 |

## 4. Clustering comparison
| branch | method | resolution | n_neighbors | mean n_clusters |
|---|---|---|---|---|
| A2_trim_default | leiden_tfidf | 1.0 | 30 | 14.82 |
| A3_res0.5 | leiden_tfidf | 0.5 | 30 | 10.53 |
| A3_res2.0 | leiden_tfidf | 2.0 | 30 | 22.41 |
| A3_nn15 | leiden_tfidf | 1.0 | 15 | 16.47 |
| A3_nn50 | leiden_tfidf | 1.0 | 50 | 13.82 |
| A3_libsize | leiden_libsize | 1.0 | 30 | 16.82 |

Cohort GEX concordance (n=17 GSMs):
- `ARI_gexleiden_vs_pas`: mean=0.4413, median=0.4564, std=0.145
- `AMI_gexleiden_vs_pas`: mean=0.638, median=0.6565, std=0.1118
- `ARI_celltype_vs_pas`: mean=0.4631, median=0.4627, std=0.158
- `AMI_celltype_vs_pas`: mean=0.6259, median=0.6616, std=0.1275

## 5. Switch cumulative (biological headline)
- **24 cell types** analyzed: **20 shortening** (3'UTR PDUI decreasing with stage), **4 lengthening**, 0 flat/other (mean slope -0.002738, median -0.002133)

### Recurrent switch genes (top 20)
_Genes with `|spearman| >= 0.8` in >= 3 cell types._
| gene_id | n_celltypes_trending | dominant_direction | mean_slope | mean_abs_spearman |
|---|---|---|---|---|
| ENSG00000092820 | 19 | decreasing | -0.0116 | 0.953 |
| ENSG00000124151 | 19 | decreasing | -0.0052 | 0.953 |
| ENSG00000143702 | 19 | decreasing | -0.0181 | 0.949 |
| ENSG00000152484 | 19 | decreasing | -0.0073 | 0.939 |
| ENSG00000136448 | 18 | decreasing | -0.0426 | 0.961 |
| ENSG00000068305 | 18 | decreasing | -0.0124 | 0.956 |
| ENSG00000011566 | 18 | decreasing | -0.0413 | 0.953 |
| ENSG00000102218 | 18 | decreasing | -0.0207 | 0.953 |
| ENSG00000204899 | 18 | decreasing | -0.0079 | 0.947 |
| ENSG00000150977 | 18 | decreasing | -0.0210 | 0.946 |
| ENSG00000123130 | 18 | decreasing | -0.0218 | 0.943 |
| ENSG00000149177 | 18 | decreasing | 0.0304 | 0.941 |
| ENSG00000153989 | 17 | decreasing | -0.0145 | 0.950 |
| ENSG00000172164 | 17 | decreasing | -0.0133 | 0.950 |
| ENSG00000124209 | 17 | decreasing | -0.0124 | 0.947 |
| ENSG00000136205 | 17 | decreasing | -0.0408 | 0.939 |
| ENSG00000083457 | 17 | decreasing | -0.0208 | 0.939 |
| ENSG00000105287 | 17 | decreasing | -0.0093 | 0.939 |
| ENSG00000164933 | 17 | decreasing | -0.0336 | 0.938 |
| ENSG00000206199 | 17 | decreasing | -0.0229 | 0.938 |

### Fisher significant-PAS hit counts by stage contrast
| contrast | n_celltypes | total_tests | total_sig (FDR<0.05) | mean frac sig |
|---|---|---|---|---|
| Normal_vs_StageI | 15 | 253,850 | 151,917 | 0.5834 |
| Met_vs_StageI | 18 | 247,384 | 109,164 | 0.4438 |
| Met_vs_Normal | 16 | 256,265 | 99,659 | 0.3888 |
| IVprimary_vs_StageI | 15 | 175,708 | 68,386 | 0.3836 |
| IVprimary_vs_Normal | 13 | 196,029 | 68,219 | 0.3570 |
| IVprimary_vs_Met | 13 | 104,092 | 43,788 | 0.4029 |
