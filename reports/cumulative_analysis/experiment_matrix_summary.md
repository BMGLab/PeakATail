### A. Peak-calling grid (5 runs, full 17-dataset cohort input)

| run | strategy | atlas mode | ip mode | max PAS/peak | λ window | λ method | fold-change | PAS gap | min prominence |
|---|---|---|---|---|---|---|---|---|---|
| `lg_ip_filter` | lambda_gradient | annotate | filter | 5 | 5000 | median | 2.0 | 100 | 5.0 |
| `lg_ip_off` | lambda_gradient | annotate | off | 5 | 5000 | median | 2.0 | 100 | 5.0 |
| `lp_annotate` | lambda_poisson | annotate | annotate | 5 | 5000 | median | 2.0 | 100 | 5.0 |
| `si_annotate` | sierra_iterative | annotate | annotate | 5 | 5000 | median | 2.0 | 100 | 5.0 |
| `lg_annotate` | lambda_gradient | annotate | annotate | 5 | 5000 | median | 2.0 | 100 | 5.0 |

### B. Trim / annotation-window branches

| branch | max_gene_distance | utr_multiplier | include_extended |
|---|---|---|---|
| `A2_trim_d10000_ext` | 10000 | 2.0 | yes |
| `A2_trim_d1000_ext` | 1000 | 2.0 | yes |
| `A2_trim_d2000_ext` | 2000 | 2.0 | yes |
| `A2_trim_d3000_ext` | 3000 | 2.0 | yes |
| `A2_trim_d5000_ext` | 5000 | 2.0 | yes |
| `A2_trim_default` | 5000 | 2.0 | no |
| `A2_trim_mult1.5` | 5000 | 1.5 | no |
| `A2_trim_mult3.0` | 5000 | 3.0 | no |

### C. Clustering branches

| branch | method | resolution | n_neighbors |
|---|---|---|---|
| `A3_libsize` | leiden_libsize | 1.0 | 30 |
| `A3_nn15` | leiden_tfidf | 1.0 | 15 |
| `A3_nn50` | leiden_tfidf | 1.0 | 50 |
| `A3_res0.5` | leiden_tfidf | 0.5 | 30 |
| `A3_res2.0` | leiden_tfidf | 2.0 | 30 |

### D. Cell-level filters (identical across all runs unless noted)

| filter | value |
|---|---|
| `min_read` | 1500 |
| `min_cells` | 3 |
| `min_genes` | 50 |
| `min_pas_per_cell` | 50 |
