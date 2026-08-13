# Figure manifest — corrected sweep (`RERUN_2026-08_fixed`)

Every figure below is a **standalone file** in `reports/figures/corrected/`, as
both `.png` (200 dpi) and `.svg` (vector). None exists only inside the `.tex`.
All are byte-reproducible: two runs of the generator produce identical files
(`svg.hashsalt` pinned, SVG `Date` metadata suppressed).

**Generators** (both read only on-disk tables; nothing is hardcoded, and a figure
whose input is missing is skipped rather than drawn from stale values):

| generator | figures |
|---|---|
| `reports/generate_corrected_figures.py` | 01–14 |
| `reports/generate_clustering_figures.py` | 15–19 |

**Input tables** live in `reports/cumulative_analysis/` (`tables/` from
`ema.benchmark.sweep_analysis`, `extra/` from `reports/_server_scripts/s1b,s2,s3,s4,s5,s6`).
All server scripts read the real sweep at
`/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed`. No
fixture, mock or synthetic data is used in any figure.

Axis column maps to `reports/REPORT_ANALYSIS_PLAN.md` §1–5.

---

## Plan §1 — Peak-strategy axis

| # | File | Caption | Claim it answers |
|---|---|---|---|
| **01** | `01_strategy_benchmark.png/.svg` | Peak count per strategy (a), F1 across cutoffs (b), F1@1kb against peak count (c). | **The F1 ranking is a peak-count ranking.** F1 and n_PAS are rank-identical across all five arms — so the benchmark ranks yield, not accuracy. Retires the "sierra > lambda_gradient" claim. |
| **02** | `02_atlas_saturation.png/.svg` | Precision pinned >99.6% at every cutoff (a); recall carries all the variation (b). | **Precision cannot discriminate.** 0.9964–0.9994 for *every* arm against an 18,432,135-entry reference. Complements plan §1.2 — my panel shows the saturation, the analyst's null control quantifies the excess-over-chance (+0.245, not +0.997). **Use both together.** |
| **03** | `03_count_controlled.png/.svg` | Rarefaction to n=20,605, 200 reps (a); distance-to-nearest-atlas distribution (b); pairwise containment heatmap (c). | **Count-controlled, the arms are indistinguishable and nested.** Spread 0.15 pp, median distance 0 bp; lambda_gradient is 100% contained in sierra, only 53% the reverse. Directly serves plan §1.3 (PAS-set agreement). |
| **04** | `04_ip_and_trim.png/.svg` | Internal-priming axis (a); trim window vs n_PAS (b); trim vs cluster count (c). | ip-filter removes **8.94%** of peaks (22,633→20,609); `annotate` ≡ `off` byte-identically. Trim 1 kb→10 kb adds 754 PAS; cluster count moves 0.24 across the whole grid. Serves §1 (ip axis) and §2 (trim axis). |

## Plan §2 — Trim axis

Covered by **04b–c** above. The trim axis produces no separate figure because its
effect is null: see the 0.24-cluster range and the analyst's ΔARI = 0.007 for
`max_gene_distance`. That agreement across two independent analyses is the result.

## Plan §3 — Clustering axis

| # | File | Caption | Claim it answers |
|---|---|---|---|
| **05** | `05_clustering.png/.svg` | Resolution sweep (a); n_neighbors (b); PAS-vs-GEX agreement boxplot over 17 datasets (c). | **ARI 0.463 / AMI 0.626**, superseding the withdrawn single-run ARI = 0.869. AMI ≫ ARI ⇒ concordant but more finely split. |
| **15** | `15_umap_cohort_stages.png/.svg` | Cohort UMAP, one dataset per stage; top row by PAS-Leiden cluster, bottom by GEX cell type; same cells, same embedding. | **PAS space is structured, not noise** — discrete islands with coherent cell-type colouring. The qualitative counterpart of ARI 0.463. |
| **16** | `16_umap_clustering_grid.png/.svg` | The six clustering variants on one shared dataset (GSM3516666-Normal), by cluster and by cell type. | **The embedding is fixed; only the cut changes.** Serves plan F3.3 (tfidf vs libsize UMAP grid). |
| **17** | `17_sankey_clustering_grid.png/.svg` | PAS cluster → GEX **cell type** flow, per clustering variant, on the shared dataset. | **The median PAS cluster sends 79% of its labelled cells to one cell type**; mapping is many-to-one. Purity 0.75/0.81/0.84 at res 0.5/1.0/2.0 — *not* evidence res 2.0 is better (purity rises mechanically as clusters shrink). Serves plan F3.3 (Sankey, tfidf vs libsize). |
| **19** | `19_cluster_sizes.png/.svg` | Cluster count per dataset across **all 19 experiments** (a); size distribution across clustering variants (b); fragmentation vs cluster count (c). | Median **15 clusters/dataset (range 7–29)**; smooth size distributions, no dominant-cluster-plus-fragments pattern. Feeds the analyst's fragmentation/largest-cluster-fraction metric. |
| **18** | `18_sankey_peak_strategies.png/.svg` | Same cluster → cell type flow, per peak-calling strategy, on the shared grid dataset (GSM3516675-Normal). | **Cluster → cell-type structure survives the strategy swap** — the visual form of the analyst's "peak strategy is inert, ΔARI 0.011". Serves §1.4 (downstream invariance). |

> **Caveat carried in every panel of 15–19:** cell-type calls exist for **53.0% of
> the 55,422 clustered cells** (27–87% per dataset). The unlabelled remainder is
> drawn in grey, never dropped — discarding it would make the PAS→cell-type
> mapping look far cleaner than it is.
>
> **Each Sankey panel is a single named dataset, deliberately.** `leiden` ids are
> assigned per dataset, so cluster "0" in one GSM is unrelated to cluster "0" in
> another; pooling merges unrelated clusters into one node. An earlier draft
> pooled all 17 and was wrong.

## Plan §4 — Filter-effect scenarios

**No figures yet — Phase 3 data lands tonight.** Reserved slots, to be produced by
`generate_clustering_figures.py` + `generate_corrected_figures.py` once all four
scenarios (baseline / atlas_filter / annot_filter_3utr / ip_filter) exist on the
6-GSM subset:

| # (reserved) | Intended content |
|---|---|
| **20** | PAS retained + clustering ARI/AMI vs baseline, per scenario |
| **21** | Differential hit counts and hit-set overlap vs baseline, per scenario |
| **22** | Shortening/lengthening split + slope distribution, per scenario |
| **23** | UMAP + Sankey per filter scenario (same form as 16/17) |

Report §14 already scaffolds this with **three predictions recorded in advance**,
so the result reads as a test rather than a description.

## Plan §5 — Switch biology (cell type × stage)

| # | File | Caption | Claim it answers |
|---|---|---|---|
| **06** | `06_depth_confound.png/.svg` | Reported PDUI vs zero-coverage fraction (a); stage trajectory unconditional vs coverage-conditioned (b); per-cell-type slopes under each conditioning (c). | **The "global 3′UTR shortening" is a detection-rate artifact.** Pearson **r = −0.997** (n=76, p=3.9e−83); conditioning on ≥1 read flattens the trend, ≥5 reads reverses it (4/18 down vs 15/18). **The single most important figure in the report.** |
| **07** | `07_stage_trajectory.png/.svg` | Mean PDUI by cell type × stage, as reported (a) and conditioned on ≥5 reads (b). | The left-to-right gradient in (a) is the detection-rate gradient; it is absent in (b). |
| **08** | `08_fisher_power.png/.svg` | Distribution of per-(cell type, contrast) significance rate (a); mean rate by stage contrast (b). | **43.9% of 1,233,328 Fisher tests "significant"** — pseudoreplication over reads. Rank, never threshold. |
| **12** | `12_diff_strategy.png/.svg` | Significance rate per strategy (a); hits vs tests with the all-significant diagonal (b); containment vs Jaccard (c). | **nb_multi calls 99.8% significant** (84–100% in every cell type) — a variability screen, not a differential test. The two strategies corroborate on only **24%**. Leads with containment because a 33× set-size gap pins Jaccard near zero for arithmetic reasons. |
| **13** | `13_length_strategy.png/.svg` | Zero-coverage share per metric (a); mean as-reported vs coverage-conditioned (b); stage trajectory per metric (c). | **No length metric rescues the stage signal.** classic/shannon are 96–97% zero-coverage rows; conditioning moves both (0.012→0.381, 0.570→0.688) and neither trends. **`proportion` is marked INVALID on the figure face** — 98.2% uniform `1/n_PAS` padding. |
| **09** | `09_recurrence.png/.svg` | Cell types per trending gene (a); 25 most recurrent genes (b); effect size, recurrent vs private (c). | **Recurrence tracks testability, not shared regulation** — private genes have *larger* mean \|slope\| (0.0137 vs 0.0104, MWU p=0.029). |
| **11** | `11_gene_sets.png/.svg` | Cancer-gene membership rates (a); crude vs detectability-adjusted odds ratios (b); Hallmark enrichment (c). | **No cancer-gene or pathway signal survives adjustment.** COSMIC CGC crude OR 1.99 (p=6.7e−09) → adjusted **1.22 (p=0.28)**; 0/50 Hallmark sets at FDR<0.05. |
| **10** | `10_lost_distal_elements.png/.svg` | Element density lost vs length-matched control (a); ARE density vs AT content (b); lost-segment length distribution (c). | **The miRNA/ARE de-repression mechanism is unsupported.** Seed density 0.510/kb vs 0.520/kb (p=0.058, wrong sign). **Median lost interval 10,779 bp, 80% >3 kb** — the PDUI layer pairs PAS too far apart to be a tandem 3′UTR pair. A *pipeline* finding. |
| **14** | `14_gene_walk.png/.svg` | Per-gene PAS usage by stage (EZR, DPYD, PHACTR1), read share within stage, pooled over cell types. | Replaces the two withdrawn `gene_walk_*` panels. DPYD's 148 PAS are capped to the top 12 with an explicit "other" bar carrying its remaining ~60% of reads. ⚠ **Provisional** — see status below. |

---

## Status and caveats for the assembler

1. **14 is FINAL.** The geneview pass completed; `gene_walk_tracks.tsv` covers all
   **3 genes × 24 cell types × 4 stages** (4,795 rows) and the figure is
   regenerated at full coverage and committed. EZR 17 PAS / 12,034 reads,
   DPYD 148 PAS / 93,988 reads, PHACTR1 22 PAS / 2,004 reads. Each panel shows the
   twelve best-supported sites in genomic order with the remainder pooled into a
   labelled "other" bar — which for DPYD carries ~60% of the gene's reads, so the
   bar is load-bearing rather than cosmetic.
2. **No figure is known to be wrong.** Three were wrong during development and
   are fixed: the Sankey pooled datasets with incomparable cluster ids; the UMAP
   buried cell types under the unlabelled grey; `13` presented `proportion` in a
   way that made the broken metric look like the best one.
3. **02 and 03 should be shown with the analyst's null control**, not instead of
   it. Mine show that precision saturates and that count-control collapses the
   ranking; the null control quantifies how much of the remaining signal is
   chance. Together they retire the strategy ranking; separately either can be
   read as a weaker claim.
4. **Superseded figures.** `fig_*.png` and `01–13_*.png` in `reports/figures/`
   are withdrawn for all data panels; only the schematics remain valid
   (`fig_pipeline_flow`, `fig_strategy_explained`, `fig_math_methods`,
   `fig_tfidf_pipeline`, `*_explainer`). Six data figures that appear only in the
   `.tex` carry a `[PRE-FIX RUN — illustrative only]` caption prefix.
5. **Rebuilding the PDF** is `bash reports/build_latex.sh` — it now uses a fresh
   build directory, exits non-zero on LaTeX errors, and reports page count and
   figures embedded. Last clean build: **52 pages, 18 corrected figures, 0 errors.**

---

# Cross-experiment figures (24–36) — added by the assembler

Generator: `reports/generate_cross_experiment_figures.py`, reading the TSVs in
`reports/cumulative_analysis/cross_experiment/` produced by
`ema.benchmark.cross_experiment.run_all()` against the real
`RERUN_2026-08_fixed` sweep. Same contract as 01–19: nothing hardcoded, a figure
whose input table is missing is skipped. Numbered from 24 so as not to collide
with 01–19 or the **reserved filter-effect slots 20–23**.

These answer the *cross-run* questions — whether a difference between two
branches is real — where 01–19 answer the *per-run* ones. Where both bear on a
claim, the report shows them together.

| # | File | Caption | Claim it answers |
|---|---|---|---|
| **24** | `24_atlas_null_control` | Observed precision vs a shifted (±5–50 kb, same chromosome) null and a uniform null, per arm per cutoff. | **~70% of the headline precision is free.** A displaced PAS scores 0.702 at ±50 bp. Quantifies what 02 shows qualitatively — **use both together**. |
| **25** | `25_atlas_excess_over_null` | Excess precision over the shifted null, per arm, across cutoffs. | Excess spread across arms is **0.028** — the benchmark has no power to rank strategies at any cutoff. |
| **26** | `26_pas_set_jaccard` | Pairwise PAS-set Jaccard (±50 bp). | **`lg` vs `lp` = 0.294** — two arms scoring identical precision share 29% of their sites. Complements 03c (containment) with the symmetric measure. |
| **27** | `27_yield_vs_celltype_ari` | PAS yield against ARI vs GEX cell type, one point per arm. | Doubling the catalogue does not improve cell-type recovery (0.4002–0.4116). The scatter form of "peak strategy is inert". |
| **28** | `28_trim_sensitivity` | Trim window vs cell-type ARI, cluster count, PAS yield; ARI panel scaled against the tfidf→libsize range. | **ΔARI 0.007 over a 10× window change.** Independent confirmation of 04b–c from an external criterion. Yield panel carries the reproducibility caveat. |
| **29** | `29_celltype_reassignment` | % of cells whose cell-type call changes vs the base cohort, per branch. | Translates ARI deltas into consequences: **libsize 30.8%**, res 2.0 19.3%, res 0.5 14.4%, trim branches 6–10%, nn50 4.0%. |
| **30** | `30_reannotate_reproducibility` | Table-figure: branches with identical declared trim params, their `pas_gene.tsv` md5, rows, genes, write time. | **A reproducibility defect.** Five branches declaring `{5000, 2.0, false}` produced three md5 groups that track *write time*; `annotatedpas.bed` identical throughout. The parameter-free spread exceeds the trim effect. **Trim-axis yields are withdrawn.** |
| **31** | `31_knob_ranking` | Every swept knob ranked by the full range (max−min, default included) it spans in ARI vs cell type. | **The headline positive result: tfidf→libsize = ΔARI 0.140**, three times any other knob. Uses full range, not endpoints, because resolution is non-monotone. |
| **32** | `32_resolution_ari_silhouette` | Resolution vs ARI (left axis) and silhouette (right axis), annotated with cluster counts. | **Silhouette is the wrong objective** — it rises monotonically as resolution falls while cell-type agreement peaks at the default 1.0. Same caution as Sankey purity in 17. |
| **33** | `33_pdui_celltype_stage` | Mean PDUI by cell type × stage; rows with <3 stages flagged. | 6 of 24 cell types have only two stages, where Spearman ρ is ±1 **by construction**. Those rows must not be cited. |
| **34** | `34_trend_depth_confound` | Slope on all gene–cell pairs vs slope on informative pairs only; stage-count histogram. | **9 of 24 cell types flip sign** under the detection control; "decreasing" falls 0.833→0.708. Per-stage mean PDUI vs depth: **median r = 0.9957**. Reaches 06's conclusion from a different statistic. |
| **35** | `35_fisher_pseudoreplication` | Significant PAS vs cells in the contrast, sized by median \|Δ proportion\| of the hits. | Significance tracks cell count (**ρ = 0.556, p = 1.3e−8**) while effect size runs the *other* way: 14,000 hits at \|Δ\|=0.19 vs 800 at 0.82. |
| **36** | `36_recurrent_switches` | Recurrence distribution across cell types; 20 most recurrent genes coloured by direction consistency. | **Recurrence cannot be used as a confidence filter.** Only 18.6% of 7,941 genes are private; the top 15 are significant in 24/24 with *inconsistent* direction. Agrees with 09 from the opposite side. |

## Assembler's notes

1. **Figures 20–23 remain reserved.** Do not renumber 24–36 to close the gap; the
   report's figure index and cross-references depend on these numbers.
2. **`_load()` uses `index_col=False`.** A stray trailing delimiter in a source
   TSV had promoted column 0 to the index and silently shifted every column left
   in figure 30. The fix is in the loader, so it protects every future table too.
3. **Corroboration is deliberate.** 02+24, 03+26, 04+28, 06+34, 08+35, 09+36 are
   pairs reached by independent analyses. The report cites both in each case;
   dropping either weakens the claim to a single-method result.
