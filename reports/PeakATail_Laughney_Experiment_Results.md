# PeakATail on the Laughney LUAD cohort — what each experiment asked, and what it answered

*Experiment-driven results for the corrected `RERUN_2026-08_fixed` sweep.
Every section states the question, shows the view that discriminates the arms,
and gives the verdict — including where the verdict is "no effect" or "the
metric cannot tell".*

Companion documents: `REPORT_ANALYSIS_PLAN.md` (why these views and not others),
`ema/benchmark/cross_experiment.py` (the analyses), `scripts/analyst_figures.py`
(the figures).

---

## Executive summary

Six questions were put to the data. Two returned a usable positive result, three
returned negatives that overturn claims currently in the report, and one is still
running.

| # | Question | Answer |
|---|---|---|
| 1 | Does the peak-calling strategy matter? | **Not for the biology.** `lg` and `lp` share only **29 %** of their PAS yet differ in cell-type recovery by **0.001 ARI** |
| 2 | Is "precision ≈ 1.0 vs PolyASite" evidence of accuracy? | **No.** A PAS displaced 5–50 kb at random scores **0.70**. Excess over null is **+0.245**, and it does not separate the strategies |
| 3 | Which pipeline parameter actually controls the result? | **The normalisation method.** `leiden_tfidf` beats `leiden_libsize` by **ΔARI = +0.140**; every other knob is ≤ 0.048, and `max_gene_distance` is **0.002** |
| 4 | Do 3′UTRs shorten as the tumour progresses? | **Not demonstrably.** Per-stage mean PDUI correlates with per-stage sequencing depth at **median r = 0.996**; 41 % of cell types flip sign when depth is controlled; **no cell type survives multiple-testing correction** |
| 5 | Are the differential-APA hits trustworthy? | Fisher is pseudoreplicated over reads — see §5 |
| 6 | What does each PAS filter cost downstream? | **Pending** — Phase 3 has not completed |

The two positive results (TF-IDF normalisation; PAS clustering recovering cell
identity at ARI ≈ 0.46 / AMI ≈ 0.63 against independent GEX labels) are worth
publishing. The three negatives are worth reporting honestly, because each is
currently stated the other way round.

---

## 1. Does the peak-calling strategy change the answer?

**Experiment.** `peakcall_grid.tsv`, 6-GSM subset: `lambda_gradient` (`lg`),
`lambda_poisson` (`lp`), `sierra_iterative` (`si`), plus the internal-priming
axis (`annotate` / `filter` / `off`). Everything downstream is identical.

**What discriminates.** PAS yield and width; pairwise PAS-set Jaccard; and
cell-type recovery. *(Figures F1.2, F1.3.)*

| arm | PAS | median width | genes | multi-PAS genes | ARI vs cell type |
|---|---|---|---|---|---|
| `lg_annotate` | 22,633 | 640 bp | 12,283 | 74.9 % | 0.4002 |
| `lg_ip_off` | 22,633 | 640 bp | 12,283 | 74.9 % | 0.4002 |
| `lg_ip_filter` | 20,609 | 636 bp | 12,144 | 73.7 % | 0.4052 |
| `lp_annotate` | 21,474 | **403 bp** | 12,279 | 75.1 % | 0.4005 |
| `si_annotate` | **43,035** | 522 bp | 13,291 | 79.4 % | 0.4116 |

Pairwise PAS-set agreement (Jaccard, ±50 bp):

| pair | Jaccard | interpretation |
|---|---|---|
| `lg` vs `lg_ip_off` | **1.000** | `--ip-filter-mode annotate` is a verified no-op on the PAS set |
| `lg` vs `lg_ip_filter` | 0.911 | strict subset; the filter drops 9 % of sites |
| `lg` vs `si` | 0.520 | `si` recovers **99.2 %** of `lg` and adds ~20 k more — near-superset |
| **`lg` vs `lp`** | **0.294** | each recovers only ~44 % of the other |

**Verdict.** `lambda_gradient` and `lambda_poisson` disagree about **more than
half** of the polyadenylation sites in the genome, and produce cell-type recovery
that differs by 0.0003 ARI. Strategy choice is a decision about the **PAS
catalogue**, not about the biology — and it should be presented that way instead
of as a ranked leaderboard. `sierra_iterative` is the right choice if catalogue
completeness is the goal (2× the sites, near-superset of `lg`, most multi-PAS
genes); `lambda_poisson` gives the tightest peaks (403 bp median vs 640 bp).

---

## 2. Is the PolyASite benchmark measuring accuracy?

**Experiment.** Every arm is scored against PolyASite v3.0 by
`bench_pas_vs_atlas.py`, reporting precision/recall at 50–5000 bp.

**The problem.** The reference as deployed contains **18,432,135 single-base
sites** — one every ~168 bp of genome. A ±50 bp window is therefore "atlas-positive"
for a large fraction of the genome *by construction*, and every arm scores
precision 0.996–0.999. A number that cannot go down cannot rank anything.

**What discriminates: a null control** (`atlas_null_control`). Each called PAS is
displaced by ±U(5 kb, 50 kb) on its own chromosome — staying in the same genomic
neighbourhood — and precision is recomputed. *(Figures F1.1, F1.1b.)*

`lg_annotate`, by cutoff:

| cutoff | observed | null (shifted) | null (uniform) | **excess over null** |
|---|---|---|---|---|
| 10 bp | 0.472 | 0.261 | 0.130 | **+0.211** |
| 25 bp | 0.822 | 0.514 | 0.271 | **+0.308** |
| **50 bp** | **0.947** | **0.702** | 0.409 | **+0.245** |
| 100 bp | 0.986 | 0.842 | 0.550 | +0.144 |
| 200 bp | 0.995 | 0.923 | 0.677 | +0.072 |

Excess-over-null at 50 bp across arms: `lg` +0.245, `lp` +0.243, `si` +0.221,
`lg_ip_filter` +0.249. The spread is **0.028** — smaller than the run-to-run
noise anyone would accept.

Median distance from a called PAS to the nearest atlas site is **11 bp**, versus
**24 bp** for the shifted null.

**Verdict.** There *is* real signal — PAS land twice as close to atlas sites as
displaced controls, and at ±10 bp the caller beats chance by 0.21 — but the
headline "precision ≈ 1.0" is roughly **70 % free**, and the benchmark has no
power to rank strategies. Two concrete changes follow:

1. Report precision at **±10–25 bp with the null**, never ±50 bp bare.
2. **Retire the atlas-based strategy ranking.** §1 shows arms sharing 29 % of
   PAS score identically; that is a property of the metric, not the callers.

---

## 3. Which parameter actually controls the result?

**Experiment.** 13 `ema reannotate` branches off the *same* cohort peak-calls
(`trim_cluster_grid.tsv`), sweeping `max_gene_distance`, `utr_multiplier`,
`include_extended`, `resolution`, `n_neighbors`, and clustering method —
one factor at a time, 17 datasets each.

**What discriminates.** ARI/AMI against **independent GEX cell-type labels** —
an external criterion, not an internal cluster-compactness score. *(Figure F3.1.)*

| knob | range swept | Δ ARI vs cell type | verdict |
|---|---|---|---|
| **clustering method** | `leiden_tfidf` → `leiden_libsize` | **0.140** | **dominant** |
| resolution | 0.5 → 2.0 | 0.048 | real, optimum at 1.0 |
| `include_extended` | false → true | 0.015 | minor |
| n_neighbors | 15 → 50 | 0.011 | minor |
| peak strategy | `lg`/`lp`/`si` | 0.011 | inert |
| `utr_multiplier` | 1.5 → 3.0 | 0.002 | inert |
| **`max_gene_distance`** | 1000 → 10000 (10×) | **0.002** | **inert** |

Full branch table (mean over 17 datasets):

| branch | clusters | ARI vs ref | ARI vs cell type | AMI vs cell type | silhouette |
|---|---|---|---|---|---|
| `BASE` / `A2_trim_default` | 14.8 | 1.000 | **0.4631** | **0.6259** | 0.122 |
| `A3_nn50` | 13.8 | 0.885 | 0.4568 | 0.6216 | 0.126 |
| `A2_trim_mult3.0` | 14.8 | 0.857 | 0.4533 | 0.6253 | 0.123 |
| `A2_trim_d1000_ext` | 14.9 | 0.804 | 0.4510 | 0.6266 | 0.124 |
| `A2_trim_d10000_ext` | 14.9 | 0.783 | 0.4487 | 0.6279 | 0.127 |
| `A3_nn15` | 16.5 | 0.822 | 0.4458 | 0.6224 | 0.116 |
| `A3_res0.5` | 10.5 | 0.729 | 0.4283 | 0.6159 | **0.134** |
| `A3_res2.0` | 22.4 | 0.642 | 0.4154 | 0.6078 | 0.101 |
| **`A3_libsize`** | 16.8 | **0.359** | **0.3234** | **0.4822** | n/a |

**Verdict — three results.**

**(a) TF-IDF/LSI normalisation is the pipeline's most consequential design
choice.** Replacing it with library-size normalisation costs **ΔARI = 0.140**
(0.463 → 0.323) and **ΔAMI = 0.144** against independent cell-type labels, and
produces the most fragmented partition (largest-cluster fraction 0.228). Per our
own literature survey, no published APA tool has benchmarked PAS-based clustering
against cell-type labels with ARI/AMI — this is the sweep's most publishable
number.

**(b) Resolution has a genuine optimum at the default 1.0**, and *silhouette is
the wrong objective*: silhouette rises monotonically as resolution falls (0.134 at
res 0.5 vs 0.101 at res 2.0) while cell-type agreement peaks in the middle. Show
both curves (F3.2) and use the external metric to choose.

**(c) `max_gene_distance` is inert over 1000–10000 bp.** A 10× change moves ARI by
0.002. All `_ext` branches sit at ARI-vs-reference 0.78–0.80 *regardless of
distance* — meaning `include_extended` is the real knob on that axis and
`max_gene_distance` is not. Stop sweeping it; document the default.

> ### ⚠ Reproducibility defect blocking the trim axis
>
> Five branches declaring **identical** trim parameters (`{5000, 2.0, false}`)
> produced **three different** `pas_gene.tsv` files:
>
> | branches | md5 (12) | assignments | genes | written |
> |---|---|---|---|---|
> | `B1_cohort_full`, `A3_res0.5`, `A3_res2.0`, `A3_nn15` | `312a7ae1392d` | 51,129 | 11,686 | 14:37 |
> | `A2_trim_default`, `A3_nn50` | `ae48369d2302` | 67,409 | 11,735 | 14:47 |
> | `A3_libsize` | `916bb04e4e3a` | 82,702 | 13,425 | 14:42 |
>
> `annotatedpas.bed` is identical (184,295 rows) everywhere and no rows are
> duplicated, so the PAS catalogue is fixed and **only the PAS→gene assignment
> diverges**. Clustering parameters cannot affect gene assignment, so this is not
> a parameter effect; the grouping tracks **write time**, the signature of a
> shared-state race between concurrently running `ema reannotate` branches.
>
> The parameter-free spread (51,129 ↔ 67,409) **exceeds** the spread the trim
> sweep is trying to measure (58,982 for `d5000_ext` vs 67,409 for `default`).
> **PAS-yield numbers from the trim axis are not reportable until it is re-run
> serially.** Per-cell clustering ARIs are far less exposed (the divergence is in
> feature assignment, not cell assignment) but should be re-confirmed. The
> `libsize` comparison in (a) is partly exposed — though `libsize` received *more*
> features than the branches it loses to, so the direction is unlikely to reverse.

---

## 4. Do 3′UTRs shorten as the tumour progresses?

**Experiment.** Per cell type, PDUI is computed per cell and averaged per stage
(Normal → StageI → IVprimary → Met); `ema switch trend` fits an ordered-stage
slope and Spearman ρ. 24 cell types.

This is the sweep's central biological claim, and it does not survive three
checks. *(Figures F5.1, F5.2.)*

### 4.1 Most trends are fit on too few stages

Only **17 of 24** cell types have ≥3 stages. Seven have **two** — and over two
points Spearman ρ is **±1 by construction**. The shipped
`length_trend.json` for `MESENCHYAL` reports `spearman = -1.0` on `n_stages = 2`;
that is not evidence. With three stages ρ can take only four values
(±0.5, ±1), so ρ = −1.0 at p = 0.069 is barely more informative.

### 4.2 Mean PDUI *is* sequencing depth

PDUI is 0 whenever a gene's distal PAS has no reads in that cell. In this data
**95–99 % of gene–cell pairs have no reads at all**:

| cell type | stage | mean PDUI | mean PDUI (informative) | frac uninformative | mean reads |
|---|---|---|---|---|---|
| TREG | Normal | 0.01016 | 0.409 | 0.975 | 0.107 |
| TREG | StageI | 0.00590 | 0.382 | 0.985 | 0.067 |
| TREG | IVprimary | 0.00617 | 0.399 | 0.985 | 0.067 |
| TREG | Met | 0.00462 | 0.359 | 0.987 | 0.056 |

Mean PDUI is arithmetically (1 − frac_uninformative) × mean informative PDUI, and
informative PDUI is nearly flat. So the per-stage mean tracks the **detection
rate**. Across cell types, the correlation between per-stage mean PDUI and
per-stage mean reads is:

**median r = 0.996** — 14 of 17 interpretable cell types above 0.95, 11 above 0.99
(TREG 0.998, MONOCYTIC 0.999, ENDOTHELIAL 0.9998, CYTOTOXIC_T_NK 0.9997).

**The "3′UTR shortening" curve and the sequencing-depth curve are the same curve.**

### 4.3 The direction is not robust, and nothing is significant

Refitting on informative pairs only (`slope_informative`) flips the sign in
**7 of 17** interpretable cell types (41 %) — AE2, CYTOTOXIC_T_NK, EPITHELIAL,
MACROPHAGE_M1, MONOCYTIC, PERICYTES, SEROUS_CELLS. Slopes are ~0.001–0.01 PDUI
units per stage on a 0–1 scale. Across 17 cell types, exactly **one** nominal
p < 0.05 (SEROUS_CELLS, p = 0.044) — which is what 17 tests at α = 0.05 produce by
chance, and it does not survive any correction.

**Verdict.** There is **no defensible evidence of stage-dependent 3′UTR
shortening in this cohort**. The apparent trend is a detection-rate artefact of
the classic PDUI estimator under 95–99 % dropout. This must replace the current
claim. To ask the question properly you need either an estimator that conditions
on detection (report `mean_pdui_informative` with per-stage depth as a covariate),
or depth-matched subsampling across stages — neither of which this sweep ran.

---

## 5. Are the differential-APA hits trustworthy?

**Experiment.** `ema switch diff` per cell type: `fisher` (exhaustive within-gene
screen) and `nb_multi` (omnibus LRT). `nb_pairwise` was **not run**.

**Structural caveat, restated because it governs every hit count.** Fisher's
exact test is applied to **read** counts pooled over cells, so each cell
contributes many correlated observations. Significance therefore scales with cell
count rather than effect size — classic pseudoreplication. `fisher_power_check`
and `pseudoreplication_signature` quantify this as Spearman(n_cells,
n_significant), and F5.3 plots it with the median |Δ proportion| of the hits as
point size, so the "many hits, tiny effects" pattern is visible.

Two further limits belong next to any hit list:

- Contrasts are **incomplete and unbalanced**: only 13–18 cell types have each
  of the six stage pairs, so a gene "found in Met vs Normal but not Met vs StageI"
  may simply lack the comparison.
- **`ema switch match` output is unusable.** `cluster_match.tsv` reports
  `match_confidence = 0.0102` for *every* pair, and each cluster matches
  essentially every other cluster — the `marker_overlap` matching is saturated.
  **`canonical_cluster` must not be used** as a cross-dataset identity.

The recurrent-vs-private split (`recurrent_switches`, F5.4) is the right way to
read the hit lists: a gene called in one cell type is as likely a power artefact
as biology, while genes recurring across independent cell types *with a consistent
direction* are the defensible candidates.

---

## 6. What does each PAS filter cost downstream?

**Experiment (Phase 3).** Four scenarios off identical saved peak-calls —
`baseline` (keep everything), `atlas_filter`, `annot_filter_3utr`, `ip_filter` —
each flipping exactly one filter, then re-running clustering, GEX cell typing and
the switch tests.

This is the experiment that decides whether the atlas-snapping that dropped
43–68 % of PAS in the pre-correction era mattered scientifically or only
cosmetically. The framing that matters is **not** "how many PAS were dropped" but
**"did dropping them change a conclusion"** — measured as Δ clusters, Δ ARI, % of
cells whose cell-type call changes, Δ significant hits, and **trend sign flips**.

**Status: PENDING.** The server has only `baseline` and `atlas_filter` directories,
from an aborted earlier launch (killed 18:42, and restricted to a single GSM
rather than the intended 6-GSM subset). No `B2_gex_celltyping` or `B3_switch`
outputs exist for any scenario. `filter_effect_deltas()` is implemented and wired
to emit all five delta tables the moment the re-launch completes; nothing from
this axis should enter the report before then.

---

## 7. What to change in the current report

| current claim | replace with |
|---|---|
| "precision ≈ 1.0 against PolyASite v3.0" | excess-over-null +0.245 at ±50 bp; report ±10–25 bp with the null (§2) |
| strategy ranked by atlas score | strategies are indistinguishable on the metric and share as little as 29 % of PAS (§1, §2) |
| "global 3′UTR shortening across stages" | not demonstrable; mean PDUI ≡ depth at median r = 0.996, 41 % sign flips, nothing significant (§4) |
| per-cell-type Spearman ρ from `length_trend.json` | ρ is ±1 by construction for the 7 two-stage cell types (§4.1) |
| trim-parameter yield comparisons | blocked by a concurrency race larger than the effect (§3) |
| `canonical_cluster` / `switch match` | unusable — matching is saturated (§5) |
| — *(add)* | TF-IDF beats library-size normalisation by ΔARI 0.140 (§3a) — the headline positive result |
