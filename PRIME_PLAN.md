# PRIME_PLAN — branch `peakAtail-prime`

Branch cut from `origin/develop` at **`9dfdefb`** (= Stage 1b/1c + the IP-strand fix #96 + the
performance work #97). `9dfdefb` is called **v2** throughout; it is the code that produced the
numbers currently in the manuscript.

**Goal.** Lift PeakATail's precision/recall curve. The tool already leads all de novo callers on
precision (PBMC P@100 0.7062) but its detected-gene recall (R_det@100 0.1754) sits below polyApipe
(0.199) and far below scAPAtrap (0.300). Filtering only slides along the curve; these changes are
meant to move it.

> **REVISION 2 (2026-08-21, second foundation pass).** Revision 1 of this file was written from the
> raw A1–A4 headroom reports. An adversarial verifier has since **withdrawn three of those reports'
> positive claims**, and `results/algo_headroom/VERIFY/VERIFICATION_REPORT.md` +
> `manuscript/23_algorithm_roadmap.md` now supersede A1–A4 wherever they disagree. Four changes in
> Revision 1 were wrong as written and are corrected here:
>
> 1. Revision 1 planned to default internal priming to a **scored covariate** (`--ip-filter-mode
>    score`). The verifier's finding is the opposite: **the IP veto is the single largest measured
>    lift in the whole programme (+7.7 % to +12.8 % relative recall at matched atlas precision) and
>    must stay a HARD GATE.** Every configuration that let a score override the veto looked
>    spectacular on the atlas and no better on long reads.
> 2. Revision 1 quoted A3's **+26 % / +30 % / +39 %** score gains. Those are atlas-trained and
>    atlas-scored, i.e. circular, and are **withdrawn**. The honest, verified bracket is **+6.6 %
>    relative recall with +7.0 precision points at matched call count**, or **+14.3 % relative recall
>    at matched atlas precision** — one third the advertised size.
> 3. Revision 1 listed A1's `edge100` boundary refinement as "the one positive" of the resolution
>    work. **Withdrawn**: at the same call count, simply keeping the calls with the most clip
>    molecules — a number the caller already writes to `pas_support.tsv` — beats it on independent
>    long-read precision (P@10 0.7765 vs 0.6959, P@100 0.8410 vs 0.7506).
> 4. Revision 1 put the read-geometry change first. The measurement says **one change dominates**
>    (the calibrated score) and one is free (`--ip-filter` default-on), so the order now follows
>    `23_algorithm_roadmap.md` §3.
>
> Number policy for this file: **[V]** = confirmed by the verifier; **[A*]** = measured by A1–A4 and
> *not* independently re-verified; **[19]/[24]** = from the manuscript's gate documents. Anything
> unlabelled is a design statement, not a measurement.

---

## Non-negotiables for every commit on this branch

1. **Every behavioural change is behind a CLI flag / config option.** No unconditional behaviour
   change, ever.
2. **v2 output stays reproducible byte-for-byte** through explicit flags, and through the single
   compatibility switch `--compat v2` (Change 0). `tests/test_prime_v2_compat_golden.py` asserts it
   on the committed fixture in the normal suite; `scripts/prime/identity_check.py` re-validates it
   on the real PBMC chr19+21 slice after every behavioural change.
3. **No new hard runtime dependency.** Anything model-shaped is trained OFFLINE and shipped as
   coefficients/thresholds evaluated with numpy. scikit-learn must not enter `ema`'s import path.
4. **Test suite stays green.** At the branch cut: `1,277 passed, 9 skipped, 4 xfailed, 1 xpassed`
   with exactly 2 known environment failures in `tests/test_pyproject_install.py`. The compat
   golden test adds 4, so the current baseline is **1,281 passed / 2 known failures**. Any *new*
   failure is a stop signal.
5. **Defaults are justified by a measurement, named in the commit message.** Where the measurement
   says an idea does not pay, the flag ships **defaulted OFF** and this file says so.
6. Success is defined in `manuscript/24_prime_preregistration.md`, written before any prime accuracy
   number existed. Read it — especially §3.1 (the primary criterion) and §3.3 (Rule T, threshold
   selection) — before changing a default.

Evidence base, in precedence order:
`results/algo_headroom/VERIFY/VERIFICATION_REPORT.md` > `manuscript/23_algorithm_roadmap.md` >
`results/algo_headroom/{A1_end_pileup,A2_clip_sensitivity,A3_scoring_model,A4_resolution,A5_cross_sample}/`.

## Three facts that bound everything below — read before proposing anything new

* **The ceiling.** Of the 11,606 chr19+21 detected-gene atlas sites the v2 default misses at 100 bp,
  only **21.4 %** have *any* candidate of any tier within 100 bp, and only **13.0 %** have a tier-1
  & IP-pass candidate [V §3]. **Four in five sites we miss have nothing to promote, re-score, split
  or relax into existence.** No re-ranking change can exceed that bound.
* **Do not add the gains up.** In their headline configurations the ideas below recover
  162 / 608 / 474 / 0 of those missed sites: naive sum 1,244, true union **1,040** (16 % double-
  counted). In their most aggressive configurations, naive sum 3,584, true union **2,360** = 20.3 %,
  which is already **95 % of the ceiling** [V §3]. They compete for the same small pool.
* **The dev slice flatters the caller.** chr19+21 is 27 % clip-richer (0.7254 % vs 0.5730 %
  genome-wide), 18 % tier-1-richer and **19 % easier on recall** (R_det 0.2085 vs 0.1754) [V §5].
  A negative measured on the slice is conservative; **any positive measured only on the slice is
  optimistic and must be discounted.** Say which you have, every time.

---

## Change 0 — compatibility mode and the identity harness  *(prerequisite, no behaviour change)*

**What.** Add `--compat v2` (config `compat: v2`): one switch that pins every prime option
introduced by Changes 1–5 to its v2 value, so a reviewer reaches v2 behaviour without reconstructing
a flag list.

**Flags.** `--compat {off,v2}`, default `off`.
**Default on prime.** `off` (prime defaults active).

**Status: half done.** The regression test exists and is green *before* any behavioural change:
`tests/test_prime_v2_compat_golden.py` runs the real caller over
`tests/fixtures/cellranger_pbmc_tiny.bam` on two arms (`clip_seeded` at `--seq-len 91`, the
manuscript geometry; `original` at the shipped `--seq-len 150` backstop) and asserts the SHA-256 of
every output byte against goldens generated from the **frozen v2 worktree**
`tools/pa-polya-run-9dfdefb3`. Its `_v2_settings()` function is the single place where each new
prime option is pinned to its v2 value, and `--compat v2` must be exactly equivalent to it.
Demonstrated to bite: simulating Change 1 without a flag (deleting the `span > seq_len` rule in
`read.py`) fails both arms.

**What is left.** The `--compat v2` switch itself, and extending the golden test to run each arm
twice — prime defaults and `--compat v2` — once there is a behavioural difference to see.

**Verification after every later change.**
`scripts/prime/identity_check.py <v2_run> <prime_compat_run>` must print `VERDICT: IDENTICAL` on the
PBMC chr19+21 slice. Validated on the unmodified branch: 54 files, 50 SAME + 3 SAME_NORM +
1 SAME_LOG, exit 0; and it catches a single changed character in a same-length data file, an extra
output file, a drifted setting inside `run_config.json`, and a deleted file
(`results/prime/identity_*.txt`, `results/prime/harness_validation.tsv`).

**Why first.** Every later change is only reviewable if the baseline is reproducible. This is the
one change that must never be skipped.

---

## Change 1 — `--ip-filter` default-on, and the clip-rate constant  *(free; the largest measured lift)*

**Where.** `ema/cli/config_schema.py:529` (`--ip-filter`, currently `is_flag=True`, opt-in);
`ema/countmatrix/polya.py:10` and `:1038` (the clip-rate constant and the QC estimator).

**Measured [V §7.1, `VERIFY/tables/v8_ipveto_value.tsv`].** Sweeping `tier1 & >=k molecules` with and
without the veto and interpolating recall at matched precision:

| matched precision | recall WITH veto | WITHOUT | relative gain |
|---|---:|---:|---:|
| atlas P@100 = 0.60 | 0.2033 | 0.1887 | **+7.7 %** |
| atlas P@100 = 0.7062 (the default) | 0.1754 | 0.1612 | **+8.8 %** |
| atlas P@100 = 0.8358 | 0.1419 | 0.1258 | **+12.8 %** |
| Kinnex KinP_t5 = 0.7647 | 0.2722 | 0.2313 | **+17.7 %** |
| Kinnex KinP_t5 = 0.85 | 0.2422 | 0.1971 | **+22.9 %** |

It *lifts* rather than slides because it brings in information the thresholds do not have — genomic
sequence. It is bigger than every detector change tested, combined, and it is already implemented.

**Flags.** `--ip-filter` becomes **default-on when `--genome-fasta` is supplied**; without a FASTA it
must **disable itself with a loud warning**, never silently. `--no-ip-filter` restores opt-out.
`--ip-filter-mode` keeps its v2 values `{annotate,filter}` and its v2 default `filter`.

**Default on prime.** ON with a FASTA. **v2 reachable via** `--no-ip-filter` (and `--compat v2`).
This does *not* change the manuscript's pre-registered arm, which already runs with `--ip-filter`
[19]; it changes what a user gets who does not read the flag list.

**Sub-item, same commit, no flag needed (log/doc only, but declare it).** The docstring constant
"1.152 % of CB reads" is **wrong**: the genome-wide qualifying poly(A) clip rate is **0.5730 %**
(3,195,067 / 557,564,408 accepted CB reads), independently reproduced from a from-scratch
reimplementation on chr21 at 0.003669 [V §5]. The head-sampling QC estimator
(`polya.py::estimate_clip_rate`) reports **2.2565 %** where the truth on the same library is
**0.5364 %** — 4× off, in the direction that would mask a genuinely destroyed evidence channel.
Correct the constant, fix or remove the estimator (GitHub #99 §2 is its open issue).

**Do NOT** turn the veto into a score. See Change 2 and the Revision-2 note at the top.

---

## Change 2 — a calibrated per-site score as a RE-RANKER inside the existing gates  *(the big one)*

**Where.** new `ema/countmatrix/pas_score.py`; consumed by `ema/countmatrix/paswrite.py`.

**Today.** Tier-1 sites are kept on a hard rule: `tier1 ∩ IP-pass ∩ ≥2 distinct (CB,UMI) molecules`.
**72 %** of clip-supported sites carry exactly one clip molecule and are discarded by that rule [19];
the coverage-only tier-2 (~166 k sites at P@100 0.0568 [19]) is kept but near-worthless.

**What changes: the `≥2` threshold, and ONLY that.** Tier-1 membership and the IP veto stay hard
gates in front of the score. This is not a stylistic preference — every configuration in which the
score was allowed to override the IP veto looked spectacular on the atlas and no better than the
rule on long reads [V §1.4, §7.1].

**Measured headroom [V §1.5 and §7.2]** — a gradient-boosted model trained **only on the Kinnex
long-read truth**, all 12 downstream-A features dropped, chromosome-disjoint folds, evaluated on the
PolyASite atlas it has never seen:

| arm | n | P@10 | P@25 | P@100 | R_det@100 | KinP_t5@25 | decoy@25 |
|---|---:|---:|---:|---:|---:|---:|---:|
| RULE `tier1 ∩ IP-pass ∩ ≥2 mol` (v2 default) | 46,524 | 0.5209 | 0.6389 | 0.7062 | 0.1754 | 0.7647 | 0.1300 |
| score, **at matched call count and matched expression** | 46,524 | **0.5567** | **0.6835** | **0.7766** | **0.1870** | — | **0.0954** |
| score, **at matched atlas precision** | 56,286 | 0.4934 | — | 0.7059 | **0.2005** | **0.8013** | 0.1048 |

**The honest bracket: +6.6 % relative recall at matched call count *with* +7.0 precision points, or
+14.3 % relative recall at matched precision.** Not +39 %, not +48 % — those figures are atlas-
trained and atlas-scored and are withdrawn [V §6(b)]. Two caveats that must travel with the number:
**P@10 does not improve at matched precision** (0.4934 vs 0.5209 — a re-ranker cannot move a call),
and the model was fitted and tested on **one PBMC library**; chromosome-disjoint CV tests positional
generalisation, not dataset or species transfer.

**Deliverable order is inverted on purpose — the transfer test comes BEFORE the implementation.**

1. **Transfer test first.** Refit the verifier's exact recipe (`VERIFY/code/v3_a3_decontaminate.py`:
   `HistGradientBoostingClassifier(max_iter=400, learning_rate=0.06, max_leaf_nodes=31,
   min_samples_leaf=100, l2_regularization=1.0)`, Kinnex-t5 label, all 12 downstream-A features
   dropped) on **GSE104556 mouse 1** and evaluate it unchanged on **PBMC** and **mouse 2** — the
   primary protocol of [24] §3.3. Mouse has no Kinnex truth, so the mouse-side label must be the
   mouse atlas, which makes the mouse→PBMC direction atlas-trained and therefore a generalisation
   test only, never a headline. **If it does not transfer, stop and report that as the result.**
2. **The bar the score must clear before it is worth shipping at all.** On separating genuine
   long-read termini from IP decoys, the tool's own **inverted `ip_tool_afrac` alone scores AUC
   0.7790, beating the entire A3 model's 0.7655** [V §6(c)]. A score that does not beat one
   covariate the caller already computes is not a score; it is a rename.
3. Only if it transfers: implement as **numpy constants evaluated inline** (non-negotiable 3),
   behind `--pas-score`, **defaulted OFF** until [24] §3.1 is met on all three datasets.
4. **Ship the probability as a column regardless of the default**, in `pas_support.tsv`. A
   calibrated per-site confidence is a better paper claim than a better set, and it turns `≥2` from
   a commitment into a default. Calibration against the training label is near-perfect
   (ECE 0.0038); against the *other* truth it is optimistic in the mid-range (predicted 0.6 →
   observed 0.48) [A3] — which is itself a publishable result about atlas benchmarking.

**Flags.**
* `--pas-score {none,calibrated}` (config `pas_score`) — `none` = v2.
* `--pas-score-min FLOAT` — the default operating point, chosen under **Rule T** ([24] §3.3):
  **fixed on GSE104556 mouse 1, evaluated unchanged on PBMC and mouse 2**, never by reading the
  number being reported. Kinnex truth selects nothing, ever.
* `--pas-score-model PATH` — override the shipped coefficients (offline-trained; default = shipped).
* `--polya-min-umis` keeps working and keeps its v2 default.

**Default on prime.** `calibrated` **only after** the transfer test passes and [24] §3.1 is met;
`none` until then. **v2 reachable via** `--pas-score none`.

**Also here.** Write the score into `pas_support.tsv`, **not** into column 5 of `pasbed.bed` — that
column is the molecule count every `$5>=2` filter in `scripts/benchmark_tools/` depends on, and
overwriting it silently breaks the whole benchmark suite.

**Step 0 of this change is LANDED: `--pas-features {off,on}`, default `on`** (TASK C). The 24
per-site covariates the score needs are now emitted into `pas_support.tsv` — 2 written by the caller
(`clip_positions`, `clip_span`), 15 from genomic sequence and 7 from the local candidate set. They
ride inside the internal-priming pass, so the genome is opened exactly as many times as before
(asserted, not argued), and they change no call: BED, matrix and every pre-existing sidecar column
are byte-identical, which `tests/test_prime_v2_compat_golden.py` now covers for the sidecar too.
`ip_tool_flag` / `ip_tool_afrac` / `ip_tool_arun` come from the string the veto itself tested, so
point 2 above ("a score that does not beat one covariate the caller already computes is a rename")
can be checked against the caller's own number rather than a re-derivation. See `CHANGELOG.md` and
`docs/cli/run.md` for the column reference, and `results/prime/TASK_C_pas_features.md` for the
runtime/RSS cost and the cross-check against the offline feature table.

**Do NOT build** (measured, [V], roadmap §3 Step 1.5): a second BAM pass for clip-cell counts or
end-pileup statistics (worth 0.000–0.002 held-out AUC), or a model that merely re-weights the
existing clip counts (+1.9 % to +5.8 %, inside the noise of a threshold change).

---

## Change 3 — read acceptance geometry: stop truncating and discarding reads  *(LANDED; default NOT moved)*

**STATUS: implemented, measured, and the default deliberately left at v2.** Flags are
`--read-geometry {fixed,keep,true}` (default `fixed`) and `--read-exclude-flags N` (default `0`),
not the `--read-span-mode {truncate,exact}` this plan originally proposed — a third value was added
because the ablation turned out to matter (see below). Evidence:
`results/prime/read_geometry_census_{pbmc,mouse1}_slice.tsv`,
`results/prime/taskA_read_geometry_slice.tsv`, `results/prime/taskA_read_geometry_mass.tsv`,
`results/prime/taskA_prereg_evaluation.tsv`, `results/prime/taskA_called_position_shift.tsv`,
`results/prime/identity_taskA_branch_default_vs_v2.txt`.

**Where.** `ema/countmatrix/read.py` (`read_check`).

**Today (`fixed`).** A read whose *reference* span exceeds `--seq-len` is **discarded**; a shorter
read has its end coordinate **rewritten** to `read_start + seq_len`.

**Census, re-measured on this branch (A2's genome-wide figures held up).** On the PBMC chr19+21
slice the discard removes **24.19 %** of valid-CB reads, **98.08 %** of them spliced, and costs
12,427 qualifying clip reads (**+4.65 %**); on a GSE104556 mouse1 chr18+19 slice, 19.96 % / 98.72 %
/ +1.54 %. The rewrite fabricates the 3' end of 9.77 % (PBMC, mean +13.53 bp) and 21.70 % (mouse,
mean +25.17 bp) of kept reads. Genome-wide: 13.74 %, 96 % spliced, 88.8 M reads — the slices are
spliced-richer, so slice figures over-state it.

**Two design decisions the measurement forced, both away from the obvious choice:**

1. **Acceptance is on the de-introned reference footprint, not on query length.** A query-length
   test is *not* a strict relaxation: a short alignment with a long terminal soft clip has
   span ≤ `--seq-len` (v2 keeps it) but query length > `--seq-len`. That is the shape of a poly(A)
   clip read. Caught live by `tests/test_polya_three_path_agreement.py`.
2. **N is removed from the span.** Over 7.82 M valid-CB slice reads the introns carry **10.02 Gb**,
   ~14× the real read mass over the same 105 Mb. Feeding genomic spans to the coverage state machine
   would turn it into a gene-body detector.

**Outcome against [24] §3.1** (needs ΔP@100 ≥ −0.005, ΔR_det ≥ +0.010, ΔF1 > 0 on every dataset),
default precision arm, `true` vs `fixed`:

| dataset | ΔP@100 | ΔR_det@100 | ΔF1 | precision of the ADDED calls | verdict |
|---|---:|---:|---:|---:|---|
| PBMC chr19+21 | −0.0026 | +0.0031 | +0.0030 | 0.519 (base 0.640, null 0.031) | fails (ii) |
| mouse1 chr18+19 | **−0.0292** | +0.0063 | +0.0042 | 0.368 (base 0.716, null 0.014) | fails (i) by 5.9×, and (ii) |

**FAIL ⇒ the default stays at v2**, exactly as §3.1 says. The ablation is the useful part: `keep`
carries essentially the whole precision loss (mouse −0.0283 of the −0.0292), so it is *stopping the
discard* that costs precision, not the 3'-end fix — `keep`→`true` is ΔP −0.0001 / ΔR +0.0001 (PBMC)
and −0.0010 / +0.0004 (mouse). The recovered spliced reads carry real evidence (12–26× the genic
null) but weaker evidence than what is already in hand, so they slide along the curve.

**What the flag is still for.** Raw count-matrix mass **+31.2 %** (PBMC) / **+23.5 %** (mouse),
~+15.9 % genome-wide. §3.1 is a *detection* criterion and does not see this at all. Anyone wanting
`true` as a default needs a quantification criterion first; that decision belongs to whoever owns
[24].

**Not a threat to the compatibility guarantee.** A no-flag run of the branch reproduces v2 on the
real slice: **all 50 data files byte-identical**, the only difference in the whole tree being the
two new keys in `run_config.json` / `run_manifest.json`, both recording the v2 value.

**Compute.** PBMC slice at 8 threads: wall 1.32×, peak RSS **1.53×** — through §3.2.4's 1.5× guard
rail, declared. (It was 2.21× until `ClipStream.finalize()` stopped round-tripping through Python
ints.) Mouse slice: 1.17× wall, 1.01× RSS.

**`--read-exclude-flags` (the multimapper item, measured as this plan required).** `read_check`
applies no `-F`, so a 5-way multimapper contributes 5 coverage reads and 5 matrix counts; 10.42 %
(PBMC) / 12.09 % (mouse) of valid-CB reads are secondary. `256` is a clean precision-for-recall
trade on the default arm — ΔP@100 +0.0125 / ΔR −0.0008 (PBMC), +0.0218 / −0.0031 (mouse) — i.e. a
move *along* the curve, which is what this branch is explicitly not trying to do. **Default 0.**

---

## Change 4 — clip detection: the genomic-A gate only, never the relaxation  *(default OFF)*

**Where.** `ema/countmatrix/polya.py` (`clip_site`, `ClipAccumulator`).

**RELAXING `--polya-min-clip` IS DEAD. Do not do it.** [V §2.2, roadmap row B] Lowering the
threshold costs **0–10 % recall at matched precision** and more than half of what it adds is internal
priming: the increment's truth:decoy ratio is **0.83 : 1** against 3.67 : 1 for the evidence we
already have — *worse than parity despite the truth class being 4.2× deeper than the decoy class*,
so depth-matching pushes it further down. The clip channel structurally tops out near **35 %** of the
detected atlas even at maximum sensitivity (`min_clip = 1`: 34.9 % of detected-atlas sites have ≥1
clip molecule within ±25 bp, 22.7 % have ≥2; the current criteria give 24.0 % / 15.3 %)
[A2 `T7_recovery_and_fp.tsv`, measured on the chr19+21 slice, which is 27 % clip-richer than the
genome — so the genome-wide ceiling is *lower* than this].

**What might still pay: the gate, not the relaxation.** A2 measured a genomic-A gate on the clip
call (reject a clip whose *downstream genomic* context is A-rich) worth +5.8 % at matched atlas
P@100 = 0.60 and +12.3 % at 0.70, with `k` unchanged at 6 [A2 `T11_genomicA_gate`, `T13` — **not
independently verified**]. **It is not monotone across operating points** — the same table gives
+7.1 % at P 0.40, **−0.6 % at P 0.50** and +0.8 % at 0.55 — which is a second reason to ship it OFF
and measure it here rather than believe the two good rows. Note also that it is the same information
as Change 1's IP veto applied one stage earlier, so **its gain overlaps Change 1's and the two must
never be summed** (see "Do not add the gains up"). Measure it *on top of* Change 1, not against v2.

**Flags.**
* `--polya-genomic-a-gate / --no-polya-genomic-a-gate` (new). Requires `--genome-fasta`; with no
  FASTA it must disable itself with a loud warning, never silently change behaviour.
* `--polya-a-gate-window 'UP,DOWN'`, `--polya-a-gate-frac FLOAT` — the gate's geometry, starting from
  the tool's own IP rule (6-A run or A-fraction ≥ 0.7 over transcript-relative −9..+30), which
  reproduces the shipped flag on 402,765/402,765 candidates [A3 `PROVENANCE.txt`].
* `--polya-min-clip` — **stays 6.** No prime default change.

**Default on prime.** **OFF**, pending this branch's own measurement on top of Change 1 under Rule T.
The evidence for it is [A2]-only and its gain is entangled with Change 1's. **v2 reachable via**
`--no-polya-genomic-a-gate --polya-min-clip 6`.

**Do NOT implement** the tail-composition gate: A2 measured it at ≈ 0 everywhere.

---

## Change 5 — cleavage resolution: an `inferred_cleavage` COLUMN, and nothing else

**Where.** `ema/countmatrix/paswrite.py` (the column), `ema/cli/config_schema.py:503` (a docstring
that must be retracted), `ema/countmatrix/cleavage_offset.py` (a guard).

**The only positive, and it is resolution-only [V, `VERDICT_FOR_PI.tsv`].** The median signed
distance from our call to the nearest Kinnex 3' end is **−2 bp** on both strands (n = 40,052). A −2 bp
constant improves **base-pair-exact** matching on *both* truths — atlas P@1 0.2119 → **0.2960**,
Kinnex P@1 0.1782 → **0.2642** — and moves **nothing** at W ≥ 25 (P@100 0.7064 → 0.7064). It is not
universal: on mouse the P@1 optimum is **−1 bp**, not −2 [A4], so it is dataset/chemistry-dependent,
and it costs a little at intermediate windows on the atlas (P@10 0.5209 → 0.5093) while gaining on
Kinnex (0.4979 → 0.5108) [A4].

**Therefore: emit it as an extra column, never as a move of the reported coordinate.** As a column it
changes no headline number and needs no re-verification; as a coordinate move it would change every
reported number and would need its own pre-registration ([24] §3.4).

**Retract a shipped recommendation, same commit.** `config_schema.py:503` suggests trying
`cleavage_offset = 95`. That offset was estimated for the *coverage* caller and is **destructive**
for the clip-seeded one: P@10 0.5209 → **0.0551**, P@100 0.7062 → 0.6655 [A4]. `--auto-cleavage-offset`
estimates **+95 bp** on this library, so it must stay off by default and should **refuse or warn
loudly** when combined with `--peak-strategy clip_seeded`.

**Flags.**
* `--emit-inferred-cleavage / --no-emit-inferred-cleavage` (new) — write `inferred_cleavage =
  cleavage − 2 bp` (transcript orientation) as an extra column in `pas_support.tsv`. Default **ON**
  (it moves no existing byte in `pasbed.bed`; it adds a column to the sidecar, so the compat mode
  must still be able to suppress it — that is what makes it a flag).
* `--cleavage-offset INTEGER` — already exists, default **0**, unchanged on prime.
* `--auto-cleavage-offset` — unchanged default (off) + the clip_seeded guard.

**Default on prime.** column ON, `--cleavage-offset 0`. **v2 reachable via**
`--no-emit-inferred-cleavage --cleavage-offset 0`.

---

## Explicitly NOT implemented, with the measurement that killed it

Record kept so nobody re-proposes them. All four were measured, and two got *more* negative under
verification.

* **An end-position / read-3'-end pileup detection channel.** A1 measured it end to end; the
  verifier then showed A1's null was not depth-matched (median local depth **2** molecule-ends per
  ±100 bp against 8 at missed atlas sites and 243.5 at recovered ones). After direct standardisation
  to the missed-atlas depth distribution, **a molecule-end pileup at a true missed PAS is exactly as
  likely as at a random position of the same local depth — 0.99×, 0.95×, 1.00×, 1.01× at the
  mol5 ≥ 2/5/10/20 bars — against 11.36× for the clip channel** [V §2.1]. The cause is structural:
  **96.1–96.3 %** of accepted reads carry no 3'-side soft clip at all, so their "3' end" is
  `start + seq_len`, a shifted copy of coverage [A1]. **Do not build it.**
* **Relaxing `--polya-min-clip`.** See Change 4.
* **Deconvolving merged calls.** Zero of the 11,606 missed atlas sites recovered. The scorer's
  matching is many-to-one — 3,058 slice atlas sites are recalled at 100 bp using only **2,050
  distinct calls** — and tripling the call set at (mode−1, mode, mode+1) leaves P@100 **exactly**
  unchanged at 0.7385 while *raising* P@10 [V §4]. Any gain splitting appears to show is the
  scorer's permissiveness. **Do not build it.**
* **The `edge100` 3'-boundary "resolution refinement".** Applied as A1 recommended (gate the 2,895
  slice default calls at `edge100 ≥ 0.85`, keeping 1,736) it gives Kinnex retention 0.9271 — but at
  the same call count, keeping the 1,736 calls with the **most clip molecules** gives retention
  0.9233 with far better independent-truth precision (Kinnex P@10 **0.7765 vs 0.6959**, P@100
  **0.8410 vs 0.7506**) [V §6(a)]. `edge100` is confounded with clip support (median clip molecules
  4 → 9 across the gate): it is a more expensive way to raise the molecule threshold.
* **A second BAM pass for cross-cell statistics** (distinct clip cells, top-cell share, end counts,
  pileup sharpness, end-position entropy): 0.000–0.002 held-out AUC [V].
* **Cross-sample corroboration in single-sample mode.** Real but modest and cohort-only, it fails
  [24] §3.1 against the shipped default (+7.8 % relative recall for −0.0075 precision) and cannot
  help PBMC, which is one library. It also cannot suppress residual internal priming, which
  replicates across samples by construction — the genome is the same in both mice
  [A5, unverified]. Belongs in the Stage-3 cohort path, not in the caller.

If a future agent implements any of these anyway, it ships **defaulted OFF** and cites this section.

---

## Implementation order and gates between steps

| # | change | flag(s) | prime default | gate before moving on |
|---|---|---|---|---|
| 0 | compat mode + identity/golden test | `--compat v2` | `off` | golden test green; `identity_check.py` IDENTICAL on the slice |
| 1 | `--ip-filter` default-on + clip-rate constant | `--ip-filter` / `--no-ip-filter` | ON with `--genome-fasta` | slice + full-BAM P/R unchanged for the pre-registered arm; warning fires without a FASTA |
| 2 | calibrated per-site score (re-ranker inside tier-1 ∩ IP) | `--pas-score`, `--pas-score-min`, `--pas-score-model` | OFF until the transfer test passes, then `calibrated` | **transfer test first** (mouse 1 → PBMC, mouse 2); beats inverted `ip_tool_afrac` AUC 0.7790; numpy-only at runtime; threshold under Rule T; [24] §3.1 on all three datasets |
| 3 | read acceptance geometry | `--read-span-mode` | `exact` if it clears §3.1 | full three-dataset re-run (it changes the matrix too); compute guard rail checked here |
| 4 | genomic-A gate on the clip detector | `--polya-genomic-a-gate` | **OFF** | measured *on top of* Change 1, not against v2 |
| 5 | `inferred_cleavage` column, offset retraction | `--emit-inferred-cleavage`, `--auto-cleavage-offset` guard | column ON, offset 0 | no existing byte of `pasbed.bed` moves |

After each step, without exception: `identity_check.py` on the slice in compat mode, the full test
suite, and a slice score through `scripts/prime/score_pas.sh`. **Nothing lands with a red suite or a
broken compat mode.**

The decision rule for prime's default as a whole — and the statement that **every prime number is
exploratory** with respect to the manuscript's pre-registered gates in `manuscript/13` — is in
`manuscript/24_prime_preregistration.md`. Prime cannot retroactively become "the pre-registered
default"; adoption needs its own pre-registration and its own verification pass.

---

## Repo gotchas the next agent will hit

* **`.gitignore` swallows this file.** Line 50 is `*PLAN*.md` (with `*plan*.md`, `*prompt*.md`,
  `*analysis*.md` around it — the rules that keep ad-hoc agent scratch out of the repo).
  `PRIME_PLAN.md` matches, so it was committed with `git add -f` at the PI's explicit request for
  this path. Once tracked, `.gitignore` no longer applies to it, so ordinary edits commit normally.
  If the maintainer would rather it did not sit against the repo's own rule, move it to
  `docs/prime_plan.md` — nothing references it by path.
* **`.gitignore` also swallows `*.bed`, `*.tsv`, `*.csv`, `*.mtx`, `*.h5ad`.** The only exception is
  `!tests/fixtures/*.bam`, which is why `tests/fixtures/cellranger_pbmc_tiny.bam` is the single
  tracked fixture — and why the compat golden test stores **hashes in the test file** rather than a
  checked-in expected output. If you ever need a real expected-output file, add a matching
  `!tests/fixtures/...` negation rather than `git add -f`, so it cannot be lost to `git clean`.
* **`--version` and every `run` invocation create `switch_combined/<log>` in the current working
  directory** (gitignored). Run the caller from an output directory, not from the worktree root —
  `scripts/prime/run_slice.sh` already `cd`s into the run directory.
* **Five worktrees of this package exist.** `PYTHONPATH` alone is not proof you imported the right
  one; assert `ema.__file__` at the start of every run. `scripts/prime/common.sh::assert_code` does
  it for three modules at once, and every harness run writes the assertion to `modpath.txt`.
* **`n` means two things.** `score_tool.py` reports `n_query` (lines in the BED) and `n_matched`
  (calls on a contig the reference set covers); the manuscript quotes `n_matched`. On mouse 1 the v2
  default arm is 26,263 lines but **26,255** scored. Both appear in the literature of this project;
  `scripts/prime/score_pas.sh` prints them as `n (scored / raw)`.
