# PRIME_PLAN — branch `peakAtail-prime`

Branch cut from `origin/develop` at **`9dfdefb`** (= Stage 1b/1c + the IP-strand fix #96 + the
performance work #97). `9dfdefb` is called **v2** throughout; it is the code that produced the
numbers currently in the manuscript.

**Goal.** Lift PeakATail's precision/recall curve. The tool already leads all de novo callers on
precision (PBMC P@100 0.7062) but its detected-gene recall (R_det@100 0.1754) sits below polyApipe
(0.199) and far below scAPAtrap (0.300). Filtering only slides along the curve; these five changes
are meant to move it.

---

## Non-negotiables for every commit on this branch

1. **Every behavioural change is behind a CLI flag / config option.** No unconditional behaviour
   change, ever.
2. **v2 output stays reproducible byte-for-byte** through explicit flags, and through the single
   compatibility switch `--compat v2` (Change 0). A regression test asserts this on a committed
   fixture; `scripts/prime/identity_check.py` re-validates it on the real PBMC chr19+21 slice.
3. **No new hard runtime dependency.** Anything model-shaped is trained OFFLINE and shipped as
   coefficients/thresholds evaluated with numpy. scikit-learn must not enter `ema`'s import path.
4. **Test suite stays green**: baseline on this branch is `1,277 passed, 9 skipped, 4 xfailed,
   1 xpassed` with exactly 2 known environment failures in `tests/test_pyproject_install.py`.
   Any new failure is a stop signal.
5. **Defaults are justified by a measurement, named in the commit message.** Where the measurement
   says an idea does not pay, the flag ships **defaulted OFF** and this file says so.
6. Success is defined in `manuscript/24_prime_preregistration.md`, written before any prime accuracy
   number existed. Read it before changing a default.

Evidence base: `results/algo_headroom/{A1_end_pileup,A2_clip_sensitivity,A3_scoring_model,A4_resolution}/`.
Every number quoted below is from those measurements, with its table named.

---

## Change 0 — compatibility mode and the identity harness  *(prerequisite, no behaviour change)*

**What.** Add `--compat v2` (config `compat: v2`): one switch that pins every prime option
introduced by Changes 1–5 to its v2 value, so a reviewer reaches v2 behaviour without reconstructing
a flag list. Add a regression test that runs the caller twice on a small committed fixture — once
with prime defaults, once with `--compat v2` — and asserts the `--compat v2` outputs are
byte-identical to a checked-in v2 expectation.

**Flags.** `--compat {off,v2}`, default `off`.
**Default on prime.** `off` (prime defaults active).
**Verification.** `scripts/prime/identity_check.py <v2_run> <prime_compat_run>` must print
`VERDICT: IDENTICAL` on the PBMC chr19+21 slice. Done once now on the unmodified branch (the
harness-validation run) and again after each of Changes 1–5.

**Why first.** Every later change is only reviewable if the baseline is reproducible. This is also
the only change that must never be skipped.

---

## Change 1 — read acceptance geometry: stop truncating and discarding reads

**Where.** `ema/countmatrix/read.py` (`read_check`).

**Today.** A read whose reference span exceeds `--seq-len` is **discarded**; a shorter read has its
end coordinate **rewritten** to `read_start + seq_len`. Downstream geometry is therefore partly
synthetic.

**Measured cost of today's rule.**
* `A2/tables/T0_headline_numbers.tsv`: the span rule discards **13.74 %** of valid-CB reads
  genome-wide (88,778,573 / 646,342,981); **96.0 %** of the discards are spliced alignments. It
  throws away **56,335** reads that carry a qualifying poly(A) clip — **+1.76 %** on top of the
  3,195,067 clip reads the caller keeps.
* `A1/tables/stage0_seqlen_acceptance.txt`: 4.8 % of CB reads dropped on chr21, 23.8 % on
  chr19:1–20 Mb, plus a **fabricated end coordinate on a further 8–12 %** of the reads it keeps.

**Measured benefit of removing it.** `A2/tables/T13_headroom_at_matched_precision.tsv`, variant
`k6_nospan`, recall at matched precision: atlas **+0.7 % to +1.4 %** relative (P 0.40–0.70);
Kinnex **+1.9 % to +2.1 %**. Small, but positive on **both** truths at every matched precision
tried — the only variant in that table with no negative entry.

**Flags.** `--read-span-mode {truncate,exact}` (config `read_span_mode`).
* `truncate` — v2: discard span > `--seq-len`, rewrite `read_end = read_start + seq_len`.
* `exact` — keep every accepted read, use its true aligned span.

**Default on prime.** `exact`. **v2 reachable via** `--read-span-mode truncate` (and `--compat v2`).

**Watch for.** This changes coverage, so it changes *everything* downstream — peak boundaries,
tier-2, the matrices. Land it alone, measure it alone, and expect the byte-identity check to show
differences everywhere in prime mode and none in compat mode. Also raises peak RSS and read counts;
check the compute guard rail (2× wall, 1.5× RSS) at this step, not at the end.

**Not in scope here** (measured, separate, both currently unfiltered on the coverage path):
`read_check` applies no `-F 3844`, so a 5-way multimapper contributes 5 coverage reads
(`A1` §10 item 7). Fix or flag it only with its own measurement.

---

## Change 2 — clip-detector sensitivity, gated on genomic A content

**Where.** `ema/countmatrix/polya.py` (`clip_site`, `ClipAccumulator`).

**Today.** A terminal soft clip counts as poly(A) evidence at `--polya-min-clip 6` and
`--polya-min-purity 0.8`. Only ~0.57 % of read_check-accepted CB reads qualify genome-wide
(`A2/T0`: 3,195,067 / 557,564,408). At Kinnex-confirmed isolated atlas sites, **38.4 %** of reads
terminating within ±3 bp carry a qualifying clip, but a further **13.9 %** carry a short 1–5 nt A
tail the detector rejects — and that short-tail class is **22.8×** enriched at verified PAS over
random genic points (`A2/T0`, `T12`).

**The measured trap.** Relaxing `--polya-min-clip` **alone is negative**
(`A2/T13`, recall at matched precision P 0.50, relative to k6): k5 −0.9 %, k4 −0.4 %, k3 −2.8 %,
k2 −7.5 %, k1 −28.1 %. Wrong-end specificity collapses with k (86.7× at k6 → 38.1× at k3, `A2/T0`).
**Do not relax k on its own.**

**What pays.** The same relaxation *plus* a genomic-A gate (`A2/T11_genomicA_gate`, `A2/T13`):

| variant | atlas R at matched P 0.60 | at matched P 0.70 | Kinnex R at matched P 0.60 |
|---|---|---|---|
| `cur_k6` (today) | 0.2094 | 0.1712 | 0.4343 |
| `cur_k6_ipgate` | 0.2215 (**+5.8 %**) | 0.1923 (**+12.3 %**) | 0.4591 (**+5.7 %**) |
| `k4_ipgate` | 0.2197 (+5.0 %) | 0.1885 (+10.1 %) | 0.4740 (**+9.2 %**) |
| `k3_ipgate` | 0.2170 (+3.6 %) | 0.1846 (+7.8 %) | 0.4785 (**+10.2 %**) |
| `k4` (no gate) | 0.2071 (−1.1 %) | 0.1635 (−4.5 %) | 0.4246 (−2.2 %) |
| `*_tailgate` (tail-composition gate) | ≈ 0 everywhere | ≈ 0 | ≈ 0 |

The gain is in the **gate**, not in the relaxation; the relaxation only becomes safe once the gate is
there. The tail-composition gate does nothing and is not worth implementing.

**Flags.**
* `--polya-min-clip INTEGER` — already exists; prime may change its **default** (candidate 4, to be
  decided by this branch's own measurement under Rule T, not by the table above).
* `--polya-genomic-a-gate / --no-polya-genomic-a-gate` (new) — reject a clip whose downstream
  genomic context is A-rich. Requires `--genome-fasta`; when the FASTA is absent the gate must
  **disable itself with a loud warning**, never silently change behaviour.
* `--polya-a-gate-window 'UP,DOWN'` and `--polya-a-gate-frac FLOAT` for the gate's geometry, with the
  A3-verified tool IP rule (6-A run or A-fraction ≥ 0.7 over transcript-relative −9..+30) as the
  starting point — that rule reproduces the tool's own flag on 402,765/402,765 candidates
  (`A3/PROVENANCE.txt`).

**Default on prime.** Gate **ON** when `--genome-fasta` is given, OFF otherwise; `--polya-min-clip`
default decided by measurement, defaulting to 6 (v2) until then.
**v2 reachable via** `--no-polya-genomic-a-gate --polya-min-clip 6`.

---

## Change 3 — internal priming as a scored covariate, not a hard veto

**Where.** `ema/experimental/internal_priming.py`, `ema/countmatrix/paswrite.py`.

**Today.** `--ip-filter --ip-filter-mode filter` applies the IP rule as a **hard veto**. 27 % of our
de novo sites sit within 25 bp of a long-read internal-priming decoy, and the veto is the largest
single recall sacrifice in the default.

**Measured.** In A3's chromosome-disjoint model the **binary** tool IP flag is nearly worthless
(`ip_tool` coefficient −0.029 sd-units) while the **continuous** A-fraction it is derived from is the
second-strongest sequence feature (`ip_tool_afrac` −0.918) — i.e. the information is in the degree of
A-richness, and thresholding it into a veto throws most of it away
(`A3/tables/logreg_coefficients.tsv`, `ip_veto_vs_covariate.tsv`, `expr_matched_ipveto_arms.tsv`).

**Flags.** extend `--ip-filter-mode` to `{annotate,filter,score}`.
* `filter` — v2 hard veto.
* `annotate` — v2 annotate-only.
* `score` (new) — emit the A-fraction and A-run as per-site covariates and let Change 4's score use
  them; no site is deleted for IP alone.

**Default on prime.** `score` when `--pas-score calibrated` is active, else `filter`.
**v2 reachable via** `--ip-filter-mode filter`.

**Note.** Change 3 is only meaningful together with Change 4 — a covariate with nothing to feed is a
no-op. Land it immediately before Change 4, or as its first commit.

---

## Change 4 — a calibrated per-site score in place of the hard ≥2-molecule threshold  *(the big one)*

**Where.** new `ema/countmatrix/pas_score.py`; consumed by `ema/countmatrix/paswrite.py`.

**Today.** Tier-1 sites are kept on a hard rule: `tier1 ∩ IP-pass ∩ ≥2 distinct (CB,UMI) molecules`.
**72 %** of clip-supported sites carry exactly one clip molecule and are discarded by that rule; the
coverage-only tier-2 (166,355 sites at P@100 0.0568) is kept but near-worthless.

**Measured headroom** (`A3/tables/VERDICT_headline.tsv`, chromosome-disjoint 2-fold CV, PBMC).
`R_det@100` at the v2 default's own precision (P@100 = 0.7062), and the atlas-independent Kinnex
read-out at matched truth:decoy odds:

| scorer | R_det@100 @ P 0.7062 | vs rule | Kinnex R(t20) at matched odds | vs rule |
|---|---|---|---|---|
| RULE `tier1 ∩ IP-pass ∩ ≥k mol` (v2) | 0.1754 | — | 0.2723 | — |
| soft score, **same evidence the rule uses** | 0.1801 | +2.7 % | 0.2774 | +1.9 % |
| + sequence (hexamer, downstream A) | 0.2210 | +26 % | 0.3022 | +11 % |
| + coverage + local context (**no annotation**) | 0.2286 | **+30 %** | 0.3081 | **+13 %** |
| + GTF annotation (full) | 0.2445 | +39 % | 0.3156 | +16 % |

At matched *call count* (n = 46,524) the no-annotation model reaches P@100 0.826 / R_det 0.1955
against the rule's 0.706 / 0.1754 — up **and** right (`A3/tables/headline_arms.tsv`), with the Kinnex
decoy rate **falling** 0.1300 → 0.0774. `A3/tables/single_molecule_rescue.tsv`: at matched precision
the score selects 79,484 sites (R_det 0.2609) of which **35,995 are single-molecule sites** rescued
at P@100 0.574 — the 72 % the hard rule discards is where the recall is.

**Model choice.** A3's headline numbers come from a gradient-boosted model in a `.joblib`. That
cannot ship (constraint 3). Decide on this branch, by measurement: (a) does the logistic model in
`A3/tables/logreg_coefficients.tsv` retain enough of the gain? if yes, ship its coefficients; (b) if
not, export the tree ensemble to a plain array format evaluated with numpy, or accept the logistic
loss and **say so in the manuscript**. Prefer the **no-annotation** feature set: it is within 0.008
Kinnex-R of the full model, needs no GTF at score time, and is one less circularity worry.

**Flags.**
* `--pas-score {none,calibrated}` (config `pas_score`) — `none` = v2.
* `--pas-score-min FLOAT` — the default operating point.
* `--pas-score-model PATH` — override the shipped coefficients (offline-trained; default = shipped).
* `--polya-min-umis` keeps working and keeps its v2 default.

**Default on prime.** `calibrated`, with `--pas-score-min` fixed under **Rule T** of
`manuscript/24_prime_preregistration.md` §3.3: **threshold chosen on GSE104556 mouse 1, evaluated
unchanged on PBMC and mouse 2**. Never chosen by looking at the PBMC number being reported.
**v2 reachable via** `--pas-score none`.

**Also here.** Write the per-site score into `pas_support.tsv` and column 5 of `pasbed.bed`'s
sidecar (not column 5 of `pasbed.bed` itself — that is the molecule count the launcher's
`$5>=2` filters depend on, and changing it silently breaks every downstream script in
`scripts/benchmark_tools/`).

---

## Change 5 — cleavage-site resolution: a calibrated constant offset, and NO deconvolution

**Where.** `ema/countmatrix/cleavage_offset.py`, `ema/countmatrix/polya.py` (mode selection).

**Measured, and mostly negative — read this before writing code.**

* **Deconvolution / peak splitting does not work.** Splitting every bimodal call raises PBMC P@100
  0.7261 → 0.7617 on the slice — but a control that takes the split geometry from *other* calls
  reproduces +0.0368 of that +0.0356, and a control that duplicates each call at mode−1 / mode+1
  (zero information) reproduces +0.0365 (`A4/T4`, `A4/T4b`). The apparent gain is an artefact of
  emitting more points. The oracle recall gain is R_det@25 +0.0620, but only **1.54 %** of the
  blocked sites have even one clip molecule there, so the realistic gain is **+0.0010**
  (`A4/T7_headline_headroom.tsv`). **Do not implement splitting.**
* **A constant offset is nearly flat, and the two truths disagree.** Median signed distance from
  call to nearest Kinnex 3' end is **−2 bp** on both strands (n = 40,052; `A4/T5`). Best constant
  shift: P@100 +0.0003 (flat — no window ≥ 25 bp responds to any shift in −6..+5); P@10 +0.0215 at
  **+3 bp** on the atlas but the Kinnex optimum is **−5 bp** — *opposite signs*. Only at bp-exact
  matching do the two truths agree, both on **−2 bp** (P@1 0.2119 → 0.2960 atlas,
  0.1782 → 0.2642 Kinnex).
* **`--auto-cleavage-offset` is a trap for the clip-seeded caller.** Enabled, it estimates **+95 bp**
  on this library, which costs P@10 0.5209 → **0.0551** and P@100 0.7062 → 0.6655 (`A4/T5`, `A4/T7`).
  It is off by default (`config.auto_cleavage_offset = False`) and must stay off; add a **refusal or
  loud warning** when it is combined with `--peak-strategy clip_seeded`.
* **The one positive**: an edge-boundary criterion on existing tier-1 calls gives P@10/P@100
  retention **0.937 vs 0.864** against long reads, MAD 1 bp vs 2 bp, on ~17 % of calls (480 of 2,895)
  with a calibratable ~3 bp upstream bias (`A1/tables/stage8_resolution.tsv`, `A1` §9). A resolution
  refinement on a high-confidence subset — not a recall gain.

**Flags.**
* `--cleavage-offset INTEGER` — already exists, default 0. Prime keeps **0** unless this branch's own
  measurement (under Rule T) justifies −2, and any change is reported as a resolution result only.
* `--cleavage-edge-refine / --no-cleavage-edge-refine` (new) — the A1 §9 edge-boundary refinement,
  applied only to calls that pass the boundary criterion.
* `--auto-cleavage-offset` — unchanged default (off); add the clip_seeded guard.

**Default on prime.** `--cleavage-offset 0` (= v2); `--cleavage-edge-refine` **OFF** until this branch
measures it, and OFF permanently if it does not clear its own bar. **v2 always reachable.**

---

## Explicitly NOT implemented, with the measurement that killed it

**An end-position / read-3'-end pileup detection channel.** `A1_end_pileup` measured this end to end
and the answer is no. **96.1–96.3 %** of accepted reads carry no 3'-side soft clip at all, so their
3'-end is `start + 91` — a shifted copy of the coverage the caller already models
(`A1/tables/stage0_read_end_structure.txt`). At the default's precision the best pileup-only rule
reaches R_det **0.064–0.081** against the default's **0.209** — 2.6–3.3× *less* recall — reproduced
on each chromosome separately and under the independent Kinnex truth. Calls a pileup channel would
*add* to the default top out at **36 %** precision, half the default's 73.9 %
(`A1/tables/stage5_added_call_precision.tsv`). And the one real end-position signal, the 3' boundary,
does not separate PAS from internal priming at all: median `edge100` **0.836** at Kinnex truth sites
vs **0.839** at Kinnex internal-priming decoys.

If a future agent implements it anyway, it ships **defaulted OFF** and cites this paragraph.

---

## Implementation order and gates between steps

| # | change | flag(s) | prime default | gate before moving on |
|---|---|---|---|---|
| 0 | compat mode + identity test | `--compat v2` | `off` | `identity_check.py` prints IDENTICAL on the slice; suite green |
| 1 | read acceptance geometry | `--read-span-mode` | `exact` | slice P/R measured and reported; compute guard rail checked; compat still byte-identical |
| 2 | clip detector + genomic-A gate | `--polya-genomic-a-gate`, `--polya-min-clip` | gate ON with FASTA | slice P/R; gate proven to disable itself without `--genome-fasta` |
| 3 | IP as covariate | `--ip-filter-mode score` | `score` with Change 4 | covariates present in `pas_support.tsv`; no site deleted for IP alone |
| 4 | calibrated per-site score | `--pas-score`, `--pas-score-min`, `--pas-score-model` | `calibrated` | threshold fixed on mouse 1 under Rule T; numpy-only at runtime; full-BAM run on all three datasets |
| 5 | cleavage resolution | `--cleavage-edge-refine`, `--cleavage-offset` | refine OFF, offset 0 | measured; stays OFF unless it clears its own bar |

After each step: `identity_check.py` on the slice in compat mode, the full test suite, and a slice
score through `scripts/prime/score_pas.sh`. Nothing lands with a red suite or a broken compat mode.

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
  `!tests/fixtures/*.bam` (line 22), which is why `tests/fixtures/cellranger_pbmc_tiny.bam` is the
  single tracked fixture. Change 0's byte-identity regression test needs a checked-in **expected
  output**; add a matching `!tests/fixtures/...` exception rather than sprinkling `git add -f`, so
  the fixture cannot be lost by a later `git clean`.
* **`--version` and every `run` invocation create `switch_combined/<log>` in the current working
  directory** (gitignored). Run the caller from an output directory, not from the worktree root —
  `scripts/prime/run_slice.sh` already `cd`s into the run directory.
* **Five worktrees of this package exist.** `PYTHONPATH` alone is not proof you imported the right
  one; assert `ema.__file__` at the start of every run. `scripts/prime/common.sh::assert_code` does
  it for three modules at once.
