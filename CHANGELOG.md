# Changelog

> **Next release: 0.3.0 (minor, not patch).** The `Unreleased` entries below
> change the flagless call set — the internal-priming veto now runs by default
> when a genome FASTA is present (~17% fewer calls) — so under SemVer at 0.x this
> is a minor bump, never a patch. `pyproject.toml`, `CITATION.cff` and the
> bioconda recipe are already at 0.3.0; convert these `Unreleased` headings to
> `## 0.3.0` at tag time.

## Unreleased — `nb_pairwise` dispersion floor (issue #94)

### Fixed

- **`nb_pairwise` no longer reports a q-value for a test whose dispersion hit
  the numerical floor.** The per-PAS NB dispersion alpha is clipped to
  `[1e-4, 10]`; when the estimate collapses onto the **lower** bound the GLM is
  effectively a Poisson fit, its Wald standard error is a lower bound, and the
  p-value is anti-conservative. The estimate is also taken under the *full*
  model, so sampling noise that mimics a group difference is absorbed by the
  fitted coefficient and pushes alpha down — precisely the tests that then look
  most significant. Under the 20-run label-permutation null of issue #94, 8.5 %
  of null tests hit the floor and those tests produced **67 % of all false
  q < 0.05 calls**. Such rows now get `qvalue = NaN`, so no `qvalue < fdr`
  filter can ever call them significant.

### Added

- **`dispersion_floored` column in the `nb_pairwise` per-pair output** (`True`
  when `dispersion <= 1e-4`). The row itself and its raw `pvalue` are still
  reported so the evidence stays inspectable; only the q-value is withheld.
  `NbPairwiseStrategy.test` logs a warning naming how many PAS were affected.
  The floored rows stay in the BH input so the family size *m* is unchanged —
  dropping them would have made every *other* q-value less conservative — and
  `switch diff`'s pooled BH recomputation for `within_utr`/`between_utr`
  preserves the withheld NaN.

### Changed

- **Documented that `nb_pairwise` q-values are not permutation-calibrated.**
  Withholding the floored rows removes the dominant source of anti-conservatism
  but does not by itself make the remaining q-values calibrated; that needs a
  permutation-calibrated q (or dispersion shrinkage toward a trend), which is
  still open on issue #94. `docs/strategies/diff.md` and
  `docs/cli/switch-diff.md` now say so.

## Unreleased — `switch diff --isoform-agg between_utr` row identity (issue #110)

### Fixed

- **`ema switch diff --isoform-agg between_utr` now populates `gene_id`.** Under
  that scope the unit of a row is a 3'UTR isoform, so `pas_id` reads
  `GENE::TRANSCRIPT` — which is not a key into either annotation lookup, both of
  which join on a PAS id. `gene_id` was therefore written **empty** for every
  row and the gene appeared only in `diff_group_id`. The failure was silent:
  joining `between_utr` output to `per_gene` output on `gene_id` — the obvious
  comparison, and the column was present — matched **nothing** and reported 0 %
  overlap rather than raising (the correct figure on the same data was 96 %).
  `gene_id` now carries the group's gene, i.e. the same value as
  `diff_group_id`.

### Changed

- **`between_utr` output no longer carries the `chrom`, `start`, `end` and
  `strand` columns.** A whole 3'UTR has no single cleavage position, so these
  were always blank; they are now **omitted from the TSV entirely** so that a
  join keyed on them fails with a missing-column error instead of quietly
  matching nothing. `per_gene` and `within_utr` output is unchanged — both keep
  PAS as the row unit, so their coordinates are meaningful and their `gene_id`
  was never blank. Consumers of `between_utr` TSVs that read coordinates
  positionally or by name must be updated; `ema.switch_test.long_output`
  (`findings_long`) already reads them defensively and is unaffected (it emits
  an empty `pas_uid` for such rows, exactly as it did for the blank columns).
## Unreleased — spliced 3'UTRs in `switch diff --isoform-agg`

### Fixed

- **`--isoform-agg within_utr` aborted on every run over a spliced 3'UTR, and
  `between_utr` silently double-counted one.** A spliced 3'UTR contributes one
  `three_prime_utr` record per exon, so a PAS falling inside more than one of
  them named the same `(gene, transcript)` repeatedly and was appended to that
  UTR group more than once. In `within_utr` the duplicate columns made
  `fisher`'s `int(agg1[p])` receive a Series and raise
  `TypeError: cannot convert the series to <class 'int'>`. In `between_utr`
  nothing raised — the UTR's reads were simply counted once per record, so the
  denominator was inflated (on the regression fixture: 32 reads reported where
  the honest total is 16). Members are now de-duplicated per PAS, keeping
  first-seen order.

  Reported in #109. The `NameError` that PR also describes was already fixed on
  `develop` in `7b50a5f`; this lands the remaining half.
## Unreleased — the command is now `peakatail` (`ema` deprecated)

### Breaking changes

- **The console script is `peakatail`.** The package has been named `peakatail`
  since 0.2.0 while the command it installed was still `ema`, which was a
  long-standing source of confusion in issues and in the docs.

### Deprecations

- **`ema` still works and is deprecated.** It is installed as an alias that
  behaves identically and prints a one-line notice to *stderr* — never stdout,
  so a caller parsing output is unaffected. Existing scripts, Snakemake and
  Nextflow rules keep running unchanged. It will be removed in a future release.
  Note that Click derives the program name from `argv[0]`, so `ema --version`
  and the `Usage:` line of any `--help` still say `ema` — that output is
  unchanged from before the rename, and the deprecation notice never touches
  stdout.

### Changed

- Every CLI invocation in the README, the `docs/` tree, the wizard and the
  in-code help now reads `peakatail`. 555 invocations across 106 files.
- The Click group and the argparse pre-parser both report `peakatail` in usage
  and help output.

### Not changed

- The importable package is still `ema` (`from ema.cli import main`), and every
  `ema/...` source path is unchanged. Renaming the Python package is a separate,
  larger change and is deliberately not bundled here.

## Unreleased — `switch diff` marker pre-selection double-dip (issue #94)

### Breaking changes

- **`ema switch diff --marker-top-n` now defaults to `0` (was `200`).** The
  old default pre-filtered the tested PAS to the union of the top-200 marker
  PAS per cluster, ranked by `scanpy.tl.rank_genes_groups` on the *same*
  `--cluster-key` labels the differential test then contrasts — a label
  double-dip — and additionally shrank the within-gene Fisher denominator
  (the "rest of the gene" background became the same label-selected subset).
  Under a 20-run label-permutation null on a correctly-keyed matrix:

  | Configuration | Null p < 0.05 | Runs with a q < 0.05 hit |
  |---|---|---|
  | `fisher --count-mode reads --marker-top-n 200` | 20.3 % | 20 / 20 |
  | `fisher --count-mode cells --marker-top-n 200` | 13.0 % | 19 / 20 |
  | `nb_pairwise --marker-top-n 200` | 24.7 % | 19 / 20 |
  | `fisher --count-mode cells --marker-top-n 0` | 3.0 % | 0 / 20 |

  `--marker-top-n 0` was the only FDR-controlled configuration measured, so a
  flagless `ema switch diff` is now calibrated. **A run that relied on the old
  default tested a different (smaller, label-selected) PAS set and reported
  anti-conservative q-values; re-run it at the new default.** The same applies
  to the `marker_top_n` YAML key and to `ema.switch_test.runner.run_diff`,
  whose signature default is unchanged (it is a required argument) but whose
  callers now pass `0`.

### Fixed

- **The within-gene Fisher denominator no longer changes when the tested PAS
  are restricted.** `FisherStrategy` computed each gene's "reads / cells at
  the OTHER PAS of this gene" total from the matrix it was handed, so a
  `--marker-top-n N` run silently redefined a gene total as "the
  marker-selected PAS of this gene" — the same PAS scored an `n_reads_gene` of
  1,746 restricted vs 5,289 unrestricted, and only 629 of 6,453 p-values
  agreed between a marker-on and a marker-off run of the same input.
  `run_diff` now hands `fisher` the unrestricted matrix alongside the
  restricted one (new `full_count_matrix` argument, threaded through
  `run_one_pair` / `_dispatch_pair`), so `--marker-top-n N` decides only
  *which* PAS are tested and reported and every reported p-value is identical
  to the unrestricted run's. Genes whose selected PAS number fewer than two
  are now tested too, since the gene itself still has a background.
  The remaining (unfixable-by-code) half of issue #94 is the label double-dip
  in the selection itself, which is why the default stays `0`.
  Exception: under `--isoform-agg within_utr|between_utr` the denominator *is*
  the group's own columns by design, so combining it with `--marker-top-n > 0`
  still narrows the background — that combination now warns.

### Changed

- **Any non-zero `--marker-top-n` now logs a loud warning** naming the
  double-dip and the measured null inflation, from
  `ema/switch_test/runner.py::run_diff` so the library path warns too, and
  from `ema/cli/switch_diff.py` in `--marker-top-n` vocabulary. `run_diff`
  gained `warn_marker_top_n` (default `True`) which the CLI sets to `False`,
  so a CLI user reads the warning once rather than twice. The flag remains
  available as a speed shortcut / ranking screen for large datasets — it is
  not a statistical filter.
- `--marker-top-n` gained real `--help` text (it previously showed only
  `[default: 200]`), and `docs/cli/switch-diff.md` gained a
  "Why `--marker-top-n` defaults to 0" section with the permutation-null table.
- Docs no longer promise a `markers.tsv` from a flagless run: the quickstart
  output listing, `docs/concepts/output-files.md` and
  `docs/concepts/data-flow.md` now say it appears only with
  `--marker-top-n > 0`.

### Added

- `tests/test_marker_top_n_double_dip_i94.py` — a label-permutation null
  regression on a small synthetic matrix (fixed seeds, no I/O): at the default
  settings the null p < 0.05 rate stays near nominal with zero q < 0.05 hits
  across every seed, while `--marker-top-n 200` inflates it several-fold.
  Plus denominator regressions: a hand-built sparse two-gene fixture whose
  non-selected PAS carry counts (asserting the gene totals, p-value,
  `odds_ratio`, `delta_proportion` and `log2fc` of a restricted run match the
  unrestricted run in both `--count-mode` values, with a guard test proving
  the fixture really does expose the defect), and an end-to-end `run_diff`
  comparison of a `--marker-top-n 200` run against a `--marker-top-n 0` run of
  the same input over the 400 PAS they share.
## Unreleased — switch length: degenerate proximal==distal pairs

### Fixed

- **`ema switch length --strategy classic --isoform-agg per_isoform` no
  longer emits pairs whose proximal and distal endpoints are the SAME PAS
  (`ema/quantification/strategies/classic.py`, issue #98).** The
  `transcript_id -> [(pas_id, rank), ...]` map carries **one entry per UTR
  exon a PAS intersects**, so a transcript with a single usable PAS still
  arrived at the `len(...) < 2` guard holding two or more entries; ranking
  those duplicates then picked that one PAS as *both* endpoints. Such a row
  has `proximal_pas_id == distal_pas_id`, `proximal_reads == distal_reads`
  and `pdui` exactly `0.5` by construction — no length information at all —
  and it was broadcast to every cell. Repeated PAS are now collapsed on
  `pas_id` (keeping each site's most proximal rank) *before* the endpoints
  are chosen, and a transcript left with a single **distinct** PAS is
  excluded exactly like any other `<2` PAS transcript. A `WARNING` names the
  gene and how many transcripts were dropped.

  **User-visible output change — `pdui_classic.tsv` now has FEWER ROWS.** On
  the GSE104556 testis runs reported in the issue this removes 18,116 of
  13,709,930 rows on mouse1 (7 genes; 14 `(gene, transcript)` units ×
  1,294 cells) and 15,004 rows on mouse2 (8 genes). Every removed row
  carried `pdui == 0.5`, so no real PDUI value is lost — but any tally of
  rows, of transcripts, or of the `0.5` bin will differ from a pre-fix run,
  and per-isoform PDUI histograms lose a spike at 0.5 that was an artefact,
  not biology. `--isoform-agg per_gene` is **unaffected**: it selects its
  endpoints from the deduplicated PAS coordinate table and never produced
  such a pair.

- **The proximal/distal strand guard no longer blames strand selection for a
  degenerate pair (`ema/switch_test/runner.py::_assert_pdui_strand_convention`).**
  The test was `~(proximal_start < distal_start)` on `+` (and its mirror on
  `-`), so `proximal_start == distal_start` fell into the *inverted* bucket
  and the `RuntimeError` reported strand-blind selection — the wrong cause
  for all 18,116 rows above, which had zero true inversions and zero
  different-strand rows. Equality is now its own class: a pair is
  `DEGENERATE` when the two PAS ids are equal or the two starts are equal.

  **User-visible message change** — anything grepping this `RuntimeError`
  must be updated. It now reads:

  ```
  <source>: strand/proximal-distal convention violated on N of M
  coordinate-bearing rows (G gene(s)): I INVERTED, X with proximal and
  distal on DIFFERENT strands, D DEGENERATE. ...
  ```

  where the previous text was `... rows (G gene(s)); X of those have
  proximal and distal on DIFFERENT strands.` The guard still **fails** the
  run on a degenerate row — it is a real defect, just not a strand one — and
  the message now says so.
## Unreleased — concurrent-branch safety

### Fixed

- **`ema reannotate` refuses to start when another live `ema reannotate`
  already holds the same `--out` (issue #65).** Two branches pointed at one
  output dir write the same `04_pas_gene_assignment/<ds>/pas_gene.tsv`,
  `05_annotated_matrix/<ds>/*` and `07_clustering/<ds>/clusters.h5ad` paths,
  so the surviving artifacts are an arbitrary interleaving of two different
  parameter sets — that is how five sweep-grid branches sharing a
  `branch_name` produced three different `pas_gene.tsv` row counts for
  identical declared params, grouped by write time. `reannotate_run()` now
  takes an exclusive `flock` on `<out>/.ema_reannotate.lock` before any work
  begins; if the claim fails it raises `ReannotateError`, which the CLI
  surfaces as a `click.ClickException` — **the command exits 1 having written
  nothing**, and names the pid/host/start-time holding the directory. The
  claim lives on the open file description, so the kernel drops it when the
  process exits for any reason: a crashed or killed branch leaves no stale
  lock, and `.ema_reannotate.lock` never needs deleting by hand. This is a
  **user-visible behaviour change** — a workflow that (accidentally or
  deliberately) ran two concurrent branches into one `--out` now gets a hard
  error instead of silently corrupted output; give every branch its own
  `--out`. Documented in `docs/cli/reannotate.md`.
- **Every text artifact is written atomically (`ema/outputs.py`).** The new
  `atomic_write()` context manager writes to a temp file in the **same**
  directory and `os.replace()`s it onto the final path, replacing the plain
  `open(path, "w")` / `to_csv(path)` / `write_text(path)` calls behind
  `pas_gene.tsv`, `annotatedpas.bed`, `pasbed.bed`, the annotated-matrix
  `pas_ids`/`barcodes` sidecars, `run_config.json`, `run_manifest.json` and
  the per-stage `*_stats.json`. A partially written file is therefore never
  visible at a final artifact path: a reader sees either the complete
  previous content or the complete new one, two racing writers can only
  produce one of the two complete files (never a byte-level mix), and a
  write that raises leaves the previous file untouched. File **contents** are
  unchanged — only the instant at which they become visible.
## Unreleased — PAS→gene assignment in overlapping loci

### Fixed

- **A PAS in an overlapping locus is no longer handed to the spanning /
  readthrough model (issue #99).** `find_close` called
  `bedtools closest -t first`. Where a readthrough model covers a neighbour's
  3'UTR every PAS there is distance **0** from *both* models, so `-t first`
  resolved the tie by file order — always the spanning model, because it
  starts further upstream. On the 17-sample Laughney cohort every PAS in
  chr17:7,579,636–7,582,386 was labelled `SENP3-EIF4A1` while lying in
  **CD68**'s 3'UTR, and CD68 was left with **0 assigned PAS cohort-wide**
  (same shape for PTPRCAP←CORO1B, FKBP11←AC073610.2, GPX1←RHOA, MDK←DGKZ,
  STARD10←ARAP1, ARPC1A←AC004922.1, FCER1G←NDUFS2, KRTCAP3←NRBP1). 33 % of
  the top-30 replicated gene switches carried the wrong gene name (21 % of
  the top-100, 12 % of the top-500, 7 % overall), the strongest hits being
  enriched precisely because marker genes live in such loci.

  The call is now `-t all`, and `_resolve_overlapping_genes` picks one gene
  per PAS from **what the annotation says**, not from proximity alone.
  Nearest-annotated-3'-terminus is not safe on its own: the nearest terminus
  in such a locus is very often a nested miRNA/snRNA/pseudogene or a lncRNA
  with no `three_prime_utr` record at all (MIR33B and MIR6777 inside SREBF1,
  RNU6-862P inside NCOR1, MIR1288 inside PIGL, MIR6778 inside SHMT1,
  AC016876.3 over CD68). Handing the PAS to one of those grades it TIER_3 —
  `utr_lengths` has no entry — and the default tier filter then **deletes**
  the PAS from `annotatedpas.bed` and from the count matrix, which is worse
  than a wrong label. Candidates are therefore ranked, best first:

  0. the candidate does not turn a kept PAS into a filtered-out one;
  1. the candidate has a `three_prime_utr` record — **hard rule: a gene
     without one never wins a tie against a gene with one**;
  2. the PAS lies inside the candidate's annotated 3'UTR footprint;
  3. the PAS lies inside the candidate's *unextended* gene body;
  4. distance to the candidate's annotated 3' terminus;
  5. the candidate's annotated span, then its gene ID, for determinism.

  `gtf_bed`'s 3'-end extension is now the named constant
  `gtftobed.GENE_EXTENSION_BP` (unchanged at 5000) so the annotated gene body
  and terminus can be recovered from the extended record;
  `find_close(gene_extension=…)` overrides it for a gene BED that `gtf_bed`
  did not build. A minus-strand record clamped at a contig start (where
  `gtf_bed` writes `start=1` instead of subtracting the extension) has no
  recoverable terminus, so it now gets a terminus *range* and the widest
  possible body rather than an invented point thousands of bp from the real
  one.

  Measured on GRCh38.99 chr17 with all 7,801 annotated 3'UTR termini as the
  PAS set: **7,801/7,801 assigned (was 7,716) and 97.0 % assigned to a gene
  that really has a 3'UTR record ending there (was 94.1 %)**. Both CD68 PAS
  in chr17:7,579,636–7,582,386 now go to CD68.

  **This changes gene labels.** PAS coordinates, counts and tiers are
  untouched — only the `gene_id` a contested PAS is assigned to — but
  `annotatedpas.bed`, `pas_gene.tsv` and every downstream gene-level result
  differ in overlapping loci, so switch tables computed before this fix
  should be regenerated. PAS with a single candidate gene are unaffected.
  Regression tests: `tests/test_pas_gene_overlapping_loci_i99.py`, including
  a committed slice of the real GRCh38.99 chr17 annotation
  (`tests/data/GRCh38.99_chr17_overlapping_loci.gtf`).

Issue #99's second defect — the poly(A) clip-rate warning sampling only the
head of a coordinate-sorted BAM (`ema/countmatrix/polya.py`,
`max_reads=200_000`) — is **not** addressed here: branch `peakAtail-prime`
(PR #100) already replaces that estimator with an exact per-pass count
(`--clip-rate-sampling pass`), and duplicating it would collide.

### How this composes with `--pas-gene-rescue` (PR #100 merge)

Both land in `find_close`. They do different jobs and now run in a fixed
order: **the precedence above decides WHICH gene owns a PAS; only then does
`--pas-gene-rescue inside` re-grade the PAS it was given.** The rescue is
applied to the surviving row, after `_resolve_overlapping_genes` and before
the tier filter — deliberately never to the candidate rows. Rank 0
(`_would_drop_pas`) outranks rank 1 (`_no_utr_record`), so re-grading
candidates first would let a nested miRNA's row survive the tier filter while
its protein-coding host's did not, and the miRNA would take the PAS — exactly
the hard rule this section exists to enforce, and the one whose absence
silently deleted 124 chr17 PAS in an earlier draft. `--pas-gene-rescue` keeps
its `off` default, so nothing changes unless it is asked for. Pinned by
`test_the_rescue_cannot_hand_a_pas_to_a_gene_with_no_utr_record` and
`test_the_rescue_does_not_change_which_gene_owns_a_pas` in
`tests/test_pas_gene_overlapping_loci_i99.py`.
## Unreleased — `--dynamic-threshold` no longer aborts the run (issue #101)

**`--dynamic-threshold --lambda-fold-change 2.0` — the parameter reference's
own first entry under "find more PAS" — used to abort the caller with a bare
`IndexError: list index out of range`, raised from inside a spawned
chromosome worker with no message naming a flag, a contig or a remedy.** Both
peak-calling loops computed the current peak's right edge as
`l_end = data_array[-current_threshold]`
(`ema/countmatrix/peackcalling.py`, `ema/countmatrix/peak_pipeline.py`), the
dynamic estimator sets
`current_threshold = max(floor_threshold, int(local_lambda * lambda_fold_change))`,
and nothing bounded it by `len(data_array)`. Found by the 2026-08 parameter
sweep on the PBMC chr19+21 dev slice and reproduced independently on the
GSE104556 mouse1 chr18+19 slice, so the defect is species-independent.

### Fixed

* **The look-back index is bounded by the live read window, by default.** The
  rule lives in one place, `ema.countmatrix.dynamic_threshold.resolve_l_end()`,
  and carries an identity guarantee that is tested exhaustively: for every
  in-range threshold the bounded and unbounded branches return the *same*
  element, so bounding can only change a run that would otherwise have
  aborted. A run that trips the bound now logs a census of how often it fired
  instead of dying. **No output changes for any run that did not crash** —
  `--dynamic-threshold` is off by default, `--lambda-fold-change` is read at
  exactly two places in the tree and both are inside `if dynamic_threshold:`,
  and `tests/test_dynamic_threshold_bounds.py` pins byte-identity of the BED
  and the count matrix on a dynamic run that never trips the bound.
* **`--floor-threshold` is validated.** It was an unchecked INTEGER, and
  `--floor-threshold 0` made `data_array[-0]` return `data_array[0]` — the
  OLDEST end in the window rather than the peak edge. No crash, no warning,
  a wrong peak boundary. Click now takes `IntRange(min=1)` and
  `peak_calling()` re-checks the value, so YAML and library callers (and every
  spawned tile / chromosome worker, which all re-enter through it) get the
  same refusal.

### Documentation

Two help texts that did not match the code, both proven by the same sweep:

* **`--pas-gap` does nothing on a single-BAM run.** It is consumed at exactly
  one place, `merge_pas_beds` on the multi-dataset unified path
  (`ema/main.py`), and `--pas-gap 25` and `200` are byte-identical to the
  baseline on a single BAM. Its help said "minimum gap between PAS within a
  peak"; it now says it is a multi-dataset merge parameter, and points at
  `--min-pas-spacing` for the within-peak behaviour it was mistaken for.
* **`--min-cells` and `--min-pas-per-cell` filter the AnnData, not the call
  set.** Both act in `preprocessing()`, after `pasbed.bed` is written, and
  neither changes a byte of it. Their help, `docs/cli/run.md` and the
  quickstart now say which outputs they touch.

### Note for PR #100 (`peakAtail-prime`) — RESOLVED IN THE MERGE

That branch introduced the same module behind `--dynamic-threshold-clamp`,
defaulted **off** to keep v2 byte-identical. Merging `develop` into
`peakAtail-prime` consolidated the two on THIS implementation, as this note
required: `DYNAMIC_THRESHOLD_CLAMP_DEFAULT` stays `True`, the bound is
unconditional on every production path, and no command line, YAML key or
library keyword can turn it back off.

`--dynamic-threshold-clamp` survives only as an **accepted no-op**, so that
the branch's documented v2-compat command line and every existing config keep
parsing; both of its values now mean "bounded", as does the
`dynamic_threshold_clamp=` keyword the peak loops still accept. The branch's
byte-for-byte v2 promise is unharmed, and exactly because of the identity
guarantee: every v2 run that *produced output* produces the same output here.
The only runs that differ are the ones v2 aborted, and those emitted nothing
to be identical to. `resolve_l_end(data_array, current_threshold, clamp)`
keeps its `clamp` argument for the unit tests that pin that guarantee — they
have to be able to evaluate the unbounded branch to prove the bounded one
agrees with it — and no production caller passes it.

## Unreleased (branch `peakAtail-prime`) — the dynamic-threshold crash guard

**`--dynamic-threshold` at its documented companion setting
(`--lambda-fold-change 2.0`, the parameter reference's first entry under "find
more PAS") aborts the caller with a bare `IndexError`.** Both peak-calling
loops compute the current peak's right edge as
`l_end = data_array[-current_threshold]`
(`ema/countmatrix/peackcalling.py`, `ema/countmatrix/peak_pipeline.py`), the
dynamic estimator sets
`current_threshold = max(floor_threshold, int(local_lambda * lambda_fold_change))`,
and nothing bounds it by `len(data_array)`: when the threshold outgrows the
live read window the index runs off the front of the list. Found by the
2026-08 parameter sweep (`results/paramsweep/VERDICTS.md` §5) on the PBMC
chr19+21 dev slice; the sweep's accuracy verifier reproduced the same abort on
the GSE104556 mouse1 chr18+19 slice
(`results/paramsweep/VERIFY/rerun_identity.tsv`), so the defect is
species-independent. It is a v2 (`9dfdefb`) defect, reachable **only** with
`--dynamic-threshold` (off by default): `--lambda-fold-change` is read at
exactly two places in the tree, both inside `if dynamic_threshold:`, so no
published PeakATail number is affected.

### Added

> **SUPERSEDED BY THE MERGE WITH `develop` (issue #102/#101).** The bound
> below is now **on unconditionally** and `--dynamic-threshold-clamp` is an
> accepted no-op kept for compatibility: neither value can restore the
> `IndexError`. Everything this section says about the *mechanism* still
> holds; only the default, and the ability to choose, are gone. See
> "`--dynamic-threshold` no longer aborts the run (issue #101)" above.

* **`--dynamic-threshold-clamp {off,on}`, default `off`** *(historical — see
  the note above; the flag is now a no-op)*. `off` was v2 to the
  character — the same expression, the same `IndexError`, the same abort —
  because the branch's cardinal rule is that v2 stays reproducible
  byte-for-byte and a fix that changes what a run emits would have to clear
  `manuscript/28` §2 first. On the default path the loops evaluate v2's own
  literal expression with nothing wrapped around it — same expression, same
  cost, same bare `IndexError` (which cost the sweep a whole arm before
  anyone read the source; the flag that stops it is now named in this file,
  in the issue, and in the module docstring).
  `on` bounds the look-back index to the live window (returning the oldest
  end still in it) and logs how often it fired. The single copy of the rule
  is `ema.countmatrix.dynamic_threshold.resolve_l_end()`, and its identity
  guarantee — for every in-range threshold both settings return the same
  element, so `on` can only change a run that would otherwise have aborted —
  is tested exhaustively in `tests/test_prime_dynamic_threshold_clamp.py`,
  which also reproduces the crash itself on a 14-read synthetic BAM through
  the monolithic and pipeline paths.
* The flag is threaded through all four dispatch paths (monolithic,
  per-chromosome parallel, 3-stage pipeline, tiled — including the
  single-tile short-circuit), and `V2_COMPAT_FLAGS` gains
  `--dynamic-threshold-clamp off` (14 flag/value pairs), so the documented
  compat command stays the complete list.

## Unreleased (branch `peakAtail-prime`) — the default that was a no-op

**`--ip-filter-mode`'s default was v2's `annotate`, so the branch's one
behavioural default dropped nothing.** `--ip-filter-default auto` turned the
internal-priming veto ON; `--ip-filter-mode annotate` kept every flagged site.
Two flags, two questions, and nothing checked the pair — so a prime run at its
defaults flagged 5,015 of 32,752 candidates on the PBMC chr19+21 slice, dropped
**0**, and emitted v2's exact call set (`pas 18,865`) while this file
advertised "+7.7 %–12.8 % relative recall at matched atlas precision". Every
number published for that change was measured on a run that also passed
`--ip-filter-mode filter` explicitly.

### BREAKING — the flagless call set changes

**TWO populations lose 17.09 % of their calls, and they do not have the same
escape.** The loss is the internally-primed sites (PBMC 10k v3 full BAM:
402,765 → 333,920, −68,845; `results/prime_bench/fourway_headline.txt` §3).
Nothing about either command line changes; the output does. That is the
intended behaviour of the default flip, and it is the definition of a breaking
default change.

1. **The flagless user who supplies `--genome-fasta` and no IP flag at all.**
   `--ip-filter-default auto` turns the veto ON for them and the new
   `--ip-filter-mode auto` resolves it to `filter`. Escapes:
   `--ip-filter-default off`, `--no-ip-filter`, or `--ip-filter-mode annotate`.
2. **The user who already typed `--ip-filter` and never named a mode.** In v2
   that command line meant `--ip-filter-mode annotate`: the veto ran, flagged,
   and dropped nothing. On this branch the *same* command line resolves the
   mode to `filter` and drops the same ~17 %. **`--ip-filter-default off` does
   NOT rescue this user** — `ema.main._ip_filter_decision()` gives an explicit
   `--ip-filter` precedence over the `ip_filter_default` policy, so the veto
   still runs, in `filter` mode. Their escapes are **`--ip-filter-mode
   annotate`** (exactly v2: the veto runs, the flag is recorded, nothing is
   dropped), `--no-ip-filter` (v2's call set, but the veto no longer runs so
   the IP annotation goes too), or the full 14-pair `V2_COMPAT_FLAGS` command
   line (v2 byte-for-byte, every file).

That precedence is deliberate and is not changing: `--ip-filter` means "run the
veto", `--ip-filter-default` only decides what happens when nobody said, and
`--ip-filter-mode` is the single flag that decides whether the veto DROPS.
Making an explicit `--ip-filter` imply `annotate` would rebuild the exact trap
this branch removed — two flags silently disagreeing about whether anything is
dropped. `ema.main._ip_filter_decision()` is the one copy of the rule, it is
logged, it is written to `run_config.json`, and
`tests/test_prime_ip_default_mode.py` pins it.

**With no readable `--genome-fasta` the veto cannot run**, so the flagless
*call set* is v2's: the derived BEDs (`pas.bed`, `pas_tier1.bed`,
`pas_tier2.bed`, `pasbed.bed`) and the count matrix are byte-identical to v2's.
**`pas_support.tsv` is NOT** — `--pas-features on` and
`--emit-inferred-cleavage on` are branch defaults and append sidecar columns
whatever the FASTA situation is. `--pas-features off
--emit-inferred-cleavage off` (both in `V2_COMPAT_FLAGS`) restores that file
too.

### The literal

`ema/cli/config_schema.py`, field `ip_filter_mode`:
`default="annotate"` → `default="auto"`, with
`choice=("auto", "annotate", "filter")`. `"auto"` is a **sentinel meaning "the
user did not name a mode"**, not a mode:
`ema.cli.config_schema.resolve_ip_filter_mode()` turns it into `"filter"`, and
`ema.experimental.internal_priming.filter_internal_priming` still refuses any
value outside `{annotate, filter}`, so an unresolved sentinel raises rather
than being guessed at. An explicit `--ip-filter-mode annotate` is honoured
verbatim — the D9 behaviour (an internally-primed peak is candidate
alternative-PAS signal, not noise) is one flag away, not gone.

**PRIME_PLAN.md said v2's default here was `filter`. It is `annotate`** — read
off the frozen worktree `tools/pa-polya-run-9dfdefb3`. That one wrong word is
the whole defect; the plan is corrected in place.

### v2 compatibility

`V2_COMPAT_FLAGS` gains **`--ip-filter-mode annotate`** (13 flag/value pairs).
It was previously complete only by accident: the list pins
`--ip-filter-default off`, which stops the veto running at all, so the mode
never mattered. It matters the moment anyone types `--ip-filter`, so it is
stated. The compat runs are re-verified byte-for-byte against the frozen v2
worktree on both dev slices.

### So this class of defect cannot hide again

* **The resolved mode is in the run log.** `_resolve_ip_filter()` logs
  `internal-priming veto: ON (<why>), mode=<mode> [<how the mode was decided>]`
  and still says `NOTHING IS DROPPED` when the mode is `annotate`.
* **The resolved mode is in the run record.** `run_config.json` and
  `run_manifest.json` gain an `internal_priming` block:
  `ip_filter`, `ip_filter_why`, `ip_filter_mode_requested` (the raw value,
  sentinel and all), `ip_filter_mode_resolved` (`null` when the veto does not
  run), `ip_filter_mode_why`, `genome_fasta`. The `args` block cannot show
  this — it carries the unresolved sentinel, which is the surface the defect
  hid on for a whole verification pass.
* **`tests/test_prime_ip_default_mode.py`** drives the real seam
  (`ema.main._apply_pas_filters`) with `args` seeded from the schema exactly as
  a no-flag `ema run` seeds it, plus a genome FASTA and nothing else, and
  asserts a flagged site is GONE from the BED on disk. Restoring
  `default="annotate"` makes it fail with `filtered == 0` and the flagged row
  still present.
* **`tests/test_prime_compat_flags.py::test_every_v2_option_whose_default_this_branch_moved_is_on_the_command_line`**
  diffs RunConfig **defaults** against v2's, not just field names. The old
  completeness check only asked "which fields are NEW?", which is precisely why
  a moved default on an existing v2 field went unnoticed. Both completeness
  checks used to read v2 from one absolute path on one machine behind a
  `skipif`, so **they skipped everywhere else, CI included** — the branch's
  strongest guarantee behind a green checkmark that had never run it. v2's
  defaults are now read from `tests/fixtures/v2_runconfig_defaults.json` (a
  committed dump of v2's `RunConfig`, so the checks run unconditionally) or,
  when `PEAKATAIL_V2_WORKTREE` names a real v2 checkout, live from it — and a
  third test asserts the two agree whenever both are available, so the frozen
  copy cannot go quietly stale.
* **`tests/test_prime_v2_compat_golden.py::test_branch_defaults_still_write_v2s_beds_and_matrix`**
  pins the headline property directly: at the caller, **typing nothing** writes
  v2's BED and count-matrix bytes (the veto that moves the call set lives above
  `peak_calling()`), while `pas_support.tsv` is asserted to DIFFER — so the
  "byte-identical" claim in this file can never again be written unscoped.
* `_v2_settings()` now pins `args.ip_filter_mode`, and
  `test_prime_compat_flags.py` reads pins from `args` as well as
  `variable_config`, closing the structural gap that file already declared.

### Measured, default vs default

See `results/prime_bench/` for the re-validation: the v2 compatibility pin is
still byte-for-byte v2 on both slices, and the branch default now differs from
v2 in exactly the expected way (fewer, more precise calls).

## Unreleased (branch `peakAtail-prime`) — the cleavage-offset column, the QC that stopped lying, and the seam that drops lncRNAs

Four small, separately-measured changes (TASK E). **One of them is a genuine
accuracy default (`--ip-filter` turns itself on), one replaces an estimator
that was 4.2x wrong with an exact count, one ships a column and leaves the
coordinate alone, and one is a defect I found, measured, and then shipped
DEFAULTED OFF because the measurement said it costs precision.**

### 1. `--ip-filter` runs by default when a genome FASTA is available

The internal-priming veto is the largest measured accuracy lift in this caller
— **+7.7 % to +12.8 % relative recall at matched atlas precision, +17.7 % to
+22.9 % at matched long-read precision** (`results/algo_headroom/VERIFY/tables/
v8_ipveto_value.tsv`) — bigger than every detector change tested put together,
because it is the only stage that reads genomic sequence. It was an opt-in
flag, so every benchmark arm in the manuscript ran with it and every user who
did not read the flag list did not.

`--ip-filter-default auto` (the new default) turns it on whenever a readable
`--genome-fasta` is present, and **warns loudly, naming the cost, when there is
none**. `--no-ip-filter` forces it off; `--ip-filter-default off` restores v2.

Measured, default-vs-default on two dev slices: a prime run with **no
`--ip-filter` flag** reproduces the v2 `--ip-filter` run *exactly* —
PBMC chr19+21 `pas 15,925 | tier1 8,524 | tier2 7,401 | tier1>=2mol 2,883` and
P@100 0.7392 / R_det@100 0.2080 / Kinnex t5 P@25 0.7839 on both; GSE104556
mouse 1 chr18+19 `pas 4,141 | tier1 3,006 | tier2 1,135 | tier1>=2mol 1,549` on
both. It changes nothing for anyone who already passed the flag, which is the
point: it changes what a user gets who does not.

> **CORRECTION (adversarial verification pass, 2026-08-22), then FIXED the
> same day — see "the default that was a no-op" at the top of this file.** The
> paragraph above was measured with `--ip-filter-mode filter` supplied on BOTH
> sides. At the time it was **not** what the branch default alone did:
> `--ip-filter-mode` defaulted to `annotate` and `--ip-filter-default auto`
> decides only whether the veto **runs**, not whether it **drops**.
> Re-measured on the PBMC chr19+21 slice with those branch defaults, a genome
> FASTA and no other flag: `peak_filters_stats.json` reported mode `annotate`,
> **flagged 5,015, filtered 0**, and the run wrote
> `pas 18,865 | tier1 10,886 | tier2 7,979 | tier1>=2mol 3,643` — v2's call
> set to the row, not 15,925. `--ip-filter-mode`'s default is now the sentinel
> `auto`, which resolves to `filter`, so the paragraph above IS the branch
> default; v2's `annotate` is one explicit flag away and is pinned in
> `V2_COMPAT_FLAGS`.

### 2. The poly(A) clip-rate QC: an exact count instead of a 4.2x estimate

`check_clip_rate` read the **first 200,000 CB reads** of the file. On a
coordinate-sorted BAM that is the head of the first contig, and on PBMC 10k v3
it returns **2.2565 %** where the whole-file rate on the same denominator is
**0.5364 %** — a 4.2x over-estimate, reproduced here from the tool's own code
path. Over-estimating is the dangerous direction: this estimator exists to
shout when the poly(A) evidence channel has been destroyed, and 4.2x high would
MASK exactly that.

`--clip-rate-sampling {head,strided,pass}`, **default `pass`**:

* **`pass`** — do not estimate. **Count**, during the peak-calling pass the run
  performs anyway, every read the caller accepts and every qualifying clip
  among them, and report the EXACT rate per (contig, strand). Zero extra I/O:
  `clip_site` is already called on each of those reads.
* `strided` — the "sample across the BAM" design: coordinate-uniform windows,
  every read whose start falls in one. **Shipped, and measured to be
  unreliable**: on the full PBMC BAM against the 0.5364 % truth it returns
  1.5258 % (200k reads) and 1.8068 % (1M) — more reads did not help, because
  clips are rare AND concentrated at the 3' ends of expressed genes, so the
  estimate is dominated by which windows happen to hit one. Kept because it is
  the obvious design and its failure is worth being able to reproduce.
* `head` — v2, **byte-for-byte including its log line**.

A first `strided` implementation gave each stratum an equal READ quota; it was
biased, and the bias moved with the stratum count (0.9111 % at 25 strata,
0.1853 % at 2,500, against a 0.6169 % slice truth). A second took 100 % of MT
and of every scaffold — 36.2 M reads for a 200,000-read budget. Both are
recorded in `results/prime/taskE/clip_rate_estimators.tsv`; the shipped window
sampler is unbiased in construction and still not good enough, which is why the
default is the exact count.

**The docstring constant is corrected.** `ema/countmatrix/polya.py` stated that
**1.152 %** of CB reads carry a poly(A) clip. The genome-wide truth is
**0.5730 %** on the caller's accepted-read denominator (3,195,067 / 557,564,408)
and **0.5364 %** on `check_clip_rate`'s, 0.3669 % (chr21) to 0.8179 % (chr19)
per chromosome.

### 3. `inferred_cleavage`: report the offset, do not move the coordinate

Long-read termini sit a few bp UPSTREAM of the clip-seeded caller's reported
cleavage base. TASK E asked whether that offset can be estimated per library
from the clip-anchored subset. **It cannot, and that is the result.**

* **The caller's own clip evidence says the offset is zero.** With
  `--pas-features on` every tier-1 PAS now carries `clip_offset_mean`: the
  read-weighted mean of (clip position − reported base) over its own cluster,
  in transcript orientation. Aggregated per library that is **−0.334 bp** over
  267,520 clip reads in 16,338 clusters (PBMC chr19+21; **76.10 %** of clip
  reads sit *exactly* on the reported base) and **+0.214 bp** over 124,308
  reads in 8,722 clusters (mouse 1 chr18+19). Both round to **zero**, so
  `--cleavage-offset auto` is a measured no-op on both libraries.
* **The external truths say −1 bp, and only at base-pair resolution.** On the
  PBMC default arm the exact-match optimum is **−1 bp against the atlas**
  (P@1 0.3926 → 0.4346) and **−2 bp against Kinnex long reads** (t5 P@1
  0.2844 → 0.3698); on mouse 1, −1 bp against the atlas (0.3919 → 0.4274).
  **Nothing at W >= 25 responds**: P@100 moves by 0.0028 over the whole
  −8…+8 sweep, and the two truths disagree about the optimum at every
  intermediate window (atlas P@10 wants +5, Kinnex t5 P@10 wants −4). The
  optimum is otherwise stable — −1 on both strands, both chromosomes and every
  support bin except the 1-molecule tier. Full sweeps in
  `results/prime/taskE/offset_sweep_*.tsv`.

So: **`--cleavage-offset {none,auto,<int>}`, default `none`.** The estimate is
computed and reported regardless (`01_peak_calling/cleavage_offset_stats.json`,
`run_manifest.json`), and `--emit-inferred-cleavage on` (default) appends
`inferred_cleavage` to `pas_support.tsv` — the coordinate the offset implies —
**without moving a byte of `pasbed.bed`**.

A signed offset now works: v2 accepted only positive values and silently
ignored negative ones. **The sign decides which tier moves**, and that is
measured, not stylistic: a POSITIVE offset is the coverage correction
(issue #72, ~+95 bp) and under `clip_seeded` must not touch the clip-anchored
tier — v2's `skip_supported` rule, and applying +95 there costs 46 points of
P@10 — while a NEGATIVE offset is the base-pair correction and belongs to that
tier, tier 2's own P@1 being 0.0016 and flat.

**`--auto-cleavage-offset` is now refused under `--peak-strategy clip_seeded`**
rather than warned about. It searches a 60–120 bp band for a genomic A-fraction
crest, so it *cannot* return anything but a large positive number, and the +95
it estimates is destructive here. The `"A sane data-driven constant is ~90-100
(try 95)"` advice is removed from `--cleavage-offset`'s help and from
`docs/cli/run.md`.

**Applying an offset is not free, and the measurement says so.** The shift lands
before the internal-priming veto and gene assignment, so `--cleavage-offset -1`
also changes the call set: PBMC slice 15,925 → 15,823 PAS, ΔP@100 +0.0038,
ΔR_det −0.0009, ΔF1 −0.0007, ΔKinnex t5 P@25 +0.0096, with atlas P@1
0.3926 → 0.4379. A base-pair-resolution gain, a call-set change, and no
movement at the windows the manuscript reports — one more reason `none` is the
default.

**Also measured, and reported rather than fixed:** `cluster_clip_sites` breaks
read-count ties toward the LOWEST COORDINATE on both strands, which is not
transcript-oriented. **27.6 %** of multi-position clusters on chr21 have a tied
top read count, and flipping the tie-break flips the sign of the within-cluster
asymmetry on both strands (+1.504 → −2.400 on `+`, −1.159 → +1.510 on `−`).
Changing it would move every tier-1 coordinate, so it is documented here and
left alone.

### 4. The gene-assignment seam drops PAS in lncRNAs, and no rescue is free

**31–33 % of tier-1 clip clusters never reach `pasbed.bed`.** Decomposed on the
PBMC chr19+21 slice (16,338 tier-1 clusters): 23.3 % dropped by the
internal-priming veto, **20.4 % by the gene-assignment tier gate**, 4.1 % for
having no counts left after the `--min-read` cell filter, 52.2 % kept.

The gene gate is the interesting one, and the mechanism is an **annotation gap,
not a distance judgement**. `assign_tier` can only award TIER_1/TIER_2 when the
assigned gene has an annotated 3'UTR **length**, so a PAS at distance 0 — INSIDE
its gene body — falls through to TIER_3 and is dropped whenever that gene has no
UTR record. **1,762 of the 3,338 tier-1 clusters the gate drops (52.8 %) are
inside a gene body**, hosted by 514 genes of which **469 are lncRNA**, 20 miRNA,
7 snRNA and 18 protein-coding. PeakATail cannot report a PAS in a non-coding
gene, however much clip evidence backs it. On mouse 1 chr18+19 the same gate
drops 29.8 % of tier-1 clusters, 25.7 % of them inside a gene body.

`--pas-gene-rescue {off,inside}` (+ `--pas-gene-rescue-min-mol N`) grades those
TIER_2. **It ships `off`, because it is measured to cost precision:**

| slice | arm | n | P@100 | R_det@100 | F1 | Kinnex t5 P@25 |
|---|---:|---:|---:|---:|---:|---:|
| PBMC chr19+21 | default (v2) | 2,883 | **0.7392** | 0.2080 | 0.3246 | **0.7839** |
| | `inside` | 3,590 | 0.6646 | 0.2132 | 0.3228 | 0.7253 |
| | `inside`, >=10 molecules | 3,080 | 0.7166 | 0.2102 | 0.3250 | 0.7662 |
| mouse 1 chr18+19 | default (v2) | 1,549 | **0.7650** | 0.1979 | 0.3144 | — |
| | `inside` | 1,882 | 0.7428 | 0.1994 | 0.3143 | — |
| | `inside`, >=10 molecules | 1,805 | 0.7607 | 0.1984 | 0.3147 | — |

The dropped set is genuinely worse (P@100 0.3607 and Kinnex t5 P@25 0.4866 for
the inside-gene drops, against 0.7392 / 0.7839 for what is kept), and **no
support floor makes the rescue free**: even the >=10-molecule casualties reach
only P@100 0.4138 / Kinnex t5 0.5655. `manuscript/24` §3.1 fails on (i) and
(ii) for every setting on both slices. The recall column also cannot reward it
— these PAS are, by construction, in genes the detected-gene denominator does
not contain — so the seam is documented, the flag exists, and the default does
not move.

### Flags

| flag | values | default | v2 |
|---|---|---|---|
| `--ip-filter-default` | `off`, `auto` | **`auto`** | `off` |
| `--no-ip-filter` | flag | off | (n/a; `--ip-filter` absent) |
| `--clip-rate-sampling` | `head`, `strided`, `pass` | **`pass`** | `head` |
| `--cleavage-offset` | `none`, `auto`, signed int | `none` | `0` / `none` |
| `--emit-inferred-cleavage` | `off`, `on` | **`on`** | `off` |
| `--pas-gene-rescue` | `off`, `inside` | `off` | `off` |
| `--pas-gene-rescue-min-mol` | int | `0` | `0` |

`--pas-features on` gains a **25th** column, `clip_offset_mean` (TASK C shipped
24); `--emit-inferred-cleavage on` appends `inferred_cleavage` after it. Both
groups stay header-addressed, and `--pas-features off` is still v2 exactly.

### Compatibility

All 50 data files of a PBMC chr19+21 slice run with the v2 pins
(`--read-geometry fixed --read-exclude-flags 0 --pas-features off --pas-score
none --pas-score-model prime1 --pas-score-min -1 --cleavage-offset none
--emit-inferred-cleavage off --clip-rate-sampling head --ip-filter-default off
--ip-filter-mode annotate --pas-gene-rescue off --pas-gene-rescue-min-mol 0
--dynamic-threshold-clamp off` — the canonical copy of this list is
`ema.cli.config_schema.V2_COMPAT_FLAGS`, and `tests/test_prime_compat_flags.py`
checks this paragraph against it) are **byte-identical to the v2 reference run**,
including its clip-rate log line; only `run_config.json` / `run_manifest.json`
differ, by the keys that record the new options.
`tests/test_prime_v2_compat_golden.py::_v2_settings()` pins all five new knobs.

## Unreleased (branch `peakAtail-prime`) — the calibrated per-site score

**One new capability, three flags, and it is OFF by default because the
pre-registered criterion said so.**

The measurement programme found exactly one change that *lifts* this caller's
precision/recall curve rather than sliding along it: a calibrated per-site score
used as a **re-ranker inside the existing tier-1 and internal-priming gates**,
replacing only the `>= 2 clip molecules` threshold. This release implements it,
and reports the two things that matter about it.

**It transfers.** The model is fitted on GSE104556 testis mouse 1 alone and
applied **unchanged** to PBMC 10k v3 — a different species, chemistry and
aligner — and to testis mouse 2. At matched call count it is up-and-right
against the current default on **both** held-out datasets:

| dataset | arm | n | P@10 | P@100 | R_det@100 |
|---|---:|---:|---:|---:|---:|
| PBMC 10k v3 | current default | 46,524 | 0.5209 | 0.7062 | 0.1754 |
| | score, matched n | 46,524 | **0.6020** | **0.7538** | **0.1809** |
| testis mouse 2 | current default | 26,526 | 0.5945 | 0.7572 | 0.2080 |
| | score, matched n | 26,526 | **0.6570** | **0.8043** | **0.2195** |

At matched precision that is **+11.7 % (PBMC) / +16.4 % (mouse 2) / +20.9 %
(mouse 1) relative recall**. Ranking the same candidates by molecule count
instead — the incumbent ordering — gives −0.0007 / −0.0002 on PBMC, so this is
re-ranking and not a rename.

**The default still does not move.** The adoption criterion is *default vs
default* at the threshold the model actually ships with, and it needs
`dP@100 >= -0.005`, `dR_det@100 >= +0.010` and `dF1 > 0` on all three datasets.
The shipped threshold — the calibrated decision boundary `p >= 0.50`, fixed on
mouse 1 without reading any metric — lands at a precision-first operating point
and **fails criterion (ii) on all three** (PBMC: dP **+0.0673**, dR **−0.0048**).
Two other threshold rules, also fixed on mouse 1, fail differently: matching the
incumbent's call count still gives dR +0.0058 on PBMC, and matching the
incumbent's precision on mouse 1 costs 0.047 of precision on PBMC. **No fixed
probability clears the criterion everywhere**, because 33 % of mouse candidates
are atlas-positive against 13.9 % of PBMC's: **the probability transfers as a
ranking, not as an absolute scale.** So the flag ships OFF, and the probability
ships as a column regardless.

### Added

- **`--pas-score {none,calibrated,select}`** (`pas_score`), **default `none`**
  (= previous behaviour).
  - `calibrated` appends one column, `pas_score`, to `pas_support.tsv`. It
    **adds, drops and moves no PAS**.
  - `select` additionally uses the score **in place of** the molecule-count
    threshold. **Tier-1 membership and the internal-priming veto stay hard
    gates in front of it**: the score runs at the same seam, immediately after
    the veto, on the BEDs the veto has already rewritten, so it can only ever
    *remove* a tier-1 candidate — never promote a coverage-only one, never
    rescue a vetoed one. A test pins that by giving an internally-primed
    candidate a `0.0` threshold and asserting it stays dropped.
- **`--pas-score-model NAME|PATH`** (`pas_score_model`), default `prime1` — a
  model shipped with the package, or a JSON produced offline.
- **`--pas-score-min FLOAT`** (`pas_score_min`), default `-1` = "the threshold
  the model was shipped with".
- `run_config.json` records all three under `variables`.
- `ema/countmatrix/models/pas_score_model_prime1.json` — 161 trees, 9,821 nodes,
  0.46 MB of constants. Fitted offline by `scripts/prime/taskD_fit_model.py`;
  the fit is deterministic (re-running reproduces the file byte-for-byte).

### Engineering notes for the reviewer

- **scikit-learn never enters the run-time path.** `ema/countmatrix/pas_score.py`
  evaluates node arrays with numpy; the offline exporter refuses to write a
  model whose numpy evaluation differs from scikit-learn's by more than 1e-9 on
  any training row (the shipped one agrees to **2.2e-16**), and a subprocess
  test asserts that loading a model and scoring with it imports no `sklearn`.
- **Every one of the 21 features is a `pas_support.tsv` column, verbatim**
  (`tier`, not a derived indicator), so the offline fitter and the tool are one
  computation and any `pas_score` can be recomputed from the row beside it. On a
  real run the tool's column and an independent offline recomputation agree to
  **5.0e-7** — half of the last printed digit.
- **One pass, not two.** The score rides the pass `--ip-filter` already makes and
  reuses the feature collector's rows through a lazy `ScoredFeatures` view, so
  the sidecar is still appended to once and a genome-wide run never holds two
  copies of 650 k rows.
- **`seq_ok == 0` means `NA`, not zero.** A candidate whose sequence window could
  not be read is exempt from `select` rather than silently dropped.
- **Cost.** PBMC chr19+21 slice at 8 threads with `--pas-score calibrated`:
  6 m 26.5 s / 1.169 GB against the reference 6 m 35.4 s / 1.164 GB — 1.00x wall,
  1.004x peak RSS. Genome-wide the score step is **24.4 s for PBMC's 652,665
  candidates** (18.9 s of that is the 161-tree traversal), about 1 % of that
  run's wall time, with no extra pass over the BAM, the FASTA or the BED.

### Measured on all three full BAMs

Full transcript, every threshold rule, the calibration curves and the
reproduction commands: `results/prime/TASK_D_pas_score.md`.

**A score inherits the truth that trained it — record this one.** Trained on a
curated atlas, the score buys atlas agreement and *loses* long-read agreement:
on PBMC at the shipped cut, Kinnex x3p t5 P@25 falls **0.7647 -> 0.7294** (t20
0.5584 -> 0.5061) while the internal-priming decoy rate falls **0.1300 ->
0.0650**. An otherwise identical model trained on the Kinnex long reads instead
*raises* t5 P@25 to 0.7904. On separating long-read termini from decoys the
shipped score reaches AUC 0.7185 against **0.7790 for the tool's own inverted
`ip_tool_afrac`** — i.e. it does not beat one covariate the caller already
computes on that axis, and no atlas-trained variant does. Calibration says the
same thing twice: expected calibration error 0.0210 on mouse 2 against the
training label, 0.0399 on PBMC against the same label, and **0.1300 on PBMC
against long reads**.

**The development slice reverses the sign of that read-out.** On chr19+21 the
same comparison gives Kinnex t5 **+0.0129**; on the other 22 contigs it is
**−0.0380**. Any atlas-independent claim from this branch has to be genome-wide.

## Unreleased (branch `peakAtail-prime`) — per-site scoring features

**One new flag. It appends columns to a sidecar and changes nothing else.**

The one measured-positive algorithmic change left for this caller is a
calibrated per-site score used as a **re-ranker inside the existing tier-1 and
internal-priming gates**, replacing only the `>= 2 clip molecules` threshold.
Fitted on long-read termini and evaluated on a curated atlas it had never
seen, it is worth **+6.6 % relative recall together with +7.0 precision
points at matched call count and matched expression**, or **+14.3 % relative
recall at matched precision**. That score has to be fitted offline, and to be
fitted at all it needs per-site covariates written at call time. This release
writes them; it does not score anything.

### Added

- **`--pas-features {off,on}`** (`pas_features`), **default `on`** —
  append 24 columns to `pas_support.tsv`. `off` is the previous behaviour,
  byte-for-byte.

  It **adds, drops and moves no PAS**. `pasbed.bed` stays BED6, the count
  matrix is untouched, and every pre-existing sidecar column keeps its name,
  its position and its value — columns are only ever appended. It is a flag
  and not unconditional because the sidecar's *bytes* change, and on this
  branch previous output must stay reachable exactly.

  | group | columns |
  |---|---|
  | written by the caller | `clip_positions` `clip_span` |
  | genomic sequence | `seq_ok` `ip_tool_flag` `ip_tool_afrac` `ip_tool_arun` `a_count_d18` `a_frac_d18` `a_run_d18` `a_frac_d30` `a_run_d30` `kin_ip_flag` `hex_strong` `hex_any12` `hex_n_types` `hex_best_off` `hex_strong_off` |
  | local candidate context | `d_prev_cand` `d_next_cand` `n_cand_100` `n_cand_500` `mol_500_sum` `is_local_mol_max` `mol_frac_local` |

  Every column is documented in `docs/cli/run.md`.

- `run_config.json` records `pas_features` under `variables`, so a run tree
  states which sidecar schema produced it.

- Multi-BAM runs, whose PAS ids are re-keyed at merge time and which therefore
  have no run-root `pas_support.tsv` to extend, get a standalone
  `<run>/pas_features.tsv` in the merged id space instead.

### Measured on two development slices

Full transcript and reproduction commands: `results/prime/TASK_C_pas_features.md`.

**The call set does not move.** `pas.bed` has the same md5 in every arm — v2, prime with the
features on, and prime with them off — on PBMC chr19+21 (`bb65b0e0…`, n 15,925 | tier1 8,524 |
tier2 7,401 | tier1≥2mol 2,883) and on GSE104556 mouse1 chr18+19 (`57b86628…`, n 4,141), with
and without `--ip-filter`.

**Previous output stays reachable.** `--pas-features off` against a run of the pre-branch code:
**51 / 51 data files byte-identical** on both slices, `pas_support.tsv` included. The only
differences in either tree are the keys `read_geometry`, `read_exclude_flags` and `pas_features`
now recorded under `variables` in `run_config.json` / `run_manifest.json`.

**Cost.** Whole-run peak RSS **+0.28 %** (PBMC) and **+0.55 %** (mouse1). The whole-run wall
delta is not usable — the `off` arm, which does strictly less work, was the slowest PBMC run in
the set, so the box's run-to-run spread (30 s on a 6-minute run) exceeds the effect. Measured
directly instead, one arm per process over the real caller BEDs and the real genome:

| slice | candidates | internal-priming pass alone | + all 22 seam columns, appended |
|---|---:|---:|---:|
| PBMC chr19+21 | 32,752 | 0.286 s / 17.0 MB | **1.319 s / 39.6 MB** |
| mouse1 chr18+19 | 11,444 | 0.109 s / 16.1 MB | **0.485 s / 24.2 MB** |

**+1.03 s and +22.6 MB** (PBMC), **+0.38 s and +8.1 MB** (mouse1) — 32 µs and ~0.49 kB per
candidate. Projected onto a genome-wide PBMC run's 651,957 PAS: ≈ 21 s and ≈ 320 MB held while
the seam runs.

**The values are the offline analysis's, verified twice.** Joined on
`(contig, cleavage, strand)`: every sequence column agrees **exactly** with the stored offline
feature table on all **15,925** coordinates the two runs share, and with a fresh run of that
analysis's own script over `bedtools getfasta -s` windows on **3,000** sampled sites. The
context columns agree on 95.4–99.9 %, and every disagreement is the expected consequence of the
two runs having different candidate sets (32,752 here against 19,242 in the genome-wide table on
the same two contigs) — anyone fitting a model must recompute the context block from the tool's
own output rather than reusing the offline one.

**Internal consistency.** `ip_tool_flag == 1` on exactly the 5,015 sites the internal-priming
veto dropped on the PBMC slice — the covariate and the gate agree on every site, which is the
point of taking it from the string the veto tested.

### Cost: no extra pass over anything

`clip_positions` / `clip_span` come from numbers the caller already has —
tier 1 from the single-linkage cluster's own member list, tier 2 from the same
`±--polya-window` slice its four clip counts already come from.

The 22 sequence and context columns are computed at the internal-priming
stage, **inside the pass that already walks the PAS BED with the genome
open**: one extra 71 nt slice per PAS from a contig record already in hand.
The genome is opened exactly as many times with the features on as with them
off, and that is asserted by a test rather than argued. When `--genome-fasta`
is supplied but `--ip-filter` is not, the feature scan *is* that single pass
(it drops nothing and writes no BED). With no FASTA at all the sequence
columns are `NA` — announced by a warning, never silently — and the context
columns are still real.

### Strand

Every window is transcript-relative around the cleavage base (BED `end - 1` on
`+`, BED `start` on `-`) and is reverse-complemented on `-`, the same
convention the internal-priming window was corrected to. The first base's
offset is derived from the bases *actually returned*, so clamping at a contig
start and truncation at a contig end are accounted for rather than assumed
away: a window that does not fully cover `r in [-40, +30]` reports `seq_ok 0`
and `NA`, never a silently short denominator. `seq_ok` is `NA` rather than `0`
when no genome was supplied, so "no genome" stays distinguishable from
"contig edge".

### `ip_tool_afrac` is the veto's own window

On separating genuine long-read 3' termini from internal-priming decoys, the
tool's own A-fraction over its internal-priming window — **inverted** — is a
stronger discriminator (AUC 0.7790) than the entire 49-feature model that
motivated this column set (0.7655). It is therefore taken from the very string
the veto tests, so it tracks `--ip-window-left` / `--ip-window-right` exactly
and cannot drift from the rule it summarises.

`ip_tool_flag` is emitted **in addition to** the internal-priming veto and
never as a replacement for it: every configuration in which a score was
allowed to override that veto looked excellent against a curated atlas and no
better than the plain rule against long reads. The veto stays a hard gate.

### Not emitted, on purpose

* A second BAM pass for per-cell clip counts, top-cell share, end counts at
  ±5/25/100, pileup sharpness or end-position entropy: **0.000–0.002** held-out
  AUC.
* Any molecule-end pileup statistic: after matching on local read depth, a
  molecule-end pileup at a true missed site is as likely as at a random
  position of the same depth (0.95–1.01×), against 11.4× for the clip channel.
* Anything needing a sequence window wider than `r in [-40, +30]`.

---

## Unreleased (branch `peakAtail-prime`) — read acceptance geometry

**Two new flags. Both default to the previous behaviour, and one of them
does so because the measurement said to.**

`ema/countmatrix/read.py` normalised every accepted read to exactly
`--seq-len` bp of *reference* span: a read whose span was longer was
**discarded**, a shorter one had its end rewritten to `start + seq_len`. A
spliced alignment's reference span includes its introns, so the discard fell
almost entirely on spliced reads — **before the poly(A) clip detector and
before the count matrix**.

Census on the PBMC 10k v3 chr19+21 dev slice (50,898,456 records,
`--seq-len 91`) and a GSE104556 mouse1 chr18+19 slice (11,580,142 records,
`--seq-len 98`):

| | PBMC slice | mouse1 slice |
|---|---:|---:|
| valid-CB reads reaching the rule | 48,647,964 | 11,343,946 |
| **discarded**, reference span > `--seq-len` | 11,768,752 (24.19 %) | 2,263,738 (19.96 %) |
| …spliced | 98.08 % | 98.72 % |
| **end rewritten**, span < `--seq-len` | 4,752,308 (9.77 %) | 2,462,026 (21.70 %) |
| mean fabricated 3'-end shift, downstream | 13.53 bp | 25.17 bp |
| qualifying poly(A) clip reads lost | 12,427 (+4.65 %) | 1,913 (+1.54 %) |
| secondary alignments among valid-CB reads | 10.42 % | 12.09 % |

Genome-wide the discard is 13.74 % of valid-CB reads, 96 % of them spliced —
**88.8 M reads on the full PBMC BAM**. The slices are spliced-richer than
average, so every slice figure here over-states the genome-wide effect.

`keep` and `true` recover **11,526,534** of those reads on the PBMC slice
(23.69 % of valid-CB reads, 97.94 % of what was discarded) and **2,227,952**
on the mouse slice (19.64 % / 98.42 %) — **100 % of them spliced**, by
construction: only an intron can shrink a read's footprint below its span, so
an unspliced read the old rule rejected is still rejected. Scaling onto the
genome-wide census, ≈ **87 M reads** on the full PBMC BAM.

### Added

- **`--read-geometry {fixed,keep,true}`** (`read_geometry`), **default
  `fixed`** — the previous behaviour.
  * `fixed` — discard span > `--seq-len`, pad shorter reads to
    `start + seq_len`. **Reproduces pre-branch output byte-for-byte**:
    validated on the real chr19+21 slice (all 50 data files identical to a
    reference run made from the frozen pre-branch tree) and on a committed
    fixture by `tests/test_prime_v2_compat_golden.py`.
  * `keep` — stop discarding; keep the fixed-length interval. The ablation arm.
  * `true` — additionally use the read's real aligned reference footprint:
    soft clips excluded at both ends (they are not aligned, and the terminal
    poly(A) clip is the clip detector's evidence, not coverage), deletions
    inside the span, **introns removed**. Introns must be removed: over
    7.82 M valid-CB slice reads they carry 10.02 Gb, ~14× the real read mass
    over the same 105 Mb, so admitting them turns the coverage state machine
    into a gene-body detector.

  Acceptance moves to the read's **de-introned reference footprint**, not to
  its query length. A query-length rule is *not* a strict relaxation: a short
  alignment with a long terminal soft clip has span ≤ `--seq-len` (kept
  before) but query length > `--seq-len` — the shape of a poly(A) clip read.
  Every read `fixed` accepts, `keep` and `true` accept.

- **`--read-exclude-flags N`** (`read_exclude_flags`), **default `0`** = the
  previous behaviour = no filtering. The coverage/count channel applies no
  `samtools -F`, so a read aligned to *N* places contributes *N* reads of
  coverage and *N* matrix counts. `256` drops secondary alignments.

- `run_config.json` now records both under `variables`, so a run tree states
  which geometry produced it.

### Measured — and why the defaults did not move

Default precision arm (tier-1 ∩ ≥2 molecules), against `fixed`:

| dataset | arm | n | ΔP@100 | ΔR_det@100 | ΔF1 | precision of the ADDED calls |
|---|---:|---:|---:|---:|---:|---:|
| PBMC slice | `keep` | 3,722 | −0.0026 | +0.0030 | +0.0029 | 0.519 |
| PBMC slice | `true` | 3,724 | −0.0026 | +0.0031 | +0.0030 | 0.519 |
| mouse1 slice | `keep` | 1,964 | −0.0283 | +0.0059 | +0.0039 | 0.353 |
| mouse1 slice | `true` | 1,977 | **−0.0292** | +0.0063 | +0.0042 | 0.368 |

(base P@100: PBMC 0.6401, mouse1 0.7156; null P@100 0.031 / 0.014.)

`manuscript/24` §3.1 requires ΔP@100 ≥ −0.005, ΔR_det ≥ +0.010 and ΔF1 > 0 on
every dataset. **`true` fails (ii) on PBMC and fails (i) on mouse1 by 5.9×**,
so per the pre-registration it stays behind a non-default flag. The mouse
precision loss survives `--read-exclude-flags 256` (−0.0117) and the ablation
puts essentially all of it on *stopping the discard*, not on the 3'-end fix:
`keep`→`true` is worth ΔP −0.0001 / ΔR +0.0001 on PBMC and −0.0010 / +0.0004
on mouse1. The recovered spliced reads carry real evidence — their calls agree
with the atlas 12–26× better than the genic null — but less precisely than the
evidence already in hand, so admitting them slides along the curve instead of
lifting it.

**What `true` does buy, and what §3.1 does not measure:** raw count-matrix
mass **+31.2 %** (PBMC slice) and **+23.5 %** (mouse1 slice), ~**+15.9 %**
genome-wide from the 13.74 % genome-wide discard. That is the reason the flag
exists, and moving the default on the strength of it needs a quantification
criterion the pre-registration does not currently have. (The *filtered*
matrix and cell counts move much more — +89 % mass, 7,133 → 13,516 cells —
but that is a slice artefact: `--min-read 1500` bites hard when a cell is only
seen on 7 % of the genome.)

`--read-exclude-flags 256` is a clean precision-for-recall trade, not a lift:
on the default arm it is ΔP@100 +0.0125 / ΔR −0.0008 (PBMC) and +0.0218 /
−0.0031 (mouse1). Filtering can only slide along the curve, so it stays off.

The geometry adds calls rather than moving them: **97.6 %** (PBMC) and
**91.1 %** (mouse1) of `true`'s default-arm calls sit at *exactly* the same
coordinate as a `fixed` call, and the median signed distance from a call to
the nearest Kinnex long-read 3' end is unchanged at −2 bp.

Compute on the PBMC slice at 8 threads: wall 6:31 → 8:38 (1.32×), peak RSS
1.16 GB → 1.74 GB (1.53×, over the 1.5× guard rail of `manuscript/24`
§3.2.4 — declared, not hidden). On the mouse slice, 1.17× wall and 1.01× RSS.

## Unreleased — caller memory and CPU

Peak RSS and wall time only: **every output file is byte-identical** at the
default settings. Verified end-to-end by re-running two full benchmark arms
on the same inputs and comparing bytes against the runs they are supposed to
reproduce: the chr19+21 PBMC slice (21 files), GSE104556 mouse1 testis with
`--ip-filter` (18 files) and the full PBMC 10k v3 CellRanger BAM (18 files)
— `pasbed.bed`, `pas_support.tsv`, the raw and filtered matrices, the cell
list, `annotatedpas.bed`, `annotated_matrix.mtx`, `pas_gene.tsv` — plus
`preprocessed.h5ad` and `clusters.h5ad` (`X`, `X_lsi`, Leiden labels) on all
three.

| Run | Wall before | Wall after | Peak RSS before | Peak RSS after |
|---|---|---|---|---|
| PBMC 10k v3, chr19+21 slice (`--threads 16`) | 14 min 41 s | 6 min 11 s | 5.56 GB | 1.16 GB |
| GSE104556 mouse1 testis (`--threads 12 --ip-filter`) | 1 h 01 min | 9 min 03 s | 23.12 GB | 3.68 GB |
| PBMC 10k v3, full BAM (`--threads 16`) | 3 h 45 min 53 s | 27 min 43 s | 293.74 GB | 12.45 GB |

Not all of it is parallelism: the slice re-run with `--peak-workers 1` (the
legacy single-process caller) takes 11 min 27 s / 1.71 GB, with peak calling
at 613 s instead of 733 s and the cell-barcode filter at 25 s instead of 64 s.

### Fixed

- **TF-IDF no longer densifies the count matrix
  (`ema/clustering/strategies/leiden_tfidf.py`).** `_tfidf_signac_method1`
  called `X.toarray()` and then held four dense float64 copies at once
  (`astype`, `tf`, `tf * idf`, `log1p`), i.e. `4 × n_cells × n_PAS × 8`
  bytes. That single function was **>98 % of peak RSS on every dataset
  measured** — the formula predicts 104.3 GB / 140.7 GB / 236.8 GB /
  291.2 GB / 21.5 GB against measured peaks of 105.6 / 140.2 / 239.5 /
  293.7 / 23.1 GB — and most of the 31-minute clustering stage of the PBMC
  run. The same scalar sequence applied to the stored entries needs
  `~3 × nnz × 8` bytes and produces the **same bits**: structural zeros map
  to `log1p(0) == 0` and were dropped by the trailing `csr_matrix()` anyway,
  and the row sums / per-PAS cell counts are integer sums, exact in float64
  in any order. Checked against the previous implementation on the real
  `preprocessed.h5ad` of two runs (CSR structure, data, `X_lsi`, Leiden
  labels all identical): 0.14 s / 0.56 GB vs 3.1 s / 4.5 GB on the slice,
  0.76 s / 1.55 GB vs 17.4 s / 22.2 GB on mouse1.
- **The cell-barcode filter streams integers instead of a per-row string
  frame (`ema/matrixfilter.py`).** The vectorised reader kept an object
  `cb_str` column per non-zero (~95 B/row): 2.71 GB for the slice's 14.9 M
  rows, an estimated 20–27 GB for the ~200 M rows of the full PBMC run —
  the next peak once the TF-IDF was fixed. Files are now parsed once in
  4 M-row chunks into three integer arrays (int32 where the values fit,
  12 B/row) and both passes run on those. Row semantics are unchanged; the
  Pass-1 grouping moves from the raw token to its integer, which is the same
  partition while every token is canonical, so a non-canonical token
  (`"007"`, `"+7"`, `"7.0"`), a >3-token row, a non-integer token or an
  absurd column index now defers the whole call to `_filter_cb_legacy` —
  the byte-for-byte copy of the reference algorithm. On such pathological
  input the result is therefore the reference's, which the previous fast
  path only approximated (it truncated float-valued tokens instead of
  skipping the row).

### Added

- **Peak calling runs one worker per (contig, strand)
  (`ema/countmatrix/chrom_parallel.py`, `--peak-workers`).** It was a single
  pure-Python thread streaming the whole BAM twice — 3 h 01 min of the
  3 h 46 min PBMC run at ~143 % CPU — and `--threads` never reached it
  (it fed only the `--tiles` pool and the downstream pool). Each job calls
  `peak_calling(region=(contig, 0, length))` in a spawned worker with its own
  `BarcodeIndex` and `PeakCallingState`; the merge then walks the jobs in the
  order the sequential run emitted them (`+` strand first, contigs in BAM
  header order, then `-`), renumbers `pas_id` across both strands (continuing
  from the previous BAM's last id in a multi-dataset run, as the legacy loop's
  never-reset `Peak.pasnumber` does), and
  rebuilds the shared barcode index by appending each job's local `cb.tsv` in
  local first-write order — which reproduces the sequential singleton's
  column assignment exactly. Matrices and support sidecars are re-keyed
  through those two maps.
  - Default when the BAM has an index; without one (or with
    `--peak-workers 1`) the legacy single-process caller runs unchanged and
    logs why. Worker count comes from `ResourceManager` (2.5 GB per worker;
    measured 1.5 GB on PBMC chr19 (+)), so `--threads` now governs peak
    calling, and BGZF threads are budgeted per worker so
    `workers × bam_threads` stays inside the ceiling.
  - Wall time is bounded by the largest contig, not by the BAM: 253 s for
    130 jobs on mouse1 (12 workers), and 674 s for the 252 jobs of the full PBMC BAM (16 workers, largest
    worker 3.0 GB) against 3 h 01 min of sequential streaming.
  - Per-job timings, PAS/cell counts and peak RSS land in
    `peakcalling/<id>_<n>.peak_jobs.json`.
- **`--peak-workers`** (YAML `peak_workers`).

### Changed

- `read_check` tests the strand — one flag bit — before the `get_tag(CB)`
  lookup, the `reference_end` CIGAR walk and the composite build; each
  strand pass rejects roughly half the BAM on exactly that test. Every
  rejection returns the same sentinel, so this reorder cannot change a
  result (pinned against a verbatim copy of the previous implementation over
  the full grid of read shapes).
- The `"<RG>_<CB>"` composite is interned per distinct pair, so every
  structure holding it shares one string object per cell — tracemalloc found
  102 MB of duplicate composites alive at a single chromosome flush of the
  slice.
- Region-mode `peak_calling` honours `bam_threads` instead of forcing
  `threads=1`.
- `ema/main.py` passes `default_threshold` / `merge_len` from
  `variable_config` explicitly on the single-process path. They were
  previously left to `peak_calling`'s default arguments, which froze the
  values `variable_config` held when `peackcalling.py` was imported — before
  the CLI/YAML bridge ran — so `--default-threshold` / `--merge-len` were
  silently ignored there (the `--tiles` path always read them live). At the
  defaults (5 / 100) nothing changes; a non-default value is now honoured on
  every path.

### Notes

- `--tiles` remains unusable for a real run and is untouched by this work:
  its workers crash under spawn because `barcode_tag` / `cb_len` / `seqlen`
  are never propagated, it builds jobs for every header contig, and
  `merge_tiles` sorts contigs as strings and unions barcodes per strand.
  The new path deliberately passes every config value in the job spec.
- Still single-process after peak calling: `make_dataframe`'s
  `mmread` → CSC → CSR → CSC chain and `preprocessing`'s transpose are the
  remaining multi-GB copies on a whole-genome run.

## Unreleased — read-level poly(A) evidence

### Added
- **`ema/countmatrix/polya.py`: PeakATail now reads poly(A) evidence off the
  reads.** `clip_site()` is a strand-aware terminal soft-clip A/T detector
  returning the inferred cleavage coordinate. Measured on the pbmc_10k_v3
  BAM: 1.152% of CB-bearing reads carry a qualifying clip, the wrong-end
  control fires on 0.0125% (92x specificity), and 73.9% of clip sites fall
  within 100 bp of a PolyASite 2.0 site (13.9% for arbitrary read 3' ends).
- **BED column 5 of `pasbed.bed` now carries per-PAS clip-read support.** It
  was hardcoded to `0`; every downstream reader already names it `score` and
  ignores it, and `annotatedpas.bed` inherits it. `--polya-evidence off`
  restores the previous byte-identical output.
- **New `clip_seeded` peak strategy — the accuracy fix.** Clip-site clusters
  (single-linkage, `--polya-seed-window`, read-weighted modal position) are
  PRIMARY PAS candidates; coverage peaks that overlap no cluster are emitted
  as an explicit SECOND TIER with `score == 0`. Seeding rather than filtering
  is required because ~60% of the clip evidence lies outside every coverage
  peak window.

  Measured end-to-end on chr19+chr21 of pbmc_10k_v3 and scored with
  `scripts/benchmark_tools/score_tool.py` (detected-gene-restricted
  PolyASite 2.0 atlas, both arms taken through the identical gene-assignment
  stage the shipped benchmark reports at):

  | arm | n | P@100 | R@100 | F1@100 |
  |---|---:|---:|---:|---:|
  | shipped caller (`lambda_gradient`) | 12,928 | 0.1610 | 0.1956 | **0.1766** |
  | `clip_seeded`, both tiers | 19,354 | 0.2371 | 0.3782 | **0.2914** |
  | `clip_seeded`, clip-supported tier only | 11,318 | 0.3591 | 0.3379 | **0.3482** |
  | `clip_seeded`, coverage-only tier | 8,036 | 0.0652 | 0.0404 | 0.0499 |

  Precision and recall both roughly double. The clip-supported tier's recall
  (0.3379) sits at the measured ceiling for this slice — 33.56% of its
  detected-gene atlas sites carry any clip read within 100 bp — so that tier
  extracts essentially all the recall this evidence type can provide, and the
  coverage-only tier is what carries the rest. **Always report the two tiers
  separately, with their n.**
- New flags alongside the `ip_*` block: `--polya-evidence`, `--polya-mode`
  (`annotate` default / `filter` / `require`), `--polya-min-clip`,
  `--polya-min-purity`, `--polya-window`, `--polya-seed-window`,
  `--polya-min-umis` (formerly `--polya-min-reads`, still accepted),
  `--polya-clip-filter`. `filter` and `require` are enforced at the existing
  `_apply_pas_filters()` seam and drop exactly the coverage-only tier.
- **`pas_support.tsv` — per-PAS poly(A) support sidecar.** Written next to
  every caller BED (`<bed>.support.tsv`) and, for single-BAM runs, merged to
  the run root as `pas_support.tsv`, keyed by the BED's own `pas_id`:
  `clip_reads` (raw), `clip_umis` (== BED column 5), `clip_reads_f3844`,
  `clip_umis_f3844`, `window_reads` (the reads counted into that PAS's
  matrix row) and `tier`. `pasbed.bed` stays plain BED6 — no parser sees a
  new column and the name field is untouched. All three peak-calling paths
  write it (the tile merge re-keys it with the same pasnumber remapping the
  BED gets).
- **Startup warning when the observed clip rate is below 0.3% of CB reads** —
  a pipeline that trims poly(A) before alignment destroys this evidence
  channel, and the caller now says so instead of silently emitting an
  unsupported call set.

### Fixed
- **Internal-priming filter now tests the correct side on the `-` strand.**
  `ema/experimental/internal_priming.py` applied the `+` genomic window
  `[pos-10, pos+30)` to both strands, so on `-` it scanned 30 bp *upstream* /
  10 bp downstream of the cleavage site in transcript orientation — mostly
  the wrong side for an oligo-dT priming artifact. The window is now defined
  in transcript orientation (`--ip-window-left` upstream, `--ip-window-right`
  downstream, both strands) and mirrored in genomic coordinates on `-`
  (`[pos-30, pos+10)`, reverse-complemented before the A-run / A-fraction
  scan). `+`-strand calls are byte-identical; flags, defaults and the
  `--ip-*` CLI are unchanged. New pure helpers `ip_window()`,
  `call_internal_priming()` and `check_internal_priming()` are unit-tested
  in `tests/test_internal_priming_strand.py` (implanted runs on a C/G-only
  genome, both strands, fraction rule, contig edges, end-to-end through
  `filter_internal_priming` and `apply_filters`).

  Measured post hoc on the Stage-2 final run (4efeb12) by re-applying both
  rules to the emitted `pasbed.bed` without re-running the caller: `+`
  strand, 0 changes on every dataset. On `-`, the corrected rule flags
  5,599 / 160,464 (3.5 %) of the PBMC sites that had passed the old filter
  (testis mouse1 1,505 / 35,487 = 4.2 %, mouse2 1,302 / 35,946 = 3.6 %),
  and of the 39,244 `-` sites the old rule removed from the PBMC candidate
  set, 11,217 would be kept by the corrected rule. On the pre-registered
  precision default (PBMC, tier-1, ≥2 molecules; PolyASite 2.0 within
  100 bp) removing the 386 newly-flagged sites moves P@100 0.7167 → 0.7184;
  the full corrected filter (newly flagged removed AND wrongly removed sites
  restored, 46,524 vs 44,394 sites) gives P@100 0.7062, R@100 0.1754 (was
  0.1707). Every pre-registered gate still passes; benchmark conclusions
  are unchanged. Any use of the per-site `internal_priming` flag from a
  pre-fix run should re-run the filter (`ema reannotate --genome-fasta`).
- **BED column 5 is now DISTINCT MOLECULES, the unit the flag always gated
  on.** `--polya-min-reads` documented "distinct molecules" and
  `ClipSeeder.flush()` did gate on them, but the score column wrote the raw
  clip-READ count — PCR duplicates and secondary alignments included. On the
  Stage-2 PBMC run that made 18,864 of the 79,751 "≥2-read" tier-1 sites
  (23.7%) a single molecule counted twice, and 47.7% of the exactly-2-read
  sites. The gate is unchanged, so **no PAS moves and none is added or
  dropped at the default `--polya-min-umis 1`**: only column 5 changes value
  (proven on the chr19+21 PBMC slice — columns 1–4 and 6 byte-identical to
  the previous commit, all four region BEDs). Raw clip reads remain
  available in `pas_support.tsv`.
  - `--polya-min-reads` → **`--polya-min-umis`**; the old spelling still
    works (CLI and YAML) and logs a deprecation warning.
  - A molecule is `(cell barcode, UMI)`. `read_check` keys cells by
    `RG + "_" + barcode`, so molecules deliberately drop the read-group
    prefix (`molecule_cb`): a BAM with one read group per lane would
    otherwise count one molecule once per lane (+4.7% on the PBMC chr21 (+)
    slice: 13,100 vs 12,483 molecules over identical clusters). Matrix cell
    identity is unchanged.
  - Reads with no `UB` tag still count as one molecule each.
- **Alignment filtering on the clip-evidence channel (`--polya-clip-filter`,
  `clip_read_ok` == samtools `-F 3844`).** `read_check` applies no
  secondary / supplementary / duplicate / qcfail / unmapped filter (it is
  the coverage path's contract; widening it would move every peak), so every
  clip read is now classified and BOTH counts are reported per PAS.
  `--polya-clip-filter f3844` makes the filtered reads stop being evidence
  (score, gate and tier tag all use the `-F 3844` molecule count).
  **The default stays `none`, deliberately**: UMI de-duplication already
  makes PCR duplicates uncountable (duplicates share their `(CB, UMI)` key),
  while dropping flagged reads outright also drops every cluster whose
  evidence is entirely such alignments — measured on PBMC as −8.3% of
  chr19 (+) and −27.0% of chr21 (+) tier-1 clusters (STARsolo emits
  multimappers as secondary; CellRanger flags ~50% duplicates). That is a
  call-set change and must be scored, not slipped into a default; the
  `*_f3844` sidecar columns let it be measured post hoc without a re-run.
- **The tier-1 clip fallback no longer counts a read the neighbouring
  cluster already counted.** A cluster whose candidate partition and
  midpoint-clipped cleavage window are both empty falls back to its own clip
  reads; those reads' `end1` could lie past the midpoint, inside the
  neighbour's window, and were counted twice (verified at `13:43135248+`,
  pas#23403 on mouse1: 37 clip reads with `end1` 43,135,267–43,135,294 past
  its clipped window's 43,135,262 upper bound; the row summed to 371 where
  the shipped candidate + remainder is 334). The fallback is now
  midpoint-clipped and candidate-excluded exactly like the window, so
  "no read end is ever counted twice" is true for every row; a fallback row
  can legitimately end up empty (its reads belong to the neighbour) and the
  new `empty_tier1_rows` stat counts those. On mouse1 this touched 5,250
  rows and ~8k of 140.7M counts (~0.006%).
- **`clip_seeded` tier-1 PAS are now counted from the reads that pile up at
  the cleavage site, not from their poly(A)-clipped reads.** The first
  full-data run of the strategy (GSE104556 mouse1 testis) surfaced a ~10x
  quantification regression: every tier-1 row of the count matrix held only
  the cluster's clip reads (~1% of the reads), so the annotated matrix carried
  4.1M counts against 41.9M from the shipped caller and the `min_read` cell
  filter lost 23% of the real STARsolo cells. The BED coordinates and tier
  tags were — and remain — byte-identical; only the matrix changes:
  - a cluster that suppresses a coverage candidate takes that candidate's
    accumulated `cb_positions` counts, restricted to its partition of the
    candidate (several clusters split a candidate at the midpoints between
    their anchors — the same `partition_peak_region` rule `find_pas` applies
    to multi-PAS peaks — each share counted with `reconstruct_cb_dict`
    semantics on a compact `_CandidateSlice` of the candidate's own
    `cb_positions`, which is what every positional coverage strategy's
    `get_cb_dict_for_pas` delegates to; a non-positional strategy such as
    `original` is detected by the shares not summing to the candidate and
    the candidate then goes whole to the nearest cluster), so the shipped
    mass is conserved exactly;
  - every cluster additionally counts the accepted read ends in its cleavage
    window (`--polya-count-window`, default `auto,25` ==
    `[site - seq_len, site + 25]` in transcript orientation) that belong to
    no coverage candidate, clipped at the midpoint to neighbouring clusters,
    so nothing is counted twice. For a cluster outside every coverage peak
    this is its whole count; inside a peak it recovers the reads the
    coverage caller counted into `cb_positions` but never put in a candidate
    (the streaming loop drops the ends still in `data_array` when a peak
    closes, so every candidate interval stops ~`seq_len` bp short of the
    pile's 3' end — on `+` those are exactly the reads at the cleavage site);
  - the clip-read count stays in BED column 5 as the support annotation.
  Tier-2 rows are untouched. Implemented once in `ClipSeeder.flush()` and
  fed by all three peak-calling paths (`ClipStream` records every accepted
  read end as two int32 arrays per chromosome — 8 bytes/read, released at
  each chromosome flush); `tests/test_polya_three_path_agreement.py` now pins
  per-PAS counts as well as coordinates, and
  `tests/test_polya_tier1_counts.py` pins the accounting on a synthetic BAM
  (a cluster inside a 200-read peak carries the peak, not ~2; an isolated
  cluster with 30 in-window reads + 3 clip reads counts 33; tier-2 rows are
  byte-for-byte the coverage strategy's). Each strand pass logs a mass
  accounting (`clip_seeded counting: {...}`) so a run can show where every
  count came from.

  Measured on GSE104556 mouse1 (`--seq-len 98`, 12 threads, 55 min wall):
  every caller BED (`posbed.bed`, `negbed.bed`, `01_peak_calling/*`,
  `annotatedpas.bed`) is byte-identical to the Stage-2 run. Raw matrices:
  140.7M counts vs 132.9M shipped (1.059x) — the suppressed candidates'
  partitions plus the tier-2 rows reproduce the shipped per-strand totals
  to the read (`+` 58,221,832 + 5,131,499 = 63,353,331; `-` 64,246,477 +
  5,330,874 = 69,577,351), the +7.8M is the cleavage-window remainder.
  `min_read` keeps 10,339 cells including 1,294/1,294 STARsolo cells (the
  clip-only run kept 995). Gene-assigned mass on the shipped run's cells:
  51.8M vs 49.2M shipped when the shipped `filterdmatrix.mtx` is keyed
  correctly (1.053x; per-gene Spearman 0.991, median per-gene ratio 1.039).
  The shipped `annotated_matrix.mtx` itself is mis-keyed (bug 0a: 226 of
  45,921 rows aligned) and sums to 31.8M counts — the often-quoted
  41,944,353 is that figure plus its MatrixMarket size line. Post-filter
  `pasbed.bed` grows from 58,848 to 85,993 PAS because 27k tier-1 rows now
  carry enough counts to pass `min_cells`; 120 one-to-two-read clusters
  drop out because their counts moved to cells below `min_read`.

  **Measured on GSE104556 mouse1** (`--seq-len 98`, 12 threads, 54 min wall,
  26.9 GB peak RSS — unchanged from the previous run's 55 min / 26.9 GB;
  `results/benchmark_tools/gse104556/peakatail_clipseeded_v3/mouse1`):
  - every caller BED (`posbed.bed`, `negbed.bed`, `01_peak_calling/*`,
    `peakcalling/*`, `raw/pas.bed`) is identical in columns 1–4 and 6 to the
    previous run; only column 5 changes (`annotatedpas.bed` is byte-identical
    — it does not carry the score column);
  - clip reads 2,181,334 → 2,091,279 molecules (0.959×) over 148,715 tier-1
    rows; 13.3% (+) / 13.2% (−) of tier-1 rows change value; **3,269 of the
    62,159 "≥2 clip read" sites (5.3%) are a single molecule** (the same
    statistic is 23.7% on PBMC, where the duplication rate is far higher);
  - counting is otherwise unchanged: `reads_from_candidates`,
    `reads_from_window`, `reads_tier2`, `tier1`, `tier2` and
    `suppressed_candidates` are identical to the read; the raw caller
    matrices lose exactly the double-counted fallback reads
    (140,744,000 → 140,743,255, −745 counts, −0.0005%) and the annotated
    matrix 52,703,320 → 52,703,105 (−215, −0.0004%, 85,816 vs 85,993 PAS —
    the emptied fallback rows fall below `min_cells`);
    562 of the 5,250 fallback rows (98 on `+`, all 464 on `−`) were entirely
    reads a neighbouring cluster had already counted;
  - `min_read` keeps the same 10,339 cells, including 1,294/1,294 STARsolo
    cells;
  - scoring is unmoved: tier-1 P@100 0.4992, R_det 0.2991, F1 0.3741 —
    identical to 4 dp before and after (both tiers 0.4093/0.3231/0.3611).
  - **post-hoc sensitivity, labelled as such**: gating tier-1 at ≥2
    MOLECULES gives n 30,346, P@100 0.6923, R_det 0.2164, F1 0.3298. On the
    chr19+21 PBMC slice the same comparison, at equal nominal thresholds,
    separates the two units: ≥2 reads n 4,749 / P 0.5519 / R 0.2442 /
    F1 0.3386 versus ≥2 molecules n 3,664 / P 0.6395 / R 0.2245 / F1 0.3323
    (≥1 is unchanged at n 11,318 / 0.3591 / 0.3379 / 0.3482). A read
    threshold and a molecule threshold are not the same operating point.

### Notes
- The evidence is computed in all three peak-calling paths (monolithic,
  `--pipeline`, `--tiles`) from the `AlignedSegment` at the call site;
  `read_check`'s 5-tuple return is deliberately unchanged.
  `tests/test_polya_three_path_agreement.py` pins that the three paths
  produce identical BED output and identical per-PAS counts.
## Unreleased

### Added
- **`--cleavage-offset` — data-driven 3' cleavage-site offset correction
  (issue #72).** Called peak 3' ends stop ~90–105 nt short of the true
  cleavage site because 10x R2 coverage runs out before the poly(A) junction
  (AATAAA density peaks +75 nt downstream of the peak end; genomic A-fraction
  crests at +98 nt). Tight-cutoff benchmarks then punish the offset, not the
  calls. `--cleavage-offset N` (default `0` = legacy/off) shifts each reported
  PAS 3' end downstream by `N` bp (strand-aware, clamped at 0) after peak
  calling, so `pasbed.bed`, `annotatedpas.bed`, gene assignment, atlas
  matching, and the benchmark harness all use the inferred cleavage position.
  A sane fallback constant is `95`. The per-run data-driven estimator
  (`ema/countmatrix/cleavage_offset.py::estimate_cleavage_offset`) is stubbed
  with a clear `TODO(issue #72)`; it currently returns the constant.

## 0.2.0 (2026-05-10) — Product CLI

### Breaking changes
- The `ema_switch`, `ema_merge`, `ema_parse_gtf` console scripts are removed.
  Use `ema switch {diff,length,match}`, `ema merge`, `ema parse-gtf` instead.
- `ema --pdui-*` and `ema --diff-method` flags are dropped (moved to `ema switch length` / `ema switch diff` in 0.1.b).
- Output directories now have a timestamp suffix (e.g. `emaout_2026-05-10_181523/`).
- **BREAKING (bug fix): YAML `min_read` / `min_cells` / `min_genes` /
  `min_pas_per_cell` now actually take effect.** Previously the legacy
  argparse shim only bridged `min_pas_per_cell` from YAML — the other three
  filter thresholds silently kept their dataclass defaults (2000, 3, 50)
  regardless of YAML values. The atlas baseline in
  `reports/baseline_apa_completeness/` was regenerated to reflect this
  fix: per-dataset clusters now have shape (393, 5755) instead of
  (256, 5467) when `min_read: 1500` is set in
  `test_run/atlas_full.yaml`. Peak counts, BAM column counts, and
  atlas-snapped unified peaks remain byte-identical to the prior
  baseline — only downstream cluster sizes change because more cells
  survive the looser min_read threshold.
- **BREAKING (bug fix): `run_merge` now actually writes to `--output`.**
  Previously the function called `merge(..., threads=...)` without
  forwarding the `output` argument, so pysam wrote to the default
  `Aligned.sortedByCoord.merged.out.bam` and the user's `--output` flag
  was silently ignored.
- **BREAKING (bug fix): `--cluster-method` / `--resolution` / `--n-pcs` /
  `--random-seed` / `--external-clusters` now take effect.** Previously
  the pipeline called `clustering(adata=adata)` with no kwargs in both
  the single-sample and multi-sample paths, so all five flags were
  silently dropped. Default behavior is unchanged (defaults match the
  old hard-coded values), but explicitly setting any of these now
  affects the output.
- **BREAKING (bug fix): strategy hyperparameters `--lambda-method` /
  `--max-pas` / `--smoothing-window` / `--min-prominence` now take
  effect.** Previously these were declared in argparse but never
  forwarded into `get_strategy(name, **kwargs)`. Now passed through
  using `inspect.signature` filtering so each strategy receives only
  the kwargs it accepts.

### Loud no-op warnings
- `--ip-filter` / `--genome-fasta` / `--annot-filter` / `--ip-a-stretch` /
  `--benchmark` / `--validate-db` are accepted by `ema run` but emit a
  WARNING when set, telling the user the flag is currently a no-op
  (the underlying `apply_filters` and `run_benchmark` integrations are
  not wired into the new pipeline body). No silent fallbacks.

### Added
- Click-based CLI with subcommand structure
- Interactive wizard launched by bare `ema`
- Rich logging + per-run timestamped log file
- Nested progress bars, spawn-worker safe
- pyproject.toml + uv.lock for both pip and uv users
- `ema run --list-strategies` discovery
- New `--pas-gap` CLI flag (was YAML-only)
- New common options on every subcommand: `-v/-vv/-q`, `--log-level`,
  `--no-log-file`, `--no-progress`

### Changed
- Package PyPI name: `peakatail` (was `ema`); import path stays `ema`
- Build backend: setuptools via pyproject.toml (was setup.py)
