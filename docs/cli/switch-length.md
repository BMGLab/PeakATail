# `peakatail switch length`

`peakatail switch length` quantifies 3' UTR shortening and lengthening patterns by
computing a per-cell, per-gene score from the PAS count matrix in one or more
`clusters.h5ad` files. Three strategies are available: `classic` (2-PAS PDUI),
`proportion` (full per-PAS proportion vector), and `shannon` (entropy of PAS
usage). Each strategy writes a separate TSV file with a different score column.

When `--output` is left at its default and the `--h5ad` files come from a
recognisable `peakatail_runs/<run>/` directory, output is routed **inside** the
originating run directory as `peakatail_runs/<run>/switch_length_<timestamp>/`.
Source: `ema/cli/common.py::resolve_subcommand_output_dir`.

!!! note "When to use it"
    - You want a per-cell PDUI score (fraction of reads at the distal PAS)
      for each gene to quantify 3' UTR lengthening.
    - You want the full per-PAS proportion vector to feed into downstream
      differential proportion analyses.
    - You want to measure how dispersed PAS usage is within a gene using
      Shannon entropy.

!!! warning "When NOT to use it"
    - You want pairwise statistical tests between clusters. Use `peakatail switch diff`.
    - You are using `--isoform-agg=per_isoform` without providing `--gtf`. The
      command falls back to `per_gene` and emits a warning.
    - Your h5ad `var` has no `gene_id` column (older files from before the
      annotation refactor). PDUI output will be empty; re-run `peakatail run` first.

## Quick example

```bash
# Classic PDUI (default strategy)
uv run peakatail switch length \
  --h5ad peakatail_runs/emaout_2026-05-11_120000/per_dataset/sample1/clusters.h5ad \
  --strategy classic \
  --isoform-agg per_gene \
  --pdui-pseudocount 1.0

# Shannon entropy with isoform-level aggregation
uv run peakatail switch length \
  --h5ad peakatail_runs/emaout_2026-05-11_120000/per_dataset/sample1/clusters.h5ad \
  --strategy shannon \
  --isoform-agg per_isoform \
  --gtf /path/to/gencode.v44.annotation.gtf
```

What lands on disk after the first command:

- `peakatail_runs/emaout_.../switch_length_<ts>/pdui_classic.tsv` — long-format PDUI scores.
- `peakatail_runs/emaout_.../switch_length_<ts>/peakatail_<ts>.log` — run log.
- Heatmap and UMAP overlay figures in `switch_length_<ts>/figures/` (when plotting is enabled).

## Full `--help` output

```text
Usage: peakatail switch length [OPTIONS]

  3'UTR shortening / lengthening quantification (PDUI variants).

Options:
  --list-strategies               Print available length strategies and exit.
  --threads INTEGER               Max parallel workers (auto-detected if not
                                  set). Respected by ResourceManager as an
                                  absolute ceiling.
  -v, --verbose                   Increase verbosity. -v = DEBUG for ema.*;
                                  -vv = DEBUG everywhere.
  -q, --quiet                     WARNING and up only. Overrides --verbose.
  --log-level TEXT                Explicit logger level (DEBUG/INFO/WARNING
                                  /ERROR) or `logger.name=LEVEL` (repeatable:
                                  comma-separated).
  --no-log-file                   Don't write peakatail_<ts>.log next to the
                                  outputs.
  --no-progress                   Suppress Rich progress bars.
  -c, --config PATH               YAML config; CLI flags override individual
                                  keys.
  -o, --output PATH               Output directory (timestamp suffix added
                                  automatically).
  --plot-engine TEXT              Engines: 'matplotlib' (default), 'plotly',
                                  'both', 'none', or comma list.
  --plot-format TEXT              Restrict output formats. Default 'all' =
                                  png+svg+html as appropriate.
  --no-plots                      Disable all plotting (alias for --plot-
                                  engine none).
  -i, --h5ad PATH                 [required]
  --gtf PATH                      Required when --isoform-agg=per_isoform.
  --cluster-pairs TEXT
  --cluster-key TEXT              [default: leiden]
  -s, --strategy TEXT             [default: classic]
  --isoform-agg [per_gene|per_isoform]
                                  Aggregation level — must match strategy
                                  vocabulary (per_gene collapses isoforms,
                                  per_isoform keeps them).  [default:
                                  per_gene]
  --isoform-collapse [none|mean|majority]
                                  How to collapse multiple isoforms when
                                  --isoform-agg=per_gene and the strategy
                                  tracks isoforms internally.  [default: none]
  --pdui-pseudocount FLOAT        Pseudocount added to counts before
                                  PDUI/entropy computation. Default 0.0
                                  (original behaviour). Use 1.0 to avoid NaN
                                  on zero-count cells.  [default: 0.0]
  --help                          Show this message and exit.
```

## Flags

### Inputs

| Flag | Type | Default | Description |
|---|---|---|---|
| `--h5ad` / `-i` | PATH (repeatable) | — | One or more `clusters.h5ad` files from `peakatail run`. Results are computed per h5ad independently; only the last result is returned to the viz hooks. Required. |
| `--gtf` | PATH | — | Ensembl/GENCODE GTF file. Required when `--isoform-agg=per_isoform`. When absent and `per_isoform` is requested, the runner logs a warning and falls back to `per_gene`. |
| `--cluster-key` | TEXT | `leiden` | The `adata.obs` column holding cluster labels. Used to add a `cluster` column to the augmented output TSV. |
| `--cluster-pairs` | TEXT | — | Currently unused by the length analysis. The flag is accepted but any value logs a warning and has no effect. Reserved for future per-pair PDUI comparisons. |

### Strategy options

| Flag | Type | Default | Description |
|---|---|---|---|
| `--strategy` / `-s` | TEXT | `classic` | PDUI quantification method. Options: `classic`, `proportion`, `shannon`. Run `peakatail switch length --list-strategies` to see all registered names. |
| `--isoform-agg` | CHOICE | `per_gene` | Aggregation level. `per_gene` collapses all isoforms of a gene and selects proximal/distal by genomic rank. `per_isoform` computes the score independently per transcript using UTR structure from the GTF. Use `per_gene` for speed; `per_isoform` for isoform-resolution results. |
| `--isoform-collapse` | CHOICE | `none` | How to collapse isoform-level scores when `--isoform-agg=per_gene`. `none` leaves them separate; `mean` averages across isoforms; `majority` takes the dominant value. Not used by the `classic` or `shannon` strategies; relevant for `proportion`. |
| `--pdui-pseudocount` | FLOAT | 0.0 | Pseudocount added to each per-cell PAS count before computing PDUI, proportions, or entropy. The default `0.0` preserves the original behaviour exactly. Set to `1.0` to avoid `NaN` in the output for cells with zero reads at a gene. Note that any non-zero pseudocount shifts entropy toward uniformity. |

### Output

| Flag | Type | Default | Description |
|---|---|---|---|
| `--output` / `-o` | PATH | `switch_out` | Output directory base name. Auto-routed inside the originating run dir when inputs are from a `peakatail_runs/` path. |

### Performance

| Flag | Type | Default | Description |
|---|---|---|---|
| `--threads` | INT | auto | Absolute thread ceiling passed to `ResourceManager`. The length runner uses a `threading` joblib backend (not `loky`) to avoid pickling the large count matrix into each worker. Source: `ema/switch_test/runner.py::run_length` lines 581–588. |

## Strategies and output files

### `classic` — 2-PAS PDUI

**Output file:** `pdui_classic.tsv`

Computes the fraction of reads at the distal PAS relative to proximal + distal:

```
PDUI = distal_reads / (proximal_reads + distal_reads + pseudocount)
```

For genes with more than 2 PAS, proximal is the first and distal is the last
in transcription order (strand-aware). Genes with fewer than 2 PAS are
excluded.

With `--isoform-agg=per_isoform` the same rule is applied per transcript, and
it counts **distinct** PAS: transcripts with fewer than 2 **distinct** PAS are
excluded. The same PAS can be listed more than once for a transcript — one
entry per UTR exon it intersects — and those repeats are collapsed on `pas_id`
before the endpoints are chosen. A transcript left with a single usable PAS
would otherwise be emitted as a degenerate `proximal == distal` pair whose
`pdui` is `0.5` by construction, carrying no length information; it is dropped
instead, and a `WARNING` names the gene and the number of transcripts dropped
(issue #98). Source: `ema/quantification/strategies/classic.py`.

**Columns:**

| Column | Type | Description |
|---|---|---|
| `gene_id` | str | Gene identifier. |
| `transcript_id` | str | `"_gene_"` sentinel when `--isoform-agg=per_gene`; actual transcript ID when `per_isoform`. |
| `proximal_pas_id` | int | PAS ID of the proximal site used. Never equal to `distal_pas_id` — degenerate pairs are dropped, not written (see above). |
| `distal_pas_id` | int | PAS ID of the distal site used. |
| `cell` | str | Cell barcode. |
| `pdui` | float64 | PDUI in [0, 1]. `NaN` when `proximal_reads + distal_reads = 0` and `pseudocount = 0.0`. |
| `proximal_reads` | float64 | Raw read count at the proximal PAS for this cell. |
| `distal_reads` | float64 | Raw read count at the distal PAS for this cell. |
| `total_reads` | float64 | Sum of proximal + distal reads (denominator before pseudocount). |
| `cluster` | str | Cluster label from `--cluster-key` (added by augment helper). |
| `chrom`, `start`, `end`, `strand` | str/int | Genomic coordinates from `pasbed.bed` (empty if not found). |

### `proportion` — full per-PAS proportion vector

**Output file:** `proportion.tsv`

For each (gene, cell), computes the fraction of reads falling at each PAS so
that all proportions for a gene in a cell sum to 1.0. Unlike `classic`, this
preserves information about all PAS, not just the proximal/distal pair.
Source: `ema/quantification/strategies/proportion.py`.

**Columns:**

| Column | Type | Description |
|---|---|---|
| `gene_id` | str | Gene identifier. |
| `transcript_id` | str | `"_gene_"` or transcript ID (matches `--isoform-agg`). |
| `pas_id` | int | PAS identifier. |
| `rank` | int | Proximal-to-distal rank (0 = most proximal). |
| `cell` | str | Cell barcode. |
| `proportion` | float64 | Fraction of gene reads at this PAS. Sums to 1.0 per (gene, cell). `NaN` when total is 0. |
| `reads_at_pas` | float64 | Raw read count at this PAS (before pseudocount). |
| `total_reads_gene` | float64 | Total reads across all PAS of this gene in this cell (per_gene). Or `total_reads_transcript` when `per_isoform`. |
| `cluster` | str | Cluster label. |
| `chrom`, `start`, `end`, `strand` | str/int | Coordinates from `pasbed.bed`. |

### `shannon` — entropy of PAS usage

**Output file:** `entropy_shannon.tsv`

Computes the Shannon entropy (in bits) of the PAS read distribution for each
(gene, cell):

```
H = -Σ p_i * log2(p_i)
H_norm = H / log2(N)    (N = number of PAS with p_i > 0)
```

`H = 0` when all reads go to one PAS (completely focused usage). `H = log2(N)`
for uniform distribution across N PAS (maximally dispersed usage).
`H_norm` normalises to [0, 1]. Source: `ema/quantification/strategies/shannon.py`.

**Columns:**

| Column | Type | Description |
|---|---|---|
| `gene_id` | str | Gene identifier. |
| `transcript_id` | str | `"_gene_"` or transcript ID. |
| `pas_ids` | str | Semicolon-separated PAS IDs included in this gene/transcript. |
| `cell` | str | Cell barcode. |
| `entropy` | float64 | Shannon entropy H in bits. `NaN` when total is 0. |
| `normalized_entropy` | float64 | H / log2(N). `NaN` for single-PAS genes (log2(1) = 0). |
| `n_pas` | int | Number of PAS in the gene/transcript. |
| `total_reads_gene` | float64 | Total reads at all PAS for this cell (per_gene). Or `total_reads_transcript` when per_isoform. |
| `cluster` | str | Cluster label. |

## How it relates to other commands

- **[`peakatail run`](run.md)** — produces the `clusters.h5ad` inputs. The `gene_id` column in `adata.var` is required for all strategies; it is written by `peakatail run` when `--gtf` is provided.
- **[`peakatail switch geneview`](switch-geneview.md)** — draws per-cluster PAS coverage tracks. It overlays differential results via `--diff-tsv`; it does not consume the length TSVs.
- **[`peakatail switch diff`](switch-diff.md)** — complementary pairwise test; combine with `switch length` to characterise both significance and magnitude of APA changes.

## See also

- Strategy pages in [`../strategies/`](../strategies/) — mathematical details of PDUI computation.
- Tutorial in [`../tutorials/`](../tutorials/) — length quantification walkthrough.

## Flags added in the strand/counts correctness fix (2026-08)

| Flag | Meaning |
|---|---|
| `--pasbed PATH` | PAS coordinate BED used for strand-aware proximal/distal ordering. **Required** when the h5ad lacks embedded coordinates; the command now RAISES instead of silently falling back to a plus-strand convention. |
| `--counts-layer NAME` | AnnData layer holding raw integer counts (default `counts`). PDUI is computed on this layer, never on `.X`, which clustering overwrites with TF-IDF weights. |
| `--allow-non-count-matrix` | Escape hatch: proceed even if the selected matrix is non-integral. Off by default; without it a non-count matrix is an error, not a warning. |

**Breaking change:** runs that previously "worked" by silently using TF-IDF-normalised `.X` and
plus-strand ordering now fail fast with an actionable message. This is deliberate: those runs produced
inverted PDUI for 100% of minus-strand genes.

## Degenerate `proximal == distal` pairs are dropped (2026-09, issue #98)

Before this fix, `--isoform-agg=per_isoform` could select the **same PAS** as
both the proximal and the distal endpoint of a transcript, because a PAS is
listed once per UTR exon it intersects and those repeats were ranked as if
they were separate sites. The resulting rows have
`proximal_pas_id == distal_pas_id`, `proximal_reads == distal_reads` and
`pdui` exactly `0.5` — a value forced by the arithmetic, not measured from the
data.

**What changed for you:**

- `pdui_classic.tsv` written with `--isoform-agg=per_isoform` now has **fewer
  rows** than a pre-fix run. On the GSE104556 testis runs reported in the
  issue, 18,116 of 13,709,930 rows (7 genes) disappear on mouse1 and 15,004
  (8 genes) on mouse2. Every dropped row had `pdui == 0.5`, so no real PDUI
  value is lost — but row counts, transcript counts and the `0.5` bin of a
  per-isoform PDUI histogram will not match a pre-fix run.
- `--isoform-agg=per_gene` output is unchanged: it never produced such a pair.
- The strand/proximal-distal guard now reports three separate counts —
  `I INVERTED, X with proximal and distal on DIFFERENT strands, D DEGENERATE`
  — instead of folding `proximal_start == distal_start` into the inverted
  count and blaming strand-blind selection for it. A `DEGENERATE` row is still
  a hard error, but it points at single-PAS / duplicate UTR-exon mapping, not
  at strand.
