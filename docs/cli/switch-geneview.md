# `peakatail switch geneview`

`peakatail switch geneview` renders gene-track panels showing per-cluster PAS
coverage and proportions for a set of genes. For each gene, it builds a panel
with depth-normalised per-cluster PAS read bars, optional isoform structure
(from a GTF), and optional differential-result overlays (from `--diff-tsv`).
The output is one figure file per gene, a `meta.json` manifest, and a
`figures_INDEX.md` index.

Genes are selected from the union of:

1. Explicit `--gene-id` arguments (normalised to uppercase).
2. Top-N genes auto-ranked from `--diff-tsv` results by volcano score.

At least one of `--diff-tsv` or `--gene-id` must be provided.

When `--output` is left at its default and the `--h5ad` file comes from a
recognisable `peakatail_runs/<run>/` directory, output is routed **inside** the
originating run directory as `peakatail_runs/<run>/switch_geneview_<timestamp>/`.
Source: `ema/cli/common.py::resolve_subcommand_output_dir`.

!!! note "When to use it"
    - You have run `peakatail switch diff` and want to visualise the top differentially
      used genes as per-cluster PAS track plots.
    - You have a specific gene of interest and want to see how PAS usage differs
      across all clusters.
    - You want to overlay PDUI scores from `peakatail switch length` on the per-cluster
      coverage bars to connect quantitative scores to the raw data.

!!! warning "When NOT to use it"
    - You have not provided either `--diff-tsv` or `--gene-id`. The command raises
      a `UsageError` and exits without writing any files.
    - The gene identifier you supply does not appear in `adata.var` or `pasbed.bed`.
      The command logs a warning and skips that gene (no figure is written).
    - The `--h5ad` file has no `gene_id` column in `var`. Gene-to-PAS mapping will
      fail and no panels will be rendered.

## Quick example

```bash
# Auto-pick top 10 genes from a diff result
uv run peakatail switch geneview \
  --diff-tsv peakatail_runs/emaout_.../switch_diff_<ts>/differential/fisher_0_vs_1.tsv \
  --pasbed peakatail_runs/emaout_.../per_dataset/sample1/pasbed.bed \
  --h5ad peakatail_runs/emaout_.../per_dataset/sample1/clusters.h5ad \
  --top-genes 10

# Explicit gene with isoform overlay and PDUI score
uv run peakatail switch geneview \
  --gene-id ACTB \
  --gene-id MYC \
  --diff-tsv peakatail_runs/emaout_.../switch_diff_<ts>/differential/fisher_0_vs_1.tsv \
  --pasbed peakatail_runs/emaout_.../per_dataset/sample1/pasbed.bed \
  --h5ad peakatail_runs/emaout_.../per_dataset/sample1/clusters.h5ad \
  --gtf /path/to/gencode.v44.annotation.gtf \
  --plot-engine plotly
```

What lands on disk after the first command:

- `peakatail_runs/emaout_.../switch_geneview_<ts>/figures/gene_<GENEID>.png` — one PNG per rendered gene.
- `peakatail_runs/emaout_.../switch_geneview_<ts>/figures/meta.json` — manifest listing all rendered figures.
- `peakatail_runs/emaout_.../switch_geneview_<ts>/figures/figures_INDEX.md` — Markdown index of all figure files.
- `peakatail_runs/emaout_.../switch_geneview_<ts>/peakatail_<ts>.log` — run log.

## Full `--help` output

```text
Usage: peakatail switch geneview [OPTIONS]

  Gene-track visualisation: per-cluster PAS coverage and proportions.

  Renders one panel per gene showing:
  - Per-cluster PAS read coverage (depth-normalised)
  - Per-cluster within-gene PAS proportions
  - Gene isoform structure (when --gtf is supplied)

  Genes are selected from the union of --gene-id and the top-N genes
  ranked by volcano score from --diff-tsv results.

Options:
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
  --diff-tsv PATH                 Differential TSV from `peakatail switch diff`
                                  (repeatable). Used to auto-select top-N genes
                                  by volcano score. Optional when --gene-id is
                                  supplied.
                                  switch length`. When present, overlays
                                  strategy scores on per-cluster bars.
  --pasbed PATH                   pasbed.bed — PAS coordinates (BED6 format).
                                  [required]
  -i, --h5ad PATH                 clusters.h5ad with leiden cluster labels and
                                  gene_id in var.  [required]
  --gtf PATH                      Ensembl GTF for isoform structure (optional).
  --top-genes INTEGER             Number of top genes to auto-pick from --diff-
                                  tsv results.  [default: 10]
  --gene-id TEXT                  Explicit gene ID(s) to render (repeatable).
                                  Union-ed with --top-genes list.
  --cluster-key TEXT              obs column carrying cluster labels.
                                  [default: leiden]
  --help                          Show this message and exit.
```

## Flags

### Inputs

| Flag | Type | Default | Description |
|---|---|---|---|
| `--h5ad` / `-i` | PATH | — | Single `clusters.h5ad` file from `peakatail run`. Must have `gene_id` in `adata.var` (written when `--gtf` was provided to `peakatail run`). The Leiden cluster labels (`adata.obs[cluster_key]`) drive the per-cluster display. Required. |
| `--pasbed` | PATH | — | PAS BED file in BED6 format (chrom, start, end, pas_id, score, strand). Loaded to retrieve genomic coordinates for PAS track positioning. Required. |
| `--diff-tsv` | PATH (repeatable) | — | One or more differential TSV files from `peakatail switch diff`. Each TSV is parsed for `cluster1`, `cluster2`, `gene_id`, and volcano columns. Multiple TSVs for the same cluster pair are concatenated. Optional when `--gene-id` is provided; required otherwise. |
| `--gtf` | PATH | — | Ensembl/GENCODE GTF file. When provided, isoform structure is loaded via `ema.viz._gene_track_helpers.load_isoforms_for_gene` and drawn as a transcript model track below the PAS bars. Optional; failures are logged as warnings and do not abort the run. |

### Gene selection

| Flag | Type | Default | Description |
|---|---|---|---|
| `--top-genes` | INT | 10 | Number of top genes to auto-select from `--diff-tsv` results, ranked by volcano score via `ema.viz._gene_track_helpers.rank_top_genes`. Only used when `--diff-tsv` is provided. |
| `--gene-id` | TEXT (repeatable) | — | Explicit gene ID(s) to render regardless of differential ranking. Normalised to uppercase and stripped of whitespace. Repeated to specify multiple genes: `--gene-id ACTB --gene-id MYC`. Added to the gene list before auto-ranked genes; the union is deduplicated while preserving order. |

### Display

| Flag | Type | Default | Description |
|---|---|---|---|
| `--cluster-key` | TEXT | `leiden` | `adata.obs` column carrying cluster labels. Change this if `peakatail run` was called with `--external-clusters` or a custom labelling. |

### Output

| Flag | Type | Default | Description |
|---|---|---|---|
| `--output` / `-o` | PATH | `switch_geneview` | Output directory base name. Auto-routed inside the originating run dir when `--h5ad` is from a `peakatail_runs/` path. |
| `--plot-engine` | TEXT | `matplotlib` | Rendering engine. `matplotlib` writes PNG (and SVG with `--plot-format svg`). `plotly` writes interactive HTML. `both` writes all formats. The renderer auto-discovers engines via `ema.viz.render_all`. |
| `--plot-format` | TEXT | `all` | Restrict output formats. `png`, `svg`, `html`, or comma list. |
| `--no-plots` | FLAG | off | Disable all figure output. The gene list and manifest are still computed, but no figure files are written. |

## Gene selection logic

The command implements the following selection pipeline (source:
`ema/cli/switch_geneview.py` lines 161–218):

1. Parse `--gene-id` arguments: normalise to uppercase, strip whitespace.
2. If `--diff-tsv` is provided, call `rank_top_genes(pair_results, n=top_genes)` on the accumulated diff DataFrames. This ranks genes by volcano score across all pairs.
3. Build the final gene list as `explicit_genes + auto_genes`, deduplicated while preserving order (explicit genes appear first).
4. For each gene, call `build_gene_panel` and `render_all`. Genes not found in the AnnData or pasbed are skipped with a warning.

## Output files

All figure files are written to `<out_dir>/figures/`.

**`figures/gene_<GENE_ID>.<ext>`**

One file per gene per engine/format combination. The base name is always
`gene_<GENE_ID>` where `<GENE_ID>` is the gene identifier as normalised (uppercase).
Extensions: `.png` (matplotlib), `.svg` (matplotlib with svg format), `.html`
(plotly). Source: `ema/cli/switch_geneview.py` line 269.

**`figures/meta.json`**

Written by `ema.viz._meta.write_figures_index`. JSON manifest of all figure
files in the directory with their paths and the originating command
(`"peakatail switch geneview"`). Failures to write this file are caught and logged
as warnings (the figures are not affected).

**`figures/figures_INDEX.md`**

Markdown index of figure files in the `figures/` directory, generated alongside
`meta.json`. Suitable for inclusion in a report or MkDocs documentation site.

## How it relates to other commands

- **[`peakatail run`](run.md)** — produces `clusters.h5ad` and `pasbed.bed`. Gene annotation in `adata.var["gene_id"]` is required.
- **[`peakatail switch diff`](switch-diff.md)** — produces the `differential/*.tsv` files consumed via `--diff-tsv` for auto gene ranking.
- **[`peakatail switch length`](switch-length.md)** — produces the PDUI / proportion / entropy TSVs. `geneview` does not read them; it overlays differential results supplied with `--diff-tsv`.

## See also

- Strategy pages in [`../strategies/`](../strategies/) — details on the volcano ranking and gene-track rendering.
- Tutorial in [`../tutorials/`](../tutorials/) — step-by-step gene-track visualisation guide.
