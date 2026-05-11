# `ema merge`

Merge multiple BAM files into a single sorted + indexed BAM. Thin convenience wrapper around `samtools merge` + `samtools sort` + `samtools index` with PeakATail's logging plumbing — useful as a pre-step when a sample is split across lanes / sequencing runs.

!!! note "When to use it"
    - You have technical replicates of the same library that need to be merged before peak calling.
    - You want one combined BAM for `ema run`'s `merge_strategy: before` flow but you'd rather do the merge separately so you can reuse the merged BAM across runs.

!!! warning "When NOT to use it"
    - **Different libraries / different cell-barcode spaces** — never merge BAMs whose CB-tag space differs (e.g. two distinct 10x runs); use multi-dataset mode in `ema run` instead.
    - **You're already running `ema run` with `merge_strategy: before`** — that path does the merge internally; running `ema merge` first would be redundant.

## Quick example

```bash
uv run ema merge \
    -i sample1_lane1.bam \
    -i sample1_lane2.bam \
    -o sample1_merged.bam \
    --threads 4
```

Produces:

- `sample1_merged.bam` — coord-sorted merged BAM
- `sample1_merged.bam.bai` — BAM index (samtools auto-builds)

## Full `--help` output

```text
Usage: ema merge [OPTIONS]

  Merge multiple BAM files into one sorted+indexed BAM.

Options:
  --threads INT                Max parallel workers passed to samtools merge / sort.
  -v, --verbose                Increase verbosity. -v = DEBUG for ema.*.
  -q, --quiet                  WARNING and up only.
  --log-level TEXT             Explicit logger level.
  --no-log-file                Don't write peakatail_<ts>.log next to the output.
  --no-progress                Suppress Rich progress bars.
  -c, --config FILE            YAML config (rarely used for merge).
  --output, -o PATH            Output BAM path (will be sorted+indexed). Default: merged.bam
  --bam-files, -i PATH         BAM files to merge. Repeat -i for each. [required]
  -h, --help                   Show this message and exit.
```

## Flags

| Flag | Type | Default | Description |
|---|---|---|---|
| `--bam-files` / `-i` | PATH (repeatable, required) | — | Input BAMs to merge. Pass `-i` once per file. |
| `--output` / `-o` | PATH (file) | `merged.bam` | Output BAM path. Coord-sorted + indexed automatically. |
| `--threads` | INT | auto (≥4) | Forwarded to `samtools merge --threads`. |

Plus the common logging flags from [`common_options`](run.md#common-flags) (`-v`, `-q`, `--log-level`, `--no-log-file`, `--no-progress`, `--config`).

## Output files

- **`<output>.bam`** — coord-sorted merged BAM
- **`<output>.bam.bai`** — samtools-generated index

## How it relates to other commands

`ema merge` is a one-shot helper; its output BAM can be fed straight into `ema run` as a `datasets[].bams[]` entry. See [`ema run`](run.md) for the main pipeline.

## See also

- [`ema run`](run.md) — main pipeline entry point (uses `merge_strategy: before` internally when configured).
- [Quickstart tutorial](../tutorials/02-quickstart.md)
