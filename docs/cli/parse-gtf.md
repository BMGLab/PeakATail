# `peakatail parse-gtf`

Pre-parse and cache an Ensembl-style GTF so subsequent `peakatail run` / `peakatail switch length --isoform-agg per_isoform` invocations skip the expensive multi-minute GTF parse and load from disk in milliseconds.

The cache lives at `~/.cache/peakatail/gtf/<gtf-hash>/` by default. Both `peakatail run` and `peakatail switch length` look it up automatically — `parse-gtf` is just a way to warm it explicitly (e.g. as part of an environment-setup script, or before a Docker build).

!!! note "When to use it"
    - First-time setup for a fresh environment — warm the cache so the first real run doesn't pay the parse cost.
    - CI / Docker images that bake the cache into the layer so containers start instantly.
    - You changed GTFs (e.g. Ensembl version bump) and want to confirm the new file parses cleanly before the pipeline tries to use it.

!!! warning "When NOT to use it"
    - You're running a one-off pipeline once — `peakatail run` will populate the cache anyway on first invocation. Pre-warming has no benefit.

## Quick examples

Pre-parse one GTF:

```bash
uv run peakatail parse-gtf -g /data/gtfs/Homo_sapiens.GRCh38.99.gtf
```

List existing cache entries:

```bash
uv run peakatail parse-gtf --show-cache
```

Force re-parse on cache hit (e.g. after a code change):

```bash
uv run peakatail parse-gtf -g /data/gtfs/Homo_sapiens.GRCh38.99.gtf --force
```

## Full `--help` output

```text
Usage: peakatail parse-gtf [OPTIONS]

  Pre-warm the GTF cache so subsequent runs hit instantly.

Options:
  --threads INT                Max parallel workers (parse is multiprocessed).
  -v, --verbose                Increase verbosity.
  -q, --quiet                  WARNING and up only.
  --log-level TEXT             Explicit logger level.
  --no-progress                Suppress Rich progress bars.
  -c, --config FILE            YAML config.
  -g, --gtf PATH               GTF file to pre-parse (required unless --show-cache).
  --cache-dir PATH             Override the global ~/.cache/peakatail/gtf/ location.
  --force                      Re-parse even on cache hit.
  --show-cache                 List cache entries with sizes/age and exit.
  -h, --help                   Show this message and exit.
```

## Flags

| Flag | Type | Default | Description |
|---|---|---|---|
| `--gtf` / `-g` | PATH | — | GTF to pre-parse. Required unless `--show-cache` is used. |
| `--cache-dir` | PATH (dir) | `~/.cache/peakatail/gtf/` | Override the cache location (e.g. for a CI runner with limited home-dir size). |
| `--force` | flag | False | Re-parse and overwrite the cache entry even on hit. |
| `--show-cache` | flag | False | List cached GTFs with size + age, then exit. |
| `--threads` | INT | auto | Parallel workers for the parse (`parse_isoform_utrs` uses up to `min(--threads, 4)` chromosome workers). |

Plus the common logging flags (`-v`, `-q`, `--log-level`, `--no-progress`, `--config`). No `--output` flag — output is the cache directory.

## Output files

The parse writes to `<cache-dir>/<gtf-content-hash>/`:

- **`<hash>/transcripts.parquet`** — per-transcript exon + 3'UTR table
- **`<hash>/source.txt`** — original GTF path + hash for provenance
- **`<hash>/meta.json`** — parse timestamp, peakatail version, transcript count

Future invocations of `peakatail run` / `peakatail switch length` compute the same content-hash and skip directly to loading the parquet.

## How it relates to other commands

The cache is consumed transparently by:

- `peakatail run` — when annotating peaks against transcripts.
- `peakatail switch length --isoform-agg per_isoform` — when mapping PAS to isoform UTRs.
- `peakatail switch geneview --gtf <gtf>` — when rendering the isoform structure track.

No manual flag is required to use a cached entry; PeakATail computes the hash and looks it up automatically.

## See also

- [`peakatail run`](run.md)
- [`peakatail switch length`](switch-length.md) — per-isoform mode
- [`peakatail switch geneview`](switch-geneview.md) — isoform overlay
