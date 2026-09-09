# `peakatail wizard` (and bare `ema`)

An interactive `questionary`-based wizard that walks you through assembling a config for any of PeakATail's main commands, saves it to `.peakatail/last_run.yaml`, and dispatches into the same code path as the explicit CLI invocations.

`ema` with **no subcommand** also drops you into this wizard — handy for first-time users who haven't memorised the flag set yet.

!!! note "When to use it"
    - First exposure to PeakATail — you don't know the flag names yet.
    - Quickly assembling a config to save and re-use as `--config last_run.yaml`.
    - Teaching workshops — the prompts double as a discoverable parameter inventory.

!!! warning "When NOT to use it"
    - Scripted / unattended runs (CI, Docker, batch jobs) — `questionary` requires a TTY. Use the explicit subcommands with `--config` instead.

## Quick example

```bash
# Bare ema also routes here
uv run ema

# Or explicitly
uv run peakatail wizard
```

The wizard first asks **which mode** you want:

```
? What do you want to do?
  > Run full pipeline (peak calling -> cluster)
    Differential APA test (peakatail switch diff)
    3'UTR shortening / lengthening (peakatail switch length)
    Cross-dataset cluster matching (peakatail switch match)
    Merge BAMs (peakatail merge)
    Pre-warm GTF cache (peakatail parse-gtf)
```

Then sequences mode-specific prompts (paths to BAMs, GTF, strategy, FDR, output dir, etc.) and finally asks:

```
? Save this config to .peakatail/last_run.yaml? (Y/n)
? Proceed to run now? (Y/n)
```

If you decline the run, the YAML is still on disk so you can launch it later via:

```bash
uv run peakatail run --config .peakatail/last_run.yaml
```

## Full `--help` output

```text
Usage: peakatail wizard [OPTIONS]

  Launch the interactive wizard explicitly.

Options:
  -h, --help  Show this message and exit.
```

No flags — the wizard takes everything via prompts.

## Modes dispatched

| Wizard option | Dispatches to |
|---|---|
| Run full pipeline | [`peakatail run`](run.md) |
| Differential APA test | [`peakatail switch diff`](switch-diff.md) |
| 3'UTR shortening/lengthening | [`peakatail switch length`](switch-length.md) |
| Cross-dataset cluster matching | [`peakatail switch match`](switch-match.md) |
| Merge BAMs | [`peakatail merge`](merge.md) |
| Pre-warm GTF cache | [`peakatail parse-gtf`](parse-gtf.md) |

## Output files

`.peakatail/last_run.yaml` — the configuration you assembled (only written if you answer "yes" to the save prompt). Re-usable via any of the matching subcommands' `--config` flag.

After dispatch, the chosen subcommand writes its own outputs to `peakatail_runs/<name>_<ts>/` (see each subcommand's page).

## See also

- All the subcommands the wizard dispatches into — linked in the table above.
- [Quickstart tutorial](../tutorials/02-quickstart.md) — covers the same ground in static command form.
