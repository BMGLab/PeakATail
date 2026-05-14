# Post-Detection PAS-Summit Merger

A strategy-agnostic merge step that runs **after** every
`strategy.find_pas(peak)` call, collapsing adjacent PAS regions that are too
close together or share too shallow a valley to be biologically distinct.

## Why this step exists

The multi-PAS [peak-calling strategies](peak-calling.md)
(`lambda_gradient`, `sierra_iterative`) produce one or more summit positions
per peak, then split the peak's genomic span among them by midpoint cut. This
split has two failure modes that no upstream knob catches:

1. **Twin summits.** On low-coverage regions, the smoothed gradient can show
   multiple shallow oscillations within a few base pairs. Each one clears the
   per-summit prominence threshold individually and gets emitted, producing
   pairs of PAS regions that differ by 1 bp at the partition boundary. No
   single read can span across both — they are not biologically
   distinguishable.
2. **Shoulder peaks.** Two summits 20–80 bp apart with only a shallow dip
   between them are emitted as separate PAS even when the valley sits at the
   local Poisson background. These are one peak with a noisy double-top.

Neither `--merge-len` (read-cluster grouping, upstream of any strategy) nor
`--pas-gap` (cross-dataset BED merging, downstream of the full pipeline)
addresses within-peak summit spacing. The post-detection merger fills that
gap.

## Two-tier design

The merger lives in
[`ema/countmatrix/peackcalling.py`](https://github.com/BMGLab/PeakATail/blob/develop/ema/countmatrix/peackcalling.py)
and is wired into all three live execution paths (monolithic, pipeline,
tile). Strategy-agnostic by construction — it operates on the
`[(start, end), …]` contract that every strategy already returns.

### Tier 1 — Distance gate

```
if next_region.start - prev_region.end < min_pas_spacing:
    merge them
```

`min_pas_spacing` defaults to **`-1`** which triggers auto-detection of the
median aligned read length from the BAM header (first 1,000 mapped reads).
Two summits closer than one aligned read length cannot be resolved by any
short-read library — a single read straddles both, so they must be the same
PAS by construction.

Typical auto-detected values:

| Library | Auto-detected `min_pas_spacing` |
|---|---|
| 10x Chromium v2 | ~98 bp |
| 10x Chromium v3 | ~91 bp |
| Smart-seq2 | ~150 bp |
| Long-read (Iso-Seq) | very large — Tier 2 dominates |

Set `--min-pas-spacing 0` to disable Tier 1 (e.g. when you trust the
strategy's own distance handling).

### Tier 2 — Lambda-anchored valley

If the gap passes Tier 1, the valley between the two regions is compared
against a threshold supplied by the **strategy itself**:

```python
class PeakFinderStrategy:
    def valley_threshold(self, peak, user_min_pas_prominence):
        return user_min_pas_prominence   # base: static fallback

class LambdaGradientStrategy(PeakFinderStrategy):
    def valley_threshold(self, peak, _user):
        return compute_lambda(           # override: dynamic, per-peak
            [h for _, h in peak.peak_list],
            method=self.lambda_method,
            floor_fraction=self.floor_fraction,
        )
```

- **Lambda strategies** (`lambda_poisson`, `lambda_gradient`) override to
  return their own `compute_lambda(heights)` — the same background estimate
  they use for the region significance gate. Fully dynamic; the user's
  `--min-pas-prominence` is **ignored** for these strategies.
- **Non-lambda strategies** (`original`, `sierra_iterative`) fall through to
  the base implementation which returns `min_pas_prominence` unchanged. The
  default value is `5.0` coverage units.

The pair is merged when `valley_height < threshold`. Set
`--min-pas-prominence -1` to disable Tier 2 regardless of strategy.

## Reads in the merge gap are absorbed automatically

The pipeline's `peak.cb_positions` dictionary is populated densely across
the entire parent peak — every read end position contributes, not only
summit positions. The next pipeline call after merging,
`strategy.get_cb_dict_for_pas(peak, merged_start, merged_end)`, routes
through `reconstruct_cb_dict` which slices `cb_positions` by **inclusive
range**. Reads that fell in the gap between the two pre-merge regions are
therefore added to the survivor's barcode counts without any explicit
summing logic. The merger only needs to widen the genomic span; read
absorption is a free consequence of how reconstruction already works.

## Configuration

| CLI flag | YAML key | Default | What 0/disabled means |
|---|---|---|---|
| `--min-pas-spacing N` | `min_pas_spacing: N` | `-1` (auto = median read length) | `0` disables Tier 1 |
| `--min-pas-prominence X` | `min_pas_prominence: X` | `5.0` | `-1` (or any negative) disables Tier 2 |

Combinations:

```bash
# Defaults: both tiers active (recommended)
ema run --config example.yaml

# Distance only — for strategies that already have good valley logic
ema run --config example.yaml --min-pas-prominence -1

# Prominence only — for libraries where you don't trust the read-length
# auto-detect, but still want lambda-driven valley collapse
ema run --config example.yaml --min-pas-spacing 0

# Baseline (no merge at all) — exact pre-merger behaviour, for reproducibility
ema run --config example.yaml --min-pas-spacing 0 --min-pas-prominence -1
```

## Wiring across execution paths

The merger runs in all three live peak-calling code paths:

- **Monolithic** (`peackcalling.py`, default sequential): 5 `find_pas` call
  sites, each followed by a one-line `merge_close_or_low_prominence(…)`.
- **Pipeline** (`peak_pipeline.py`, `--pipeline`): merger applied inside the
  writer subprocess after every `find_pas`. Parameters passed at spawn time.
- **Tile** (`tile_runner.py`, `--tiles`): each tile worker calls the
  monolithic `peak_calling()` and inherits the merger automatically. Sentinel
  resolution happens once at pool startup so tile workers don't re-scan the
  BAM header.

Equivalence between the monolithic and pipeline paths is verified on the
chr22 test BAM with `lambda_gradient`: both produce bit-identical 252-PAS
output.

## Effect on the reference v9 dataset

On SRR8325947 with `lambda_gradient` and default merger settings:

| Metric | Before merger | After merger | Δ |
|---|---|---|---|
| Total PAS | 15,948 | 15,568 | −380 |
| Same-strand pairs < 30 bp | 371 | **0** | −371 |
| Filtered cells | 752 | 752 | 0 |
| Leiden clusters | 11 | 11 | 0 |
| Genes (classic length) | 1,925 | 1,898 | −27 |
| Genes (proportion length) | 5,339 | 5,407 | +68 |
| Genes (shannon entropy) | 5,339 | 5,407 | +68 |

Counter-intuitively, the gene count for proportion/shannon length
**increases** — when twin PAS whose narrow split-boundaries straddled a gene
boundary are merged into one wider region, the merged region's overlap with
a single gene becomes unambiguous and the gene becomes eligible for
length-shift analysis.

## See also

- [Peak-Calling Strategies](peak-calling.md) — how PAS summits are produced upstream
- [`ema run`](../cli/run.md) — `--min-pas-spacing` and `--min-pas-prominence` flags
- [PeakATail Technical Report](https://github.com/BMGLab/PeakATail/blob/develop/reports/PeakATail_Technical_Report.pdf) — Section 5 covers the merger in depth
