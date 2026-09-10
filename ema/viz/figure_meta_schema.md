# Figure Provenance Sidecar Schema

Every figure written by a `VizStrategy` produces a companion
`<stem>.meta.json` next to it (e.g. `umap_ds1.meta.json` beside
`umap_ds1.png`).

## Fixed keys (all figures)

| Key | Type | Description |
|---|---|---|
| `figure_name` | string | Filename stem (no extension) |
| `viz_strategy` | string | Strategy class `name` attribute |
| `generated_at` | ISO 8601 UTC | Timestamp of render |
| `peakatail_version` | string | Package version or `"unknown"` |

## Plot-type-specific keys

### `umap_*` (matplotlib / plotly / scanpy)
| Key | Type |
|---|---|
| `dataset_id` | string |
| `n_cells` | int |
| `n_genes_used` | int |
| `color_key` | string or null (`"leiden"` or `null`) |

### `cluster_sizes_*`
| Key | Type |
|---|---|
| `dataset_id` | string |
| `cluster_key` | string (always `"leiden"`) |
| `n_clusters` | int |
| `n_observations` | int |

### `pdui_distribution_*` (matplotlib / plotly / scanpy)
| Key | Type |
|---|---|
| `cluster_key` | string |
| `score_key` | string (e.g. `"mean_pdui"`) |
| `n_clusters` | int |
| `n_observations` | int |

### `length_shifts_*`
| Key | Type |
|---|---|
| `cluster_key` | string |
| `n_genes_shown` | int (capped at 50) |
| `n_cluster_pairs` | int |
| `cluster_pairs` | list[string] |

### `diff_agreement_*`
| Key | Type |
|---|---|
| `sig_set_names` | list[string] |
| `sig_set_sizes` | dict[string, int] |

### `volcano_*`
| Key | Type |
|---|---|
| `fdr` | float |
| `log2fc_thresh` | float |
| `n_tested` | int |
| `n_significant` | int |
| `n_up` | int |
| `n_down` | int |

### `peak_qc_*`
| Key | Type |
|---|---|
| `n_chromosomes` | int |
| `n_peaks_total` | int |
| `n_cells` | int |

### `pas_overlap_*`
| Key | Type |
|---|---|
| `n_datasets` | int |
| `dataset_ids` | list[string] |
| `dataset_sizes` | dict[string, int] |

### `atlas_snap_diag_*`
| Key | Type |
|---|---|
| `n_snapped` | int |
| `n_unsnapped` | int |
| `snap_rate` | float or null |
| `n_snap_distances` | int |

### `cluster_match_sankey_*`
| Key | Type |
|---|---|
| `n_datasets` | int |
| `dataset_ids` | list[string] |
| `n_canonical_clusters` | int |
| `n_match_rows` | int |

### `match_confidence_*`
| Key | Type |
|---|---|
| `n_source_rows` | int |
| `n_canonical_clusters` | int |
| `datasets` | list[string] |

### `tile_timing_*`
| Key | Type |
|---|---|
| `n_datasets` | int |
| `dataset_ids` | list[string] |
| `n_tiles_total` | int |

### `resource_timeline_*`
| Key | Type |
|---|---|
| `n_samples` | int |
| `elapsed_range_s` | [float, float] |
| `n_stage_annotations` | int |
