#!/usr/bin/env python
"""Build the experiment matrix: every run executed, with the parameters it used.

The reports describe results per axis but never state, in one place, which runs
were actually executed and under which resolved parameters.  This reconstructs
that table from the runs themselves rather than from the sweep definition, so it
records what the pipeline *did*, not what it was asked to do -- the two can
differ when a CLI default overrides a grid cell.

Sources, in precedence order:

1. ``run_manifest.json`` -> ``resolved_config`` for each run directory.  This is
   authoritative: it is written by the run and contains the fully resolved
   ``args``, ``filters`` and ``variables``.
2. ``experiments/laughney/*.tsv`` -- the grid definitions, used only to recover
   the human-readable axis labels and notes, and to flag any grid row that has no
   corresponding run (i.e. a planned experiment that did not execute).

Emits ``experiment_matrix.csv`` (one row per run, all resolved parameters) plus
``experiment_matrix_summary.md`` (the compact table for the report).

Usage
-----
    python reports/build_experiment_matrix.py [--manifests DIR] [--grids DIR] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
MANIFESTS = HERE / "cumulative_analysis" / "manifests"
GRIDS = HERE.parent / "experiments" / "laughney"
OUT = HERE / "cumulative_analysis"

# Parameters worth surfacing per run, grouped by the axis they belong to.
PEAK_ARGS = [
    "strategy", "atlas_mode", "atlas_distance", "ip_filter", "ip_filter_mode",
    "annot_filter", "annot_filter_mode", "default_threshold", "lambda_window",
    "lambda_method", "lambda_fold_change", "max_pas", "smoothing_window",
    "min_prominence", "floor_threshold", "pas_gap", "min_pas_spacing",
    "min_pas_prominence", "merge_len", "dynamic_threshold",
]
ANNOT_ARGS = ["max_gene_distance", "utr_multiplier", "include_extended"]
CLUSTER_ARGS = ["cluster_method", "resolution", "n_neighbors", "n_comps"]
FILTER_KEYS = ["min_read", "min_cells", "min_genes", "min_pas_per_cell"]

AXIS = {  # run -> which axis the sweep was varying
    "grid_lg_annotate": "peak-calling strategy (baseline)",
    "grid_lp_annotate": "peak-calling strategy",
    "grid_si_annotate": "peak-calling strategy",
    "grid_lg_ip_filter": "internal-priming mode",
    "grid_lg_ip_off": "internal-priming mode",
    "B1_cohort_full": "full cohort (downstream analyses)",
}


def load_manifests(directory: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(directory.glob("*.json")):
        try:
            man = json.loads(path.read_text())
        except json.JSONDecodeError as exc:
            print(f"  WARN unreadable {path.name}: {exc}")
            continue
        rc = man.get("resolved_config", {})
        args = rc.get("args", {}) or {}
        filters = rc.get("filters", {}) or {}
        dirs = rc.get("directories", {}) or {}
        datasets = dirs.get("datasets") or man.get("datasets") or []

        row = {
            "run": path.stem,
            "axis": AXIS.get(path.stem, "trim / clustering re-annotation"),
            "run_id": man.get("run_id"),
            "timestamp": man.get("timestamp"),
            "n_datasets": len(datasets) if isinstance(datasets, list) else None,
            "n_bams": sum(len(d.get("bams", [])) for d in datasets
                          if isinstance(d, dict)) if isinstance(datasets, list) else None,
            "merge_strategy": next((d.get("merge_strategy") for d in datasets
                                    if isinstance(d, dict)), None) if datasets else None,
            "gtf": Path(dirs["gtf_dir"]).name if dirs.get("gtf_dir") else None,
            "atlas": Path(dirs["atlas"]).name if dirs.get("atlas") else None,
        }
        for k in PEAK_ARGS + ANNOT_ARGS + CLUSTER_ARGS:
            if k in args:
                row[k] = args[k]
        # `ip_filter_mode` stays "annotate" even when the filter is disabled, so
        # the mode alone does not identify the arm; `ip_filter` is the real switch.
        if "ip_filter" in args:
            row["ip_mode_effective"] = (
                args.get("ip_filter_mode") if args.get("ip_filter") else "off")
        if "annot_filter" in args:
            row["annot_mode_effective"] = (
                args.get("annot_filter_mode") if args.get("annot_filter") else "off")
        row["param_source"] = "run_manifest"
        for k in FILTER_KEYS:
            row[k] = filters.get(k)
        # entity counts, when the run recorded them
        for k, v in (man.get("entity_counts") or {}).items():
            row[f"count_{k}"] = v
        rows.append(row)
    return pd.DataFrame(rows)


def load_grids(directory: Path) -> dict[str, pd.DataFrame]:
    out = {}
    for name in ("peakcall_grid", "trim_cluster_grid", "sweep_grid", "contrasts", "stage_groups"):
        p = directory / f"{name}.tsv"
        if p.exists():
            out[name] = pd.read_csv(p, sep="\t")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifests", type=Path, default=MANIFESTS)
    ap.add_argument("--grids", type=Path, default=GRIDS)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    df = load_manifests(args.manifests)
    if df.empty:
        print(f"no manifests under {args.manifests}")
        return 1
    grids = load_grids(args.grids)
    print(f"[manifests] {len(df)} runs; [grids] {', '.join(grids) or 'none'}")

    # A2_*/A3_* are re-annotation branches of a base run: their run_manifest.json
    # carries the BASE run's args (max_gene_distance 5000, resolution 1.0, ...),
    # not the branch's. Taking the manifest at face value would report every trim
    # branch as identical. For those runs the grid TSV is authoritative, so the
    # varied columns are overwritten from it and `param_source` records the swap.
    tc = grids.get("trim_cluster_grid")
    if tc is not None:
        tc = tc.set_index("branch_name")
        varied = ["max_gene_distance", "utr_multiplier", "include_extended",
                  "cluster_method", "resolution", "n_neighbors"]
        for i, r in df.iterrows():
            branch = r["run"].replace("reannotate_", "")
            if branch not in tc.index:
                continue
            for col in varied:
                if col in tc.columns:
                    v = tc.loc[branch, col]
                    if isinstance(v, str) and v.lower() in {"true", "false"}:
                        v = v.lower() == "true"
                    df.at[i, col] = v
            df.at[i, "param_source"] = "trim_cluster_grid.tsv (re-annotation branch)"
        n_swap = int((df["param_source"] != "run_manifest").sum())
        print(f"[provenance] {n_swap} runs took varied params from the grid, not the manifest")

    df = df.sort_values(["axis", "run"]).reset_index(drop=True)
    df.to_csv(args.out / "experiment_matrix.csv", index=False)
    print(f"wrote {args.out / 'experiment_matrix.csv'}  ({df.shape[0]} runs x {df.shape[1]} fields)")

    # cross-check: grid rows that never produced a run
    missing: list[str] = []
    for gname, key in (("peakcall_grid", "run_name"), ("trim_cluster_grid", "branch_name"),
                       ("sweep_grid", "branch_name")):
        g = grids.get(gname)
        if g is None or key not in g:
            continue
        executed = {r.split("_", 1)[-1] if r.startswith(("grid_", "reannotate_")) else r
                    for r in df["run"]}
        for name in g[key]:
            if name not in executed:
                missing.append(f"{gname}:{name}")
    if missing:
        print(f"[gap] {len(missing)} grid rows with no run: {', '.join(missing[:12])}"
              + (" ..." if len(missing) > 12 else ""))

    # ---------------- compact markdown for the report ----------------
    def fmt(v) -> str:
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return "—"
        if isinstance(v, bool):
            return "yes" if v else "no"
        return str(v)

    lines: list[str] = []
    lines.append("### A. Peak-calling grid (5 runs, full 17-dataset cohort input)\n")
    lines.append("| run | strategy | atlas mode | ip mode | max PAS/peak | λ window | λ method | "
                 "fold-change | PAS gap | min prominence |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for _, r in df[df["run"].str.startswith("grid_")].iterrows():
        lines.append(
            f"| `{r['run'].replace('grid_', '')}` | {fmt(r.get('strategy'))} | "
            f"{fmt(r.get('atlas_mode'))} | {fmt(r.get('ip_mode_effective'))} | "
            f"{fmt(r.get('max_pas'))} | {fmt(r.get('lambda_window'))} | "
            f"{fmt(r.get('lambda_method'))} | {fmt(r.get('lambda_fold_change'))} | "
            f"{fmt(r.get('pas_gap'))} | {fmt(r.get('min_pas_prominence'))} |")

    trim = df[df["run"].str.contains("A2_trim")]
    if not trim.empty:
        lines.append("\n### B. Trim / annotation-window branches\n")
        lines.append("| branch | max_gene_distance | utr_multiplier | include_extended |")
        lines.append("|---|---|---|---|")
        for _, r in trim.iterrows():
            lines.append(
                f"| `{r['run'].replace('reannotate_', '')}` | {fmt(r.get('max_gene_distance'))} | "
                f"{fmt(r.get('utr_multiplier'))} | {fmt(r.get('include_extended'))} |")

    clus = df[df["run"].str.contains("A3_")]
    if not clus.empty:
        lines.append("\n### C. Clustering branches\n")
        lines.append("| branch | method | resolution | n_neighbors |")
        lines.append("|---|---|---|---|")
        for _, r in clus.iterrows():
            lines.append(
                f"| `{r['run'].replace('reannotate_', '')}` | {fmt(r.get('cluster_method'))} | "
                f"{fmt(r.get('resolution'))} | {fmt(r.get('n_neighbors'))} |")

    lines.append("\n### D. Cell-level filters (identical across all runs unless noted)\n")
    lines.append("| filter | value |")
    lines.append("|---|---|")
    for k in FILTER_KEYS:
        vals = df[k].dropna().unique()
        lines.append(f"| `{k}` | {fmt(vals[0]) if len(vals) == 1 else ' / '.join(map(fmt, vals))} |")

    (args.out / "experiment_matrix_summary.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out / 'experiment_matrix_summary.md'}")
    print("\n".join(lines[:12]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
