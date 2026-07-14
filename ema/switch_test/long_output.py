"""E5: normalized long outputs for `ema switch` — findings_long / length_long.

Reshapes the per-pair `switch diff` wide TSVs (and `switch length` outputs) into
one tidy long table conforming to peakatail_contract.models.FindingRow /
LengthRow, keyed by the content-addressed ``pas_uid`` + ``canonical_cluster``,
with an explicit ``direction`` (never NA-silent — bug D8 / E5).

ID minting reproduces the frozen peakatail-contract ``ids`` formulas exactly
(kept inline so the engine has no hard dependency on the hub repo):
  * pas_uid       = chrom:pos:strand   (pos = end-1 on +, start on -)  [E1]
  * finding_uid   = arm:strategy:resolved_label:var_name:pas_uid

DIRECTION (conservative + honest): the `switch diff` wide TSV does not carry
per-PAS proximal/distal polarity, so shorten-vs-lengthen cannot be determined
from a diff row alone. We therefore emit:
  * "flat"          — not significant (qvalue >= fdr) with a finite delta
  * "undetermined"  — significant change, or missing stats, where polarity is
                      not yet computable
We NEVER emit "shorten"/"lengthen" from diff without a polarity signal (that
requires isoform-rank wiring + real-data validation — flagged for follow-up).
This satisfies the contract's exhaustive Direction enum without fabricating
biology.
"""
from __future__ import annotations

import logging
import math

import pandas as pd

from ema.datasets.pas_merge import pas_uid_of

log = logging.getLogger(__name__)

# FindingRow column order (peakatail_contract.models.FindingRow).
FINDING_LONG_COLUMNS = [
    "finding_uid", "pas_uid", "gene_id", "canonical_cluster", "comparison_cluster",
    "celltype", "strategy", "arm", "direction", "utr_class",
    "qvalue", "pvalue", "delta_proportion", "log2fc", "odds_ratio",
    "n_cells", "n_reads", "n_cells_subject", "n_cells_comparison",
    "n_reads_subject", "n_reads_comparison",
]

_SEP = ":"


def _sanitize(value) -> str:
    """Make a value safe as a separator-free finding_uid field."""
    return str(value).replace(_SEP, "_")


def finding_uid(arm: str, strategy: str, resolved_label: str, var_name: str, pas_uid: str) -> str:
    """Mint finding_uid exactly as peakatail_contract.ids.finding_uid."""
    return _SEP.join([
        _sanitize(arm), _sanitize(strategy), _sanitize(resolved_label),
        _sanitize(var_name), str(pas_uid),
    ])


def _pas_uid_for_row(chrom, start, end, strand) -> str:
    """Content-addressed pas_uid, or "" when coordinates are missing."""
    if chrom in (None, "", "nan") or strand not in ("+", "-"):
        return ""
    try:
        int(start)
        int(end)
    except (TypeError, ValueError):
        return ""
    return pas_uid_of(str(chrom), start, end, strand)


def _finite(x):
    try:
        f = float(x)
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def classify_direction(qvalue, delta_proportion, fdr: float = 0.05) -> str:
    """Conservative, non-NA direction (see module docstring)."""
    q = _finite(qvalue)
    d = _finite(delta_proportion)
    if q is None or d is None:
        return "undetermined"
    if q >= fdr:
        return "flat"
    return "undetermined"  # significant, but shorten/lengthen polarity unknown


def _canon(label, canonical_map) -> str:
    if canonical_map and label in canonical_map:
        return str(canonical_map[label])
    return str(label)


def findings_long(
    pair_dfs: dict,
    *,
    strategy: str,
    arm: str,
    canonical_map: dict | None = None,
    fdr: float = 0.05,
) -> pd.DataFrame:
    """Build a FindingRow-conformant long DataFrame from per-pair diff DataFrames.

    ``pair_dfs`` maps ``(cluster1, cluster2) -> augmented diff DataFrame`` (the
    columns ``_augment_diff_df`` produces: pas_id, gene_id, chrom, start, end,
    strand, cluster1, cluster2, + per-strategy stats). One output row per
    PAS x comparison.
    """
    rows: list[dict] = []
    for (c1, c2), df in pair_dfs.items():
        if df is None or len(df) == 0:
            continue
        canon_subj = _canon(c1, canonical_map)
        canon_comp = _canon(c2, canonical_map)
        for r in df.to_dict("records"):
            pas_id = str(r.get("pas_id", ""))
            puid = _pas_uid_for_row(r.get("chrom"), r.get("start"), r.get("end"), r.get("strand"))
            direction = classify_direction(r.get("qvalue"), r.get("delta_proportion"), fdr)
            n_reads_1 = _finite(r.get("n_reads_pas_cluster1"))
            n_reads_2 = _finite(r.get("n_reads_pas_cluster2"))
            n_reads = None
            if n_reads_1 is not None or n_reads_2 is not None:
                n_reads = int((n_reads_1 or 0) + (n_reads_2 or 0))
            rows.append({
                "finding_uid": finding_uid(arm, strategy, canon_subj, pas_id, puid),
                "pas_uid": puid,
                "gene_id": str(r.get("gene_id", "") or ""),
                "canonical_cluster": canon_subj,
                "comparison_cluster": canon_comp,
                "celltype": None,  # gated on A1 (ema celltype); never a raw dir name
                "strategy": strategy,
                "arm": arm,
                "direction": direction,
                "utr_class": None,
                "qvalue": _finite(r.get("qvalue")),
                "pvalue": _finite(r.get("pvalue")),
                "delta_proportion": _finite(r.get("delta_proportion")),
                "log2fc": _finite(r.get("log2fc")),
                "odds_ratio": _finite(r.get("odds_ratio")),
                "n_cells": int(_finite(r.get("n_cells")) or 0),
                "n_reads": n_reads if n_reads is not None else 0,
                "n_cells_subject": _as_int(r.get("n_cells_cluster1")),
                "n_cells_comparison": _as_int(r.get("n_cells_cluster2")),
                "n_reads_subject": _as_int(r.get("n_reads_pas_cluster1")),
                "n_reads_comparison": _as_int(r.get("n_reads_pas_cluster2")),
            })
    return pd.DataFrame(rows, columns=FINDING_LONG_COLUMNS)


def _as_int(x):
    f = _finite(x)
    return int(f) if f is not None else None


# LengthRow column order (peakatail_contract.models.LengthRow).
LENGTH_LONG_COLUMNS = [
    "strategy", "gene_id", "transcript_id", "cell_uid", "canonical_cluster",
    "value", "pas_uid", "rank", "direction",
]


def cell_uid(dataset_id: str, barcode: str) -> str:
    """Mint cell_uid as peakatail_contract.ids.cell_uid (``dataset_id:barcode``).

    If ``barcode`` already looks namespaced (contains the separator), it is
    returned unchanged so we never double-namespace.
    """
    b = str(barcode)
    if _SEP in b:
        return b
    return f"{_sanitize(dataset_id)}{_SEP}{b}"


def length_long(
    df: pd.DataFrame,
    *,
    strategy: str,
    value_col: str,
    dataset_id: str = "",
    canonical_map: dict | None = None,
    cell_col: str = "cell",
    cluster_col: str = "cluster",
    gene_col: str = "gene_id",
    transcript_col: str | None = "transcript_id",
    pas_col: str | None = None,
    rank_col: str | None = None,
    pas_uid_map: dict | None = None,
) -> pd.DataFrame:
    """Reshape a `switch length` output into a LengthRow-conformant long table.

    ``strategy`` in {classic, proportion, shannon}; ``value_col`` is the native
    value column (pdui / proportion / entropy). ``pas_col``/``rank_col`` are only
    passed for the per-PAS ``proportion`` strategy; ``pas_uid_map`` maps a pas_id
    to its content-addressed pas_uid (from ``unified/pas_uid.tsv`` / the ledger).
    The ``_gene_`` isoform sentinel is normalized to ``None`` (contract: "no
    specific isoform" is an explicit null, not a fabricated transcript).

    NOTE: not auto-wired into run_length yet — the per-strategy column names and
    cell→dataset namespacing need confirmation against real length outputs. This
    is the tested LengthRow producer for hub-team; wiring is flagged as follow-up.
    """
    rows: list[dict] = []
    for r in df.to_dict("records"):
        transcript = None
        if transcript_col and transcript_col in r:
            t = r.get(transcript_col)
            if t not in (None, "", "_gene_", "nan"):
                transcript = str(t)
        pas_id = str(r.get(pas_col)) if (pas_col and pas_col in r) else None
        puid = None
        if pas_id is not None and pas_uid_map:
            puid = pas_uid_map.get(pas_id)
        rank = _as_int(r.get(rank_col)) if (rank_col and rank_col in r) else None
        rows.append({
            "strategy": strategy,
            "gene_id": str(r.get(gene_col, "") or ""),
            "transcript_id": transcript,
            "cell_uid": cell_uid(dataset_id, r.get(cell_col, "")),
            "canonical_cluster": _canon(r.get(cluster_col, ""), canonical_map),
            "value": _finite(r.get(value_col)),
            "pas_uid": puid,
            "rank": rank,
            # Direction only for proportion (per-PAS) rows; left None otherwise
            # (contract: classic/shannon are gene-level trends without a per-row
            # direction). Polarity itself is a follow-up (needs isoform rank).
            "direction": None,
        })
    return pd.DataFrame(rows, columns=LENGTH_LONG_COLUMNS)


def write_long_table(df: pd.DataFrame, base_path: str) -> str:
    """Write ``df`` as parquet (contract format); fall back to TSV if no parquet
    engine is installed. Returns the actual path written.

    The engine env may lack pyarrow/fastparquet; in that case a ``.tsv`` sibling
    is written and a warning is logged so the data is never lost. Install pyarrow
    in the pipeline env for true ``.parquet`` output.
    """
    parquet_path = base_path if base_path.endswith(".parquet") else base_path + ".parquet"
    try:
        df.to_parquet(parquet_path, index=False)
        return parquet_path
    except Exception as exc:  # missing pyarrow/fastparquet, etc.
        tsv_path = parquet_path[: -len(".parquet")] + ".tsv"
        df.to_csv(tsv_path, sep="\t", index=False)
        log.warning(
            "parquet engine unavailable (%s); wrote %s instead of %s. Install "
            "pyarrow in the pipeline env for contract parquet output.",
            exc, tsv_path, parquet_path,
        )
        return tsv_path
