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

# FindingRow column order (peakatail_contract.models.FindingRow) + one extra
# informational column (direction_basis) the contract ignores on validate().
FINDING_LONG_COLUMNS = [
    "finding_uid", "pas_uid", "gene_id", "canonical_cluster", "comparison_cluster",
    "celltype", "strategy", "arm", "direction", "utr_class",
    "qvalue", "pvalue", "delta_proportion", "log2fc", "odds_ratio",
    "n_cells", "n_reads", "n_cells_subject", "n_cells_comparison",
    "n_reads_subject", "n_reads_comparison", "direction_basis",
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


def _summit_pos(start, end, strand):
    """Strand-aware 3' summit position (contract ids.pas_summit_pos)."""
    try:
        s, e = int(start), int(end)
    except (TypeError, ValueError):
        return None
    if strand == "+":
        return e - 1
    if strand == "-":
        return s
    return None


def structural_direction_by_gene(rows, fdr: float = 0.05, eps: float = 1e-9) -> dict:
    """Deterministic, geometry-derived 3'UTR direction per gene (D8).

    Rank each gene's PAS by 3' position along the UTR (strand-aware): the
    *distal* PAS is the one farthest in the transcription direction — largest
    coordinate on ``+``, smallest on ``-``. The gene's direction is the sign of
    the change in **distal** PAS usage between the two compared groups:
    Δdistal > 0 → ``lengthen``, < 0 → ``shorten``, ≈ 0 or not significant →
    ``flat``, and ``undetermined`` when the geometry is unresolvable (single-PAS
    gene, missing coordinates, or no usage delta).

    Δdistal is taken from the distal PAS's ``delta_proportion`` (proportion delta
    of subject − comparison group) when present, else its ``log2fc``. This is
    definitional geometry, NOT fitted biology — mark rows ``direction_basis=
    "structural"``; the *magnitude/significance calibration* still wants a real
    no-atlas re-run.

    Returns ``{gene_id: direction}``.
    """
    by_gene: dict[str, list[dict]] = {}
    for r in rows:
        gid = str(r.get("gene_id", "") or "")
        if not gid:
            continue
        by_gene.setdefault(gid, []).append(r)

    out: dict[str, str] = {}
    for gid, grows in by_gene.items():
        # PAS with resolvable geometry.
        located = []
        for r in grows:
            pos = _summit_pos(r.get("start"), r.get("end"), r.get("strand"))
            if pos is not None and r.get("strand") in ("+", "-"):
                located.append((pos, r))
        if len(located) < 2:
            # Single-PAS / unlocated gene has no distal-vs-proximal contrast — no
            # structural call is possible; omit it so the caller falls back to the
            # significance-based classifier for these rows.
            continue
        strand = located[0][1].get("strand")
        # Distal = farthest in transcription direction.
        distal = (max if strand == "+" else min)(located, key=lambda t: t[0])[1]
        d = _finite(distal.get("delta_proportion"))
        if d is None:
            d = _finite(distal.get("log2fc"))
        q = _finite(distal.get("qvalue"))
        if d is None:
            out[gid] = "undetermined"
        elif q is not None and q >= fdr:
            out[gid] = "flat"
        elif abs(d) <= eps:
            out[gid] = "flat"
        elif d > 0:
            out[gid] = "lengthen"   # distal usage up → longer 3'UTR
        else:
            out[gid] = "shorten"    # distal usage down → shorter 3'UTR
    return out


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
        recs = df.to_dict("records")
        # D8: deterministic per-gene structural direction for THIS comparison.
        gene_dir = structural_direction_by_gene(recs, fdr=fdr)
        for r in recs:
            pas_id = str(r.get("pas_id", ""))
            puid = _pas_uid_for_row(r.get("chrom"), r.get("start"), r.get("end"), r.get("strand"))
            gid = str(r.get("gene_id", "") or "")
            # Structural direction (geometry) where the gene had ≥2 located PAS;
            # else fall back to the conservative significance-only classifier.
            _sdir = gene_dir.get(gid)
            if _sdir is not None:
                direction = _sdir
                direction_basis = "structural"
            else:
                direction = classify_direction(
                    r.get("qvalue"), r.get("delta_proportion"), fdr
                )
                direction_basis = "significance"
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
                # How the direction was derived: "structural" (geometry, ≥2
                # located PAS) or "significance" (single-PAS fallback). Geometry
                # polarity is not yet calibrated on a real no-atlas re-run.
                "direction_basis": direction_basis,
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
