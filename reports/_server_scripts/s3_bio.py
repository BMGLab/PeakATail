#!/usr/bin/env python
"""Biological layer: recurrent vs cell-type-private switches, and whether the
lost distal 3'UTR segment carries AU-rich / miRNA-seed regulatory elements.

Four outputs, all written under ``analysis_extra/``:

1. ``gene_symbols.tsv``            -- ENSG -> HGNC symbol from the exact GTF the
   run used (Ensembl GRCh38.99), so no symbol in the report is guessed.
2. ``celltype_gene_trend_matrix.csv`` -- gene x cell-type stage slope, from the
   small per-gene trend TSVs.  Drives the recurrent-vs-private decomposition.
3. ``gene_interval_map.tsv``       -- per gene, the proximal and distal PAS used
   by the classic PDUI layer, and the "lost distal segment": the stretch of
   3'UTR present only in the long isoform, which shortening removes.
4. ``lost_distal_element_scan.csv`` -- sequence composition of that lost segment
   (and of a length-matched proximal control) for every gene: AU-rich element
   counts and miRNA seed-match counts.

The ARE measure is sequence-defined and needs no external database: the class-II
ARE is the ``ATTTA`` pentamer, and the high-confidence functional form is the
``WTTTATTTAW`` nonamer (Zubiaga 1995).  The miRNA panel is small, explicit, and
computed from the mature sequences written below -- it is a seed-match scan,
NOT TargetScan: it does not model site conservation, accessibility, or context,
so it bounds an upper limit on candidate sites rather than predicting repression.
"""
from __future__ import annotations

import gzip
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

SWEEP = Path("/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed")
B1 = SWEEP / "runs" / "B1_cohort_full"
TREND = B1 / "B3_switch" / "trend"
LENGTH = B1 / "B3_switch" / "length"
GTF = Path("/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf")
FASTA = Path("/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
OUT = SWEEP / "analysis_extra"
OUT.mkdir(parents=True, exist_ok=True)

STAGE_ORDER = ["Normal", "StageI", "IVprimary", "Met"]

# ---------------------------------------------------------------------------
# 1. ENSG -> symbol from the run's own GTF
# ---------------------------------------------------------------------------
sym_path = OUT / "gene_symbols.tsv"
if not sym_path.exists():
    pat_id = re.compile(r'gene_id "([^"]+)"')
    pat_nm = re.compile(r'gene_name "([^"]+)"')
    seen: dict[str, str] = {}
    opener = gzip.open if str(GTF).endswith(".gz") else open
    with opener(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            mi, mn = pat_id.search(f[8]), pat_nm.search(f[8])
            if mi:
                seen[mi.group(1)] = mn.group(1) if mn else ""
    pd.Series(seen, name="symbol").rename_axis("gene_id").to_csv(sym_path, sep="\t")
    print("wrote gene_symbols.tsv", len(seen))
symbols = pd.read_csv(sym_path, sep="\t").set_index("gene_id")["symbol"]

# ---------------------------------------------------------------------------
# 2. gene x celltype trend matrix -> recurrent vs private
# ---------------------------------------------------------------------------
frames = []
for ct_dir in sorted(TREND.iterdir()):
    tsv = ct_dir / "length_trend_by_gene.tsv"
    if not tsv.exists():
        continue
    df = pd.read_csv(tsv, sep="\t")
    df["celltype"] = ct_dir.name
    frames.append(df)
trend = pd.concat(frames, ignore_index=True)
trend["symbol"] = trend["gene_id"].map(symbols)
trend.to_csv(OUT / "celltype_gene_trend_matrix.csv", index=False)
print("wrote celltype_gene_trend_matrix.csv", trend.shape, trend.celltype.nunique(), "celltypes")

# recurrence decomposition: a gene "trends" in a cell type when |spearman| >= 0.8
TREND_CUT = 0.8
t = trend[trend["spearman"].abs() >= TREND_CUT].copy()
t["is_dec"] = t["slope"] < 0
rec = (
    t.groupby("gene_id")
    .agg(
        n_celltypes_trending=("celltype", "nunique"),
        n_decreasing=("is_dec", "sum"),
        mean_slope=("slope", "mean"),
        mean_abs_spearman=("spearman", lambda s: s.abs().mean()),
    )
    .reset_index()
)
rec["n_increasing"] = rec["n_celltypes_trending"] - rec["n_decreasing"]
rec["dominant_direction"] = np.where(rec["n_decreasing"] >= rec["n_increasing"], "decreasing", "increasing")
rec["direction_consistency"] = rec[["n_decreasing", "n_increasing"]].max(axis=1) / rec["n_celltypes_trending"]
# how many cell types even tested this gene -> recurrence must be normalised by it
tested = trend.groupby("gene_id")["celltype"].nunique().rename("n_celltypes_tested")
rec = rec.merge(tested, on="gene_id")
rec["recurrence_frac"] = rec["n_celltypes_trending"] / rec["n_celltypes_tested"]
rec["symbol"] = rec["gene_id"].map(symbols)
rec = rec.sort_values(["n_celltypes_trending", "mean_abs_spearman"], ascending=False)
rec.to_csv(OUT / "recurrence_decomposition.csv", index=False)
print("wrote recurrence_decomposition.csv", rec.shape)

# ---------------------------------------------------------------------------
# 3. gene -> proximal/distal PAS interval, from the classic PDUI tables
# ---------------------------------------------------------------------------
iv_path = OUT / "gene_interval_map.tsv"
if not iv_path.exists():
    parts = []
    for ct_dir in sorted(LENGTH.iterdir()):
        tsv = ct_dir / "classic" / "pdui_classic.tsv"
        if not tsv.exists():
            continue
        # unique (gene, proximal, distal) triples only -- streamed, never loaded whole
        cmd = (
            f"awk -F'\\t' 'NR>1 {{print $1\"\\t\"$11\"\\t\"$12\"\\t\"$13\"\\t\"$14\"\\t\""
            f"$15\"\\t\"$16\"\\t\"$17\"\\t\"$18}}' {tsv} | sort -u"
        )
        res = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, check=True)
        if not res.stdout.strip():
            continue
        df = pd.DataFrame(
            [l.split("\t") for l in res.stdout.strip().split("\n")],
            columns=[
                "gene_id", "prox_chrom", "prox_start", "prox_end", "prox_strand",
                "dist_chrom", "dist_start", "dist_end", "dist_strand",
            ],
        )
        parts.append(df)
        print("  intervals", ct_dir.name, len(df))
    iv = pd.concat(parts, ignore_index=True).drop_duplicates(subset=["gene_id"])
    for c in ["prox_start", "prox_end", "dist_start", "dist_end"]:
        iv[c] = iv[c].astype(int)
    # lost distal segment = the stretch between the two PAS, present only in the
    # long isoform.  Orientation matters: on '-' the distal PAS is at lower coord.
    plus = iv["prox_strand"] == "+"
    iv["lost_start"] = np.where(plus, iv["prox_end"], iv["dist_end"])
    iv["lost_end"] = np.where(plus, iv["dist_start"], iv["prox_start"])
    iv["lost_len"] = iv["lost_end"] - iv["lost_start"]
    # proximal control = equally long stretch immediately upstream of proximal PAS
    iv["ctrl_start"] = np.where(plus, iv["prox_start"] - iv["lost_len"], iv["prox_end"])
    iv["ctrl_end"] = np.where(plus, iv["prox_start"], iv["prox_end"] + iv["lost_len"])
    iv = iv[(iv["lost_len"] > 50) & (iv["lost_len"] < 50000) & (iv["ctrl_start"] > 0)]
    iv.to_csv(iv_path, sep="\t", index=False)
    print("wrote gene_interval_map.tsv", iv.shape)
iv = pd.read_csv(iv_path, sep="\t", dtype={"prox_chrom": str, "dist_chrom": str})

# ---------------------------------------------------------------------------
# 4. sequence scan of the lost distal segment vs proximal control
# ---------------------------------------------------------------------------
# Mature miRNA sequences (miRBase, 5'->3'); the seed is positions 2-8.  A target
# 7mer-m8 site is the reverse complement of that seed.  Written out explicitly so
# every seed in the report is auditable rather than asserted.
MIRNAS = {
    "miR-34a-5p":  "UGGCAGUGUCUUAGCUGGUUGU",
    "let-7a-5p":   "UGAGGUAGUAGGUUGUAUAGUU",
    "miR-21-5p":   "UAGCUUAUCAGACUGAUGUUGA",
    "miR-221-3p":  "AGCUACAUUGUCUGCUGGGUUUC",
    "miR-17-5p":   "CAAAGUGCUUACAGUGCAGGUAG",
    "miR-155-5p":  "UUAAUGCUAAUCGUGAUAGGGGU",
}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def seed_site(mature: str) -> str:
    """7mer-m8 target site = reverse complement of mature positions 2-8."""
    seed = mature.replace("U", "T")[1:8]
    return "".join(COMP[b] for b in reversed(seed))


SEED_SITES = {name: seed_site(seq) for name, seq in MIRNAS.items()}
print("seed sites:", SEED_SITES)

ARE_PENT = "ATTTA"
ARE_NONA = re.compile(r"(?=([AT]TTTATTTA[AT]))")


def fetch(regions: list[str]) -> list[str]:
    """Batch samtools faidx; returns uppercase sequences in input order."""
    if not regions:
        return []
    out: list[str] = []
    CH = 500
    for i in range(0, len(regions), CH):
        chunk = regions[i : i + CH]
        res = subprocess.run(
            ["samtools", "faidx", str(FASTA), *chunk], capture_output=True, text=True
        )
        cur: list[str] = []
        for line in res.stdout.splitlines():
            if line.startswith(">"):
                if cur:
                    out.append("".join(cur).upper())
                cur = []
            else:
                cur.append(line)
        if cur:
            out.append("".join(cur).upper())
    return out


def count_all(seq: str, sub: str) -> int:
    """Overlapping substring count."""
    n = start = 0
    while (idx := seq.find(sub, start)) != -1:
        n += 1
        start = idx + 1
    return n


rows = []
CHUNK = 2000
for lo in range(0, len(iv), CHUNK):
    part = iv.iloc[lo : lo + CHUNK]
    lost_r = [f"{r.prox_chrom}:{int(r.lost_start)+1}-{int(r.lost_end)}" for r in part.itertuples()]
    ctrl_r = [f"{r.prox_chrom}:{int(r.ctrl_start)+1}-{int(r.ctrl_end)}" for r in part.itertuples()]
    lost_s, ctrl_s = fetch(lost_r), fetch(ctrl_r)
    if len(lost_s) != len(part) or len(ctrl_s) != len(part):
        print("  WARN length mismatch", lo, len(lost_s), len(ctrl_s), len(part))
        continue
    for r, ls, cs in zip(part.itertuples(), lost_s, ctrl_s):
        if r.prox_strand == "-":  # report the mRNA-sense strand
            ls = "".join(COMP.get(b, "N") for b in reversed(ls))
            cs = "".join(COMP.get(b, "N") for b in reversed(cs))
        rec_row = {
            "gene_id": r.gene_id,
            "strand": r.prox_strand,
            "lost_len": len(ls),
            "ctrl_len": len(cs),
            "lost_at_frac": (ls.count("A") + ls.count("T")) / max(len(ls), 1),
            "ctrl_at_frac": (cs.count("A") + cs.count("T")) / max(len(cs), 1),
            "lost_are_pent": count_all(ls, ARE_PENT),
            "ctrl_are_pent": count_all(cs, ARE_PENT),
            "lost_are_nona": len(ARE_NONA.findall(ls)),
            "ctrl_are_nona": len(ARE_NONA.findall(cs)),
        }
        for name, site in SEED_SITES.items():
            rec_row[f"lost_{name}"] = count_all(ls, site)
            rec_row[f"ctrl_{name}"] = count_all(cs, site)
        rows.append(rec_row)
    print("  scanned", lo + len(part), "/", len(iv))

scan = pd.DataFrame(rows)
scan["symbol"] = scan["gene_id"].map(symbols)
scan["lost_seed_total"] = scan[[f"lost_{n}" for n in SEED_SITES]].sum(axis=1)
scan["ctrl_seed_total"] = scan[[f"ctrl_{n}" for n in SEED_SITES]].sum(axis=1)
scan.to_csv(OUT / "lost_distal_element_scan.csv", index=False)
print("wrote lost_distal_element_scan.csv", scan.shape)
print("DONE s3")
