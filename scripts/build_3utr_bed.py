#!/usr/bin/env python3
"""Extract a merged, strand-aware 3'UTR-only BED from a GTF.

Used by the FILTER_EFFECT experiment's ``annot_filter_3utr`` scenario:
``peakatail run --annot-filter --annotation-bed <this output>`` keeps only PAS
that fall inside an annotated transcript 3'UTR, which is a narrower (and
biologically more specific) region set than the default ``--annot-filter``
fallback (``ema/annotate/gtftobed.py``'s gene-body BED, source_type="gene").

Reads every ``three_prime_utr`` feature row (GTF column 3) directly from the
GTF -- no dependency on ema's gene/isoform caches -- converts to 0-based
BED, then merges overlapping/bookended intervals per-strand with bedtools
so overlapping transcript UTRs collapse into one interval (keeps the filter
fast and avoids double-counting).

Usage
-----
    python scripts/build_3utr_bed.py --gtf <path> --out <out.bed>

Output BED is 6-column: chrom, start(0-based), end, name("3UTR"), score(0),
strand. Sorted by chrom,start (required by downstream bedtools intersect).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def extract_three_prime_utrs(gtf_path: Path, raw_out: Path) -> int:
    """Stream three_prime_utr rows from *gtf_path* into a raw (unmerged) BED.

    Returns the number of rows written. GTF is 1-based closed; BED is
    0-based half-open, so start is decremented by one.
    """
    n = 0
    with open(gtf_path) as fh, open(raw_out, "w") as out:
        for line in fh:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or fields[2] != "three_prime_utr":
                continue
            chrom, start, end, strand = fields[0], fields[3], fields[4], fields[6]
            out.write(f"{chrom}\t{int(start) - 1}\t{end}\t3UTR\t0\t{strand}\n")
            n += 1
    return n


def sort_and_merge(raw_bed: Path, out_bed: Path) -> None:
    """Sort by chrom,start then bedtools-merge per strand, write BED6."""
    sorted_bed = raw_bed.with_suffix(".sorted.bed")
    subprocess.run(
        f"sort -k1,1 -k2,2n {raw_bed} > {sorted_bed}", shell=True, check=True
    )
    with open(out_bed, "w") as out:
        for strand in ("+", "-"):
            strand_bed = raw_bed.with_suffix(f".{'plus' if strand == '+' else 'minus'}.bed")
            subprocess.run(
                f"awk -F'\\t' '$6==\"{strand}\"' {sorted_bed} > {strand_bed}",
                shell=True, check=True,
            )
            if strand_bed.stat().st_size == 0:
                continue
            merged = subprocess.run(
                ["bedtools", "merge", "-s", "-c", "4,5,6", "-o", "distinct,distinct,distinct",
                 "-i", str(strand_bed)],
                check=True, capture_output=True, text=True,
            ).stdout
            for line in merged.splitlines():
                f = line.split("\t")
                # bedtools merge -c 4,5,6 -o distinct may collapse score to
                # "0" already; force score=0, name=3UTR, strand preserved.
                out.write(f"{f[0]}\t{f[1]}\t{f[2]}\t3UTR\t0\t{strand}\n")
    subprocess.run(f"sort -k1,1 -k2,2n -o {out_bed} {out_bed}", shell=True, check=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", type=Path, required=True, help="Input GTF (e.g. Ensembl GRCh38.99)")
    ap.add_argument("--out", type=Path, required=True, help="Output merged 3'UTR BED path")
    ap.add_argument("--tmpdir", type=Path, default=None,
                     help="Scratch dir for intermediates (default: system tmp; "
                          "MUST be set on the BioLab node -- see run_experiment.sh)")
    args = ap.parse_args(argv)

    if not args.gtf.exists():
        sys.exit(f"GTF not found: {args.gtf}")
    args.out.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=str(args.tmpdir) if args.tmpdir else None) as td:
        raw = Path(td) / "raw_3utr.bed"
        n_raw = extract_three_prime_utrs(args.gtf, raw)
        if n_raw == 0:
            sys.exit(f"No three_prime_utr rows found in {args.gtf} -- wrong GTF/feature naming?")
        sort_and_merge(raw, args.out)

    n_merged = sum(1 for _ in open(args.out))
    print(f"three_prime_utr rows in GTF : {n_raw}")
    print(f"merged 3'UTR intervals      : {n_merged}")
    print(f"wrote                       : {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
