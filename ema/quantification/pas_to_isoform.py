"""PAS-to-isoform mapping for PeakATail.

For each PAS in the unified BED file, finds every transcript whose
3'UTR contains that PAS (using ``bedtools intersect -wa -wb -s``) and
computes:

- ``transcript_pos``: strand-aware position within the transcript's UTR
  (sum of upstream UTR-exon lengths + offset within the current exon).
- ``rank``: 0-based proximal-to-distal rank of the PAS within the
  transcript's PAS list (0 = most proximal / closest to stop codon).
- ``total_pas_in_transcript``: total number of PAS assigned to that
  transcript (used by proportion/shannon strategies).

Output format::

    {
        pas_id (int): [
            (gene_id, transcript_id, transcript_pos, rank, total_pas_in_transcript),
            ...
        ]
    }

The heavy genomic intersection is delegated to ``bedtools intersect``
(a C binary) via ``subprocess``.  Post-processing is done in pandas
with vectorized operations.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd


# ---------------------------------------------------------------------------
# BED builder for isoform UTRs
# ---------------------------------------------------------------------------

def _build_utr_bed(isoform_utrs: dict[str, dict[str, list[tuple]]]) -> str:
    """Serialize the isoform UTR map to BED6+2 text.

    Columns: chrom, start, end, gene_id, transcript_id, strand, exon_rank.

    Returns the BED text as a string (written to a temp file by the caller).
    """
    lines: list[str] = []
    for gene_id, t_map in isoform_utrs.items():
        for transcript_id, exons in t_map.items():
            for chrom, start, end, strand, exon_rank in exons:
                lines.append(
                    f"{chrom}\t{start}\t{end}\t{gene_id}\t{transcript_id}\t{strand}\t{exon_rank}\n"
                )
    return "".join(lines)


# ---------------------------------------------------------------------------
# Bedtools intersection
# ---------------------------------------------------------------------------

def _run_bedtools_intersect(pasbed_path: Path, utr_bed_text: str) -> pd.DataFrame:
    """Run ``bedtools intersect -wa -wb -s`` and return the joined DataFrame.

    Columns returned:
        pas_chrom, pas_start, pas_end, pas_id, pas_score, pas_strand,
        utr_chrom, utr_start, utr_end, gene_id, transcript_id,
        utr_strand, exon_rank

    Args:
        pasbed_path: Path to the PAS BED file (BED6 with integer pas_id in col 4).
        utr_bed_text: BED text string for all UTR exons.

    Returns:
        DataFrame of intersections (may be empty).

    Raises:
        RuntimeError: If bedtools returns a non-zero exit code.
    """
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bed", delete=False, prefix="isoform_utr_"
    ) as tmp:
        tmp.write(utr_bed_text)
        utr_bed_path = tmp.name

    try:
        cmd = [
            "bedtools", "intersect",
            "-a", str(pasbed_path),
            "-b", utr_bed_path,
            "-wa", "-wb",
            "-s",
        ]
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode not in (0, 1):
            raise RuntimeError(
                f"bedtools intersect failed (rc={proc.returncode}): {proc.stderr[:500]}"
            )

        if not proc.stdout.strip():
            return pd.DataFrame(columns=[
                "pas_chrom", "pas_start", "pas_end", "pas_id", "pas_score",
                "pas_strand", "utr_chrom", "utr_start", "utr_end",
                "gene_id", "transcript_id", "utr_strand", "exon_rank",
            ])

        from io import StringIO
        df = pd.read_csv(
            StringIO(proc.stdout),
            sep="\t",
            header=None,
            names=[
                "pas_chrom", "pas_start", "pas_end", "pas_id", "pas_score",
                "pas_strand", "utr_chrom", "utr_start", "utr_end",
                "gene_id", "transcript_id", "utr_strand", "exon_rank",
            ],
            dtype={
                "pas_start": int, "pas_end": int, "pas_id": int,
                "utr_start": int, "utr_end": int, "exon_rank": int,
            },
        )
        return df
    finally:
        Path(utr_bed_path).unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Transcript position calculation
# ---------------------------------------------------------------------------

def _compute_transcript_pos(row: Any, exon_index: dict) -> int:
    """Compute position within transcript coordinates.

    Transcript position = sum of lengths of all UTR exons that are
    *upstream* (5'-ward) of the current exon, plus the offset of the PAS
    within the current exon.

    Args:
        row: A pandas Series from the merged intersection DataFrame.
        exon_index: {(gene_id, transcript_id): [(chrom, start, end, strand, rank), ...]}

    Returns:
        Integer position in transcript-coordinate space (0-based from 5' end
        of the UTR).
    """
    key = (row["gene_id"], row["transcript_id"])
    exons = exon_index.get(key, [])
    if not exons:
        return 0

    strand = row["pas_strand"]
    pas_mid = (row["pas_start"] + row["pas_end"]) // 2

    # Exons are already sorted in transcription order (5'->3') with exon_rank
    exons_by_rank = sorted(exons, key=lambda e: e[4])

    upstream_len = 0
    for chrom, ex_start, ex_end, ex_strand, ex_rank in exons_by_rank:
        if ex_start <= pas_mid < ex_end:
            # PAS falls in this exon
            if strand == "+":
                offset = pas_mid - ex_start
            else:
                offset = ex_end - pas_mid - 1
            return upstream_len + max(0, offset)
        elif (
            (strand == "+" and ex_end <= pas_mid)
            or (strand == "-" and ex_start > pas_mid)
        ):
            upstream_len += ex_end - ex_start
        # Exons downstream of PAS are skipped

    return upstream_len


# ---------------------------------------------------------------------------
# Rank assignment
# ---------------------------------------------------------------------------

def _assign_ranks(df: pd.DataFrame) -> pd.DataFrame:
    """Assign proximal-to-distal ranks per transcript.

    For a given (gene_id, transcript_id), all PAS that intersect the UTR
    are ranked 0 (proximal) .. N-1 (distal) based on transcript_pos.

    Args:
        df: DataFrame with columns [pas_id, gene_id, transcript_id,
            transcript_pos].

    Returns:
        DataFrame with ``rank`` and ``total_pas_in_transcript`` columns added.
    """
    df = df.copy()
    df["rank"] = (
        df.groupby(["gene_id", "transcript_id"])["transcript_pos"]
        .rank(method="first", ascending=True)
        .astype(int) - 1
    )
    df["total_pas_in_transcript"] = df.groupby(["gene_id", "transcript_id"])[
        "pas_id"
    ].transform("count")
    return df


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def map_pas_to_isoforms(
    pasbed_path: Path,
    isoform_utrs: dict[str, dict[str, list[tuple]]],
) -> dict[int, list[tuple[str, str, int, int, int]]]:
    """Map each PAS to every isoform whose UTR contains it.

    Uses ``bedtools intersect -wa -wb -s`` for the genomic intersection
    (C-native, fast) and pandas vectorized operations for post-processing.

    Args:
        pasbed_path: Path to the PAS BED file.  Column 4 (0-based) must
            contain integer PAS identifiers.
        isoform_utrs: Output of :func:`parse_isoform_utrs`::

                {gene_id: {transcript_id: [(chrom, start, end, strand, exon_rank), ...]}}

    Returns:
        Dict keyed by ``pas_id`` (int)::

            {
                pas_id: [
                    (gene_id, transcript_id, transcript_pos, rank, total_pas_in_transcript),
                    ...
                ]
            }

        PAS that fall outside every known UTR are absent from the dict.

    Raises:
        FileNotFoundError: If *pasbed_path* does not exist.
        RuntimeError: If bedtools is not on PATH or fails.
    """
    pasbed_path = Path(pasbed_path)
    if not pasbed_path.exists():
        raise FileNotFoundError(f"PAS BED file not found: {pasbed_path}")

    if not isoform_utrs:
        return {}

    # Build a flat exon index for O(1) lookup during transcript-pos computation
    exon_index: dict[tuple, list] = {}
    for gene_id, t_map in isoform_utrs.items():
        for transcript_id, exons in t_map.items():
            exon_index[(gene_id, transcript_id)] = exons

    # Serialize UTR exons to BED text
    utr_bed_text = _build_utr_bed(isoform_utrs)

    # Genomic intersection via bedtools
    intersect_df = _run_bedtools_intersect(pasbed_path, utr_bed_text)

    if intersect_df.empty:
        return {}

    # Compute transcript-coordinate position (vectorized via apply)
    intersect_df["transcript_pos"] = intersect_df.apply(
        _compute_transcript_pos, axis=1, exon_index=exon_index
    )

    # Assign ranks
    intersect_df = _assign_ranks(intersect_df)

    # Build output dict
    result: dict[int, list[tuple[str, str, int, int, int]]] = {}
    for row in intersect_df.itertuples(index=False):
        entry = (
            row.gene_id,
            row.transcript_id,
            int(row.transcript_pos),
            int(row.rank),
            int(row.total_pas_in_transcript),
        )
        pas_id = int(row.pas_id)
        result.setdefault(pas_id, []).append(entry)

    return result
