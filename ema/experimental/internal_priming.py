"""Internal priming filter for PeakATail.

Filters peaks near genomic A-rich stretches that cause false positive
poly(A) signals from internal priming during reverse transcription.

Reference implementations:
- SAPAS: >=6 A's genomic check
- scraps: A-content analysis
- scAPAtrap: 6 continuous A's within -140 to +10 bp
- scPAISO: 6-mer AAAAAA within -5/+20 bp

This filter is OPTIONAL and configurable via CLI:
  --internal-priming-filter     Enable the filter
  --genome-fasta PATH           Path to genome FASTA (required if filter enabled)
  --ip-window-left INT          Left window from PAS (default: -10)
  --ip-window-right INT         Right window from PAS (default: +30)
  --ip-a-stretch INT            Minimum consecutive A's to flag (default: 6)
  --ip-a-fraction FLOAT         Alternative: fraction of A's in window (default: 0.7)
"""

import logging
from typing import List, Tuple, Optional

logger = logging.getLogger("peakatail.filters.internal_priming")


def filter_internal_priming(bed_path: str, genome_fasta: str,
                            output_path: str,
                            window_left: int = 10,
                            window_right: int = 30,
                            a_stretch: int = 6,
                            a_fraction: float = 0.7) -> dict:
    """Filter peaks near genomic A-rich stretches.

    Reads a BED file, checks the genomic sequence around each PAS for
    A-rich stretches that indicate internal priming artifacts. Writes
    filtered peaks to output and returns statistics.

    Args:
        bed_path: Input BED file with PAS peaks.
        genome_fasta: Path to genome FASTA file (must be indexed with .fai).
        output_path: Path to write filtered BED output.
        window_left: bp upstream of PAS to check (default 10).
        window_right: bp downstream of PAS to check (default 30).
        a_stretch: Minimum consecutive A's to flag as internal priming (default 6).
        a_fraction: Alternative: flag if A-fraction in window exceeds this (default 0.7).

    Returns:
        Dict with filtering statistics:
        - total: total peaks processed
        - passed: peaks that passed the filter
        - filtered: peaks removed as internal priming
        - filtered_fraction: fraction removed
    """
    try:
        from pyfaidx import Fasta
    except ImportError:
        logger.error("pyfaidx not installed. Run: pip install pyfaidx")
        logger.error("Skipping internal priming filter.")
        # Copy input to output unchanged
        import shutil
        shutil.copy2(bed_path, output_path)
        return {"total": 0, "passed": 0, "filtered": 0, "filtered_fraction": 0,
                "error": "pyfaidx not installed"}

    genome = Fasta(genome_fasta)
    a_pattern = "A" * a_stretch

    total = 0
    passed = 0
    filtered = 0

    with open(bed_path) as infile, open(output_path, 'w') as outfile:
        for line in infile:
            total += 1
            parts = line.strip().split('\t')
            if len(parts) < 6:
                outfile.write(line)
                passed += 1
                continue

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            strand = parts[5]

            # PAS position depends on strand
            if strand == '+':
                pas_pos = end  # 3' end for positive strand
            else:
                pas_pos = start  # 3' end for negative strand

            # Extract genomic sequence around PAS
            try:
                seq_start = max(0, pas_pos - window_left)
                seq_end = pas_pos + window_right
                seq = str(genome[chrom][seq_start:seq_end]).upper()
            except (KeyError, ValueError):
                # Chromosome not in FASTA or out of range — keep the peak
                outfile.write(line)
                passed += 1
                continue

            # Check for consecutive A stretch
            is_internal_priming = False

            if strand == '+':
                # For positive strand, check for A's downstream
                if a_pattern in seq:
                    is_internal_priming = True
            else:
                # For negative strand, check for T's (complementary)
                t_pattern = "T" * a_stretch
                if t_pattern in seq:
                    is_internal_priming = True

            # Also check A-fraction in the window
            if not is_internal_priming and len(seq) > 0:
                if strand == '+':
                    a_count = seq.count('A')
                else:
                    a_count = seq.count('T')
                if a_count / len(seq) >= a_fraction:
                    is_internal_priming = True

            if is_internal_priming:
                filtered += 1
                logger.debug(f"Internal priming: {chrom}:{pas_pos} strand={strand} seq={seq[:20]}...")
            else:
                outfile.write(line)
                passed += 1

    filtered_fraction = filtered / total if total > 0 else 0

    stats = {
        "total": total,
        "passed": passed,
        "filtered": filtered,
        "filtered_fraction": round(filtered_fraction, 4)
    }

    logger.info(f"Internal priming filter: {total} total, {passed} passed, "
                f"{filtered} filtered ({filtered_fraction:.1%})")

    return stats
