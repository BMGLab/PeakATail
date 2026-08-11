"""Configurable post-processing filters for PeakATail peaks.

All filters are OPTIONAL and can be enabled/disabled via CLI flags.
This allows users to choose between discovery mode (all peaks) and
publication mode (filtered, high-confidence peaks).

Usage:
    from ema.experimental.peak_filters import apply_filters
    stats = apply_filters(input_bed, output_bed, config)
"""

import os
import logging
from typing import Optional

logger = logging.getLogger("peakatail.filters")


def apply_filters(input_bed: str, output_bed: str,
                  genome_fasta: Optional[str] = None,
                  annotation_bed: Optional[str] = None,
                  enable_internal_priming: bool = False,
                  enable_annotation_filter: bool = False,
                  ip_window_left: int = 10,
                  ip_window_right: int = 30,
                  ip_a_stretch: int = 6,
                  ip_a_fraction: float = 0.7,
                  ip_mode: str = "annotate") -> dict:
    """Apply all enabled filters sequentially to a BED file.

    Filters are applied in order:
    1. Internal priming filter (if enabled + genome FASTA provided)
    2. Annotation filter (if enabled + annotation BED provided)

    Each filter writes to a temp file, which becomes input for the next.
    The final result is written to output_bed.

    Args:
        input_bed: Path to input BED file.
        output_bed: Path to write filtered output.
        genome_fasta: Path to genome FASTA (required for internal priming filter).
        annotation_bed: Path to gene annotation BED (required for annotation filter).
        enable_internal_priming: Enable internal priming filter.
        enable_annotation_filter: Enable annotation region filter.
        ip_window_left: Left window for internal priming check.
        ip_window_right: Right window for internal priming check.
        ip_a_stretch: Minimum consecutive A's for internal priming.
        ip_a_fraction: Minimum A-fraction for internal priming.
        ip_mode: ``"annotate"`` (default, keep + flag) or ``"filter"`` (drop)
            -- forwarded to :func:`ema.experimental.internal_priming.
            filter_internal_priming`. Only affects the internal-priming
            filter; the annotation-region filter (``enable_annotation_filter``)
            is a distinct "keep only PAS overlapping a gene region" concept
            and always drops non-overlapping peaks, unaffected by this flag.

    Returns:
        Dict with statistics from each applied filter. When the internal
        priming filter runs, ``all_stats["internal_priming_flags"]`` is the
        ``{pas_id: bool}`` map from :func:`filter_internal_priming`.
    """
    all_stats = {"filters_applied": []}
    current_input = input_bed

    # Filter 1: Internal priming
    if enable_internal_priming:
        if genome_fasta and os.path.exists(genome_fasta):
            from ema.experimental.internal_priming import filter_internal_priming
            ip_output = output_bed + ".ip_tmp"
            ip_stats = filter_internal_priming(
                current_input, genome_fasta, ip_output,
                window_left=ip_window_left,
                window_right=ip_window_right,
                a_stretch=ip_a_stretch,
                a_fraction=ip_a_fraction,
                mode=ip_mode,
            )
            all_stats["internal_priming"] = ip_stats
            all_stats["internal_priming_flags"] = ip_stats.get("flags", {})
            all_stats["filters_applied"].append("internal_priming")
            current_input = ip_output
            logger.info(f"Internal priming (mode={ip_mode}): {ip_stats['filtered']}/{ip_stats['total']} "
                       f"removed, {ip_stats.get('flagged', 0)}/{ip_stats['total']} flagged")
        else:
            logger.warning("Internal priming filter enabled but no genome FASTA provided. Skipping.")

    # Filter 2: Annotation filter (keep only peaks in annotated regions)
    if enable_annotation_filter:
        if annotation_bed and os.path.exists(annotation_bed):
            annot_output = output_bed + ".annot_tmp"
            annot_stats = _filter_by_annotation(current_input, annotation_bed, annot_output)
            all_stats["annotation"] = annot_stats
            all_stats["filters_applied"].append("annotation")
            current_input = annot_output
            logger.info(f"Annotation filter: {annot_stats['filtered']}/{annot_stats['total']} "
                       f"removed ({annot_stats['filtered_fraction']:.1%})")
        else:
            logger.warning("Annotation filter enabled but no annotation BED provided. Skipping.")

    # Copy final result to output
    if current_input != input_bed:
        import shutil
        shutil.move(current_input, output_bed)
        # Clean up any remaining temp files
        for suffix in [".ip_tmp", ".annot_tmp"]:
            tmp = output_bed + suffix
            if os.path.exists(tmp):
                os.remove(tmp)
    else:
        # No filters applied — copy input to output
        import shutil
        shutil.copy2(input_bed, output_bed)

    # Count final peaks
    final_count = sum(1 for _ in open(output_bed))
    input_count = sum(1 for _ in open(input_bed))
    all_stats["input_peaks"] = input_count
    all_stats["output_peaks"] = final_count
    all_stats["total_removed"] = input_count - final_count
    all_stats["total_removed_fraction"] = round(
        (input_count - final_count) / max(input_count, 1), 4
    )

    return all_stats


def _filter_by_annotation(bed_path: str, annotation_bed: str,
                          output_path: str) -> dict:
    """Keep only peaks that overlap annotated genomic regions.

    Uses pybedtools intersect to keep peaks overlapping 3'UTR, terminal
    exons, or other annotated features.

    Args:
        bed_path: Input BED file.
        annotation_bed: Annotation BED file (from GTF).
        output_path: Output filtered BED.

    Returns:
        Dict with filtering statistics.
    """
    from pybedtools import BedTool

    peaks = BedTool(bed_path)
    annotation = BedTool(annotation_bed)

    total = peaks.count()

    # Keep peaks that overlap ANY annotated region
    filtered = peaks.intersect(annotation, u=True)
    filtered.saveas(output_path)

    passed = sum(1 for _ in open(output_path))
    removed = total - passed

    return {
        "total": total,
        "passed": passed,
        "filtered": removed,
        "filtered_fraction": round(removed / max(total, 1), 4)
    }
