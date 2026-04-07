import os
import logging
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger("peakatail.validation.db")

# Default cutoffs tuned for APA analysis.
APA_DISTANCE_CUTOFFS: List[int] = [50, 100, 200, 500, 1000, 5000]

# Precision targets per cutoff window.
# Rationale: if our PAS are real, precision should rise as the matching window
# widens because almost every true PAS will lie within a few kb of a catalogued
# site.  These targets encode that expectation explicitly.
PRECISION_TARGETS: Dict[int, float] = {
    50: 0.70,
    200: 0.80,
    1000: 0.90,
    5000: 0.95,
}


def validate_against_database(peakatail_bed: str, reference_bed: str,
                               windows: Optional[List[int]] = None,
                               name: str = "reference") -> Dict:
    """Compare PeakATail peaks against a reference PAS database.

    Uses pybedtools for BED intersection at multiple distance cutoffs.
    Precision is the primary metric: "Are the PAS we found real?"

    Args:
        peakatail_bed: Path to PeakATail output BED file.
        reference_bed: Path to reference database BED file (PolyASite/PolyA_DB).
        windows: Distance cutoffs in bp.
            Defaults to APA_DISTANCE_CUTOFFS [50, 100, 200, 500, 1000, 5000].
        name: Name of the reference database for logging.

    Returns:
        Dict with validation results at each cutoff, precision-first ordering.
    """
    from pybedtools import BedTool

    if windows is None:
        windows = APA_DISTANCE_CUTOFFS

    if not os.path.exists(peakatail_bed):
        logger.error(f"PeakATail BED file not found: {peakatail_bed}")
        return {"error": "PeakATail BED not found"}

    if not os.path.exists(reference_bed):
        logger.error(f"Reference BED file not found: {reference_bed}")
        return {"error": "Reference BED not found"}

    predicted = BedTool(peakatail_bed)
    reference = BedTool(reference_bed)

    n_predicted = predicted.count()
    n_reference = reference.count()

    logger.info(
        f"Validating against {name}: {n_predicted} predicted, {n_reference} reference PAS"
    )

    results = {
        "database": name,
        "n_predicted": n_predicted,
        "n_reference": n_reference,
        "cutoffs": {}
    }

    for window in sorted(windows):
        # Precision (PRIMARY): how many of our peaks match a known PAS?
        overlap_pred = predicted.window(reference, w=window)
        matched_pred_ids: set = set()
        for feature in overlap_pred:
            fields = str(feature).strip().split('\t')
            if len(fields) >= 4:
                matched_pred_ids.add(fields[3])

        precision = len(matched_pred_ids) / n_predicted if n_predicted > 0 else 0.0

        # Recall: how many known PAS did we find?  Kept for completeness but
        # not used as a primary quality signal.
        overlap_ref = reference.window(predicted, w=window)
        matched_ref_ids: set = set()
        for feature in overlap_ref:
            fields = str(feature).strip().split('\t')
            if len(fields) >= 4:
                matched_ref_ids.add(fields[3])

        recall = len(matched_ref_ids) / n_reference if n_reference > 0 else 0.0

        # F1 retained for completeness.
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

        # Potentially novel PAS: predicted peaks not near any reference entry.
        # These may represent real undiscovered PAS not yet in the database.
        novel = n_predicted - len(matched_pred_ids)
        novel_fraction = novel / n_predicted if n_predicted > 0 else 0.0

        cutoff_result = {
            "precision": round(precision, 4),
            "recall": round(recall, 4),
            "f1": round(f1, 4),
            "matched_predicted": len(matched_pred_ids),
            "matched_reference": len(matched_ref_ids),
            "potentially_novel_pas": novel,
            "novel_fraction": round(novel_fraction, 4),
        }

        results["cutoffs"][window] = cutoff_result
        logger.info(
            f"  @{window}bp: P={precision:.3f} R={recall:.3f} F1={f1:.3f} "
            f"potentially_novel={novel} ({novel_fraction:.1%})"
        )

    return results


def classify_peak_regions(peakatail_bed: str, annotation_bed: str) -> Dict[str, int]:
    """Classify peaks by genomic region (3'UTR, exon, intron, intergenic).

    Args:
        peakatail_bed: Path to PeakATail output BED.
        annotation_bed: Path to genomic annotation BED (from GTF).

    Returns:
        Dict with counts per region type.
    """
    from pybedtools import BedTool

    peaks = BedTool(peakatail_bed)
    annotation = BedTool(annotation_bed)

    # Intersect peaks with annotation
    intersected = peaks.intersect(annotation, wo=True)

    region_counts: Dict[str, int] = {}
    for feature in intersected:
        fields = str(feature).strip().split('\t')
        # The annotation field typically contains region type
        if len(fields) >= 10:
            region_type = fields[9] if len(fields) > 9 else "unknown"
            region_counts[region_type] = region_counts.get(region_type, 0) + 1

    # Count intergenic (peaks not overlapping any annotation)
    non_intersected = peaks.intersect(annotation, v=True)
    n_intergenic = non_intersected.count()
    region_counts["intergenic"] = n_intergenic

    return region_counts


def generate_validation_report(results: Dict, output_path: str) -> str:
    """Generate a text validation report with precision as the primary metric.

    Precision targets by distance window (encoding the expectation that real
    PAS should be increasingly well-covered as the window widens):

    - @50bp:   P >= 70%  (PASS/FAIL gate)
    - @200bp:  P >= 80%
    - @1000bp: P >= 90%
    - @5000bp: P >= 95%

    "Potentially novel PAS" (predicted peaks not near any reference) are
    reported positively — they may represent real undiscovered sites.

    Args:
        results: Output from validate_against_database().
        output_path: Path to write the report.

    Returns:
        Path to the written report file.
    """
    lines = []
    lines.append("=" * 70)
    lines.append("PeakATail Database Validation Report")
    lines.append(f"Reference: {results.get('database', 'unknown')}")
    lines.append("=" * 70)
    lines.append(f"\nPredicted PAS : {results['n_predicted']}")
    lines.append(f"Reference PAS : {results['n_reference']}")
    lines.append("")
    lines.append("Primary metric: Precision — 'Are the PAS we found real?'")
    lines.append("")

    header = (
        f"{'Cutoff':>8}  {'Precision':>10}  {'Recall':>8}  {'F1':>6}  "
        f"{'Novel PAS':>10}  {'Novel%':>7}  {'Target':>10}  {'Status':>12}"
    )
    lines.append(header)
    lines.append("-" * 80)

    for cutoff, metrics in sorted(results.get("cutoffs", {}).items()):
        p = metrics["precision"]
        r = metrics["recall"]
        f1 = metrics["f1"]
        novel = metrics["potentially_novel_pas"]
        novel_frac = metrics["novel_fraction"]

        target = PRECISION_TARGETS.get(cutoff)
        if target is not None:
            target_str = f"P>={target:.0%}"
            status = "PASS" if p >= target else "BELOW TARGET"
        else:
            target_str = ""
            status = ""

        lines.append(
            f"{cutoff:>6}bp  {p:>10.3f}  {r:>8.3f}  {f1:>6.3f}  "
            f"{novel:>10}  {novel_frac:>6.1%}  {target_str:>10}  {status:>12}"
        )

    lines.append("")
    lines.append("Precision targets (encoding expected PAS quality profile):")
    for cutoff, target in sorted(PRECISION_TARGETS.items()):
        lines.append(f"  @{cutoff}bp  : P >= {target:.0%}")
    lines.append("")
    lines.append(
        "Note: 'Novel PAS' = predicted peaks with no match in reference database.\n"
        "      These may represent real, undiscovered polyadenylation sites."
    )
    lines.append("=" * 70)

    report_text = "\n".join(lines)
    with open(output_path, 'w') as f:
        f.write(report_text)

    print(report_text)
    return output_path
