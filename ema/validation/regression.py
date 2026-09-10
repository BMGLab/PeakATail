import os
import json
import logging
from typing import Dict, Optional

logger = logging.getLogger("peakatail.validation.regression")


class RegressionResult:
    """Holds results of regression comparison against golden reference."""

    def __init__(self):
        self.checks = []
        self.reports = []

    def add_check(self, name: str, passed: bool, message: str):
        self.checks.append((name, passed, message))
        level = logging.INFO if passed else logging.WARNING
        status = "PASS" if passed else "FAIL"
        logger.log(level, f"[{status}] {name}: {message}")

    def add_report(self, name: str, message: str):
        """Non-pass/fail observation for logging."""
        self.reports.append((name, message))
        logger.info(f"[REPORT] {name}: {message}")

    @property
    def passed(self) -> bool:
        return all(p for _, p, _ in self.checks)

    def summary(self) -> str:
        lines = ["Regression Check Summary:"]
        passed = sum(1 for _, p, _ in self.checks if p)
        lines.append(f"  {passed}/{len(self.checks)} checks passed")
        for name, p, msg in self.checks:
            status = "PASS" if p else "FAIL"
            lines.append(f"  [{status}] {name}: {msg}")
        if self.reports:
            lines.append("\nObservations:")
            for name, msg in self.reports:
                lines.append(f"  {name}: {msg}")
        return "\n".join(lines)


def capture_golden(bed_path: str, mtx_path: str, output_dir: str) -> str:
    """Capture current output as golden reference for future comparison.

    Args:
        bed_path: Path to BED output file
        mtx_path: Path to MTX output file
        output_dir: Directory to save golden reference

    Returns:
        Path to stats.json
    """
    os.makedirs(output_dir, exist_ok=True)

    stats = {"peaks": {}, "matrix": {}}

    # Parse BED
    peaks_per_chrom = {}
    strand_plus = 0
    strand_minus = 0
    widths = []
    peak_ids = set()

    if os.path.exists(bed_path):
        with open(bed_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    chrom = parts[0]
                    width = int(parts[2]) - int(parts[1])
                    widths.append(width)
                    peaks_per_chrom[chrom] = peaks_per_chrom.get(chrom, 0) + 1
                    peak_ids.add(parts[3])
                    if parts[5] == '+':
                        strand_plus += 1
                    else:
                        strand_minus += 1

        # Copy BED file
        import shutil
        shutil.copy2(bed_path, os.path.join(output_dir, os.path.basename(bed_path)))

    total_peaks = len(peak_ids)
    stats["peaks"] = {
        "total": total_peaks,
        "per_chromosome": peaks_per_chrom,
        "chromosomes": sorted(peaks_per_chrom.keys()),
        "strand_plus": strand_plus,
        "strand_minus": strand_minus,
        "strand_ratio": strand_plus / max(total_peaks, 1),
        "mean_width": sum(widths) / max(len(widths), 1) if widths else 0
    }

    # Parse MTX
    mtx_rows = set()
    mtx_entries = 0
    if os.path.exists(mtx_path):
        with open(mtx_path) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('%'):
                    parts = line.split()
                    if len(parts) >= 3:
                        mtx_rows.add(parts[0])
                        mtx_entries += 1

        import shutil
        shutil.copy2(mtx_path, os.path.join(output_dir, os.path.basename(mtx_path)))

    stats["matrix"] = {
        "unique_peaks": len(mtx_rows),
        "total_entries": mtx_entries
    }

    stats_path = os.path.join(output_dir, "stats.json")
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=2)

    logger.info(f"Golden reference captured: {total_peaks} peaks, {mtx_entries} MTX entries")
    logger.info(f"Saved to: {output_dir}")

    return stats_path


def compare_to_golden(bed_path: str, mtx_path: str,
                      golden_dir: str) -> RegressionResult:
    """Compare current output against golden reference using similarity thresholds.

    NOT an exact match — checks that output is SIMILAR ENOUGH to detect
    catastrophic failures while allowing expected algorithmic changes.

    Args:
        bed_path: Path to current BED output
        mtx_path: Path to current MTX output
        golden_dir: Path to golden reference directory (from capture_golden)

    Returns:
        RegressionResult with check results
    """
    result = RegressionResult()

    # Load golden stats
    stats_path = os.path.join(golden_dir, "stats.json")
    if not os.path.exists(stats_path):
        result.add_check("golden_exists", False, f"Golden stats not found: {stats_path}")
        return result

    with open(stats_path) as f:
        golden = json.load(f)

    golden_peaks = golden["peaks"]
    golden_total = golden_peaks["total"]

    # Parse current output
    current_peaks = 0
    current_chroms = {}
    current_plus = 0
    current_minus = 0
    current_widths = []

    if os.path.exists(bed_path):
        with open(bed_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    current_peaks += 1
                    chrom = parts[0]
                    current_chroms[chrom] = current_chroms.get(chrom, 0) + 1
                    current_widths.append(int(parts[2]) - int(parts[1]))
                    if parts[5] == '+':
                        current_plus += 1
                    else:
                        current_minus += 1

    # --- PASS/FAIL CHECKS (similarity thresholds) ---

    # Peak count within ±50%
    if golden_total > 0:
        ratio = current_peaks / golden_total
        result.add_check("peak_count_similar",
                        0.5 <= ratio <= 1.5,
                        f"Current: {current_peaks}, Golden: {golden_total}, Ratio: {ratio:.2f}")
    else:
        result.add_check("peak_count_similar",
                        current_peaks >= 0,
                        f"Golden had 0 peaks, current has {current_peaks}")

    # All golden chromosomes present
    golden_chroms = set(golden_peaks.get("chromosomes", []))
    current_chrom_set = set(current_chroms.keys())
    missing_chroms = golden_chroms - current_chrom_set
    result.add_check("chromosomes_present",
                    len(missing_chroms) == 0,
                    f"Missing: {missing_chroms}" if missing_chroms else "All chromosomes present")

    # No chromosome lost >90% of peaks
    chrom_ok = True
    bad_chroms = []
    for chrom, golden_count in golden_peaks.get("per_chromosome", {}).items():
        current_count = current_chroms.get(chrom, 0)
        if golden_count > 0 and current_count < golden_count * 0.1:
            chrom_ok = False
            bad_chroms.append(f"{chrom}: {current_count}/{golden_count}")
    result.add_check("no_chromosome_devastated", chrom_ok,
                    f"Severely reduced: {bad_chroms}" if bad_chroms else "No chromosome lost >90% peaks")

    # Strand ratio within ±20%
    current_total = current_plus + current_minus
    if current_total > 0:
        current_ratio = current_plus / current_total
        golden_ratio = golden_peaks.get("strand_ratio", 0.5)
        result.add_check("strand_ratio_stable",
                        abs(current_ratio - golden_ratio) <= 0.2,
                        f"Current: {current_ratio:.2f}, Golden: {golden_ratio:.2f}")

    # MTX not empty
    mtx_entries = 0
    if os.path.exists(mtx_path):
        with open(mtx_path) as f:
            for line in f:
                if line.strip() and not line.startswith('%'):
                    mtx_entries += 1
    result.add_check("mtx_not_empty", mtx_entries > 0,
                    f"MTX entries: {mtx_entries}")

    # Mean peak width within ±100%
    if current_widths:
        current_mean_w = sum(current_widths) / len(current_widths)
        golden_mean_w = golden_peaks.get("mean_width", 100)
        if golden_mean_w > 0:
            w_ratio = current_mean_w / golden_mean_w
            result.add_check("peak_width_reasonable",
                            0.25 <= w_ratio <= 4.0,
                            f"Current mean: {current_mean_w:.0f}bp, Golden: {golden_mean_w:.0f}bp")

    # --- REPORT (always logged, not pass/fail) ---

    result.add_report("total_peaks",
                     f"Current: {current_peaks}, Golden: {golden_total}")

    # Per-chromosome comparison
    for chrom in sorted(golden_chroms | current_chrom_set):
        g = golden_peaks.get("per_chromosome", {}).get(chrom, 0)
        c = current_chroms.get(chrom, 0)
        diff = c - g
        result.add_report(f"chrom_{chrom}", f"Golden: {g}, Current: {c}, Diff: {diff:+d}")

    logger.info(result.summary())
    return result
