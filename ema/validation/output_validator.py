import os
import logging
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger("peakatail.validation")


class ValidationResult:
    """Holds results of output validation checks."""

    def __init__(self):
        self.checks = []  # list of (name, passed, message)

    def add(self, name: str, passed: bool, message: str = ""):
        self.checks.append((name, passed, message))
        level = logging.INFO if passed else logging.ERROR
        status = "PASS" if passed else "FAIL"
        logger.log(level, f"  [{status}] {name}: {message}")

    @property
    def passed(self) -> bool:
        return all(p for _, p, _ in self.checks)

    @property
    def failed_checks(self) -> list:
        return [(n, m) for n, p, m in self.checks if not p]

    def summary(self) -> str:
        total = len(self.checks)
        passed = sum(1 for _, p, _ in self.checks if p)
        lines = [f"Validation: {passed}/{total} checks passed"]
        for name, p, msg in self.checks:
            status = "PASS" if p else "FAIL"
            lines.append(f"  [{status}] {name}: {msg}")
        return "\n".join(lines)


def validate_output(bed_path: str, mtx_path: str,
                    expected_chromosomes: Optional[List[str]] = None) -> ValidationResult:
    """Run all validation checks on peak calling output.

    Args:
        bed_path: Path to output BED file
        mtx_path: Path to output MTX (MatrixMarket) file
        expected_chromosomes: Optional list of expected chromosome names

    Returns:
        ValidationResult with all check results
    """
    result = ValidationResult()
    logger.info(f"Validating output: {bed_path}, {mtx_path}")

    # --- BED file checks ---
    bed_lines = _read_bed(bed_path, result)
    if bed_lines is None:
        return result

    _check_bed_not_empty(bed_lines, result)
    _check_bed_columns(bed_lines, result)
    _check_bed_coordinates(bed_lines, result)
    _check_bed_strand(bed_lines, result)
    _check_bed_unique_ids(bed_lines, result)
    _check_bed_peak_count(bed_lines, result)
    _check_bed_width_range(bed_lines, result)
    _check_bed_strand_ratio(bed_lines, result)

    if expected_chromosomes:
        _check_bed_chromosomes(bed_lines, expected_chromosomes, result)

    # --- MTX file checks ---
    mtx_lines = _read_mtx(mtx_path, result)
    if mtx_lines is not None:
        _check_mtx_not_empty(mtx_lines, result)
        _check_mtx_values(mtx_lines, result)
        _check_mtx_bed_consistency(bed_lines, mtx_lines, result)

    logger.info(result.summary())
    return result


def _read_bed(path: str, result: ValidationResult):
    """Read BED file and return parsed lines."""
    if not os.path.exists(path):
        result.add("bed_exists", False, f"BED file not found: {path}")
        return None

    lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                lines.append(line.split("\t"))
    result.add("bed_readable", True, f"Read {len(lines)} lines from BED")
    return lines


def _read_mtx(path: str, result: ValidationResult):
    """Read MatrixMarket file lines (skip header)."""
    if not os.path.exists(path):
        result.add("mtx_exists", False, f"MTX file not found: {path}")
        return None

    lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("%"):
                lines.append(line.split())
    result.add("mtx_readable", True, f"Read {len(lines)} entries from MTX")
    return lines


def _check_bed_not_empty(lines, result):
    result.add("bed_not_empty", len(lines) > 0,
               f"Peak count: {len(lines)}")


def _check_bed_columns(lines, result):
    bad = [i for i, l in enumerate(lines) if len(l) < 6]
    result.add("bed_6_columns", len(bad) == 0,
               f"{len(bad)} lines with <6 columns" if bad else "All lines have 6+ columns")


def _check_bed_coordinates(lines, result):
    bad = []
    for i, l in enumerate(lines):
        try:
            start, end = int(l[1]), int(l[2])
            if start >= end:
                bad.append(i)
        except (ValueError, IndexError):
            bad.append(i)
    result.add("bed_valid_coords", len(bad) == 0,
               f"{len(bad)} lines with start >= end" if bad else "All coordinates valid (start < end)")


def _check_bed_strand(lines, result):
    bad = [i for i, l in enumerate(lines) if len(l) > 5 and l[5] not in ("+", "-")]
    result.add("bed_valid_strand", len(bad) == 0,
               f"{len(bad)} lines with invalid strand" if bad else "All strands are + or -")


def _check_bed_unique_ids(lines, result):
    ids = [l[3] for l in lines if len(l) > 3]
    unique = len(set(ids))
    result.add("bed_unique_ids", unique == len(ids),
               f"{len(ids)} IDs, {unique} unique" if unique != len(ids) else f"All {len(ids)} peak IDs unique")


def _check_bed_peak_count(lines, result):
    n = len(lines)
    # Sanity: should not have absurdly many peaks (> 1M likely a bug)
    result.add("bed_reasonable_count", 0 < n < 1_000_000,
               f"Peak count: {n}")


def _check_bed_width_range(lines, result):
    widths = []
    for l in lines:
        try:
            widths.append(int(l[2]) - int(l[1]))
        except (ValueError, IndexError):
            pass
    if widths:
        mean_w = sum(widths) / len(widths)
        result.add("bed_reasonable_width", 1 <= mean_w <= 50000,
                    f"Mean peak width: {mean_w:.0f} bp")
    else:
        result.add("bed_reasonable_width", False, "Could not compute peak widths")


def _check_bed_strand_ratio(lines, result):
    plus = sum(1 for l in lines if len(l) > 5 and l[5] == "+")
    minus = sum(1 for l in lines if len(l) > 5 and l[5] == "-")
    total = plus + minus
    if total > 0:
        ratio = plus / total
        # Expect roughly balanced (0.2 - 0.8 range)
        result.add("bed_strand_balance", 0.1 <= ratio <= 0.9,
                    f"Strand ratio: {plus}+ / {minus}- ({ratio:.2f})")
    else:
        result.add("bed_strand_balance", False, "No strand information")


def _check_bed_chromosomes(lines, expected, result):
    found = set(l[0] for l in lines if l)
    missing = set(expected) - found
    result.add("bed_chromosomes_present",
               len(missing) <= len(expected) * 0.1,  # allow 10% missing
               f"Missing chromosomes: {missing}" if missing else "All expected chromosomes present")


def _check_mtx_not_empty(lines, result):
    result.add("mtx_not_empty", len(lines) > 0,
               f"MTX entries: {len(lines)}")


def _check_mtx_values(lines, result):
    bad = []
    for i, l in enumerate(lines):
        try:
            vals = [int(x) for x in l]
            if any(v < 0 for v in vals):
                bad.append(i)
        except ValueError:
            bad.append(i)
    result.add("mtx_valid_values", len(bad) == 0,
               f"{len(bad)} lines with invalid/negative values" if bad else "All MTX values valid")


def _check_mtx_bed_consistency(bed_lines, mtx_lines, result):
    bed_ids = set(l[3] for l in bed_lines if len(l) > 3)
    mtx_rows = set(l[0] for l in mtx_lines if l)
    # Every peak in BED should have at least 1 entry in MTX
    missing = bed_ids - mtx_rows
    result.add("mtx_bed_consistent",
               len(missing) <= len(bed_ids) * 0.05,  # allow 5% missing
               f"{len(missing)}/{len(bed_ids)} BED peaks missing from MTX" if missing
               else "All BED peaks have MTX entries")
