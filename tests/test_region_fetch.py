"""Tests for Phase 2: region-fetch path in peak_calling().

Verifies that peak_calling(region=(chrom, start, end)) produces the same
read set as the full-BAM path would see for that region, and that the BAI
index-check guard raises an actionable error when the index is missing.

The test BAM used here is tests/chr22.bam (already present in the repo and
indexed as tests/chr22.bam.bai).
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import pysam
import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

TESTS_DIR = Path(__file__).parent
TEST_BAM = TESTS_DIR / "SRR8325947_Aligned.sortedByCoord.out.bam"
TEST_BAI = TESTS_DIR / "SRR8325947_Aligned.sortedByCoord.out.bam.bai"
TILE_BAM = Path("/home/user/PeakATail/test_run/chr22.bam")
TILE_BAI = Path("/home/user/PeakATail/test_run/chr22.bam.bai")

# Use the small tile BAM (always present and indexed) for the region tests.
# Fall back to the large BAM only if the tile BAM is absent.
if TILE_BAM.exists() and TILE_BAI.exists():
    REGION_BAM = TILE_BAM
elif TEST_BAM.exists() and TEST_BAI.exists():
    REGION_BAM = TEST_BAM
else:
    REGION_BAM = None  # type: ignore[assignment]

# A small region on chr22 with a reasonable read density.
REGION_CHROM = "22"
REGION_START = 16_000_000
REGION_END = 17_000_000


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _collect_fetch_read_ids(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
) -> set[str]:
    """Return query_name + start + end strings for reads in the region."""
    ids: set[str] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            # Use a composite key that is unique enough for testing
            ids.add(f"{read.query_name}:{read.reference_start}:{read.reference_end}")
    return ids


def _collect_full_scan_read_ids(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
) -> set[str]:
    """Full-BAM iteration, filter to [chrom, start, end) by hand."""
    ids: set[str] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if (
                read.reference_name == chrom
                and read.reference_start is not None
                and read.reference_end is not None
                and read.reference_start >= start
                and read.reference_start < end
            ):
                ids.add(
                    f"{read.query_name}:{read.reference_start}:{read.reference_end}"
                )
    return ids


# ---------------------------------------------------------------------------
# Region fetch read-set equivalence
# ---------------------------------------------------------------------------


@pytest.mark.skipif(REGION_BAM is None, reason="No indexed test BAM available")
def test_region_fetch_returns_same_reads_as_full_scan() -> None:
    """pysam fetch should return the same reads as a full-scan filtered to region.

    This validates the pysam contract that peak_calling(region=...) relies on:
    bamfile.fetch(chrom, start, end) must yield exactly the reads whose
    reference_start is in [start, end).

    Note: pysam.fetch overlaps are based on ANY overlap with [start, end), but
    reads starting before 'start' can be returned. We compare only reads whose
    reference_start >= start, which matches what the peak-calling algorithm
    processes in the tile path.
    """
    fetch_ids = _collect_fetch_read_ids(str(REGION_BAM), REGION_CHROM, REGION_START, REGION_END)
    # fetch returns reads overlapping the region, including those starting before
    # REGION_START. Filter to the core window for comparison.
    with pysam.AlignmentFile(str(REGION_BAM), "rb") as bam:
        fetch_ids_core: set[str] = set()
        for read in bam.fetch(REGION_CHROM, REGION_START, REGION_END):
            if (
                read.reference_start is not None
                and read.reference_start >= REGION_START
                and read.reference_start < REGION_END
            ):
                fetch_ids_core.add(
                    f"{read.query_name}:{read.reference_start}:{read.reference_end}"
                )

    full_scan_ids = _collect_full_scan_read_ids(
        str(REGION_BAM), REGION_CHROM, REGION_START, REGION_END
    )
    assert fetch_ids_core == full_scan_ids, (
        f"Region fetch missed {len(full_scan_ids - fetch_ids_core)} reads "
        f"and had {len(fetch_ids_core - full_scan_ids)} extra reads vs full scan"
    )


@pytest.mark.skipif(REGION_BAM is None, reason="No indexed test BAM available")
def test_region_fetch_read_count_nonzero() -> None:
    """Sanity check: the test region must contain at least one read."""
    ids = _collect_fetch_read_ids(str(REGION_BAM), REGION_CHROM, REGION_START, REGION_END)
    assert len(ids) > 0, (
        f"Expected reads in {REGION_CHROM}:{REGION_START}-{REGION_END} "
        f"but got none. Check the test BAM."
    )


# ---------------------------------------------------------------------------
# peak_calling() region parameter — BAI guard
# ---------------------------------------------------------------------------


def test_peak_calling_region_raises_when_bai_missing() -> None:
    """peak_calling(region=...) must raise IOError when .bai is absent."""
    # Patch sys.argv so ema.config's argparse does not fail under pytest
    with patch.object(sys, "argv", ["ema"]):
        from ema.countmatrix.peackcalling import peak_calling

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a fake BAM file (content irrelevant — we check before opening)
        fake_bam = os.path.join(tmpdir, "fake.bam")
        Path(fake_bam).write_bytes(b"FAKEBAM")
        # Deliberately do NOT create fake.bam.bai

        with pytest.raises(IOError, match="samtools index"):
            with patch.object(sys, "argv", ["ema"]):
                peak_calling(
                    direction=False,
                    bedfilepath=os.path.join(tmpdir, "out.bed"),
                    matrixpath=os.path.join(tmpdir, "out.mtx"),
                    bamfile_dir=fake_bam,
                    region=("chr1", 0, 1000),
                )


def test_peak_calling_no_region_skips_bai_check() -> None:
    """Without region=, the BAI guard must NOT run even if .bai is absent.

    This preserves backward compatibility: existing callers that open a
    non-indexed BAM for full-scan iteration must not be broken.
    The test uses a dummy BAM path; we expect pysam to raise on open, NOT
    our IOError, so we just confirm we don't get the IOError from our guard.
    """
    with patch.object(sys, "argv", ["ema"]):
        from ema.countmatrix.peackcalling import peak_calling

    with tempfile.TemporaryDirectory() as tmpdir:
        fake_bam = os.path.join(tmpdir, "fake.bam")
        Path(fake_bam).write_bytes(b"FAKEBAM")

        # Without region= the BAI check is skipped; pysam will raise on open.
        # We care only that the raised exception is NOT our custom IOError.
        with pytest.raises(Exception) as exc_info:
            with patch.object(sys, "argv", ["ema"]):
                peak_calling(
                    direction=False,
                    bedfilepath=os.path.join(tmpdir, "out.bed"),
                    matrixpath=os.path.join(tmpdir, "out.mtx"),
                    bamfile_dir=fake_bam,
                    # region intentionally omitted
                )
        # Our guard IOError contains "samtools index"; pysam errors do not.
        assert "samtools index" not in str(exc_info.value), (
            "BAI guard must not fire when region= is not set"
        )


# ---------------------------------------------------------------------------
# peak_calling() region path — produces output files
# ---------------------------------------------------------------------------


@pytest.mark.skipif(REGION_BAM is None, reason="No indexed test BAM available")
def test_peak_calling_region_creates_output_files() -> None:
    """peak_calling(region=...) must create bed and mtx output files.

    Smoke test: we don't validate peak semantics here, only that the function
    runs to completion without error and writes non-empty-or-existing files.

    sys.argv is patched with the minimum required CLI arguments so that
    variable_config (parsed at module-import time via argparse) has valid
    values for barcode_tag, cb_len, and seqlen.
    """
    _fake_argv = [
        "ema",
        "--bamDir", str(REGION_BAM),
        "--sequenceLen", "150",
        "--CellBarcodeLen", "16",
        "--BarcodeTag", "CB",
    ]

    with patch.object(sys, "argv", _fake_argv):
        import importlib
        import ema.config as _cfg_mod
        import ema.countmatrix.read as _read_mod
        import ema.countmatrix.peackcalling as _pc_mod
        importlib.reload(_cfg_mod)
        importlib.reload(_read_mod)
        importlib.reload(_pc_mod)
        peak_calling = _pc_mod.peak_calling

    with tempfile.TemporaryDirectory() as tmpdir:
        bed_out = os.path.join(tmpdir, "peaks.bed")
        mtx_out = os.path.join(tmpdir, "peaks.mtx")

        peak_calling(
            direction=False,
            bedfilepath=bed_out,
            matrixpath=mtx_out,
            bamfile_dir=str(REGION_BAM),
            region=(REGION_CHROM, REGION_START, REGION_END),
        )

        assert os.path.exists(bed_out), "BED output file must be created"
        assert os.path.exists(mtx_out), "MTX output file must be created"
