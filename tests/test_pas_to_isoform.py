"""Tests for ema.quantification.pas_to_isoform.

Covers:
- Plus-strand rank assignment (proximal rank < distal rank)
- Minus-strand rank assignment
- Multi-isoform PAS: each PAS maps to all overlapping transcripts
- Transcript position calculation (single UTR exon)
- Transcript position calculation (split UTR exon)
- PAS outside UTR region → absent from output
- total_pas_in_transcript correctness
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd

from ema.quantification.pas_to_isoform import (
    map_pas_to_isoforms,
    rank_pas_by_genomic_position,
)


# ---------------------------------------------------------------------------
# Helpers to build synthetic data
# ---------------------------------------------------------------------------

def _make_bed6(records: list[tuple]) -> str:
    """Build BED6 text from (chrom, start, end, pas_id, score, strand) tuples."""
    return "".join(
        f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\n"
        for r in records
    )


def write_bed_file(tmp_path: Path, records: list[tuple]) -> Path:
    p = tmp_path / "pas.bed"
    p.write_text(_make_bed6(records))
    return p


# ---------------------------------------------------------------------------
# Synthetic isoform UTR maps
# ---------------------------------------------------------------------------

# Format: {gene_id: {transcript_id: [(chrom, start, end, strand, exon_rank), ...]}}
# Coordinates are 0-based half-open (BED convention).


def simple_plus_isoform_utrs() -> dict:
    """One gene, one transcript, two UTR exons, + strand.

    UTR exon 1: chr1:1000-1200  (rank 0, proximal)
    UTR exon 2: chr1:1500-1800  (rank 1, distal)
    """
    return {
        "GENE_A": {
            "TRANSCRIPT_A1": [
                ("chr1", 1000, 1200, "+", 0),
                ("chr1", 1500, 1800, "+", 1),
            ]
        }
    }


def simple_minus_isoform_utrs() -> dict:
    """One gene, one transcript, two UTR exons, - strand.

    On - strand, 3'UTR is towards lower coords.
    UTR exon at lower coord → rank 0 (proximal).
    UTR exon at higher coord → rank 1 (distal — farther from stop codon on -).

    Actually for - strand in our parser: exons sorted descending by start.
    So rank 0 = highest start coord (closest to CDS on minus strand = proximal).
         rank 1 = lower start coord (farthest from CDS = distal).
    """
    return {
        "GENE_B": {
            "TRANSCRIPT_B1": [
                ("chr1", 5000, 5300, "-", 0),   # rank 0: proximal (high coord)
                ("chr1", 4600, 4900, "-", 1),   # rank 1: distal  (low coord)
            ]
        }
    }


def multi_isoform_utrs() -> dict:
    """One gene, two transcripts that share one PAS.

    TRANSCRIPT_C1: chr2:2000-2500 (rank 0 only — one PAS per transcript)
    TRANSCRIPT_C2: chr2:2000-2500 (rank 0) + chr2:3000-3500 (rank 1)

    PAS at chr2:2200 should appear in BOTH transcripts.
    PAS at chr2:3200 should appear only in TRANSCRIPT_C2.
    """
    return {
        "GENE_C": {
            "TRANSCRIPT_C1": [
                ("chr2", 2000, 2500, "+", 0),
            ],
            "TRANSCRIPT_C2": [
                ("chr2", 2000, 2500, "+", 0),
                ("chr2", 3000, 3500, "+", 1),
            ],
        }
    }


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_path_local(tmp_path):
    return tmp_path


# ---------------------------------------------------------------------------
# Plus-strand tests
# ---------------------------------------------------------------------------

class TestPlusStrand:
    def test_single_pas_in_first_exon(self, tmp_path):
        """PAS inside exon rank 0 → transcript_pos = offset within exon."""
        isoform_utrs = simple_plus_isoform_utrs()
        # PAS at chr1:1100-1102 → mid = 1101, inside exon [1000, 1200)
        # transcript_pos = 1101 - 1000 = 101
        bed = write_bed_file(tmp_path, [
            ("chr1", 1100, 1102, 1, 0, "+"),
        ])
        result = map_pas_to_isoforms(bed, isoform_utrs)
        assert 1 in result
        entries = result[1]
        assert len(entries) == 1
        gene_id, tid, tpos, rank, total = entries[0]
        assert gene_id == "GENE_A"
        assert tid == "TRANSCRIPT_A1"
        assert tpos == 101  # mid=1101, offset within [1000,1200)

    def test_single_pas_in_second_exon(self, tmp_path):
        """PAS inside exon rank 1 → transcript_pos = exon0_len + offset."""
        isoform_utrs = simple_plus_isoform_utrs()
        # PAS at chr1:1600-1602 → mid = 1601, inside exon [1500, 1800)
        # transcript_pos = len(exon0) + (1601 - 1500) = 200 + 101 = 301
        bed = write_bed_file(tmp_path, [
            ("chr1", 1600, 1602, 2, 0, "+"),
        ])
        result = map_pas_to_isoforms(bed, isoform_utrs)
        assert 2 in result
        gene_id, tid, tpos, rank, total = result[2][0]
        assert tpos == 301  # 200 (exon0) + 101 (offset in exon1)

    def test_rank_proximal_before_distal(self, tmp_path):
        """PAS closer to CDS has lower rank than PAS farther from CDS."""
        isoform_utrs = simple_plus_isoform_utrs()
        # pas 10 in exon 0 (proximal) → rank 0
        # pas 11 in exon 1 (distal)   → rank 1
        bed = write_bed_file(tmp_path, [
            ("chr1", 1050, 1052, 10, 0, "+"),  # mid=1051, in exon 0
            ("chr1", 1550, 1552, 11, 0, "+"),  # mid=1551, in exon 1
        ])
        result = map_pas_to_isoforms(bed, isoform_utrs)
        rank_prox = result[10][0][3]
        rank_dist = result[11][0][3]
        assert rank_prox < rank_dist
        assert rank_prox == 0
        assert rank_dist == 1

    def test_total_pas_count(self, tmp_path):
        """total_pas_in_transcript reflects total PAS in that transcript."""
        isoform_utrs = simple_plus_isoform_utrs()
        bed = write_bed_file(tmp_path, [
            ("chr1", 1050, 1052, 10, 0, "+"),
            ("chr1", 1550, 1552, 11, 0, "+"),
        ])
        result = map_pas_to_isoforms(bed, isoform_utrs)
        total_10 = result[10][0][4]
        total_11 = result[11][0][4]
        assert total_10 == 2
        assert total_11 == 2


# ---------------------------------------------------------------------------
# Minus-strand tests
# ---------------------------------------------------------------------------

class TestMinusStrand:
    def test_proximal_has_lower_rank_minus(self, tmp_path):
        """On minus strand, PAS with higher genomic coord is proximal (rank 0)."""
        isoform_utrs = simple_minus_isoform_utrs()
        # pas 20: mid=5100, in exon rank 0 [5000, 5300) → proximal
        # pas 21: mid=4700, in exon rank 1 [4600, 4900) → distal
        bed = write_bed_file(tmp_path, [
            ("chr1", 5098, 5102, 20, 0, "-"),
            ("chr1", 4698, 4702, 21, 0, "-"),
        ])
        result = map_pas_to_isoforms(bed, isoform_utrs)
        assert 20 in result, "PAS 20 should intersect UTR"
        assert 21 in result, "PAS 21 should intersect UTR"
        rank20 = result[20][0][3]
        rank21 = result[21][0][3]
        assert rank20 == 0
        assert rank21 == 1
        assert rank20 < rank21

    def test_minus_transcript_pos_first_exon(self, tmp_path):
        """On minus strand, transcript_pos in rank-0 exon = offset from 3' end."""
        isoform_utrs = simple_minus_isoform_utrs()
        # pas at mid=5100: in exon [5000, 5300), strand="-"
        # offset = ex_end - mid - 1 = 5300 - 5100 - 1 = 199
        bed = write_bed_file(tmp_path, [
            ("chr1", 5098, 5102, 30, 0, "-"),
        ])
        result = map_pas_to_isoforms(bed, isoform_utrs)
        tpos = result[30][0][2]
        assert tpos == 199


# ---------------------------------------------------------------------------
# Multi-isoform tests
# ---------------------------------------------------------------------------

class TestMultiIsoform:
    def test_shared_pas_maps_to_both_transcripts(self, tmp_path):
        """PAS at chr2:2200 overlaps both isoforms of GENE_C."""
        isoform_utrs = multi_isoform_utrs()
        bed = write_bed_file(tmp_path, [
            ("chr2", 2198, 2202, 40, 0, "+"),  # shared PAS
        ])
        result = map_pas_to_isoforms(bed, isoform_utrs)
        assert 40 in result
        entries = result[40]
        tids = {e[1] for e in entries}
        assert "TRANSCRIPT_C1" in tids
        assert "TRANSCRIPT_C2" in tids
        assert len(entries) == 2

    def test_distal_pas_only_in_longer_transcript(self, tmp_path):
        """PAS at chr2:3200 is in TRANSCRIPT_C2 only (TRANSCRIPT_C1 ends at 2500)."""
        isoform_utrs = multi_isoform_utrs()
        bed = write_bed_file(tmp_path, [
            ("chr2", 3198, 3202, 41, 0, "+"),
        ])
        result = map_pas_to_isoforms(bed, isoform_utrs)
        assert 41 in result
        tids = {e[1] for e in result[41]}
        assert "TRANSCRIPT_C1" not in tids
        assert "TRANSCRIPT_C2" in tids

    def test_rank_in_longer_transcript(self, tmp_path):
        """In TRANSCRIPT_C2, shared PAS rank=0 and distal PAS rank=1."""
        isoform_utrs = multi_isoform_utrs()
        bed = write_bed_file(tmp_path, [
            ("chr2", 2198, 2202, 50, 0, "+"),  # proximal
            ("chr2", 3198, 3202, 51, 0, "+"),  # distal
        ])
        result = map_pas_to_isoforms(bed, isoform_utrs)
        # Within TRANSCRIPT_C2
        c2_entries_50 = [e for e in result[50] if e[1] == "TRANSCRIPT_C2"]
        c2_entries_51 = [e for e in result[51] if e[1] == "TRANSCRIPT_C2"]
        assert c2_entries_50, "pas 50 should map to TRANSCRIPT_C2"
        assert c2_entries_51, "pas 51 should map to TRANSCRIPT_C2"
        rank_50 = c2_entries_50[0][3]
        rank_51 = c2_entries_51[0][3]
        assert rank_50 < rank_51


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_pas_outside_utr_absent(self, tmp_path):
        """PAS that does not overlap any UTR must not appear in the result."""
        isoform_utrs = simple_plus_isoform_utrs()
        bed = write_bed_file(tmp_path, [
            ("chr1", 500, 502, 99, 0, "+"),  # before any UTR
        ])
        result = map_pas_to_isoforms(bed, isoform_utrs)
        assert 99 not in result

    def test_empty_isoform_map_returns_empty(self, tmp_path):
        bed = write_bed_file(tmp_path, [
            ("chr1", 1000, 1002, 1, 0, "+"),
        ])
        result = map_pas_to_isoforms(bed, {})
        assert result == {}

    def test_nonexistent_pasbed_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            map_pas_to_isoforms(tmp_path / "nonexistent.bed", {})


# ---------------------------------------------------------------------------
# rank_pas_by_genomic_position — used by the run_length per_gene synthesis
# path, which has no transcript/UTR model, only raw pasbed coordinates.
#
# Regression coverage for the previously-hardcoded (gene_id, "_gene_", 0, 1,
# 1) sentinel: every PAS of every gene reported rank=1/total=1, which is
# indistinguishable from a real single-PAS gene and silently broke
# classic's per_gene proximal/distal selection and the proportion
# strategy's rank column.
# ---------------------------------------------------------------------------

class TestRankPasByGenomicPosition:
    def test_plus_strand_ascending_coordinate_is_distal(self):
        """+ strand: distal PAS = HIGHER genomic coordinate -> rank N."""
        df = pd.DataFrame({
            "pas_id": [10, 11, 12],
            "gene_id": ["G1", "G1", "G1"],
            "start": [500, 100, 900],   # unsorted on purpose
            "strand": ["+", "+", "+"],
        })
        out = rank_pas_by_genomic_position(df).set_index("pas_id")
        assert out.loc[11, "rank"] == 1   # start=100 -> most proximal
        assert out.loc[10, "rank"] == 2   # start=500 -> middle
        assert out.loc[12, "rank"] == 3   # start=900 -> most distal
        assert (out["total"] == 3).all()

    def test_minus_strand_descending_coordinate_is_distal(self):
        """- strand: distal PAS = LOWER genomic coordinate -> rank N."""
        df = pd.DataFrame({
            "pas_id": [20, 21, 22],
            "gene_id": ["G2", "G2", "G2"],
            "start": [500, 100, 900],
            "strand": ["-", "-", "-"],
        })
        out = rank_pas_by_genomic_position(df).set_index("pas_id")
        assert out.loc[22, "rank"] == 1   # start=900 -> most proximal (- strand)
        assert out.loc[20, "rank"] == 2   # start=500 -> middle
        assert out.loc[21, "rank"] == 3   # start=100 -> most distal (- strand)
        assert (out["total"] == 3).all()

    def test_1_based_1_to_n(self):
        """Ranks within a gene are exactly {1, ..., N} — 1-based, unlike
        the bedtools/UTR path's 0-based ``_assign_ranks`` (documented
        divergence; the two paths never populate the same map)."""
        df = pd.DataFrame({
            "pas_id": [1, 2, 3, 4],
            "gene_id": ["G1"] * 4,
            "start": [10, 40, 20, 30],
            "strand": ["+"] * 4,
        })
        out = rank_pas_by_genomic_position(df)
        assert sorted(out["rank"].tolist()) == [1, 2, 3, 4]

    def test_single_pas_gene_rank_1_total_1(self):
        df = pd.DataFrame({
            "pas_id": [1], "gene_id": ["G1"], "start": [100], "strand": ["+"],
        })
        out = rank_pas_by_genomic_position(df)
        assert out.iloc[0]["rank"] == 1
        assert out.iloc[0]["total"] == 1

    def test_multiple_genes_ranked_independently(self):
        df = pd.DataFrame({
            "pas_id": [1, 2, 3, 4],
            "gene_id": ["G1", "G1", "G2", "G2"],
            "start": [100, 200, 500, 400],
            "strand": ["+", "+", "-", "-"],
        })
        out = rank_pas_by_genomic_position(df).set_index("pas_id")
        assert out.loc[1, "rank"] == 1 and out.loc[2, "rank"] == 2
        assert out.loc[1, "total"] == 2 and out.loc[2, "total"] == 2
        # G2 is '-' strand: start=500 (higher) is proximal -> rank 1
        assert out.loc[3, "rank"] == 1 and out.loc[4, "rank"] == 2
        assert out.loc[3, "total"] == 2 and out.loc[4, "total"] == 2

    def test_missing_coordinates_do_not_raise_and_sort_last(self):
        """A PAS with no known genomic coordinate must not crash the
        ranking — it should still get a valid rank (sorted after the
        coordinate-known PAS of its gene) rather than raising."""
        df = pd.DataFrame({
            "pas_id": [1, 2, 3],
            "gene_id": ["G1", "G1", "G1"],
            "start": [100.0, None, 50.0],
            "strand": ["+", None, "+"],
        })
        out = rank_pas_by_genomic_position(df).set_index("pas_id")
        assert out.loc[3, "rank"] == 1   # start=50, known coord, most proximal
        assert out.loc[1, "rank"] == 2   # start=100, known coord
        assert out.loc[2, "rank"] == 3   # unknown coord -> sorts last
        assert (out["total"] == 3).all()

    def test_empty_input(self):
        df = pd.DataFrame({"pas_id": [], "gene_id": [], "start": [], "strand": []})
        out = rank_pas_by_genomic_position(df)
        assert out.empty
        assert "rank" in out.columns and "total" in out.columns
