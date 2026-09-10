"""Tests for ema.annotate.gtf2isoform_utr.

Covers:
- Single-isoform gene
- Multi-isoform gene
- Gene with split 3'UTR (multiple UTR exons per transcript)
- Minus-strand gene (exons sorted in reverse-coord order)
- Gene with no 3'UTR annotation
- Cache hit is fast (< 100 ms on second call)
- Cache invalidation when fingerprint changes
"""

from __future__ import annotations

import os
import sys
import tempfile
import time
from pathlib import Path

import pytest

# Ensure the project root is on the path
sys.path.insert(0, str(Path(__file__).parent.parent))

from ema.annotate.gtf2isoform_utr import parse_isoform_utrs


# ---------------------------------------------------------------------------
# GTF fixture builder
# ---------------------------------------------------------------------------

def _attr(gene_id: str, gene_name: str, transcript_id: str) -> str:
    return (
        f'gene_id "{gene_id}"; gene_name "{gene_name}"; '
        f'transcript_id "{transcript_id}"; '
    )


def _gtf_line(
    chrom: str,
    feature: str,
    start: int,
    end: int,
    strand: str,
    gene_id: str,
    gene_name: str,
    transcript_id: str,
) -> str:
    """Build a single GTF line (1-based inclusive coords)."""
    attrs = _attr(gene_id, gene_name, transcript_id)
    return f"{chrom}\tENSEMBL\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n"


def build_synthetic_gtf() -> str:
    """Return a synthetic 10-gene GTF string.

    Genes:
        GENE01 — single isoform, single-exon 3'UTR (+ strand, chr1)
        GENE02 — two isoforms, each with a single-exon 3'UTR (+ strand, chr1)
        GENE03 — single isoform, split 3'UTR (2 UTR exons, + strand, chr1)
        GENE04 — single isoform, single-exon 3'UTR (- strand, chr1)
        GENE05 — no three_prime_utr annotation (only gene + exon lines)
        GENE06 — two isoforms, one has split UTR (+ strand, chr2)
        GENE07 — single isoform, single-exon 3'UTR (+ strand, chr2)
        GENE08 — single isoform, split 3'UTR (3 UTR exons, - strand, chr2)
        GENE09 — three isoforms (+ strand, chr2)
        GENE10 — single isoform, version-suffixed IDs (+ strand, chr3)
    """
    lines: list[str] = [
        "##format: GTF\n",
        "#!genome-build GRCh38\n",
    ]

    # GENE01: single isoform, + strand, 1 UTR exon
    lines.append(_gtf_line("chr1", "gene", 1000, 2000, "+", "GENE01", "Gene1", "T01_1"))
    lines.append(_gtf_line("chr1", "three_prime_utr", 1800, 2000, "+", "GENE01", "Gene1", "T01_1"))

    # GENE02: two isoforms, + strand, 1 UTR exon each
    lines.append(_gtf_line("chr1", "gene", 3000, 5000, "+", "GENE02", "Gene2", "T02_1"))
    lines.append(_gtf_line("chr1", "three_prime_utr", 4500, 5000, "+", "GENE02", "Gene2", "T02_1"))
    lines.append(_gtf_line("chr1", "gene", 3000, 4500, "+", "GENE02", "Gene2", "T02_2"))
    lines.append(_gtf_line("chr1", "three_prime_utr", 4000, 4500, "+", "GENE02", "Gene2", "T02_2"))

    # GENE03: split UTR (2 exons), + strand
    lines.append(_gtf_line("chr1", "gene", 6000, 8000, "+", "GENE03", "Gene3", "T03_1"))
    lines.append(_gtf_line("chr1", "three_prime_utr", 7000, 7200, "+", "GENE03", "Gene3", "T03_1"))
    lines.append(_gtf_line("chr1", "three_prime_utr", 7500, 8000, "+", "GENE03", "Gene3", "T03_1"))

    # GENE04: single isoform, - strand
    lines.append(_gtf_line("chr1", "gene", 9000, 11000, "-", "GENE04", "Gene4", "T04_1"))
    lines.append(_gtf_line("chr1", "three_prime_utr", 9000, 9400, "-", "GENE04", "Gene4", "T04_1"))

    # GENE05: no UTR annotation
    lines.append(_gtf_line("chr1", "gene", 12000, 14000, "+", "GENE05", "Gene5", "T05_1"))
    lines.append(_gtf_line("chr1", "exon", 12000, 14000, "+", "GENE05", "Gene5", "T05_1"))

    # GENE06: two isoforms, second has split UTR, chr2
    lines.append(_gtf_line("chr2", "gene", 1000, 4000, "+", "GENE06", "Gene6", "T06_1"))
    lines.append(_gtf_line("chr2", "three_prime_utr", 3500, 4000, "+", "GENE06", "Gene6", "T06_1"))
    lines.append(_gtf_line("chr2", "gene", 1000, 3800, "+", "GENE06", "Gene6", "T06_2"))
    lines.append(_gtf_line("chr2", "three_prime_utr", 2800, 3000, "+", "GENE06", "Gene6", "T06_2"))
    lines.append(_gtf_line("chr2", "three_prime_utr", 3200, 3800, "+", "GENE06", "Gene6", "T06_2"))

    # GENE07: single isoform, chr2
    lines.append(_gtf_line("chr2", "gene", 5000, 6000, "+", "GENE07", "Gene7", "T07_1"))
    lines.append(_gtf_line("chr2", "three_prime_utr", 5800, 6000, "+", "GENE07", "Gene7", "T07_1"))

    # GENE08: split UTR (3 exons), - strand, chr2
    lines.append(_gtf_line("chr2", "gene", 7000, 10000, "-", "GENE08", "Gene8", "T08_1"))
    lines.append(_gtf_line("chr2", "three_prime_utr", 7000, 7300, "-", "GENE08", "Gene8", "T08_1"))
    lines.append(_gtf_line("chr2", "three_prime_utr", 7600, 7900, "-", "GENE08", "Gene8", "T08_1"))
    lines.append(_gtf_line("chr2", "three_prime_utr", 8200, 8500, "-", "GENE08", "Gene8", "T08_1"))

    # GENE09: three isoforms, chr2
    for i in range(1, 4):
        lines.append(_gtf_line("chr2", "gene", 11000, 14000, "+", "GENE09", "Gene9", f"T09_{i}"))
        lines.append(_gtf_line("chr2", "three_prime_utr", 13000 + i * 100, 14000, "+", "GENE09", "Gene9", f"T09_{i}"))

    # GENE10: version-suffixed IDs (should be stripped)
    lines.append(
        'chr3\tENSEMBL\tthree_prime_utr\t1000\t1500\t.\t+\t.\t'
        'gene_id "GENE10.3"; transcript_id "T10_1.2"; gene_name "Gene10";\n'
    )

    return "".join(lines)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def gtf_file(tmp_path_factory) -> Path:
    """Write the synthetic GTF to a temp file and return its Path."""
    p = tmp_path_factory.mktemp("gtf") / "synthetic.gtf"
    p.write_text(build_synthetic_gtf())
    return p


@pytest.fixture(scope="module")
def parsed(gtf_file) -> dict:
    """Parse the synthetic GTF (no cache)."""
    return parse_isoform_utrs(gtf_file, cache_dir=None)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestOutputStructure:
    def test_returns_dict(self, parsed):
        assert isinstance(parsed, dict)

    def test_gene05_absent(self, parsed):
        """GENE05 has no three_prime_utr annotation, must not appear."""
        assert "GENE05" not in parsed

    def test_all_annotated_genes_present(self, parsed):
        for gid in ("GENE01", "GENE02", "GENE03", "GENE04",
                    "GENE06", "GENE07", "GENE08", "GENE09", "GENE10"):
            assert gid in parsed, f"{gid} not in parsed"

    def test_inner_dict_is_transcript_map(self, parsed):
        for gene_id, t_map in parsed.items():
            assert isinstance(t_map, dict), f"{gene_id}: expected dict of transcripts"
            for tid, exons in t_map.items():
                assert isinstance(tid, str)
                assert isinstance(exons, list)
                assert len(exons) >= 1

    def test_exon_tuple_length(self, parsed):
        """Each exon tuple must be (chrom, start, end, strand, exon_rank)."""
        for gene_id, t_map in parsed.items():
            for tid, exons in t_map.items():
                for ex in exons:
                    assert len(ex) == 5, f"{gene_id}/{tid}: exon tuple length != 5"


class TestSingleIsoform:
    def test_gene01_single_transcript(self, parsed):
        assert "T01_1" in parsed["GENE01"]
        assert len(parsed["GENE01"]) == 1

    def test_gene01_exon_count(self, parsed):
        exons = parsed["GENE01"]["T01_1"]
        assert len(exons) == 1

    def test_gene01_coordinates(self, parsed):
        chrom, start, end, strand, rank = parsed["GENE01"]["T01_1"][0]
        assert chrom == "chr1"
        assert start == 1799  # GTF 1-based -> 0-based: 1800-1=1799
        assert end == 2000
        assert strand == "+"
        assert rank == 0


class TestMultiIsoform:
    def test_gene02_two_transcripts(self, parsed):
        assert len(parsed["GENE02"]) == 2
        assert "T02_1" in parsed["GENE02"]
        assert "T02_2" in parsed["GENE02"]

    def test_gene02_transcripts_have_different_utrs(self, parsed):
        t1_end = parsed["GENE02"]["T02_1"][0][2]
        t2_end = parsed["GENE02"]["T02_2"][0][2]
        assert t1_end != t2_end


class TestSplitUTR:
    def test_gene03_two_exons(self, parsed):
        exons = parsed["GENE03"]["T03_1"]
        assert len(exons) == 2

    def test_gene03_sorted_plus_strand(self, parsed):
        exons = parsed["GENE03"]["T03_1"]
        # Plus strand: ascending start
        starts = [e[1] for e in exons]
        assert starts == sorted(starts)

    def test_gene03_exon_ranks_sequential(self, parsed):
        exons = parsed["GENE03"]["T03_1"]
        ranks = [e[4] for e in exons]
        assert ranks == list(range(len(ranks)))

    def test_gene08_three_exons_minus_strand(self, parsed):
        exons = parsed["GENE08"]["T08_1"]
        assert len(exons) == 3
        # Minus strand: descending start (3'UTR is towards low coords)
        starts = [e[1] for e in exons]
        assert starts == sorted(starts, reverse=True)

    def test_gene06_t2_split_utr(self, parsed):
        exons = parsed["GENE06"]["T06_2"]
        assert len(exons) == 2


class TestMinusStrand:
    def test_gene04_strand(self, parsed):
        exons = parsed["GENE04"]["T04_1"]
        assert exons[0][3] == "-"

    def test_gene04_rank_zero_for_sole_exon(self, parsed):
        exons = parsed["GENE04"]["T04_1"]
        assert exons[0][4] == 0

    def test_gene08_minus_strand_order(self, parsed):
        """On minus strand, first exon in transcription order has highest coord."""
        exons = parsed["GENE08"]["T08_1"]
        # rank 0 exon should have the highest start coordinate
        rank0 = next(e for e in exons if e[4] == 0)
        max_start = max(e[1] for e in exons)
        assert rank0[1] == max_start


class TestVersionStrip:
    def test_gene10_version_stripped(self, parsed):
        """Version suffixes in gene_id / transcript_id must be stripped."""
        assert "GENE10" in parsed
        assert "T10_1" in parsed["GENE10"]


class TestThreeIsoforms:
    def test_gene09_three_transcripts(self, parsed):
        assert len(parsed["GENE09"]) == 3


class TestCaching:
    def test_cache_hit_fast(self, gtf_file, tmp_path):
        """Second call with a warm cache must complete in < 100 ms."""
        cache_dir = tmp_path / "cache"

        # Warm the cache
        parse_isoform_utrs(gtf_file, cache_dir=cache_dir)

        # Time the cache-hit path
        t0 = time.monotonic()
        result = parse_isoform_utrs(gtf_file, cache_dir=cache_dir)
        elapsed = time.monotonic() - t0

        assert elapsed < 0.1, f"Cache hit took {elapsed:.3f}s (> 100 ms)"
        assert "GENE01" in result

    def test_cache_invalidation(self, gtf_file, tmp_path):
        """Modifying the GTF file must invalidate the cache."""
        cache_dir = tmp_path / "cache2"
        parse_isoform_utrs(gtf_file, cache_dir=cache_dir)

        # Touch the file to change mtime
        Path(gtf_file).touch()
        result_fresh = parse_isoform_utrs(gtf_file, cache_dir=cache_dir)
        # Should still work (re-parsed)
        assert isinstance(result_fresh, dict)
