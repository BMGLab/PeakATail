"""Tests for the gene-model region loader in ema.viz._gene_track_helpers.

Covers:
- Protein-coding transcript: 5'UTR + CDS + 3'UTR all extracted, has_typed_regions True
- Minus-strand transcript with a split 3'UTR (two UTR exons)
- Non-coding transcript (exon lines only, no CDS/UTR): has_typed_regions False
- Alternate UTR feature-type spelling ("five_prime_UTR" / "three_prime_UTR")
- Multiple transcripts of the same gene resolved independently
- Missing/unreadable GTF returns []
- Back-compat: load_isoforms_for_gene keeps working unchanged
- pas_region_labels(): geometric PAS<->UTR overlap fallback
- pas_distance_table(): falls back to region-overlap UTR labels when
  panel.pas_isoforms wasn't supplied
- format_count(): human-readable k/M count formatting
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from ema.viz._gene_track_helpers import (
    GenePanel,
    TranscriptRegions,
    format_count,
    load_isoform_regions_for_gene,
    load_isoforms_for_gene,
    pas_distance_table,
    pas_region_labels,
)


def _minimal_panel(**overrides) -> GenePanel:
    """Build a minimal-but-valid GenePanel for helper-function tests."""
    defaults = dict(
        gene_id="GENE01",
        chrom="1",
        start=1000,
        end=5000,
        strand="+",
        pas_ids=[1, 2],
        pas_positions=[1100, 4900],
        pas_starts=[1050, 4850],
        pas_ends=[1150, 4950],
        clusters=["A"],
        n_cells_per_cluster=[10],
        reads=np.array([[5.0, 3.0]]),
        reads_per_cell=np.array([[0.5, 0.3]]),
        proportions=np.array([[0.625, 0.375]]),
    )
    defaults.update(overrides)
    return GenePanel(**defaults)


def _attr(gene_id: str, transcript_id: str) -> str:
    return f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; '


def _line(feature: str, start: int, end: int, strand: str, gene_id: str, transcript_id: str,
          chrom: str = "chr1") -> str:
    attrs = _attr(gene_id, transcript_id)
    return f"{chrom}\tENSEMBL\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n"


def _write_gtf(lines: list[str]) -> Path:
    fh = tempfile.NamedTemporaryFile(mode="w", suffix=".gtf", delete=False)
    fh.writelines(lines)
    fh.close()
    return Path(fh.name)


def test_protein_coding_transcript_full_gene_model():
    """5'UTR + CDS + 3'UTR all land in the right bucket, exons too."""
    lines = [
        _line("gene", 1000, 5000, "+", "GENE01", "T1"),
        _line("exon", 1000, 2000, "+", "GENE01", "T1"),
        _line("exon", 4000, 5000, "+", "GENE01", "T1"),
        _line("five_prime_utr", 1000, 1200, "+", "GENE01", "T1"),
        _line("CDS", 1200, 2000, "+", "GENE01", "T1"),
        _line("CDS", 4000, 4800, "+", "GENE01", "T1"),
        _line("three_prime_utr", 4800, 5000, "+", "GENE01", "T1"),
    ]
    gtf = _write_gtf(lines)
    try:
        regions = load_isoform_regions_for_gene(gtf, "GENE01")
    finally:
        gtf.unlink()

    assert len(regions) == 1
    rec = regions[0]
    assert isinstance(rec, TranscriptRegions)
    assert rec.transcript_id == "T1"
    assert rec.exons == [(999, 2000), (3999, 5000)]
    assert rec.utr5 == [(999, 1200)]
    assert rec.cds == [(1199, 2000), (3999, 4800)]
    assert rec.utr3 == [(4799, 5000)]
    assert rec.has_typed_regions is True


def test_minus_strand_split_3utr():
    """Two three_prime_utr exons on a minus-strand transcript both land in utr3."""
    lines = [
        _line("gene", 1000, 5000, "-", "GENE04", "T4"),
        _line("exon", 1000, 1500, "-", "GENE04", "T4"),
        _line("exon", 4500, 5000, "-", "GENE04", "T4"),
        _line("three_prime_utr", 1000, 1300, "-", "GENE04", "T4"),
        _line("three_prime_utr", 4500, 4600, "-", "GENE04", "T4"),
        _line("CDS", 1300, 1500, "-", "GENE04", "T4"),
        _line("CDS", 4600, 5000, "-", "GENE04", "T4"),
    ]
    gtf = _write_gtf(lines)
    try:
        regions = load_isoform_regions_for_gene(gtf, "GENE04")
    finally:
        gtf.unlink()

    assert len(regions) == 1
    rec = regions[0]
    assert len(rec.utr3) == 2
    assert rec.utr3 == sorted(rec.utr3)
    assert rec.has_typed_regions is True


def test_noncoding_transcript_falls_back_to_exon_only():
    """A transcript with only exon lines has no typed regions."""
    lines = [
        _line("gene", 6000, 8000, "+", "GENE05", "T5"),
        _line("exon", 6000, 6500, "+", "GENE05", "T5"),
        _line("exon", 7500, 8000, "+", "GENE05", "T5"),
    ]
    gtf = _write_gtf(lines)
    try:
        regions = load_isoform_regions_for_gene(gtf, "GENE05")
    finally:
        gtf.unlink()

    assert len(regions) == 1
    rec = regions[0]
    assert rec.exons == [(5999, 6500), (7499, 8000)]
    assert rec.cds == []
    assert rec.utr5 == []
    assert rec.utr3 == []
    assert rec.has_typed_regions is False


def test_alternate_utr_feature_spelling():
    """'five_prime_UTR' / 'three_prime_UTR' (upper-case UTR) are recognised too."""
    lines = [
        _line("gene", 100, 900, "+", "GENE06", "T6"),
        _line("exon", 100, 900, "+", "GENE06", "T6"),
        _line("five_prime_UTR", 100, 200, "+", "GENE06", "T6"),
        _line("CDS", 200, 800, "+", "GENE06", "T6"),
        _line("three_prime_UTR", 800, 900, "+", "GENE06", "T6"),
    ]
    gtf = _write_gtf(lines)
    try:
        regions = load_isoform_regions_for_gene(gtf, "GENE06")
    finally:
        gtf.unlink()

    assert len(regions) == 1
    rec = regions[0]
    assert rec.utr5 == [(99, 200)]
    assert rec.utr3 == [(799, 900)]
    assert rec.has_typed_regions is True


def test_multiple_transcripts_resolved_independently():
    lines = [
        _line("gene", 100, 2000, "+", "GENE02", "T2_1"),
        _line("exon", 100, 600, "+", "GENE02", "T2_1"),
        _line("CDS", 300, 600, "+", "GENE02", "T2_1"),
        _line("five_prime_utr", 100, 300, "+", "GENE02", "T2_1"),
        _line("three_prime_utr", 600, 600, "+", "GENE02", "T2_1"),  # zero-width edge case OK
        _line("exon", 1200, 2000, "+", "GENE02", "T2_2"),
        _line("CDS", 1200, 1800, "+", "GENE02", "T2_2"),
        _line("three_prime_utr", 1800, 2000, "+", "GENE02", "T2_2"),
    ]
    gtf = _write_gtf(lines)
    try:
        regions = load_isoform_regions_for_gene(gtf, "GENE02")
    finally:
        gtf.unlink()

    by_tid = {r.transcript_id: r for r in regions}
    assert set(by_tid) == {"T2_1", "T2_2"}
    assert by_tid["T2_1"].utr5 == [(99, 300)]
    assert by_tid["T2_2"].utr3 == [(1799, 2000)]
    assert by_tid["T2_2"].utr5 == []


def test_missing_gtf_returns_empty_list():
    regions = load_isoform_regions_for_gene(Path("/nonexistent/path/to.gtf"), "GENE01")
    assert regions == []


def test_gene_not_in_gtf_returns_empty_list():
    lines = [
        _line("gene", 1000, 5000, "+", "GENE01", "T1"),
        _line("exon", 1000, 5000, "+", "GENE01", "T1"),
    ]
    gtf = _write_gtf(lines)
    try:
        regions = load_isoform_regions_for_gene(gtf, "GENE_NOT_PRESENT")
    finally:
        gtf.unlink()
    assert regions == []


def test_back_compat_load_isoforms_for_gene_unchanged():
    """The legacy exon-only loader keeps working exactly as before."""
    lines = [
        _line("gene", 1000, 5000, "+", "GENE01", "T1"),
        _line("exon", 1000, 2000, "+", "GENE01", "T1"),
        _line("exon", 4000, 5000, "+", "GENE01", "T1"),
        _line("CDS", 1200, 2000, "+", "GENE01", "T1"),
    ]
    gtf = _write_gtf(lines)
    try:
        isoforms = load_isoforms_for_gene(gtf, "GENE01")
    finally:
        gtf.unlink()

    assert isoforms == [("T1", [(999, 2000), (3999, 5000)])]


# ---------------------------------------------------------------------------
# pas_region_labels / pas_distance_table fallback
# ---------------------------------------------------------------------------

def test_pas_region_labels_finds_3utr_and_5utr_overlap():
    """PAS 1 overlaps T1's 3'UTR; PAS 2 overlaps T1's 5'UTR."""
    regions = TranscriptRegions(
        transcript_id="T1",
        exons=[(1000, 1200), (4800, 5000)],
        cds=[(1200, 1800)],
        utr5=[(1000, 1200)],
        utr3=[(4800, 5000)],
    )
    panel = _minimal_panel(
        pas_ids=[1, 2],
        pas_positions=[4900, 1050],
        pas_starts=[4850, 1000],
        pas_ends=[4950, 1100],
        isoform_regions=[regions],
    )
    labels = pas_region_labels(panel)
    assert labels[1] == "T1 (3'UTR)"
    assert labels[2] == "T1 (5'UTR)"


def test_pas_region_labels_no_overlap_omitted():
    regions = TranscriptRegions(transcript_id="T1", utr3=[(4800, 5000)])
    panel = _minimal_panel(
        pas_ids=[1],
        pas_positions=[2500],
        pas_starts=[2450],
        pas_ends=[2550],
        isoform_regions=[regions],
    )
    assert pas_region_labels(panel) == {}


def test_pas_region_labels_empty_without_isoform_regions():
    panel = _minimal_panel()  # isoform_regions defaults to []
    assert pas_region_labels(panel) == {}


def test_pas_distance_table_falls_back_to_region_overlap():
    """No pas_isoforms supplied -> utr_transcripts column is still populated
    from isoform_regions geometry instead of being absent."""
    regions = TranscriptRegions(transcript_id="T1", utr3=[(4800, 5000)])
    panel = _minimal_panel(isoform_regions=[regions])
    assert not panel.pas_isoforms

    df = pas_distance_table(panel)
    assert "utr_transcripts" in df.columns
    # PAS 2 (4850-4950) overlaps T1's 3'UTR (4800-5000); PAS 1 (1050-1150) doesn't.
    row = df[df["pas_id"] == 2].iloc[0]
    assert row["utr_transcripts"] == "T1 (3'UTR)"
    other = df[df["pas_id"] == 1].iloc[0]
    assert other["utr_transcripts"] == ""


def test_pas_distance_table_prefers_explicit_pas_isoforms_over_fallback():
    """When panel.pas_isoforms IS supplied, it wins over the region fallback
    (the per-isoform quantification is authoritative)."""
    regions = TranscriptRegions(transcript_id="T1", utr3=[(4800, 5000)])
    panel = _minimal_panel(
        isoform_regions=[regions],
        pas_isoforms={1: ["T_explicit"]},
    )
    df = pas_distance_table(panel)
    row = df[df["pas_id"] == 1].iloc[0]
    assert row["utr_transcripts"] == "T_explicit"


def test_pas_distance_table_no_columns_without_any_source():
    panel = _minimal_panel()
    df = pas_distance_table(panel)
    assert "utr_transcripts" not in df.columns


# ---------------------------------------------------------------------------
# format_count
# ---------------------------------------------------------------------------

def test_format_count_thousands_and_millions():
    assert format_count(12345) == "12.3k"
    assert format_count(1234567) == "1.23M"


def test_format_count_small_values():
    assert format_count(235.8) == "235.8"
    assert format_count(0) == "0.0"


def test_format_count_sub_one_keeps_precision():
    # reads/cell is frequently < 1; must not collapse to "0.0".
    assert format_count(0.0234) == "0.023"


def test_format_count_non_finite_and_none():
    assert format_count(None) == "—"
    assert format_count(float("nan")) == "—"
    assert format_count(float("inf")) == "—"


def test_format_count_negative():
    assert format_count(-12345) == "-12.3k"


# ---------------------------------------------------------------------------
# Gene ranking must not treat a withheld q-value as maximal significance.
# ---------------------------------------------------------------------------

def test_withheld_qvalue_does_not_rank_first(monkeypatch) -> None:
    """A q withheld for a dispersion-floored test (issue #94) scored 300.

    `q.where(q > 0, 1e-300)` swallowed NaN too (NaN > 0 is False), so the
    subsequent `.fillna(1.0)` never fired and -log10(1e-300) = 300 put the
    untrustworthy gene above every genuine hit in the auto gene panels and
    `switch geneview --top-genes`.
    """
    import numpy as np
    import pandas as pd

    from ema.viz._gene_track_helpers import rank_top_genes

    df = pd.DataFrame({
        "gene_id": ["G_withheld", "G_big_fc", "G_ok"],
        "qvalue": [np.nan, 0.9, 0.5],
        "log2fc": [1.0, 5.0, 2.0],
    })
    ranked = rank_top_genes({("A", "B"): df}, n=3)
    assert ranked, "fixture produced no ranking"
    assert ranked[0] != "G_withheld", (
        f"a gene whose q-value was WITHHELD ranked first ({ranked!r}); a "
        "missing q must contribute no evidence, not infinite evidence."
    )
    assert ranked[-1] == "G_withheld", (
        f"expected the withheld gene to score lowest, got {ranked!r}"
    )
