"""Regression test for issue #99 -- PAS -> gene assignment in overlapping loci.

Where a readthrough/spanning gene model covers a neighbour's 3'UTR, every PAS
in that UTR is distance 0 from BOTH models.  ``bedtools closest -t first``
then hands them all to whichever model sorts first -- always the spanning one,
because it starts further upstream -- and the real gene ends up with zero
assigned PAS.  On the 17-sample Laughney cohort that is how CD68 (inside the
SENP3-EIF4A1 readthrough) got 0 PAS cohort-wide, and 33% of the top-30
replicated "gene switches" carried the wrong gene name.

The fixtures below are synthetic but keep the real topology: a compact gene
whose 3'UTR hosts the PAS, wholly contained in a much larger model, on each
strand.  The correct answer is unambiguous -- the PAS sits 42 bp / 200 bp from
the compact gene's annotated 3' end and >7 kb from the spanning model's.
"""
import pytest

pybedtools = pytest.importorskip("pybedtools")

from ema.annotate.find_close import TIER_1, find_close
from ema.annotate.gtftobed import gtf_bed


# --- the CD68 / SENP3-EIF4A1 topology, minus strand -------------------------
NEG_SPAN_GENE, NEG_SPAN_START, NEG_SPAN_END = "ENSG00000000001", 7_572_644, 7_590_881
NEG_REAL_GENE, NEG_REAL_START, NEG_REAL_END = "ENSG00000000002", 7_579_708, 7_582_507
NEG_PAS_POS = 7_579_750  # inside the compact gene's 3'UTR, 42 bp from its end

# --- the same shape on the plus strand --------------------------------------
POS_SPAN_GENE, POS_SPAN_START, POS_SPAN_END = "ENSG00000000003", 990_000, 1_020_000
POS_REAL_GENE, POS_REAL_START, POS_REAL_END = "ENSG00000000004", 1_000_000, 1_003_000
POS_PAS_POS = 1_002_800  # inside the compact gene's 3'UTR, 200 bp from its end

CHROM = "17"


def _gtf_row(start, end, strand, feature, gene_id, gene_name):
    attrs = (f'gene_id "{gene_id}"; gene_name "{gene_name}"; '
             f'gene_biotype "protein_coding";')
    return f"{CHROM}\tema_test\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n"


@pytest.fixture
def overlapping_locus(tmp_path):
    """Build gene_end.bed + utr_lengths from a synthetic overlapping-loci GTF."""
    gtf = tmp_path / "test.gtf"
    with open(gtf, "w") as fh:
        # Minus strand: the compact gene lives inside the spanning model.
        fh.write(_gtf_row(NEG_SPAN_START, NEG_SPAN_END, "-", "gene",
                          NEG_SPAN_GENE, "SPAN_NEG"))
        fh.write(_gtf_row(NEG_SPAN_START, NEG_SPAN_START + 400, "-",
                          "three_prime_utr", NEG_SPAN_GENE, "SPAN_NEG"))
        fh.write(_gtf_row(NEG_REAL_START, NEG_REAL_END, "-", "gene",
                          NEG_REAL_GENE, "REAL_NEG"))
        fh.write(_gtf_row(NEG_REAL_START, NEG_REAL_START + 492, "-",
                          "three_prime_utr", NEG_REAL_GENE, "REAL_NEG"))
        # Plus strand: same shape, mirrored.
        fh.write(_gtf_row(POS_SPAN_START, POS_SPAN_END, "+", "gene",
                          POS_SPAN_GENE, "SPAN_POS"))
        fh.write(_gtf_row(POS_SPAN_END - 1000, POS_SPAN_END, "+",
                          "three_prime_utr", POS_SPAN_GENE, "SPAN_POS"))
        fh.write(_gtf_row(POS_REAL_START, POS_REAL_END, "+", "gene",
                          POS_REAL_GENE, "REAL_POS"))
        fh.write(_gtf_row(POS_REAL_END - 500, POS_REAL_END, "+",
                          "three_prime_utr", POS_REAL_GENE, "REAL_POS"))

    endbed = tmp_path / "gene_end.bed"
    utr_lengths = gtf_bed(
        endbeddir=str(endbed),
        gtfdir=str(gtf),
        featuresdir=str(tmp_path / "features.tsv"),
        utr_lengths_dir=str(tmp_path / "utr_lengths.tsv"),
    )
    return endbed, utr_lengths


def _run_find_close(tmp_path, endbed, utr_lengths, pos_rows, neg_rows):
    posbed = tmp_path / "pos.bed"
    negbed = tmp_path / "neg.bed"
    posbed.write_text("".join(pos_rows))
    negbed.write_text("".join(neg_rows))

    return find_close(
        posbed_dir=str(posbed),
        negbed_dir=str(negbed),
        genomebed_dir=str(endbed),
        annotatedbed_dir=str(tmp_path / "annotated.bed"),
        mergebed=str(tmp_path / "merged.bed"),
        utr_lengths=utr_lengths,
    )


def test_pas_in_a_spanned_gene_utr_is_assigned_to_that_gene(tmp_path,
                                                            overlapping_locus):
    """The PAS goes to the gene whose 3'UTR holds it, not the spanning model."""
    endbed, utr_lengths = overlapping_locus
    genes = _run_find_close(
        tmp_path, endbed, utr_lengths,
        pos_rows=[f"{CHROM}\t{POS_PAS_POS}\t{POS_PAS_POS + 1}\tpas_pos\t10\t+\n"],
        neg_rows=[f"{CHROM}\t{NEG_PAS_POS}\t{NEG_PAS_POS + 1}\tpas_neg\t10\t-\n"],
    )

    assert set(genes.index) == {"pas_pos", "pas_neg"}
    assert genes.loc["pas_neg", "gene_id"] == NEG_REAL_GENE, (
        "PAS in the compact gene's 3'UTR was handed to the spanning model")
    assert genes.loc["pas_pos", "gene_id"] == POS_REAL_GENE, (
        "PAS in the compact gene's 3'UTR was handed to the spanning model")
    assert genes.loc["pas_neg", "tier"] == TIER_1
    assert genes.loc["pas_pos", "tier"] == TIER_1


def test_spanned_gene_is_not_starved_across_many_pas(tmp_path,
                                                     overlapping_locus):
    """The cohort-wide symptom: the real gene must not end up with zero PAS."""
    endbed, utr_lengths = overlapping_locus
    neg_rows = [
        f"{CHROM}\t{NEG_REAL_START + offset}\t{NEG_REAL_START + offset + 1}"
        f"\tpas_neg_{offset}\t10\t-\n"
        for offset in (0, 40, 120, 300, 480)
    ]
    genes = _run_find_close(tmp_path, endbed, utr_lengths,
                            pos_rows=[], neg_rows=neg_rows)

    assigned = genes["gene_id"].tolist()
    assert len(assigned) == len(neg_rows)
    assert assigned.count(NEG_REAL_GENE) == len(neg_rows)
    assert NEG_SPAN_GENE not in assigned


def test_pas_only_the_spanning_model_reaches_still_gets_it(tmp_path,
                                                           overlapping_locus):
    """No tie, no change: a PAS at the spanning model's own 3' end keeps it."""
    endbed, utr_lengths = overlapping_locus
    genes = _run_find_close(
        tmp_path, endbed, utr_lengths,
        pos_rows=[],
        neg_rows=[f"{CHROM}\t{NEG_SPAN_START + 100}\t{NEG_SPAN_START + 101}"
                  f"\tpas_span\t10\t-\n"],
    )
    assert genes.loc["pas_span", "gene_id"] == NEG_SPAN_GENE


def test_one_row_per_pas(tmp_path, overlapping_locus):
    """closest -t all must never leak a duplicate PAS into the mapping frame."""
    endbed, utr_lengths = overlapping_locus
    genes = _run_find_close(
        tmp_path, endbed, utr_lengths,
        pos_rows=[f"{CHROM}\t{POS_PAS_POS}\t{POS_PAS_POS + 1}\tpas_pos\t10\t+\n"],
        neg_rows=[f"{CHROM}\t{NEG_PAS_POS}\t{NEG_PAS_POS + 1}\tpas_neg\t10\t-\n"],
    )
    assert genes.index.is_unique
