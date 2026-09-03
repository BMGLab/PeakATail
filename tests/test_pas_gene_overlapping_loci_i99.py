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
import os

import pandas as pd
import pytest

pybedtools = pytest.importorskip("pybedtools")

from ema.annotate import find_close as find_close_mod
from ema.annotate.find_close import TIER_1, _gene_geometry, find_close
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


# ---------------------------------------------------------------------------
# Real GRCh38.99 annotation (issue #99 blockers 1 and 2)
#
# ``tests/data/GRCh38.99_chr17_overlapping_loci.gtf`` is a verbatim slice of
# Homo_sapiens.GRCh38.99 -- the ``gene`` and ``three_prime_utr`` records of
# four real chr17 loci:
#
#   * CD68 (7,579,491-7,582,111) buried under the SENP3-EIF4A1 readthrough,
#     EIF4A1 and the AC016876.3 lncRNA -- the locus the issue is titled around;
#   * RNU6-862P inside NCOR1, MIR1288 inside PIGL, MIR33B/MIR6777 inside
#     SREBF1, MIR6778 inside SHMT1 -- nested small-RNA genes that carry no
#     ``three_prime_utr`` record at all and so cannot be graded above TIER_3.
#
# The PAS set is built the way the reviewer built it: every annotated 3'UTR
# terminus in the slice.  Each of those PAS has a correct answer by
# construction -- the gene whose 3'UTR record ends there.
# ---------------------------------------------------------------------------

REAL_GTF = os.path.join(os.path.dirname(__file__), "data",
                        "GRCh38.99_chr17_overlapping_loci.gtf")

CD68 = "ENSG00000129226"
SENP3_EIF4A1 = "ENSG00000277957"
AC016876_3 = "ENSG00000264772"          # lncRNA spanning CD68's 3'UTR
# The window the issue names: every PAS in it belongs to CD68.
CD68_WINDOW = (7_579_636, 7_582_386)

# Nested small-RNA genes with no three_prime_utr record, and their hosts.
NESTED_THIEVES = {
    "ENSG00000199674": "ENSG00000141027",   # RNU6-862P  <- NCOR1
    "ENSG00000221355": "ENSG00000108474",   # MIR1288    <- PIGL
    "ENSG00000274111": "ENSG00000072310",   # MIR6777    <- SREBF1
    "ENSG00000207839": "ENSG00000072310",   # MIR33B     <- SREBF1
    "ENSG00000283772": "ENSG00000176974",   # MIR6778    <- SHMT1
}


def _utr_terminus_pas(gtf_path):
    """Every annotated 3'UTR terminus in a GTF, as pos/neg PAS BED rows."""
    pos, neg, seen = [], [], set()
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "three_prime_utr":
                continue
            chrom, strand = fields[0], fields[6]
            terminus = int(fields[4]) if strand == "+" else int(fields[3])
            if (terminus, strand) in seen:
                continue
            seen.add((terminus, strand))
            row = (f"{chrom}\t{terminus - 1}\t{terminus}"
                   f"\tpas_{strand}_{terminus}\t10\t{strand}\n")
            (pos if strand == "+" else neg).append(row)
    return pos, neg


@pytest.fixture
def real_chr17_loci(tmp_path):
    """gene_end.bed + utr_lengths + PAS set from the real chr17 slice."""
    endbed = tmp_path / "gene_end.bed"
    utr_lengths = gtf_bed(
        endbeddir=str(endbed),
        gtfdir=REAL_GTF,
        featuresdir=str(tmp_path / "features.tsv"),
        utr_lengths_dir=str(tmp_path / "utr_lengths.tsv"),
    )
    pos_rows, neg_rows = _utr_terminus_pas(REAL_GTF)
    assert pos_rows and neg_rows, "fixture must exercise both strands"
    return endbed, utr_lengths, pos_rows, neg_rows


def _pas_position(pas_id):
    return int(pas_id.rsplit("_", 1)[1])


def test_cd68_gets_its_pas_on_the_real_grch38_annotation(tmp_path,
                                                         real_chr17_loci):
    """Blocker 2: the locus the issue is named after must actually be fixed.

    Both annotated 3'UTR termini of CD68 (7,581,622 and 7,582,111) fall in the
    window the issue names.  On develop both were labelled SENP3-EIF4A1; with
    nearest-terminus-wins the AC016876.3 lncRNA took 7,581,622 instead (and,
    having no 3'UTR record, got it dropped).  They are CD68's.
    """
    endbed, utr_lengths, pos_rows, neg_rows = real_chr17_loci
    genes = _run_find_close(tmp_path, endbed, utr_lengths, pos_rows, neg_rows)

    in_window = [p for p in genes.index
                 if CD68_WINDOW[0] <= _pas_position(p) <= CD68_WINDOW[1]]
    assert len(in_window) == 2, (
        f"expected both CD68 3'UTR termini in {CD68_WINDOW}, got {in_window}")
    for pas_id in in_window:
        assert genes.loc[pas_id, "gene_id"] == CD68, (
            f"{pas_id} went to {genes.loc[pas_id, 'gene_id']}, not CD68")
        assert genes.loc[pas_id, "tier"] == TIER_1

    assert SENP3_EIF4A1 not in set(genes.loc[in_window, "gene_id"])
    assert AC016876_3 not in set(genes.loc[in_window, "gene_id"])


def test_no_pas_is_lost_against_first_tie_wins(tmp_path, real_chr17_loci,
                                               monkeypatch):
    """Blocker 1: resolving a tie must never delete a PAS.

    The baseline is develop's behaviour: ``bedtools closest -t all`` emits ties
    in B-file order, so keeping the first row per PAS is exactly what
    ``-t first`` did.  Every PAS develop assigned must still be assigned, and
    the total must not fall -- a tie-break that hands a PAS to a gene with no
    ``three_prime_utr`` record grades it TIER_3 and the default tier filter
    then deletes it from annotatedpas.bed and from the count matrix.
    """
    endbed, utr_lengths, pos_rows, neg_rows = real_chr17_loci

    def _first_tie_wins(frame, *_args, **_kwargs):
        return frame.drop_duplicates(
            subset=find_close_mod.PAS_KEY_COLS, keep="first"), 0, 0

    baseline_dir = tmp_path / "baseline"
    fixed_dir = tmp_path / "fixed"
    baseline_dir.mkdir()
    fixed_dir.mkdir()

    monkeypatch.setattr(find_close_mod, "_resolve_overlapping_genes",
                        _first_tie_wins)
    baseline = _run_find_close(baseline_dir, endbed, utr_lengths,
                               pos_rows, neg_rows)
    monkeypatch.undo()
    fixed = _run_find_close(fixed_dir, endbed, utr_lengths,
                            pos_rows, neg_rows)

    lost = sorted(set(baseline.index) - set(fixed.index))
    assert not lost, f"{len(lost)} PAS silently dropped by the tie-break: {lost[:10]}"
    assert len(fixed) >= len(baseline)


def test_nested_small_rna_never_steals_its_hosts_pas(tmp_path,
                                                     real_chr17_loci):
    """A gene with no 3'UTR record must never win a tie against one that has.

    MIR33B/MIR6777 inside SREBF1, RNU6-862P inside NCOR1, MIR1288 inside PIGL
    and MIR6778 inside SHMT1 all rank nearest by 3'-terminus distance for PAS
    that belong to their protein-coding hosts.  None of them may take one.
    """
    endbed, utr_lengths, pos_rows, neg_rows = real_chr17_loci
    genes = _run_find_close(tmp_path, endbed, utr_lengths, pos_rows, neg_rows)

    for thief in NESTED_THIEVES:
        assert utr_lengths.get(thief, 0) == 0, "fixture assumption: no 3'UTR"
        assert thief not in set(genes["gene_id"])

    # and the hosts keep their own annotated termini
    for host in set(NESTED_THIEVES.values()):
        host_pas = genes.index[genes["gene_id"] == host]
        assert len(host_pas) > 0, f"host {host} was starved of every PAS"


def test_every_real_utr_terminus_keeps_a_gene(tmp_path, real_chr17_loci):
    """Each PAS here IS an annotated 3'UTR end, so none may go unassigned."""
    endbed, utr_lengths, pos_rows, neg_rows = real_chr17_loci
    genes = _run_find_close(tmp_path, endbed, utr_lengths, pos_rows, neg_rows)
    assert len(genes) == len(pos_rows) + len(neg_rows)


# ---------------------------------------------------------------------------
# Contig-start clamp on the minus strand
#
# gtf_bed extends every record's 3' end by 5 kb, EXCEPT a minus-strand gene
# that starts within 5 kb of the contig start: that one is clamped to 1 rather
# than shifted, so the annotated 3' terminus is not recoverable from the BED.
# Re-adding the extension unconditionally invents a terminus that can be
# thousands of bp from the real one (a gene at 300-4,000 is written as 1-4,000,
# and 1 + 5,000 clipped to the record end lands on 4,000 -- the gene's 5' end,
# 3.7 kb from its actual 3' terminus).  Real genes sit like this on the
# mitochondrion and on short scaffolds.
# ---------------------------------------------------------------------------

TRUE_START, TRUE_END = 300, 4000          # gtf_bed writes this as 1-4000
CLAMPED_ROW = [
    "clamp_contig", 349, 350, "pas_clamped", 10, "-",      # PAS columns
    "clamp_contig", 1, TRUE_END, "ENSG00000000010", "OWNER", "-",  # gene columns
    0,                                                     # bedtools distance
]
UNCLAMPED_ROW = [
    "clamp_contig", 20_349, 20_350, "pas_far", 10, "-",
    "clamp_contig", 15_000, 40_000, "ENSG00000000011", "RIVAL", "-",
    0,
]


def test_contig_start_clamp_yields_a_terminus_interval_not_a_wrong_point():
    """A clamped record gets a terminus *range* that brackets the real one.

    The old point estimate put this gene's 3' terminus at 4,000 -- its 5' end
    -- so a PAS in its own 3'UTR looked 3.7 kb away and any rival could take
    it.  The range says what is actually known: the terminus is somewhere in
    [record start, record start + extension], clipped to the record.
    """
    frame = pd.DataFrame([CLAMPED_ROW])
    lo, hi, body_start, body_end, minus = _gene_geometry(frame, 5000)

    assert bool(minus.iloc[0])
    assert lo.iloc[0] < hi.iloc[0], "clamped record must not get a point terminus"
    assert lo.iloc[0] <= TRUE_START <= hi.iloc[0], (
        "the real annotated terminus must lie inside the recovered range")
    # and the body must not exclude the gene's own 3'UTR
    assert body_start.iloc[0] <= TRUE_START <= body_end.iloc[0]
    assert body_start.iloc[0] <= 350 <= body_end.iloc[0]


def test_unclamped_minus_record_still_gets_an_exact_terminus():
    """Away from a contig start the extension is exactly recoverable."""
    frame = pd.DataFrame([UNCLAMPED_ROW])
    lo, hi, body_start, body_end, _ = _gene_geometry(frame, 5000)

    assert lo.iloc[0] == hi.iloc[0] == 20_000   # 15,000 + 5,000
    assert body_start.iloc[0] == 20_000
    assert body_end.iloc[0] == 40_000
