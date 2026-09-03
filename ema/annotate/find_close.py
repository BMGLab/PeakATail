import logging

from ema.annotate.gtftobed import GENE_EXTENSION_BP
from ema.config import directory_config
import pybedtools
import pandas as pd

log = logging.getLogger(__name__)


# Confidence tier constants
TIER_1 = "TIER_1"  # within known UTR
TIER_2 = "TIER_2"  # within UTR x multiplier (possible UTR extension)
TIER_3 = "TIER_3"  # within max_distance (possible novel distant PAS)
INTERGENIC = "INTERGENIC"  # beyond max_distance

# bedtools closest output schema (PAS BED 6 cols + gene BED 6 cols + distance):
# 0: pas_chro, 1: pas_start, 2: pas_end, 3: pas_id, 4: pas_score, 5: pas_strand
# 6: gene_chro, 7: gene_start, 8: gene_end, 9: gene_id, 10: gene_symbol, 11: gene_strand
# 12: distance
PAS_CHRO_COL = 0
PAS_START_COL = 1
PAS_END_COL = 2
PAS_ID_COL = 3
PAS_STRAND_COL = 5
GENE_START_COL = 7
GENE_END_COL = 8
GENE_ID_COL = 9
GENE_SYMBOL_COL = 10
GENE_STRAND_COL = 11
DISTANCE_COL = 12  # last column from closest with D="b"

# The columns that identify one PAS.  A PAS may appear on several rows of
# ``closest -t all`` output (one per tied gene); it must appear on exactly one
# row of the result.
PAS_KEY_COLS = [PAS_CHRO_COL, PAS_START_COL, PAS_END_COL, PAS_ID_COL, PAS_STRAND_COL]


def assign_tier(distance, utr_length, utr_multiplier=2.0, max_distance=5000):
    """Assign a confidence tier based on distance and UTR length.

    Args:
        distance: Absolute distance from PAS to nearest gene (bp).
        utr_length: Known 3'UTR length for the gene (0 if unknown).
        utr_multiplier: Multiplier for UTR length to define TIER_2 boundary.
        max_distance: Maximum distance for any assignment (TIER_3 boundary).

    Returns:
        str: One of TIER_1, TIER_2, TIER_3, INTERGENIC.
    """
    abs_dist = abs(distance)

    if utr_length > 0:
        if abs_dist <= utr_length:
            return TIER_1
        if abs_dist <= utr_length * utr_multiplier:
            return TIER_2

    if abs_dist <= max_distance:
        return TIER_3

    return INTERGENIC


def _three_prime_offset(frame, gene_extension):
    """Distance from each PAS to its candidate gene's annotated 3' terminus.

    ``gtf_bed`` pushes every gene record's 3' end out by ``gene_extension`` bp,
    so the record end is NOT the annotated gene end.  Undo that shift (clamped
    inside the record, which is what the ``-`` strand near a contig start gets)
    and measure from there.

    Args:
        frame: bedtools ``closest`` output frame (positional schema above).
        gene_extension: bp ``gtf_bed`` added to each record's 3' end.

    Returns:
        pd.Series: absolute bp between the PAS midpoint and the gene 3' end.
    """
    gene_start = frame.iloc[:, GENE_START_COL].astype("int64")
    gene_end = frame.iloc[:, GENE_END_COL].astype("int64")
    strand = frame.iloc[:, GENE_STRAND_COL].astype(str)

    plus_terminus = (gene_end - gene_extension).clip(lower=gene_start)
    minus_terminus = (gene_start + gene_extension).clip(upper=gene_end)
    terminus = plus_terminus.where(strand != "-", minus_terminus)

    pas_mid = (frame.iloc[:, PAS_START_COL].astype("int64")
               + frame.iloc[:, PAS_END_COL].astype("int64")) // 2
    return (pas_mid - terminus).abs()


def _resolve_overlapping_genes(frame, gene_extension):
    """Pick one gene per PAS out of ``closest -t all`` ties (issue #99).

    In an overlapping locus a readthrough/spanning model covers its
    neighbour's 3'UTR, so a PAS there is distance 0 from BOTH models and
    bedtools' own ``-t first`` hands it to whichever sorts first -- the
    spanning one, because it starts further upstream.  That is how every PAS
    in CD68's 3'UTR ended up labelled SENP3-EIF4A1, leaving CD68 with zero
    assigned PAS cohort-wide.

    A PAS is a cleavage site, so among genes tied at the minimal distance the
    owner is the one whose annotated 3' end the PAS sits nearest; the span of
    the model and then the gene ID break any remaining tie deterministically.

    Args:
        frame: bedtools ``closest`` output frame, no-hit rows already dropped.
        gene_extension: bp ``gtf_bed`` added to each record's 3' end.

    Returns:
        tuple: (one row per PAS, number of PAS that had more than one
        candidate gene).
    """
    work = frame.copy()
    work["_abs_distance"] = work.iloc[:, DISTANCE_COL].astype("int64").abs()
    work["_three_prime_offset"] = _three_prime_offset(work, gene_extension)
    work["_gene_span"] = (work.iloc[:, GENE_END_COL].astype("int64")
                          - work.iloc[:, GENE_START_COL].astype("int64"))
    work["_gene_id"] = work.iloc[:, GENE_ID_COL].astype(str)

    rank_cols = ["_abs_distance", "_three_prime_offset", "_gene_span", "_gene_id"]
    work = work.sort_values(by=PAS_KEY_COLS + rank_cols, kind="mergesort")

    n_candidates = len(work)
    work = work.drop_duplicates(subset=PAS_KEY_COLS, keep="first")
    n_ambiguous = n_candidates - len(work)

    return work.drop(columns=rank_cols), n_ambiguous


def find_close(posbed_dir=None,
               negbed_dir=None,
               genomebed_dir=None,
               annotatedbed_dir=None,
               mergebed=None,
               utr_lengths=None,
               max_distance=5000,
               utr_multiplier=2.0,
               include_extended=False,
               gene_extension=None,
               ) -> pd.DataFrame:
    """Find the closest gene for each PAS and annotate with confidence tier.

    Uses bedtools closest (not intersect) to find the nearest gene for each
    PAS, then assigns confidence tiers based on distance and known 3'UTR
    length for adaptive thresholding.

    Args:
        posbed_dir: Path to positive-strand PAS BED file.
        negbed_dir: Path to negative-strand PAS BED file.
        genomebed_dir: Path to gene endpoint BED file.
        annotatedbed_dir: Path to write annotated PAS BED.
        mergebed: Path to write merged PAS BED.
        utr_lengths: Dict mapping gene_id -> max UTR length (bp). If None,
            all PAS use max_distance as the threshold.
        max_distance: Maximum distance (bp) for gene assignment.
        utr_multiplier: Multiplier for UTR length to define TIER_2.
        include_extended: If True, also keep TIER_3 PAS (default: TIER_1 + TIER_2 only).
        gene_extension: bp that ``gtf_bed`` appended to the 3' end of every
            record in ``genomebed_dir``; used to recover the annotated 3'
            terminus when resolving overlapping-loci ties. ``None`` (default)
            uses ``gtftobed.GENE_EXTENSION_BP``. Pass ``0`` for a gene BED that
            was not built by ``gtf_bed`` and so carries no extension.

    Returns:
        pd.DataFrame: DataFrame with PAS IDs as index and gene_id as values,
            plus a 'tier' column indicating confidence level.
    """
    if utr_lengths is None:
        utr_lengths = {}

    if gene_extension is None:
        gene_extension = GENE_EXTENSION_BP

    # Resolve config-dependent defaults at call time (set_directory_config
    # changes after import don't reach function-default values otherwise).
    if posbed_dir is None:
        posbed_dir = directory_config.posbed
    if negbed_dir is None:
        negbed_dir = directory_config.negbed
    if genomebed_dir is None:
        genomebed_dir = directory_config.endbed
    if annotatedbed_dir is None:
        annotatedbed_dir = directory_config.annotatedbed
    if mergebed is None:
        mergebed = directory_config.pasbed

    posbed = pybedtools.BedTool(posbed_dir)
    negbed = pybedtools.BedTool(negbed_dir)
    genomebed = pybedtools.BedTool(genomebed_dir).sort()

    pasbed = negbed.cat(posbed, postmerge=False).sort()
    pasbed.saveas(mergebed)

    # Use closest instead of intersect to find nearest gene for each PAS
    # D="b" = report distance relative to B (gene)
    # t="all" = report EVERY gene at the minimal distance; bedtools' own
    #   t="first" resolves overlapping loci by file order, which hands the PAS
    #   to the spanning/readthrough model (issue #99). The tie is resolved by
    #   _resolve_overlapping_genes below instead.
    # s=True = require same strand
    annotated = pasbed.closest(genomebed, s=True, D="b", t="all")
    annotated.saveas(annotatedbed_dir)

    annotated_frame = pd.read_csv(annotatedbed_dir, delimiter="\t", header=None)

    # Filter: remove rows where no gene was found (gene_id column == ".")
    annotated_frame = annotated_frame[annotated_frame.iloc[:, GENE_ID_COL] != "."]

    if annotated_frame.empty:
        return pd.DataFrame(columns=["gene_id", "tier"])

    # Collapse the closest -t all ties back to one gene per PAS (issue #99).
    annotated_frame, n_ambiguous = _resolve_overlapping_genes(
        annotated_frame, gene_extension)
    if n_ambiguous:
        log.info("overlapping loci: %d candidate gene(s) beyond the first "
                 "dropped across %d PAS; kept the gene whose annotated 3' end "
                 "the PAS lies nearest", n_ambiguous, len(annotated_frame))

    # Assign confidence tiers based on distance and UTR length
    def _get_tier(row):
        gene_id = row.iloc[GENE_ID_COL]
        distance = row.iloc[DISTANCE_COL]
        utr_len = utr_lengths.get(gene_id, 0)
        return assign_tier(distance, utr_len, utr_multiplier, max_distance)

    annotated_frame["tier"] = annotated_frame.apply(_get_tier, axis=1)

    # Filter by tier: keep TIER_1 + TIER_2 by default, optionally TIER_3
    keep_tiers = {TIER_1, TIER_2}
    if include_extended:
        keep_tiers.add(TIER_3)
    annotated_frame = annotated_frame[annotated_frame["tier"].isin(keep_tiers)]

    # Save the full annotated BED (dropping internal columns for cleanliness)
    # Keep: pas_chro, pas_start, pas_end, pas_id, gene_id, gene_symbol, strand, distance, tier
    output_cols = [0, 1, 2, PAS_ID_COL, GENE_ID_COL, GENE_SYMBOL_COL, 5, DISTANCE_COL, "tier"]
    annotated_out = annotated_frame[output_cols].copy()
    annotated_out.to_csv(annotatedbed_dir, sep="\t", header=False, index=False)

    # Build the PAS-to-gene mapping frame
    genes_frame = annotated_frame[[PAS_ID_COL, GENE_ID_COL, "tier"]].copy()
    genes_frame.columns = ["pas_id", "gene_id", "tier"]
    genes_frame = genes_frame.sort_values(by="pas_id", kind="quicksort")
    genes_frame = genes_frame.set_index("pas_id")

    return genes_frame


if __name__ == "__main__":
    find_close()
