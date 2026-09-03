import logging

import numpy as np
import pandas as pd
import pybedtools

from ema.annotate.gtftobed import GENE_EXTENSION_BP
from ema.config import directory_config

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


def assign_tiers(distances, utr_lengths, utr_multiplier=2.0, max_distance=5000):
    """Vectorised :func:`assign_tier` over two aligned Series.

    Args:
        distances: Signed PAS->gene distances (bp).
        utr_lengths: Known 3'UTR length per row (0 where unknown).
        utr_multiplier: Multiplier for UTR length to define TIER_2 boundary.
        max_distance: Maximum distance for any assignment (TIER_3 boundary).

    Returns:
        pd.Series: tier label per row.
    """
    abs_dist = distances.abs()
    has_utr = utr_lengths > 0
    conditions = [
        has_utr & (abs_dist <= utr_lengths),
        has_utr & (abs_dist <= utr_lengths * utr_multiplier),
        abs_dist <= max_distance,
    ]
    return pd.Series(
        np.select(conditions, [TIER_1, TIER_2, TIER_3], default=INTERGENIC),
        index=distances.index,
    )


def _gene_geometry(frame, gene_extension):
    """Recover each candidate gene's annotated geometry from the extended BED.

    ``gtf_bed`` pushes every gene record's 3' end out by ``gene_extension`` bp,
    so a record's coordinates are NOT the annotated gene's.  Undo that shift.

    The one lossy case is a ``-`` strand gene whose start is within
    ``gene_extension`` bp of the contig start: ``gtf_bed`` clamps it to 1
    instead of subtracting, so the annotated 3' terminus is unrecoverable --
    all that is known is that it lies somewhere in
    ``[gene_start, gene_start + gene_extension]``.  Those rows therefore get a
    terminus *interval* rather than a point, and the widest possible gene body,
    so the clamp can never make a candidate look falsely bad.

    Args:
        frame: bedtools ``closest`` output frame (positional schema above).
        gene_extension: bp ``gtf_bed`` added to each record's 3' end.

    Returns:
        tuple: ``(terminus_lo, terminus_hi, body_start, body_end, minus)`` --
        the first four int64 Series, the last a boolean strand mask.
        ``terminus_lo == terminus_hi`` except on contig-start-clamped rows.
    """
    gene_start = frame.iloc[:, GENE_START_COL].astype("int64")
    gene_end = frame.iloc[:, GENE_END_COL].astype("int64")
    minus = frame.iloc[:, GENE_STRAND_COL].astype(str).eq("-")

    clamped = minus & gene_start.le(1)

    plus_terminus = (gene_end - gene_extension).clip(lower=gene_start)
    minus_terminus_hi = (gene_start + gene_extension).clip(upper=gene_end)
    minus_terminus_lo = minus_terminus_hi.where(~clamped, gene_start)

    terminus_lo = plus_terminus.where(~minus, minus_terminus_lo)
    terminus_hi = plus_terminus.where(~minus, minus_terminus_hi)
    body_start = gene_start.where(~minus, terminus_lo)
    body_end = plus_terminus.where(~minus, gene_end)

    return terminus_lo, terminus_hi, body_start, body_end, minus


def _resolve_overlapping_genes(frame, gene_extension, utr_lengths, keep_tiers):
    """Pick one gene per PAS out of ``closest -t all`` ties (issue #99).

    In an overlapping locus a readthrough/spanning model covers its
    neighbour's 3'UTR, so a PAS there is distance 0 from BOTH models and
    bedtools' own ``-t first`` hands it to whichever sorts first -- the
    spanning one, because it starts further upstream.  That is how every PAS
    in CD68's 3'UTR ended up labelled SENP3-EIF4A1, leaving CD68 with zero
    assigned PAS cohort-wide.

    Nearest-annotated-3'-terminus is NOT a safe rule on its own: the nearest
    terminus in such a locus is very often a nested miRNA/snRNA/pseudogene or
    a lncRNA that carries no ``three_prime_utr`` record at all (MIR33B inside
    SREBF1, RNU6-862P inside NCOR1, AC016876.3 over CD68).  Handing the PAS to
    one of those genes grades it TIER_3 -- because ``utr_lengths`` has no entry
    -- and the default tier filter then deletes the PAS from
    ``annotatedpas.bed`` and from the count matrix.  Silently losing a PAS is
    worse than mislabelling it, so the candidates are ranked by what the
    annotation actually says, in this order (all ascending, best first):

    0. the candidate does not turn a kept PAS into a filtered-out one;
    1. the candidate has a ``three_prime_utr`` record (hard rule: a gene
       without one never wins a tie against a gene with one);
    2. the PAS lies inside the candidate's annotated 3'UTR footprint;
    3. the PAS lies inside the candidate's UNEXTENDED gene body;
    4. distance to the candidate's annotated 3' terminus;
    5. the candidate's annotated span, then its gene ID, for determinism.

    The 3'UTR footprint is the ``utr_lengths[gene]`` bp immediately 5' of the
    annotated 3' terminus -- the longest annotated 3'UTR of the gene, anchored
    at the gene end.  ``utr_lengths`` is the only UTR geometry ``find_close``
    is given, and it is the same number ``assign_tier`` already trusts.

    Args:
        frame: bedtools ``closest`` output frame, no-hit rows dropped, with a
            per-candidate ``"tier"`` column already assigned.
        gene_extension: bp ``gtf_bed`` added to each record's 3' end.
        utr_lengths: Dict mapping gene_id -> max annotated 3'UTR length (bp).
        keep_tiers: The tiers ``find_close`` will keep; a candidate whose tier
            is outside it would delete the PAS.

    Returns:
        tuple: ``(one row per PAS, n_contested, n_dropped_candidates)`` where
        ``n_contested`` is the number of PAS that had more than one candidate
        gene at the minimal distance.
    """
    work = frame.copy()

    gene_id = work.iloc[:, GENE_ID_COL].astype(str)
    utr_len = gene_id.map(utr_lengths).fillna(0).astype("int64")
    terminus_lo, terminus_hi, body_start, body_end, minus = _gene_geometry(
        work, gene_extension)

    pas_mid = (work.iloc[:, PAS_START_COL].astype("int64")
               + work.iloc[:, PAS_END_COL].astype("int64")) // 2

    # Point-to-interval distance; the interval is a point unless the record was
    # clamped at a contig start, where the smallest consistent offset is used.
    terminus_offset = ((terminus_lo - pas_mid).clip(lower=0)
                       + (pas_mid - terminus_hi).clip(lower=0))

    utr_lo = (terminus_lo - utr_len).where(~minus, terminus_lo)
    utr_hi = terminus_hi.where(~minus, terminus_hi + utr_len)
    in_utr = utr_len.gt(0) & pas_mid.ge(utr_lo) & pas_mid.le(utr_hi)
    in_body = pas_mid.ge(body_start) & pas_mid.le(body_end)

    work["_abs_distance"] = work.iloc[:, DISTANCE_COL].astype("int64").abs()
    work["_would_drop_pas"] = (~work["tier"].isin(keep_tiers)).astype("int8")
    work["_no_utr_record"] = (~utr_len.gt(0)).astype("int8")
    work["_not_in_utr"] = (~in_utr).astype("int8")
    work["_not_in_body"] = (~in_body).astype("int8")
    work["_terminus_offset"] = terminus_offset
    work["_gene_span"] = body_end - body_start
    work["_gene_id"] = gene_id

    rank_cols = ["_abs_distance", "_would_drop_pas", "_no_utr_record",
                 "_not_in_utr", "_not_in_body", "_terminus_offset",
                 "_gene_span", "_gene_id"]
    work = work.sort_values(by=PAS_KEY_COLS + rank_cols, kind="mergesort")

    contested = work.duplicated(subset=PAS_KEY_COLS, keep=False)
    n_candidates = len(work)
    work = work.drop_duplicates(subset=PAS_KEY_COLS, keep="first")
    n_contested = int(contested.loc[work.index].sum())

    return work.drop(columns=rank_cols), n_contested, n_candidates - len(work)


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
            record in ``genomebed_dir``; used to recover the annotated gene
            body and 3' terminus when resolving overlapping-loci ties.
            ``None`` (default) uses ``gtftobed.GENE_EXTENSION_BP``. Pass ``0``
            for a gene BED that was not built by ``gtf_bed`` and so carries no
            extension.

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

    # Tier every CANDIDATE, not just the winner: the tie-break in
    # _resolve_overlapping_genes needs to know which candidates would survive
    # the tier filter, so that resolving a tie can never delete a PAS.
    keep_tiers = {TIER_1, TIER_2}
    if include_extended:
        keep_tiers.add(TIER_3)

    candidate_utr_len = (annotated_frame.iloc[:, GENE_ID_COL].astype(str)
                         .map(utr_lengths).fillna(0).astype("int64"))
    annotated_frame["tier"] = assign_tiers(
        annotated_frame.iloc[:, DISTANCE_COL].astype("int64"),
        candidate_utr_len, utr_multiplier, max_distance)

    # Collapse the closest -t all ties back to one gene per PAS (issue #99).
    annotated_frame, n_contested, n_dropped = _resolve_overlapping_genes(
        annotated_frame, gene_extension, utr_lengths, keep_tiers)
    if n_contested:
        log.info("overlapping loci: %d of %d PAS had more than one candidate "
                 "gene at the minimal distance (%d rival candidate rows "
                 "discarded); kept the candidate the annotation supports",
                 n_contested, len(annotated_frame), n_dropped)

    # Filter by tier: keep TIER_1 + TIER_2 by default, optionally TIER_3
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
