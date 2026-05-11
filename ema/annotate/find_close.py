from ema.config import directory_config
import pybedtools
import pandas as pd


# Confidence tier constants
TIER_1 = "TIER_1"  # within known UTR
TIER_2 = "TIER_2"  # within UTR x multiplier (possible UTR extension)
TIER_3 = "TIER_3"  # within max_distance (possible novel distant PAS)
INTERGENIC = "INTERGENIC"  # beyond max_distance


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


def find_close(posbed_dir=None,
               negbed_dir=None,
               genomebed_dir=None,
               annotatedbed_dir=None,
               mergebed=None,
               utr_lengths=None,
               max_distance=5000,
               utr_multiplier=2.0,
               include_extended=False,
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

    Returns:
        pd.DataFrame: DataFrame with PAS IDs as index and gene_id as values,
            plus a 'tier' column indicating confidence level.
    """
    if utr_lengths is None:
        utr_lengths = {}

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
    # t="first" = report only the first closest match (breaks ties)
    # s=True = require same strand
    annotated = pasbed.closest(genomebed, s=True, D="b", t="first")
    annotated.saveas(annotatedbed_dir)

    annotated_frame = pd.read_csv(annotatedbed_dir, delimiter="\t", header=None)

    # bedtools closest output schema (PAS BED 6 cols + gene BED 6 cols + distance):
    # 0: pas_chro, 1: pas_start, 2: pas_end, 3: pas_id, 4: pas_score, 5: pas_strand
    # 6: gene_chro, 7: gene_start, 8: gene_end, 9: gene_id, 10: gene_symbol, 11: gene_strand
    # 12: distance
    PAS_ID_COL = 3
    GENE_ID_COL = 9
    GENE_SYMBOL_COL = 10
    DISTANCE_COL = 12  # last column from closest with D="b"

    # Filter: remove rows where no gene was found (gene_id column == ".")
    annotated_frame = annotated_frame[annotated_frame.iloc[:, GENE_ID_COL] != "."]

    if annotated_frame.empty:
        return pd.DataFrame(columns=["gene_id", "tier"])

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
