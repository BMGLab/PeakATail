import numpy as np
from typing import Dict, List, Optional, Tuple
from pybedtools import BedTool


# Default cutoffs tuned for APA analysis: exact nucleotide through gene-region scale.
APA_DISTANCE_CUTOFFS: List[int] = [50, 100, 200, 500, 1000, 5000]


def compute_metrics(predicted_bed: str, reference_bed: str,
                    distance_cutoffs: Optional[List[int]] = None) -> Dict:
    """Compute precision/recall/F1 at multiple distance cutoffs.

    Args:
        predicted_bed: Path to predicted PAS BED file.
        reference_bed: Path to reference PAS BED (PolyASite/PolyA_DB).
        distance_cutoffs: List of distance windows in bp.
            Defaults to APA_DISTANCE_CUTOFFS [50, 100, 200, 500, 1000, 5000].

    Returns:
        Dict with metrics at each cutoff plus summary statistics.
    """
    if distance_cutoffs is None:
        distance_cutoffs = APA_DISTANCE_CUTOFFS

    predicted = BedTool(predicted_bed)
    reference = BedTool(reference_bed)

    n_predicted = predicted.count()
    n_reference = reference.count()

    results = {
        "n_predicted": n_predicted,
        "n_reference": n_reference,
        "cutoffs": {}
    }

    for cutoff in distance_cutoffs:
        # Precision: fraction of predicted PAS within cutoff of a reference PAS
        hits = predicted.window(reference, w=cutoff)
        n_matched_predicted = len(set(str(f).split('\t')[3] for f in hits))
        precision = n_matched_predicted / n_predicted if n_predicted > 0 else 0

        # Recall: fraction of reference PAS matched by predicted
        hits_rev = reference.window(predicted, w=cutoff)
        n_matched_reference = len(set(str(f).split('\t')[3] for f in hits_rev))
        recall = n_matched_reference / n_reference if n_reference > 0 else 0

        # F1
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

        results["cutoffs"][cutoff] = {
            "precision": round(precision, 4),
            "recall": round(recall, 4),
            "f1": round(f1, 4),
            "matched_predicted": n_matched_predicted,
            "matched_reference": n_matched_reference
        }

    return results


def compute_precision_focused_metrics(predicted_bed: str, reference_bed: str,
                                      distance_cutoffs: Optional[List[int]] = None) -> Dict:
    """Compute precision-focused metrics for APA analysis.

    Primary question: "Are the PAS we found real?"

    Precision is the primary metric at each cutoff. In place of recall, this
    function reports *gene coverage*: for each expressed gene represented in the
    predicted BED (identified by chromosome region), the fraction that has at
    least one predicted PAS matched to a known reference entry within the given
    distance window.

    Args:
        predicted_bed: Path to predicted PAS BED file. Name field (column 4)
            should encode the gene or peak identifier; the chromosome (column 1)
            is used to group predictions by gene region.
        reference_bed: Path to reference PAS BED (PolyASite/PolyA_DB).
        distance_cutoffs: List of distance windows in bp.
            Defaults to APA_DISTANCE_CUTOFFS [50, 100, 200, 500, 1000, 5000].

    Returns:
        Dict structured as::

            {
                "n_predicted": int,
                "n_reference": int,
                "n_expressed_genes": int,   # unique chromosomes/regions in predicted BED
                "cutoffs": {
                    <cutoff_bp>: {
                        "precision": float,         # PRIMARY metric
                        "matched_predicted": int,
                        "novel_pas": int,           # predicted not near any reference
                        "novel_fraction": float,
                        "gene_coverage": float,     # fraction of expressed genes
                                                    # with >= 1 matched PAS
                        "genes_with_match": int,
                    },
                    ...
                }
            }
    """
    if distance_cutoffs is None:
        distance_cutoffs = APA_DISTANCE_CUTOFFS

    predicted = BedTool(predicted_bed)
    reference = BedTool(reference_bed)

    n_predicted = predicted.count()
    n_reference = reference.count()

    # Build a mapping of gene region -> set of predicted peak name IDs.
    # Gene region is represented by the chromosome field (column 0) as a
    # coarse proxy.  Callers with richer annotation can use a pre-grouped BED.
    gene_to_peaks: Dict[str, set] = {}
    for feature in predicted:
        fields = str(feature).strip().split('\t')
        chrom = fields[0]
        peak_id = fields[3] if len(fields) >= 4 else f"{chrom}:{fields[1]}-{fields[2]}"
        gene_to_peaks.setdefault(chrom, set()).add(peak_id)

    n_expressed_genes = len(gene_to_peaks)

    results: Dict = {
        "n_predicted": n_predicted,
        "n_reference": n_reference,
        "n_expressed_genes": n_expressed_genes,
        "cutoffs": {}
    }

    for cutoff in distance_cutoffs:
        # Precision: how many of our peaks are within cutoff of a known PAS?
        hits = predicted.window(reference, w=cutoff)
        matched_ids: set = set()
        matched_chroms: set = set()
        for feature in hits:
            fields = str(feature).strip().split('\t')
            peak_id = fields[3] if len(fields) >= 4 else fields[0]
            matched_ids.add(peak_id)
            matched_chroms.add(fields[0])

        n_matched = len(matched_ids)
        precision = n_matched / n_predicted if n_predicted > 0 else 0.0
        novel = n_predicted - n_matched
        novel_fraction = novel / n_predicted if n_predicted > 0 else 0.0

        # Gene coverage: fraction of expressed genes with >= 1 matched PAS.
        genes_with_match = len(matched_chroms & set(gene_to_peaks.keys()))
        gene_coverage = genes_with_match / n_expressed_genes if n_expressed_genes > 0 else 0.0

        results["cutoffs"][cutoff] = {
            "precision": round(precision, 4),
            "matched_predicted": n_matched,
            "novel_pas": novel,
            "novel_fraction": round(novel_fraction, 4),
            "gene_coverage": round(gene_coverage, 4),
            "genes_with_match": genes_with_match,
        }

    return results


def compute_distance_distribution(predicted_bed: str, reference_bed: str) -> List[int]:
    """Compute distance from each predicted PAS to nearest reference PAS.

    Args:
        predicted_bed: Path to predicted PAS BED file.
        reference_bed: Path to reference PAS BED.

    Returns:
        List of distances (in bp) suitable for histogram plotting.
    """
    predicted = BedTool(predicted_bed)
    reference = BedTool(reference_bed)

    closest = predicted.closest(reference, d=True)
    distances = []
    for feature in closest:
        fields = str(feature).strip().split('\t')
        try:
            dist = abs(int(fields[-1]))
            if dist >= 0:  # -1 means no match found
                distances.append(dist)
        except (ValueError, IndexError):
            pass

    return distances


def compute_overlap_sets(bed_files: Dict[str, str], window: int = 50) -> Dict[str, set]:
    """Compute PAS position sets for Venn/UpSet diagram.

    Args:
        bed_files: Mapping of {strategy_name: bed_file_path}.
        window: Merge window in bp for considering two PAS as the same site.

    Returns:
        Mapping of {strategy_name: set of (chrom, rounded_pos) tuples}.
    """
    sets: Dict[str, set] = {}
    for name, path in bed_files.items():
        positions: set = set()
        with open(path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    mid = (start + end) // 2
                    # Round to nearest window to allow fuzzy matching
                    rounded = (mid // window) * window
                    positions.add((chrom, rounded))
        sets[name] = positions
    return sets
