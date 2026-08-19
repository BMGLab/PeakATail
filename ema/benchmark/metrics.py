import numpy as np
from typing import Dict, List, Optional, Tuple
from pybedtools import BedTool


# Default cutoffs tuned for APA analysis: exact nucleotide through gene-region scale.
APA_DISTANCE_CUTOFFS: List[int] = [50, 100, 200, 500, 1000, 5000]


def compute_metrics(predicted_bed: str, reference_bed: str,
                    distance_cutoffs: Optional[List[int]] = None,
                    point_strand: bool = False,
                    restricted_reference_bed: Optional[str] = None,
                    null_shuffle: bool = False,
                    null_n_seeds: int = 3,
                    null_base_seed: int = 0,
                    null_gene_body_bed: Optional[str] = None) -> Dict:
    """Compute precision/recall/F1 at multiple distance cutoffs.

    The legacy interval-based, strand-agnostic metrics (top-level keys
    ``n_predicted``/``n_reference``/``cutoffs``) are always computed and stay
    byte-compatible. Optional point-mode strand-matched scoring and a
    shuffled-null baseline are added ALONGSIDE them under new blocks.

    Args:
        predicted_bed: Path to predicted PAS BED file.
        reference_bed: Path to reference PAS BED (PolyASite/PolyA_DB).
        distance_cutoffs: List of distance windows in bp.
            Defaults to APA_DISTANCE_CUTOFFS [50, 100, 200, 500, 1000, 5000].
        point_strand: If True, add a ``point_strand`` block: every interval
            reduced to its strand-aware 3'-most base (BED end exclusive:
            '+' -> end-1, '-' -> start), matched same-strand only. See
            :func:`compute_point_strand_metrics`.
        restricted_reference_bed: Optional cohort-restricted reference BED;
            adds ``recall_restricted`` inside the ``point_strand`` block
            (implies ``point_strand=True``).
        null_shuffle: If True, add a ``null`` block: width/chromosome/strand-
            preserving random placement of the predicted intervals, scored
            identically point-mode. See :func:`compute_null_baseline`.
        null_n_seeds: Number of shuffle seeds (default 3).
        null_base_seed: First shuffle seed (deterministic given this).
        null_gene_body_bed: Optional gene-body BED constraining where
            shuffled intervals may land.

    Returns:
        Dict with legacy metrics at each cutoff plus summary statistics, and
        — when requested — ``point_strand`` and ``null`` blocks.
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

    if point_strand or restricted_reference_bed is not None:
        results["point_strand"] = compute_point_strand_metrics(
            predicted_bed, reference_bed,
            distance_cutoffs=distance_cutoffs,
            restricted_reference_bed=restricted_reference_bed)

    if null_shuffle:
        results["null"] = compute_null_baseline(
            predicted_bed, reference_bed,
            distance_cutoffs=distance_cutoffs,
            n_seeds=null_n_seeds, base_seed=null_base_seed,
            gene_body_bed=null_gene_body_bed)

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


# ---------------------------------------------------------------------------
# Point-mode, strand-matched scoring (new)
#
# The legacy metrics above treat predictions and references as intervals and
# match them strand-agnostically via ``bedtools window``. The functions below
# instead reduce every entry to its strand-aware 3'-most base (the cleavage
# site under BED's exclusive-end convention) and only match entries on the
# same chromosome AND strand. They are pure Python/numpy — no bedtools
# dependency — and their output lives in new JSON blocks ("point_strand",
# "null") alongside the untouched legacy keys.
# ---------------------------------------------------------------------------

_VALID_STRANDS = ("+", "-")

#: How a stranded interval is reduced to a single base. BED end is EXCLUSIVE,
#: so the 3'-most covered base is ``end - 1`` on '+' and ``start`` on '-'.
POINT_DEFINITION = (
    "strand-aware 3'-most base; BED end exclusive: '+' -> end-1, '-' -> start"
)


def _load_stranded_intervals(bed_path: str) -> Tuple[List[Tuple[str, int, int, str]], int]:
    """Parse a BED file into stranded interval records.

    Entries whose strand field (column 6) is missing or not in ``{+, -}``
    are excluded and counted, so callers can report the exclusion honestly.

    Returns:
        (records, n_excluded) where records is a list of
        ``(chrom, start, end, strand)`` tuples with strand in ``{+, -}``.
    """
    records: List[Tuple[str, int, int, str]] = []
    n_excluded = 0
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split('\t')
            if len(fields) < 3:
                continue
            strand = fields[5] if len(fields) >= 6 else None
            if strand not in _VALID_STRANDS:
                n_excluded += 1
                continue
            records.append((fields[0], int(fields[1]), int(fields[2]), strand))
    return records, n_excluded


def _intervals_to_points(
    records: List[Tuple[str, int, int, str]]
) -> Dict[Tuple[str, str], np.ndarray]:
    """Reduce stranded intervals to their strand-aware 3'-most base.

    BED end is exclusive: '+' strand -> ``end - 1``; '-' strand -> ``start``.
    A 1-length interval ``(p, p+1)`` therefore yields point ``p`` on both
    strands.

    Returns:
        ``{(chrom, strand): sorted int64 array of point positions}``.
    """
    by_key: Dict[Tuple[str, str], List[int]] = {}
    for chrom, start, end, strand in records:
        point = end - 1 if strand == '+' else start
        by_key.setdefault((chrom, strand), []).append(point)
    return {k: np.sort(np.asarray(v, dtype=np.int64)) for k, v in by_key.items()}


def _nearest_point_distances(
    query: Dict[Tuple[str, str], np.ndarray],
    target: Dict[Tuple[str, str], np.ndarray],
) -> np.ndarray:
    """Distance from every query point to its nearest same-(chrom, strand) target.

    Query points with no target on the same chromosome AND strand get
    ``np.inf`` (they can never match at any cutoff — this is what enforces
    strand-matched matching).

    Returns:
        Float array with one entry per query point (order not meaningful).
    """
    out = []
    for key, qpos in query.items():
        tpos = target.get(key)
        if tpos is None or len(tpos) == 0:
            out.append(np.full(len(qpos), np.inf))
            continue
        idx = np.searchsorted(tpos, qpos)
        left = np.where(
            idx > 0,
            np.abs(qpos - tpos[np.maximum(idx - 1, 0)]).astype(float),
            np.inf,
        )
        right = np.where(
            idx < len(tpos),
            np.abs(tpos[np.minimum(idx, len(tpos) - 1)] - qpos).astype(float),
            np.inf,
        )
        out.append(np.minimum(left, right))
    if not out:
        return np.empty(0, dtype=float)
    return np.concatenate(out)


def _score_point_sets(
    pred_points: Dict[Tuple[str, str], np.ndarray],
    ref_points: Dict[Tuple[str, str], np.ndarray],
    distance_cutoffs: List[int],
) -> Tuple[int, int, Dict[int, Dict]]:
    """Strand-matched point scoring at each cutoff.

    A match is ``|predicted point - reference point| <= cutoff`` on the same
    chromosome and strand (boundary inclusive, mirroring ``bedtools window``).

    The recall-like quantity is deliberately named ``reference_coverage``:
    the fraction of reference points with >= 1 predicted point within the
    cutoff. It is NOT cohort recall — the reference atlas contains sites from
    tissues/conditions the assayed cohort never expresses.

    Returns:
        (n_predicted_points, n_reference_points, {cutoff: metrics dict}).
    """
    n_pred = int(sum(len(v) for v in pred_points.values()))
    n_ref = int(sum(len(v) for v in ref_points.values()))
    pred_d = _nearest_point_distances(pred_points, ref_points)
    ref_d = _nearest_point_distances(ref_points, pred_points)

    cutoffs: Dict[int, Dict] = {}
    for cutoff in distance_cutoffs:
        matched_predicted = int((pred_d <= cutoff).sum())
        matched_reference = int((ref_d <= cutoff).sum())
        precision = matched_predicted / n_pred if n_pred > 0 else 0.0
        coverage = matched_reference / n_ref if n_ref > 0 else 0.0
        cutoffs[cutoff] = {
            "precision": round(precision, 4),
            "reference_coverage": round(coverage, 4),
            "matched_predicted": matched_predicted,
            "matched_reference": matched_reference,
        }
    return n_pred, n_ref, cutoffs


def compute_point_strand_metrics(predicted_bed: str, reference_bed: str,
                                 distance_cutoffs: Optional[List[int]] = None,
                                 restricted_reference_bed: Optional[str] = None
                                 ) -> Dict:
    """Point-mode, strand-matched benchmark of predicted vs reference PAS.

    Every predicted and reference interval is reduced to its strand-aware
    3'-most base (see :data:`POINT_DEFINITION`); matching requires the same
    chromosome AND strand. Entries with strand not in ``{+, -}`` are excluded
    and counted in the output.

    Args:
        predicted_bed: Path to predicted PAS BED file (6+ columns).
        reference_bed: Path to reference PAS BED (PolyASite/PolyA_DB).
        distance_cutoffs: Distance windows in bp (boundary inclusive).
            Defaults to APA_DISTANCE_CUTOFFS.
        restricted_reference_bed: Optional cohort-restricted reference BED
            (e.g. atlas subset to genes expressed in the cohort). When given,
            each cutoff additionally reports ``recall_restricted`` — the
            fraction of restricted reference points matched — and the
            restriction is recorded under ``restricted_reference``.

    Returns:
        Dict for the ``point_strand`` JSON block. Per cutoff: ``precision``,
        ``reference_coverage`` (honest name for the recall-like quantity
        against the full reference), matched counts, and — when a restricted
        reference is supplied — ``recall_restricted``.
    """
    if distance_cutoffs is None:
        distance_cutoffs = APA_DISTANCE_CUTOFFS

    pred_records, pred_excluded = _load_stranded_intervals(predicted_bed)
    ref_records, ref_excluded = _load_stranded_intervals(reference_bed)
    pred_points = _intervals_to_points(pred_records)
    ref_points = _intervals_to_points(ref_records)

    n_pred, n_ref, cutoffs = _score_point_sets(pred_points, ref_points,
                                               distance_cutoffs)

    block: Dict = {
        "mode": "point_strand_matched",
        "point_definition": POINT_DEFINITION,
        "n_predicted_points": n_pred,
        "n_predicted_excluded_no_strand": pred_excluded,
        "n_reference_points": n_ref,
        "n_reference_excluded_no_strand": ref_excluded,
        "cutoffs": cutoffs,
    }

    if restricted_reference_bed is not None:
        r_records, r_excluded = _load_stranded_intervals(restricted_reference_bed)
        r_points = _intervals_to_points(r_records)
        n_restricted = int(sum(len(v) for v in r_points.values()))
        r_d = _nearest_point_distances(r_points, pred_points)
        block["restricted_reference"] = {
            "bed": restricted_reference_bed,
            "n_reference_points": n_restricted,
            "n_reference_excluded_no_strand": r_excluded,
        }
        for cutoff in distance_cutoffs:
            matched = int((r_d <= cutoff).sum())
            recall_restricted = matched / n_restricted if n_restricted > 0 else 0.0
            block["cutoffs"][cutoff]["recall_restricted"] = round(recall_restricted, 4)
            block["cutoffs"][cutoff]["matched_reference_restricted"] = matched

    return block


def _load_placement_intervals(bed_path: str) -> Dict[str, List[Tuple[int, int]]]:
    """Parse a gene-body BED into ``{chrom: [(start, end), ...]}``.

    Strand is ignored: gene bodies only constrain WHERE a shuffled interval
    may land, not its strand (which is preserved from the original entry).
    """
    by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split('\t')
            if len(fields) < 3:
                continue
            by_chrom.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    return by_chrom


def _shuffle_records(records: List[Tuple[str, int, int, str]],
                     rng: "np.random.Generator",
                     gene_body: Optional[Dict[str, List[Tuple[int, int]]]] = None,
                     chrom_extent: Optional[Dict[str, int]] = None
                     ) -> Tuple[List[Tuple[str, int, int, str]], int]:
    """Randomly re-place intervals, preserving width, chromosome, and strand.

    With ``gene_body`` given, each interval lands uniformly inside a gene-body
    interval on its own chromosome (candidates weighted by how many start
    positions they can host, so placement is uniform over all valid starts).
    Otherwise it lands uniformly in ``[0, chrom_extent - width]``.

    Intervals that cannot be placed (chromosome absent from the placement
    space, or wider than every candidate) are dropped and counted.
    """
    shuffled: List[Tuple[str, int, int, str]] = []
    n_unplaced = 0
    for chrom, start, end, strand in records:
        width = end - start
        if gene_body is not None:
            candidates = gene_body.get(chrom)
            if not candidates:
                n_unplaced += 1
                continue
            caps = np.array([(ge - gs) - width + 1 for gs, ge in candidates],
                            dtype=np.int64)
            caps = np.maximum(caps, 0)
            total = int(caps.sum())
            if total <= 0:
                n_unplaced += 1
                continue
            gi = int(rng.choice(len(candidates), p=caps / total))
            gs, ge = candidates[gi]
            new_start = int(gs + rng.integers(0, (ge - gs) - width + 1))
        else:
            extent = (chrom_extent or {}).get(chrom, 0)
            high = extent - width
            if high < 0:
                n_unplaced += 1
                continue
            new_start = int(rng.integers(0, high + 1))
        shuffled.append((chrom, new_start, new_start + width, strand))
    return shuffled, n_unplaced


def compute_null_baseline(predicted_bed: str, reference_bed: str,
                          distance_cutoffs: Optional[List[int]] = None,
                          n_seeds: int = 3, base_seed: int = 0,
                          gene_body_bed: Optional[str] = None) -> Dict:
    """Shuffled-null baseline, scored identically to point-mode.

    For each seed the predicted intervals are re-placed at random —
    width-, chromosome- and strand-preserving (``numpy.random.default_rng``
    with an explicit seed; pure Python, no bedtools) — optionally constrained
    to the intervals of ``gene_body_bed``. Without a gene-body BED, each
    chromosome's placement space is ``[0, max end seen in predicted+reference
    on that chromosome)``. The shuffled set is then reduced to strand-aware
    3'-most bases and scored strand-matched against the reference, exactly
    like :func:`compute_point_strand_metrics`.

    Args:
        predicted_bed: Path to predicted PAS BED file.
        reference_bed: Path to reference PAS BED.
        distance_cutoffs: Distance windows in bp. Defaults to
            APA_DISTANCE_CUTOFFS.
        n_seeds: Number of independent shuffles (seeds are
            ``base_seed .. base_seed + n_seeds - 1``).
        base_seed: First RNG seed; the run is fully deterministic given
            (base_seed, n_seeds).
        gene_body_bed: Optional BED restricting where shuffled intervals may
            land (recorded in the output).

    Returns:
        Dict for the ``null`` JSON block with ``per_seed`` results and their
        ``mean``.
    """
    if distance_cutoffs is None:
        distance_cutoffs = APA_DISTANCE_CUTOFFS

    pred_records, _ = _load_stranded_intervals(predicted_bed)
    ref_records, _ = _load_stranded_intervals(reference_bed)
    ref_points = _intervals_to_points(ref_records)

    gene_body = _load_placement_intervals(gene_body_bed) if gene_body_bed else None
    chrom_extent: Optional[Dict[str, int]] = None
    if gene_body is None:
        chrom_extent = {}
        for chrom, _start, end, _strand in pred_records + ref_records:
            chrom_extent[chrom] = max(chrom_extent.get(chrom, 0), end)

    seeds = [base_seed + i for i in range(n_seeds)]
    per_seed: List[Dict] = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        shuffled, n_unplaced = _shuffle_records(pred_records, rng,
                                                gene_body=gene_body,
                                                chrom_extent=chrom_extent)
        shuf_points = _intervals_to_points(shuffled)
        _n_pred, _n_ref, cutoffs = _score_point_sets(shuf_points, ref_points,
                                                     distance_cutoffs)
        per_seed.append({
            "seed": seed,
            "n_shuffled": len(shuffled),
            "n_unplaced": n_unplaced,
            "cutoffs": cutoffs,
        })

    mean_cutoffs: Dict[int, Dict] = {}
    for cutoff in distance_cutoffs:
        mean_cutoffs[cutoff] = {
            "precision": round(
                float(np.mean([s["cutoffs"][cutoff]["precision"] for s in per_seed])), 4),
            "reference_coverage": round(
                float(np.mean([s["cutoffs"][cutoff]["reference_coverage"]
                               for s in per_seed])), 4),
        }

    return {
        "mode": "shuffled_null_point_strand",
        "placement": ("gene_body_constrained" if gene_body_bed
                      else "chromosome_extent"),
        "gene_body_bed": gene_body_bed,
        "n_seeds": n_seeds,
        "seeds": seeds,
        "per_seed": per_seed,
        "mean": {"cutoffs": mean_cutoffs},
    }
