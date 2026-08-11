"""Pure annotation overlay of the unified PAS set against a reference atlas.

``ema/datasets/atlas_snap.py::snap_beds_to_atlas`` is a SNAPPING function --
even in its ``mode="annotate"`` form it decides PAS *identity* (matched PAS
are snapped onto the atlas coordinate, unmatched PAS get their own 1:1 id per
input row). Run on a PER-DATASET basis before any cross-dataset merge, that
means a genuinely novel PAS observed in two datasets at the same locus is
never merged -- it becomes two fragmented "unmatched" features instead of the
one unified feature ``ema.datasets.pas_merge.merge_pas_beds`` would have
produced. That fragmentation is exactly what a scientist hunting ALTERNATIVE
polyadenylation does NOT want: it splits real signal into noise.

The fix (D9 follow-up / "atlas as pure overlay"): the pipeline's default
``atlas_mode="annotate"`` no longer routes through ``snap_beds_to_atlas`` at
all. Instead it builds the unified PAS set EXACTLY like a no-atlas run --
``merge_pas_beds`` proximity-merges ALL called PAS across datasets into
unified ids, atlas plays NO role in merging or id assignment -- and this
module's :func:`annotate_pas_against_atlas` is run AFTER that merge, purely
to LABEL each already-final unified PAS with whether it sits near a known
atlas site. It never adds, drops, renumbers, or moves a PAS.

``mode="filter"`` (opt-in, unchanged) still goes through
``snap_beds_to_atlas``, where atlas snapping and dropping IS the point.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import subprocess


def annotate_pas_against_atlas(
    unified_pasbed_path: str | Path,
    atlas_bed_path: str | Path,
    distance: int,
    output_dir: str | Path,
    *,
    atlas_id_col: int = 3,
) -> tuple[Path, Path]:
    """Overlay atlas match/distance onto an ALREADY-FINAL unified PAS bed.

    This is a read-only annotation pass: it never rewrites
    ``unified_pasbed_path``, never changes a PAS id, and never drops a row.
    Every ``unified_pas_id`` present in the input bed gets exactly one row in
    the output ``atlas_status.tsv``, whether or not it matched.

    Args:
        unified_pasbed_path: BED6 (chrom, start, end, unified_pas_id, score,
            strand) -- the POST-MERGE unified PAS set, e.g. the output of
            :func:`ema.datasets.pas_merge.merge_pas_beds`. Column 4 is
            expected to be the id used everywhere downstream (MTX rows,
            ``clusters.h5ad`` var_names, ...).
        atlas_bed_path: Reference atlas BED. Must have >=6 columns
            (chrom, start, end, atlas_id, score, strand); extra trailing
            columns (e.g. PolyASite v3.0's ``stringency``/``annotation``)
            are tolerated and ignored -- see ``atlas_id_col`` below for the
            optional hook to carry one through.
        distance: Maximum summit-to-atlas distance (bp) to count as a match.
        output_dir: Directory to write the two sidecars into (typically the
            run's ``unified/`` directory, same convention as
            ``ema.datasets.atlas_snap.snap_beds_to_atlas``).
        atlas_id_col: 0-based column index of the atlas's own id/name field
            (default 3, i.e. BED col 4). Not currently carried into the
            output -- the must-have signal is atlas_match/atlas_distance_bp,
            not the specific atlas feature id -- but is exposed as a
            parameter (rather than hardcoded) so a caller wanting to carry
            atlas-specific columns (v3.0 stringency, annotation, ...) through
            has a stable index to hook additional column extraction onto
            without having to rediscover the closest-output column layout.

    Returns:
        ``(status_path, stats_path)``:
            ``status_path`` -- ``atlas_status.tsv`` with header
            ``unified_pas_id\tatlas_match\tatlas_distance_bp``, one row per
            PAS id in the input bed. ``atlas_match`` is the Python-bool
            string ``"True"``/``"False"``. ``atlas_distance_bp`` is an
            integer string (the real distance, even when it exceeds
            ``distance``) or ``""`` when bedtools found no atlas feature on
            that contig/strand at all.
            ``stats_path`` -- ``atlas_stats.json`` with
            ``{"mode": "annotate", "n_atlas_matched", "n_atlas_unmatched",
            "atlas_match_rate"}``.

    Raises:
        FileNotFoundError: If ``unified_pasbed_path`` or ``atlas_bed_path``
            does not exist.
        subprocess.CalledProcessError: If ``sort``/``bedtools closest`` fail.
    """
    unified_pasbed_path = Path(unified_pasbed_path)
    atlas_bed_path = Path(atlas_bed_path)
    if not unified_pasbed_path.exists():
        raise FileNotFoundError(f"unified_pasbed_path not found: {unified_pasbed_path}")
    if not atlas_bed_path.exists():
        raise FileNotFoundError(f"atlas_bed_path not found: {atlas_bed_path}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Temp files ---
    summit_bed = output_dir / "_atlas_overlay_summit.bed"
    sorted_summit = output_dir / "_atlas_overlay_summit_sorted.bed"
    sorted_atlas = output_dir / "_atlas_overlay_atlas_sorted.bed"
    closest_raw = output_dir / "_atlas_overlay_closest_raw.txt"

    # --- Output files ---
    status_path = output_dir / "atlas_status.tsv"
    stats_path = output_dir / "atlas_stats.json"

    # Step 1: collapse each unified PAS to its 1bp 3' summit for distance
    # measurement (mirrors the D2 fix in atlas_snap.py -- bedtools closest on
    # the full peak interval scores peak WIDTH, not PAS accuracy). This is a
    # SEPARATE temp file; unified_pasbed_path itself is never touched, and
    # the unified_pas_id (col 4) is carried through unchanged so the overlay
    # can be re-keyed onto it. Vectorized (pandas/numpy) -- no per-line loop.
    _ub = pd.read_csv(
        unified_pasbed_path, sep="\t", header=None, dtype=str,
        na_filter=False, engine="c",
    )
    if _ub.shape[0] == 0 or _ub.shape[1] < 6:
        summit_bed.write_text("")
    else:
        _ub = _ub.iloc[:, :6]
        _start = pd.to_numeric(_ub.iloc[:, 1], errors="coerce")
        _end = pd.to_numeric(_ub.iloc[:, 2], errors="coerce")
        _valid = _start.notna() & _end.notna()  # drop unparseable-coordinate rows
        _ub, _start, _end = _ub[_valid], _start[_valid].astype(np.int64), _end[_valid].astype(np.int64)
        _strand = _ub.iloc[:, 5].to_numpy()
        _summit = np.where(_strand == "+", _end.to_numpy() - 1, _start.to_numpy())
        pd.DataFrame({
            0: _ub.iloc[:, 0].to_numpy(),   # chrom
            1: _summit,                     # summit start
            2: _summit + 1,                 # summit end
            3: _ub.iloc[:, 3].to_numpy(),   # unified_pas_id
            4: _ub.iloc[:, 4].to_numpy(),   # score
            5: _strand,                     # strand
        }).to_csv(summit_bed, sep="\t", header=False, index=False)

    # Step 2: sort summit bed + atlas bed.
    with open(sorted_summit, "w") as out:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", str(summit_bed)],
            stdout=out, check=True,
        )
    with open(sorted_atlas, "w") as out:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", str(atlas_bed_path)],
            stdout=out, check=True,
        )

    # Step 3: bedtools closest -s (strand-aware) -d (report distance).
    with open(closest_raw, "w") as out:
        subprocess.run(
            [
                "bedtools", "closest",
                "-s",
                "-d",
                "-a", str(sorted_summit),
                "-b", str(sorted_atlas),
            ],
            stdout=out, check=True,
        )

    # Step 4: vectorized parse -- the atlas can be up to ~18M rows, so this
    # is pandas/numpy end to end, no per-line Python loop (mirrors
    # ema/datasets/pas_merge.py::_read_mtx_triples's approach).
    #
    # Column layout of bedtools closest -a BED6 -b BED(6+) output:
    #   0..5   = summit-bed row (chrom, start, end, unified_pas_id, score, strand)
    #   6..N-2 = atlas row (chrom, start, end, atlas_id, score, strand, [extra...])
    #   N-1    = distance (always last; "-1" sentinel or "." atlas chrom = no hit)
    try:
        df = pd.read_csv(
            closest_raw, sep="\t", header=None, dtype=str, na_filter=False,
            engine="c",
        )
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()

    if df.shape[0] == 0 or df.shape[1] < 7:
        pas_ids: np.ndarray = np.empty(0, dtype=object)
        atlas_match = np.empty(0, dtype=bool)
        atlas_distance_bp = np.empty(0, dtype=object)
    else:
        pas_ids = df.iloc[:, 3].to_numpy()
        atlas_chrom = df.iloc[:, 6].to_numpy()
        dist_num = pd.to_numeric(df.iloc[:, -1], errors="coerce").to_numpy()
        has_hit = (atlas_chrom != ".") & ~np.isnan(dist_num) & (dist_num >= 0)
        atlas_match = has_hit & (dist_num <= distance)
        # Guard the int64 cast: only cast where has_hit is True (dist_num is
        # never NaN there); everywhere else the string is "".
        dist_safe = np.where(has_hit, dist_num, 0.0).astype(np.int64).astype(str)
        atlas_distance_bp = np.where(has_hit, dist_safe, "")

    status = pd.DataFrame({
        "unified_pas_id": pas_ids,
        "atlas_match": atlas_match,
        "atlas_distance_bp": atlas_distance_bp,
    })
    # bedtools closest reports ALL ties by default -- a query row equidistant
    # from two atlas features appears twice with the same distance. Collapse
    # to one row per unified_pas_id (match/distance are identical across ties).
    status = status.drop_duplicates(subset="unified_pas_id", keep="first")
    # Deterministic ordering for readability / diffability.
    status = status.sort_values(
        by="unified_pas_id",
        key=lambda s: s.map(_id_sort_key),
    )
    status.to_csv(status_path, sep="\t", index=False)

    n_atlas_matched = int(status["atlas_match"].sum())
    n_atlas_unmatched = int(len(status) - n_atlas_matched)
    n_total = n_atlas_matched + n_atlas_unmatched
    atlas_stats = {
        "mode": "annotate",
        "n_atlas_matched": n_atlas_matched,
        "n_atlas_unmatched": n_atlas_unmatched,
        "atlas_match_rate": round(n_atlas_matched / n_total, 4) if n_total else 0.0,
    }
    stats_path.write_text(json.dumps(atlas_stats, indent=2))

    # Step 5: clean up temp files.
    for tmp in (summit_bed, sorted_summit, sorted_atlas, closest_raw):
        tmp.unlink(missing_ok=True)

    return status_path, stats_path


def _id_sort_key(pas_id: str):
    """Numeric-first sort key so integer unified_pas_ids sort naturally."""
    try:
        return (0, int(pas_id), "")
    except (TypeError, ValueError):
        return (1, 0, str(pas_id))
