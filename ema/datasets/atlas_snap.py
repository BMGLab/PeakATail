from pathlib import Path
import json
import subprocess


def snap_beds_to_atlas(
    bed_paths: list[str | Path],
    dataset_ids: list[str],
    atlas_bed: str | Path,
    output_dir: str | Path,
    distance: int = 50,
    strands: list[str] | None = None,
    ledger=None,
    mode: str = "annotate",
) -> tuple[Path, Path]:
    """Snap called PAS coordinates to a reference atlas via bedtools closest.

    Args:
        bed_paths: List of BED files, one per dataset peak-calling run.
            Expected format per row: chrom, start, end, pasnumber, score, strand.
        dataset_ids: Parallel list — dataset_ids[i] owns bed_paths[i].
        atlas_bed: Path to reference atlas BED with columns:
            chrom, start, end, atlas_pas_id, score, strand.
        output_dir: Directory where output files are written.
        distance: Maximum distance in bp for a PAS to count as an atlas
            match (``atlas_match=True``). In ``"filter"`` mode this is also
            the hard drop cutoff (today's pre-existing behaviour).
        mode: ``"annotate"`` (default) or ``"filter"``.

            * ``"annotate"`` — the scientist is hunting ALTERNATIVE
              polyadenylation, so a PAS that doesn't match the atlas is
              signal, not noise: EVERY input PAS is kept. Matched PAS are
              snapped onto the atlas summit coordinate (today's behaviour);
              unmatched PAS keep their OWN (already summit-collapsed)
              coordinates and get their own unique PAS id (no cross-dataset
              merging of unmatched PAS is performed here — see module
              docstring / D6 report for the open question on whether that
              should happen).
            * ``"filter"`` — today's pre-D6+ behaviour: PAS with no atlas
              hit within ``distance`` are dropped (``dropped_at=
              "atlas_snap"`` in the ledger, if one is supplied).

            Regardless of mode, ``atlas_match``/``atlas_distance_bp`` are
            computed for every input PAS and written to two sidecars in
            ``output_dir`` (see Returns): ``atlas_status.tsv`` (per
            *unified* PAS id) and ``atlas_stats.json`` (run-level counts).

    Returns:
        Tuple of (snapped_bed_path, mapping_path) where:
            snapped_bed_path: BED file with one row per unique kept PAS
                (atlas coordinates for matched PAS, own coordinates for
                unmatched PAS in ``"annotate"`` mode).
                Columns: chrom, start, end, pas_id, score, strand.
            mapping_path: TSV with header ``dataset_id\told_pasnumber\tnew_pas_id``
                and one row per kept input PAS, where new_pas_id is the
                atlas_pas_id string (matched) or a per-PAS synthetic id
                (unmatched, annotate mode only). Many called PAS may map to
                the same MATCHED atlas PAS (many-to-one snap); unmatched PAS
                are always 1:1.

        Two additional sidecars are always written to ``output_dir``:
            ``atlas_status.tsv`` -- header ``new_pas_id\tatlas_match\t
                atlas_distance_bp``, one row per PAS id appearing in
                ``snapped_bed_path`` / the ``new_pas_id`` column of
                ``mapping_path``. ``atlas_distance_bp`` is ``""`` when no
                atlas feature was found at all (as opposed to found-but-
                too-far, which is a real integer). Callers (see
                ``ema/main.py``) read this to build the ``atlas_of`` map
                passed to ``ema.provenance.record_pas_drops``.
            ``atlas_stats.json`` -- ``{"n_atlas_matched", "n_atlas_unmatched",
                "atlas_match_rate", "mode"}``.

    Raises:
        FileNotFoundError: If any bed_path or atlas_bed does not exist.
        ValueError: If bed_paths and dataset_ids have different lengths, or
            mode is not one of ``{"annotate", "filter"}``.
        subprocess.CalledProcessError: If sort or bedtools closest fails.
    """
    if mode not in ("annotate", "filter"):
        raise ValueError(f"snap_beds_to_atlas: mode must be 'annotate' or 'filter', got {mode!r}")

    if len(bed_paths) != len(dataset_ids):
        raise ValueError(
            f"bed_paths length ({len(bed_paths)}) must match "
            f"dataset_ids length ({len(dataset_ids)})"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    atlas_bed = Path(atlas_bed)

    # --- Temp files ---
    input_bed = output_dir / "_atlas_input.bed"
    sorted_input = output_dir / "_atlas_input_sorted.bed"
    sorted_atlas = output_dir / "_atlas_sorted.bed"
    closest_raw = output_dir / "_closest_raw.txt"

    # --- Output files ---
    snapped_bed_path = output_dir / "atlas_snapped.bed"
    mapping_path = output_dir / "atlas_mapping.tsv"

    # Step 1: Build concatenated input BED with encoded names
    # col 4 = {dataset_id}::{pasnumber}
    #
    # D2 fix: `bedtools closest -s -d` below measures distance between the
    # FULL peak interval and the atlas. A wider peak is mechanically closer to
    # some atlas entry, so the reported distance scores peak WIDTH, not PAS
    # accuracy (see scripts/bench_summit_vs_atlas.py). The polyadenylation
    # site is the peak's 3' end, so we collapse each peak to a 1bp summit
    # interval (end-1 on "+", start on "-") here, before it ever reaches
    # bedtools. Everything downstream — sort, closest, distance filtering,
    # mapping, and the atlas coordinates written to the outputs — is
    # unchanged; only the distance measurement now reflects the summit.
    # B1: include the row strand in the encoded col-4 key so pos/neg PAS #N
    # within one dataset don't collide during count routing (see pas_merge).
    from ema.datasets.pas_merge import _encode_pas_key
    if strands is not None and len(strands) != len(bed_paths):
        raise ValueError(
            f"strands length ({len(strands)}) must match bed_paths ({len(bed_paths)})"
        )
    with open(input_bed, "w") as out:
        for _i, (bed_path, dataset_id) in enumerate(zip(bed_paths, dataset_ids)):
            bed_path = Path(bed_path)
            with open(bed_path) as f:
                for line in f:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) < 6:
                        continue
                    chrom, start, end, pasnumber, score, strand = (
                        parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
                    )
                    try:
                        start_i, end_i = int(start), int(end)
                    except ValueError:
                        continue  # unparseable coordinates; skip this row
                    # 3' summit: "+" -> last base (end-1); "-" -> first base (start)
                    summit = end_i - 1 if strand == "+" else start_i
                    key_strand = strand if strands is not None else None
                    encoded_name = _encode_pas_key(dataset_id, key_strand, pasnumber)
                    out.write(
                        f"{chrom}\t{summit}\t{summit + 1}\t{encoded_name}\t{score}\t{strand}\n"
                    )

    # Step 2: Sort input BED and atlas BED
    with open(sorted_input, "w") as out:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", str(input_bed)],
            stdout=out,
            check=True,
        )

    with open(sorted_atlas, "w") as out:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", str(atlas_bed)],
            stdout=out,
            check=True,
        )

    # Step 3: Run bedtools closest -s (strand-aware) -d (report distance)
    # Input is BED6 + BED6 → output has 12 cols + 1 distance = 13 cols total.
    # Col indices (0-based):
    #   0..5  = input row  (chrom, start, end, name=dataset::pasnumber, score, strand)
    #   6..11 = atlas row  (chrom, start, end, atlas_pas_id, score, strand)
    #   12    = distance   (integer; -1 means no closest found)
    with open(closest_raw, "w") as out:
        subprocess.run(
            [
                "bedtools", "closest",
                "-s",
                "-d",
                "-a", str(sorted_input),
                "-b", str(sorted_atlas),
            ],
            stdout=out,
            check=True,
        )

    # Step 4–6: Parse closest output, filter by distance, build mapping
    # Column layout (0-based, both inputs are BED6):
    #   0  = input chrom
    #   1  = input start
    #   2  = input end
    #   3  = input name  (dataset_id::pasnumber)
    #   4  = input score
    #   5  = input strand
    #   6  = atlas chrom
    #   7  = atlas start
    #   8  = atlas end
    #   9  = atlas_pas_id
    #   10 = atlas score
    #   11 = atlas strand
    #   12 = distance
    # Input is always BED6 (we wrote it). Atlas can be BED6 or wider (e.g. PolyASite
    # 2.0 has 12+ extra columns: scores, motifs, etc.). bedtools closest output:
    #   cols[0..5] = input  (6 cols)
    #   cols[6..6+N-1] = atlas (N cols, where N = atlas BED width >= 6)
    #   cols[-1]    = distance (always last)
    # So atlas-relative indices stay constant (chrom=0, start=1, end=2, name=3,
    # score=4, strand=5) BUT shifted by +6 from input.
    _IDX_INPUT_NAME = 3
    _IDX_INPUT_SCORE = 4
    _IDX_ATLAS_CHROM = 6
    _IDX_ATLAS_START = 7
    _IDX_ATLAS_END = 8
    _IDX_ATLAS_PAS_ID = 9
    _IDX_ATLAS_SCORE = 10
    _IDX_ATLAS_STRAND = 11
    _MIN_EXPECTED_COLS = 13  # BED6 input + BED6 atlas + distance

    from ema.datasets.pas_merge import _decode_pas_key
    # mapping_rows: (dataset_id, strand, old_pasnumber, key, snap_distance_bp)
    # `key` is the atlas_pas_id for matched PAS, or the (unique) input_name
    # for unmatched PAS kept in "annotate" mode.
    mapping_rows: list[tuple[str, str, str, str, int]] = []
    # atlas_pas_id -> (chrom, start, end, pas_id, score, strand) — MATCHED, atlas coords.
    atlas_hits: dict[str, tuple[str, str, str, str, str, str]] = {}
    # input_name -> (chrom, start, end, pas_id, score, strand) — UNMATCHED, own coords.
    unmatched_hits: dict[str, tuple[str, str, str, str, str, str]] = {}
    # per-PAS atlas status keyed by the SAME `key` used in mapping_rows, so it
    # can be re-keyed onto the final new_pas_id alongside atlas_hits/unmatched_hits.
    status_by_key: dict[str, tuple[bool, "int | str"]] = {}

    n_atlas_matched = 0
    n_atlas_unmatched = 0

    with open(closest_raw) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < _MIN_EXPECTED_COLS:
                continue
            # Distance is always the LAST column regardless of atlas width
            _IDX_DISTANCE = len(cols) - 1

            # E3: identity of THIS input PAS (available regardless of hit/miss).
            input_name = cols[_IDX_INPUT_NAME]
            dataset_id, in_strand, old_pasnumber = _decode_pas_key(input_name)
            in_chrom = cols[0]
            in_start, in_end = cols[1], cols[2]
            in_score = cols[_IDX_INPUT_SCORE]

            def _record(*, dropped_at: str, drop_reason: str, atlas_match, atlas_distance_bp) -> None:
                if ledger is None:
                    return
                ledger.record_pas(
                    orig_pas_key=input_name,
                    chrom=in_chrom, start=in_start, end=in_end, strand=in_strand,
                    snap_distance_bp=atlas_distance_bp,
                    last_stage="atlas_snap",
                    dropped_at=dropped_at,
                    drop_reason=drop_reason,
                    atlas_match=atlas_match,
                    atlas_distance_bp=atlas_distance_bp,
                )

            atlas_chrom = cols[_IDX_ATLAS_CHROM]
            atlas_pas_id = cols[_IDX_ATLAS_PAS_ID] if atlas_chrom != "." else None

            # bedtools reports "." for atlas chrom when no feature found at all.
            has_hit = atlas_chrom != "."
            dist: "int | None" = None
            unparseable = False
            if has_hit:
                try:
                    dist = int(cols[_IDX_DISTANCE])
                except ValueError:
                    unparseable = True
                if dist is not None and dist < 0:
                    # -1 sentinel: bedtools found no feature within range.
                    has_hit = False
                    dist = None

            if unparseable:
                # Malformed bedtools output — a data-integrity issue, not a
                # scientific filtering decision. Always dropped, both modes.
                _record(dropped_at="atlas_snap", drop_reason="unparseable closest distance",
                        atlas_match="", atlas_distance_bp="")
                continue

            atlas_match = bool(has_hit and dist is not None and dist <= distance)

            if atlas_match:
                n_atlas_matched += 1
            else:
                n_atlas_unmatched += 1

            if mode == "filter":
                if not atlas_match:
                    if not has_hit:
                        reason = "no atlas feature within range (-1)" if atlas_chrom != "." else "no atlas feature on contig/strand"
                    else:
                        reason = f"summit distance {dist}bp > atlas_distance {distance}bp"
                    _record(dropped_at="atlas_snap", drop_reason=reason,
                            atlas_match=False, atlas_distance_bp=dist if dist is not None else "")
                    continue
                # Matched -> keep (falls through to shared bookkeeping below).
                _record(dropped_at="", drop_reason="", atlas_match=True, atlas_distance_bp=dist)
                key = atlas_pas_id
                mapping_rows.append((dataset_id, in_strand, old_pasnumber, key, dist))
                status_by_key[key] = (True, dist)
                if key not in atlas_hits:
                    atlas_hits[key] = (
                        atlas_chrom, cols[_IDX_ATLAS_START], cols[_IDX_ATLAS_END],
                        key, cols[_IDX_ATLAS_SCORE], cols[_IDX_ATLAS_STRAND],
                    )
                continue

            # mode == "annotate": EVERY PAS is kept (never dropped for a
            # non-match — that's the signal the scientist is hunting for).
            _record(dropped_at="", drop_reason="", atlas_match=atlas_match,
                    atlas_distance_bp=dist if dist is not None else "")
            if atlas_match:
                key = atlas_pas_id
                mapping_rows.append((dataset_id, in_strand, old_pasnumber, key, dist))
                status_by_key[key] = (True, dist)
                if key not in atlas_hits:
                    atlas_hits[key] = (
                        atlas_chrom, cols[_IDX_ATLAS_START], cols[_IDX_ATLAS_END],
                        key, cols[_IDX_ATLAS_SCORE], cols[_IDX_ATLAS_STRAND],
                    )
            else:
                # Unmatched PAS keep their OWN (already summit-collapsed)
                # coordinates and are 1:1 (no cross-dataset merge here).
                key = input_name
                mapping_rows.append(
                    (dataset_id, in_strand, old_pasnumber, key, dist if dist is not None else -1)
                )
                status_by_key[key] = (False, dist if dist is not None else "")
                unmatched_hits[key] = (in_chrom, in_start, in_end, key, in_score, in_strand)

    # Step 6: Re-key PAS ids as sequential integers so they match MTX rows
    # for annotate.py join. Matched atlas hits AND (annotate mode only)
    # unmatched PAS share ONE sequential id space, sorted by coordinate so
    # ids stay geographically stable. Keep the original atlas_pas_id in a
    # sidecar lookup (unmatched PAS have no atlas string id to preserve).
    combined_hits = {**atlas_hits, **unmatched_hits}
    sorted_hits = sorted(
        combined_hits.values(),
        key=lambda row: (row[0], int(row[1])),
    )
    key_to_int: dict[str, str] = {}
    for i, hit in enumerate(sorted_hits, start=1):
        key_to_int[hit[3]] = str(i)

    # Step 6b: Write mapping TSV using integer new_pas_id.
    # Columns: dataset_id, strand, old_pasnumber, new_pas_id, snap_distance_bp.
    # B1: the strand column disambiguates pos/neg PAS #N. Readers are header-aware
    # (concat_matrices, parse_atlas_mapping) so this reordering is safe.
    with open(mapping_path, "w") as f:
        f.write("dataset_id\tstrand\told_pasnumber\tnew_pas_id\tsnap_distance_bp\n")
        for dataset_id, in_strand, old_pasnumber, key, snap_dist in mapping_rows:
            new_int_id = key_to_int.get(key)
            if new_int_id is None:
                continue
            f.write(
                f"{dataset_id}\t{in_strand}\t{old_pasnumber}\t{new_int_id}\t{snap_dist}\n"
            )

    # Step 6c: Sidecar lookup so atlas string IDs are preserved (matched only;
    # unmatched PAS have no atlas string id).
    lookup_path = output_dir / "atlas_pas_id_lookup.tsv"
    with open(lookup_path, "w") as f:
        f.write("integer_id\tatlas_pas_id\n")
        for atlas_id, int_id in key_to_int.items():
            if atlas_id in atlas_hits:
                f.write(f"{int_id}\t{atlas_id}\n")

    # Step 6d (D9): atlas_status.tsv — new_pas_id -> (atlas_match, atlas_distance_bp),
    # keyed on the SAME rekeyed id used everywhere else downstream (posbed,
    # count matrix, clusters.h5ad var_names). See docstring for consumers.
    status_path = output_dir / "atlas_status.tsv"
    with open(status_path, "w") as f:
        f.write("new_pas_id\tatlas_match\tatlas_distance_bp\n")
        for key, int_id in key_to_int.items():
            match, dist = status_by_key.get(key, ("", ""))
            f.write(f"{int_id}\t{match}\t{dist}\n")

    # Step 6e (D9): atlas_stats.json — run-level match/no-match counts.
    n_total = n_atlas_matched + n_atlas_unmatched
    atlas_stats = {
        "mode": mode,
        "n_atlas_matched": n_atlas_matched,
        "n_atlas_unmatched": n_atlas_unmatched,
        "atlas_match_rate": round(n_atlas_matched / n_total, 4) if n_total else 0.0,
    }
    with open(output_dir / "atlas_stats.json", "w") as f:
        json.dump(atlas_stats, f, indent=2)

    # Step 7: Write snapped BED — col 4 is the integer ID (matches MTX rows)
    from ema.datasets.pas_merge import pas_uid_of
    with open(snapped_bed_path, "w") as f:
        for chrom, start, end, pas_id, score, strand in sorted_hits:
            int_id = key_to_int[pas_id]
            f.write(f"{chrom}\t{start}\t{end}\t{int_id}\t{score}\t{strand}\n")

    # Step 7b (E1): content-addressed stable id sidecar (new_pas_id -> pas_uid).
    uid_path = output_dir / "pas_uid.tsv"
    with open(uid_path, "w") as f:
        f.write("new_pas_id\tpas_uid\n")
        for chrom, start, end, pas_id, score, strand in sorted_hits:
            int_id = key_to_int[pas_id]
            f.write(f"{int_id}\t{pas_uid_of(chrom, start, end, strand)}\n")

    # Step 8: Clean up temp files
    for tmp in (input_bed, sorted_input, sorted_atlas, closest_raw):
        tmp.unlink(missing_ok=True)

    # Step 9: Return paths
    return snapped_bed_path, mapping_path
