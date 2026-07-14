from pathlib import Path
import subprocess


def snap_beds_to_atlas(
    bed_paths: list[str | Path],
    dataset_ids: list[str],
    atlas_bed: str | Path,
    output_dir: str | Path,
    distance: int = 50,
) -> tuple[Path, Path]:
    """Snap called PAS coordinates to a reference atlas via bedtools closest.

    Args:
        bed_paths: List of BED files, one per dataset peak-calling run.
            Expected format per row: chrom, start, end, pasnumber, score, strand.
        dataset_ids: Parallel list — dataset_ids[i] owns bed_paths[i].
        atlas_bed: Path to reference atlas BED with columns:
            chrom, start, end, atlas_pas_id, score, strand.
        output_dir: Directory where output files are written.
        distance: Maximum distance in bp; called PAS farther than this are dropped.

    Returns:
        Tuple of (snapped_bed_path, mapping_path) where:
            snapped_bed_path: BED file with one row per unique atlas PAS that
                received at least one mapping, using atlas coordinates.
                Columns: chrom, start, end, atlas_pas_id, score, strand.
            mapping_path: TSV with header ``dataset_id\told_pasnumber\tnew_pas_id``
                and one row per kept input PAS, where new_pas_id is the
                atlas_pas_id string. Many called PAS may map to the same
                atlas PAS (many-to-one snap).

    Raises:
        FileNotFoundError: If any bed_path or atlas_bed does not exist.
        ValueError: If bed_paths and dataset_ids have different lengths.
        subprocess.CalledProcessError: If sort or bedtools closest fails.
    """
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
    with open(input_bed, "w") as out:
        for bed_path, dataset_id in zip(bed_paths, dataset_ids):
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
                    encoded_name = f"{dataset_id}::{pasnumber}"
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
    _IDX_ATLAS_CHROM = 6
    _IDX_ATLAS_START = 7
    _IDX_ATLAS_END = 8
    _IDX_ATLAS_PAS_ID = 9
    _IDX_ATLAS_SCORE = 10
    _IDX_ATLAS_STRAND = 11
    _MIN_EXPECTED_COLS = 13  # BED6 input + BED6 atlas + distance

    # mapping_rows: (dataset_id, old_pasnumber, atlas_pas_id, snap_distance_bp)
    mapping_rows: list[tuple[str, str, str, int]] = []
    # atlas_pas_id -> (chrom, start, end, pas_id, score, strand)
    atlas_hits: dict[str, tuple[str, str, str, str, str, str]] = {}

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

            atlas_chrom = cols[_IDX_ATLAS_CHROM]
            # bedtools reports "." for atlas chrom when no feature found
            if atlas_chrom == ".":
                continue

            try:
                dist = int(cols[_IDX_DISTANCE])
            except ValueError:
                continue

            # -1 means no feature found; filter by max distance
            if dist < 0 or dist > distance:
                continue

            input_name = cols[_IDX_INPUT_NAME]
            dataset_id, old_pasnumber = input_name.split("::", 1)
            atlas_pas_id = cols[_IDX_ATLAS_PAS_ID]

            mapping_rows.append((dataset_id, old_pasnumber, atlas_pas_id, dist))

            if atlas_pas_id not in atlas_hits:
                atlas_hits[atlas_pas_id] = (
                    atlas_chrom,
                    cols[_IDX_ATLAS_START],
                    cols[_IDX_ATLAS_END],
                    atlas_pas_id,
                    cols[_IDX_ATLAS_SCORE],
                    cols[_IDX_ATLAS_STRAND],
                )

    # Step 6: Re-key atlas PAS IDs as sequential integers so they match MTX rows
    # for annotate.py join. Keep the original atlas_pas_id in a sidecar lookup.
    sorted_atlas_hits = sorted(
        atlas_hits.values(),
        key=lambda row: (row[0], int(row[1])),
    )
    atlas_str_to_int: dict[str, str] = {}
    for i, hit in enumerate(sorted_atlas_hits, start=1):
        atlas_str_to_int[hit[3]] = str(i)

    # Step 6b: Write mapping TSV using integer new_pas_id.
    # 4th column 'snap_distance_bp' enables atlas_snap_diag histograms; existing
    # readers split on tabs and only require >=3 columns so this is additive.
    with open(mapping_path, "w") as f:
        f.write("dataset_id\told_pasnumber\tnew_pas_id\tsnap_distance_bp\n")
        for dataset_id, old_pasnumber, atlas_pas_id, snap_dist in mapping_rows:
            new_int_id = atlas_str_to_int.get(atlas_pas_id)
            if new_int_id is None:
                continue
            f.write(f"{dataset_id}\t{old_pasnumber}\t{new_int_id}\t{snap_dist}\n")

    # Step 6c: Sidecar lookup so atlas string IDs are preserved
    lookup_path = output_dir / "atlas_pas_id_lookup.tsv"
    with open(lookup_path, "w") as f:
        f.write("integer_id\tatlas_pas_id\n")
        for atlas_id, int_id in atlas_str_to_int.items():
            f.write(f"{int_id}\t{atlas_id}\n")

    # Step 7: Write snapped BED — col 4 is the integer ID (matches MTX rows)
    with open(snapped_bed_path, "w") as f:
        for chrom, start, end, pas_id, score, strand in sorted_atlas_hits:
            int_id = atlas_str_to_int[pas_id]
            f.write(f"{chrom}\t{start}\t{end}\t{int_id}\t{score}\t{strand}\n")

    # Step 8: Clean up temp files
    for tmp in (input_bed, sorted_input, sorted_atlas, closest_raw):
        tmp.unlink(missing_ok=True)

    # Step 9: Return paths
    return snapped_bed_path, mapping_path
