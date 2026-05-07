from pathlib import Path
import subprocess
import scipy.sparse as sp
import numpy as np


def merge_pas_beds(
    bed_paths: list[str | Path],
    dataset_ids: list[str],
    output_dir: str | Path,
    gap: int = 100,
) -> tuple[Path, Path]:
    """Merge per-dataset PAS BED files into a unified strand-aware coordinate set.

    Inputs:
        bed_paths: list of BED file paths (one per dataset peak-calling run)
        dataset_ids: parallel list — dataset_ids[i] owns bed_paths[i]
        output_dir: where to write outputs
        gap: bedtools merge -d distance (default 100bp)

    Outputs:
        merged_bed: chrom\\tstart\\tend\\tnew_pas_id\\tscore\\tstrand
        mapping_path: TSV with header "dataset_id\\told_pasnumber\\tnew_pas_id"
            One row per ORIGINAL pasnumber (so multiple original pasnumbers
            can map to the same new_pas_id when bedtools merge groups them).

    Returns: (merged_bed_path, mapping_path)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cat_path = output_dir / "multi_sample_cat.bed"
    sorted_path = output_dir / "multi_sample_cat_sorted.bed"
    bedtools_out = output_dir / "multi_sample_bedtools_merged.bed"
    merged_path = output_dir / "multi_sample_merged.bed"
    mapping_path = output_dir / "multi_sample_pas_mapping.tsv"

    # Step 1: Concatenate all BED files, encoding dataset_id::pasnumber in col 4.
    with open(cat_path, "w") as out:
        for bed_path, dataset_id in zip(bed_paths, dataset_ids):
            with open(bed_path) as f:
                for line in f:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) < 4:
                        continue
                    # parts[3] is the original pasnumber (integer string)
                    parts[3] = f"{dataset_id}::{parts[3]}"
                    out.write("\t".join(parts) + "\n")

    # Step 2: Sort by chrom then start position.
    with open(sorted_path, "w") as sorted_out:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", str(cat_path)],
            stdout=sorted_out,
            check=True,
        )

    # Step 3: Strand-aware bedtools merge.
    # -s: only merge features on the same strand (critical correctness requirement).
    # -c 4,5,6 -o collect,first,first: gather all encoded IDs, keep first score/strand.
    with open(bedtools_out, "w") as bt_out:
        subprocess.run(
            [
                "bedtools", "merge",
                "-s",
                "-d", str(gap),
                "-i", str(sorted_path),
                "-c", "4,5,6",
                "-o", "collapse,first,first",
            ],
            stdout=bt_out,
            check=True,
        )

    # Step 4: Assign new_pas_ids and build mapping.
    # bedtools collect produces a comma-separated list in col 4, e.g.:
    #   "ds1::123,ds1::456,ds2::789"
    merged_lines: list[str] = []
    mapping_rows: list[tuple[str, str, str]] = []  # (dataset_id, old_pasnumber, new_pas_id)
    new_id = 1

    with open(bedtools_out) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            # bedtools merge -s produces: chrom start end collected_col4 first_score first_strand
            if len(parts) < 6:
                continue
            chrom, start, end = parts[0], parts[1], parts[2]
            collected = parts[3]   # "ds1::123,ds1::456,..."
            score = parts[4]
            strand = parts[5]

            new_pas_id = f"PAS_{new_id}"
            new_id += 1

            # Write one mapping row per original entry.
            for entry in collected.split(","):
                entry = entry.strip()
                if "::" not in entry:
                    continue
                dataset_id, old_pasnumber = entry.split("::", 1)
                mapping_rows.append((dataset_id, old_pasnumber, new_pas_id))

            merged_lines.append(
                f"{chrom}\t{start}\t{end}\t{new_pas_id}\t{score}\t{strand}"
            )

    # Step 5: Write final merged BED.
    with open(merged_path, "w") as f:
        f.write("\n".join(merged_lines))
        if merged_lines:
            f.write("\n")

    # Step 6: Write mapping TSV.
    with open(mapping_path, "w") as f:
        f.write("dataset_id\told_pasnumber\tnew_pas_id\n")
        for dataset_id, old_pasnumber, new_pas_id in mapping_rows:
            f.write(f"{dataset_id}\t{old_pasnumber}\t{new_pas_id}\n")

    # Clean up temp files.
    cat_path.unlink(missing_ok=True)
    sorted_path.unlink(missing_ok=True)
    bedtools_out.unlink(missing_ok=True)

    return merged_path, mapping_path


def concat_matrices(
    mtx_paths: list[str | Path],
    dataset_ids: list[str],
    cb_paths: list[str | Path],
    mapping_path: str | Path,
    output_mtx: str | Path,
    output_cb: str | Path,
) -> tuple[Path, Path]:
    """Concatenate per-dataset MTX files into a single matrix with unified PAS rows.

    Inputs:
        mtx_paths: list of MatrixMarket files (parallel to dataset_ids and cb_paths)
        dataset_ids: dataset_ids[i] is the dataset that wrote mtx_paths[i]
        cb_paths: cb_paths[i] is the .cb.tsv (column index -> cb string) for mtx_paths[i]
        mapping_path: TSV from merge_pas_beds (or atlas_snap with same format):
            dataset_id\\told_pasnumber\\tnew_pas_id
        output_mtx, output_cb: write paths

    Output MTX:
        rows = unified PAS (sorted by new_pas_id numeric suffix or by mapping order)
        cols = all cells from all datasets concatenated
        values = original counts (multiple original rows mapping to same new_pas_id
            are summed via coo_matrix duplicate handling)

    Output CB:
        one CB per line, ordered to match the column ordering used in MTX.

    Returns: (output_mtx_path, output_cb_path)
    """
    mapping_path = Path(mapping_path)
    output_mtx = Path(output_mtx)
    output_cb = Path(output_cb)

    # Step 1: Load mapping. Key: (dataset_id, int(old_pasnumber)) -> new_pas_id.
    mapping: dict[tuple[str, int], str] = {}
    with open(mapping_path) as f:
        next(f)  # skip header: dataset_id\told_pasnumber\tnew_pas_id
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            ds_id, old_pas_str, new_pas_id = parts[0], parts[1], parts[2]
            try:
                old_pas_int = int(old_pas_str)
            except ValueError:
                # Non-integer pasnumber — store as-is with sentinel
                old_pas_int = hash(old_pas_str)
            mapping[(ds_id, old_pas_int)] = new_pas_id

    # Step 2: Collect all unique new_pas_ids and sort them.
    # Primary sort: numeric suffix of "PAS_N"; fall back to lexicographic for
    # atlas-style IDs that don't follow the "PAS_<int>" pattern.
    all_new_pas_ids: list[str] = sorted(
        set(mapping.values()),
        key=_pas_sort_key,
    )
    new_pas_id_to_row: dict[str, int] = {pid: i for i, pid in enumerate(all_new_pas_ids)}

    # Step 3–5: Iterate datasets, parse MTX files, build COO data.
    rows: list[int] = []
    cols: list[int] = []
    vals: list[int] = []
    merged_cbs: list[str] = []
    col_offset = 0

    for mtx_path, dataset_id, cb_path in zip(mtx_paths, dataset_ids, cb_paths):
        mtx_path = Path(mtx_path)
        cb_path = Path(cb_path)

        # Read CB list: column N in MTX = line N in cb.tsv (0-based).
        cb_list = [line.strip() for line in open(cb_path) if line.strip()]
        n_cols_local = len(cb_list)

        # Parse MTX. Detect format: if file starts with "%" comments, it's MatrixMarket
        # and the first non-% line is the dimension header to skip. If the file has no
        # comment lines, it's the legacy headerless COO format from peackcalling.py
        # and every line is a data triple.
        with open(mtx_path) as _peek:
            _first_char = _peek.read(1)
        has_mm_header = (_first_char == "%")

        skip_header = has_mm_header  # only skip first data line if MatrixMarket header is present
        with open(mtx_path) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith("%"):
                    continue
                if skip_header:
                    skip_header = False
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                try:
                    r, c, v = int(parts[0]), int(parts[1]), int(parts[2])
                except ValueError:
                    continue

                new_pas_id = mapping.get((dataset_id, r))
                if new_pas_id is None:
                    # This PAS didn't survive the merge (filtered out) — skip.
                    continue

                new_row = new_pas_id_to_row[new_pas_id]
                new_col = col_offset + (c - 1)  # 1-based to 0-based
                rows.append(new_row)
                cols.append(new_col)
                vals.append(v)

        col_offset += n_cols_local
        merged_cbs.extend(cb_list)

    # Steps 6–7: Build sparse matrix and sum duplicates (same PAS merged from
    # multiple original pasnumbers hitting the same cell column).
    n_rows = len(new_pas_id_to_row)
    n_cols = col_offset

    matrix = sp.coo_matrix(
        (
            np.array(vals, dtype=np.int32),
            (np.array(rows, dtype=np.int32), np.array(cols, dtype=np.int32)),
        ),
        shape=(n_rows, n_cols),
        dtype=np.int32,
    )
    # sum_duplicates is available via tocsr/tocoo round-trip when needed;
    # scipy coo_matrix does not expose sum_duplicates directly in all versions.
    matrix = matrix.tocsr().tocoo()

    # Step 8: Write MatrixMarket.
    output_mtx.parent.mkdir(parents=True, exist_ok=True)
    with open(output_mtx, "w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{n_rows} {n_cols} {matrix.nnz}\n")
        for r, c, v in zip(matrix.row, matrix.col, matrix.data):
            f.write(f"{r + 1} {c + 1} {v}\n")

    # Step 9: Write output CB file.
    output_cb.parent.mkdir(parents=True, exist_ok=True)
    with open(output_cb, "w") as f:
        for cb in merged_cbs:
            f.write(cb + "\n")

    return output_mtx, output_cb


def _pas_sort_key(pas_id: str) -> tuple[int, str]:
    """Sort key for PAS IDs.

    For "PAS_<int>" style IDs (from merge_pas_beds), sorts numerically by the
    integer suffix.  For atlas-style or other non-numeric IDs, falls back to
    lexicographic order via a high sentinel integer so they sort after numeric
    entries consistently.
    """
    parts = pas_id.rsplit("_", 1)
    if len(parts) == 2:
        try:
            return (int(parts[1]), "")
        except ValueError:
            pass
    # Non-numeric suffix — place after all numeric entries, then sort
    # lexicographically among themselves.
    return (10**18, pas_id)
