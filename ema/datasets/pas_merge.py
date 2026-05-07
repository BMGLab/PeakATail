from pathlib import Path
import subprocess
import scipy.sparse as sp
import numpy as np


def merge_pas_beds(bed_paths: list[str | Path], output_dir: Path, gap: int = 100) -> tuple[Path, Path]:
    """
    Merge PAS BED files from multiple datasets into a unified coordinate set.

    Concatenates all BED files, runs bedtools merge -d {gap}, assigns new
    sequential PAS IDs, and writes a mapping file (old_id -> new_id).

    Returns (merged_bed_path, mapping_path).
    """
    output_dir = Path(output_dir)
    cat_path = output_dir / "multi_sample_cat.bed"
    merged_path = output_dir / "multi_sample_merged.bed"
    mapping_path = output_dir / "multi_sample_pas_mapping.tsv"

    # Concatenate all BED files
    with open(cat_path, "w") as out:
        for bed in bed_paths:
            with open(bed) as f:
                out.write(f.read())

    # Sort concatenated BED
    sorted_cat = output_dir / "multi_sample_cat_sorted.bed"
    subprocess.run(
        ["sort", "-k1,1", "-k2,2n", str(cat_path)],
        stdout=open(sorted_cat, "w"),
        check=True,
    )

    # bedtools merge -d {gap} -c 4,5,6 -o first,first,first
    # Column 4 = PAS ID, 5 = score, 6 = strand — keep first on merge
    subprocess.run(
        ["bedtools", "merge", "-d", str(gap), "-i", str(sorted_cat),
         "-c", "4,5,6", "-o", "first,first,first"],
        stdout=open(merged_path, "w"),
        check=True,
    )

    # Re-number PAS IDs and build old->new mapping
    old_to_new: dict[str, str] = {}
    renumbered_lines: list[str] = []
    new_id = 1

    with open(merged_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            old_id = parts[3]
            new_pas = f"PAS_{new_id}"
            old_to_new[old_id] = new_pas
            parts[3] = new_pas
            renumbered_lines.append("\t".join(parts))
            new_id += 1

    with open(merged_path, "w") as f:
        f.write("\n".join(renumbered_lines) + "\n")

    with open(mapping_path, "w") as f:
        f.write("old_pas_id\tnew_pas_id\n")
        for old, new in old_to_new.items():
            f.write(f"{old}\t{new}\n")

    # Clean up temp files
    cat_path.unlink(missing_ok=True)
    sorted_cat.unlink(missing_ok=True)

    return merged_path, mapping_path


def concat_matrices(
    mtx_paths: list[str | Path],
    cb_list_paths: list[str | Path],
    mapping_path: str | Path,
    output_path: str | Path,
    output_cb_path: str | Path,
) -> tuple[Path, Path]:
    """
    Concatenate per-dataset MTX files into a single merged matrix.

    Rows = merged PAS (remapped via mapping_path).
    Columns = all cells across datasets (already unique via {dataset_id}_{CB} prefix).

    Returns (output_mtx_path, output_cb_path).
    """
    mapping_path = Path(mapping_path)
    output_path = Path(output_path)
    output_cb_path = Path(output_cb_path)

    # Load PAS ID mapping: old -> new (0-based index in merged set)
    old_to_new: dict[str, str] = {}
    with open(mapping_path) as f:
        next(f)  # skip header
        for line in f:
            old, new = line.rstrip("\n").split("\t")
            old_to_new[old] = new

    new_ids_sorted = sorted(old_to_new.values(), key=lambda x: int(x.split("_")[1]))
    new_id_to_row: dict[str, int] = {pid: i for i, pid in enumerate(new_ids_sorted)}
    n_rows = len(new_ids_sorted)

    all_rows: list[int] = []
    all_cols: list[int] = []
    all_vals: list[int] = []
    all_cbs: list[str] = []
    col_offset = 0

    for mtx_path, cb_path in zip(mtx_paths, cb_list_paths):
        mtx_path = Path(mtx_path)
        cb_path = Path(cb_path)

        if not mtx_path.exists() or not cb_path.exists():
            continue

        # Load CB list for this dataset
        with open(cb_path) as f:
            cbs = [line.strip() for line in f if line.strip()]
        n_cols = len(cbs)
        all_cbs.extend(cbs)

        # Parse MTX (MatrixMarket COO, 1-based)
        # Also need to map PAS IDs — load from PAS ID list alongside the MTX
        pas_id_path = mtx_path.with_suffix(".pas_ids.tsv")
        if not pas_id_path.exists():
            # Try to load from BED file with same stem
            col_offset += n_cols
            continue

        with open(pas_id_path) as f:
            mtx_pas_ids = [line.strip() for line in f if line.strip()]

        with open(mtx_path) as f:
            for line in f:
                if line.startswith("%"):
                    continue
                parts = line.split()
                if len(parts) == 3:
                    try:
                        r, c, v = int(parts[0]), int(parts[1]), int(parts[2])
                    except ValueError:
                        continue
                    old_pas = mtx_pas_ids[r - 1] if r - 1 < len(mtx_pas_ids) else None
                    if old_pas and old_pas in old_to_new:
                        new_pas = old_to_new[old_pas]
                        new_row = new_id_to_row[new_pas]
                        all_rows.append(new_row)
                        all_cols.append(col_offset + c - 1)
                        all_vals.append(v)

        col_offset += n_cols

    n_total_cols = col_offset
    merged = sp.coo_matrix(
        (all_vals, (all_rows, all_cols)),
        shape=(n_rows, n_total_cols),
        dtype=np.int32,
    )

    # Write MatrixMarket
    with open(output_path, "w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{n_rows} {n_total_cols} {merged.nnz}\n")
        coo = merged.tocoo()
        for r, c, v in zip(coo.row, coo.col, coo.data):
            f.write(f"{r + 1} {c + 1} {v}\n")

    # Write merged CB list
    with open(output_cb_path, "w") as f:
        for cb in all_cbs:
            f.write(cb + "\n")

    return output_path, output_cb_path
