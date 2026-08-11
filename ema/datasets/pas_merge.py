from pathlib import Path
import subprocess
import scipy.sparse as sp
import numpy as np
import pandas as pd


def _encode_pas_key(dataset_id: str, strand: str | None, pasnumber: str) -> str:
    """Encode a per-input PAS identity into a single BED col-4 token.

    B1: the same ``pasnumber`` is minted from 1 independently on each strand
    within a dataset, so ``(dataset_id, pasnumber)`` collides across strands and
    counts get routed to the wrong unified PAS row. Including the strand makes
    the key ``(dataset_id, strand, pasnumber)``. When ``strand`` is None the
    legacy 2-part encoding is used (backward compatible).
    """
    if strand:
        return f"{dataset_id}::{strand}::{pasnumber}"
    return f"{dataset_id}::{pasnumber}"


def _decode_pas_key(token: str) -> tuple[str, str, str]:
    """Inverse of :func:`_encode_pas_key` → ``(dataset_id, strand, pasnumber)``.

    ``strand`` is ``""`` for the legacy 2-part encoding. ``dataset_id`` may
    itself contain ``::`` (we split from the right).
    """
    parts = token.rsplit("::", 2)
    if len(parts) == 3:
        return parts[0], parts[1], parts[2]
    if len(parts) == 2:
        return parts[0], "", parts[1]
    return token, "", ""


def pas_uid_of(chrom: str, start, end, strand: str) -> str:
    """Content-addressed stable PAS id (E1): ``chrom:pos:strand``.

    ``pos`` is the strand-aware 3' end of the interval — the polyadenylation
    site itself: the last base (``end-1``) on ``+`` and the first base
    (``start``) on ``-``. This is stable across runs and re-minting of the
    run-local integer ``new_pas_id``, so cross-run joins can key on it.
    Falls back to ``chrom:end:strand`` if coordinates aren't integers.
    """
    try:
        s, e = int(start), int(end)
    except (TypeError, ValueError):
        return f"{chrom}:{end}:{strand}"
    pos = e - 1 if strand == "+" else s
    return f"{chrom}:{pos}:{strand}"


def merge_pas_beds(
    bed_paths: list[str | Path],
    dataset_ids: list[str],
    output_dir: str | Path,
    gap: int = 100,
    strands: list[str] | None = None,
) -> tuple[Path, Path]:
    """Merge per-dataset PAS BED files into a unified strand-aware coordinate set.

    Inputs:
        bed_paths: list of BED file paths (one per dataset peak-calling run)
        dataset_ids: parallel list — dataset_ids[i] owns bed_paths[i]
        output_dir: where to write outputs
        gap: bedtools merge -d distance (default 100bp)
        strands: optional parallel list — strands[i] ("+"/"-") is the strand of
            bed_paths[i]. When supplied, the count-routing key becomes
            ``(dataset_id, strand, pasnumber)`` (bug B1 fix); when a dataset
            contributes a separate pos and neg BED under the *same* dataset_id,
            this prevents pos/neg PAS #N from colliding.

    Outputs:
        merged_bed: chrom\\tstart\\tend\\tnew_pas_id\\tscore\\tstrand
        mapping_path: TSV with header
            "dataset_id\\tstrand\\told_pasnumber\\tnew_pas_id"
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

    # Step 1: Concatenate all BED files, encoding dataset_id[::strand]::pasnumber
    # in col 4 (B1: strand included when provided to avoid cross-strand key
    # collisions in count routing).
    if strands is not None and len(strands) != len(bed_paths):
        raise ValueError(
            f"strands length ({len(strands)}) must match bed_paths ({len(bed_paths)})"
        )
    with open(cat_path, "w") as out:
        for i, (bed_path, dataset_id) in enumerate(zip(bed_paths, dataset_ids)):
            strand_hint = strands[i] if strands is not None else None
            with open(bed_path) as f:
                for line in f:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) < 4:
                        continue
                    # parts[3] is the original pasnumber (integer string).
                    # Prefer the row's own strand column (parts[5]) when present,
                    # falling back to the per-file strand hint.
                    row_strand = parts[5] if len(parts) >= 6 else strand_hint
                    if strands is None:
                        row_strand = None
                    parts[3] = _encode_pas_key(dataset_id, row_strand, parts[3])
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
    # (dataset_id, strand, old_pasnumber, new_pas_id)
    mapping_rows: list[tuple[str, str, str, str]] = []
    # E1: (new_pas_id, pas_uid) content-addressed stable id sidecar.
    uid_rows: list[tuple[str, str]] = []
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

            # Use plain integer string so BED col 4 matches MTX row indices
            # (integer-keyed). annotate.py joins on these directly.
            new_pas_id = str(new_id)
            new_id += 1

            # Write one mapping row per original entry.
            for entry in collected.split(","):
                entry = entry.strip()
                if "::" not in entry:
                    continue
                dataset_id, entry_strand, old_pasnumber = _decode_pas_key(entry)
                mapping_rows.append(
                    (dataset_id, entry_strand, old_pasnumber, new_pas_id)
                )

            merged_lines.append(
                f"{chrom}\t{start}\t{end}\t{new_pas_id}\t{score}\t{strand}"
            )
            uid_rows.append((new_pas_id, pas_uid_of(chrom, start, end, strand)))

    # Step 5: Write final merged BED.
    with open(merged_path, "w") as f:
        f.write("\n".join(merged_lines))
        if merged_lines:
            f.write("\n")

    # Step 6: Write mapping TSV (with strand column — B1).
    with open(mapping_path, "w") as f:
        f.write("dataset_id\tstrand\told_pasnumber\tnew_pas_id\n")
        for dataset_id, strand, old_pasnumber, new_pas_id in mapping_rows:
            f.write(f"{dataset_id}\t{strand}\t{old_pasnumber}\t{new_pas_id}\n")

    # Step 6b (E1): content-addressed stable id sidecar.
    uid_path = output_dir / "pas_uid.tsv"
    with open(uid_path, "w") as f:
        f.write("new_pas_id\tpas_uid\n")
        for new_pas_id, pas_uid in uid_rows:
            f.write(f"{new_pas_id}\t{pas_uid}\n")

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
    strands: list[str] | None = None,
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

    # Step 1: Load mapping. Key: (dataset_id, strand, int(old_pasnumber)) ->
    # new_pas_id (B1: strand disambiguates pos/neg PAS #N within a dataset).
    # Header-aware so both the merge mapping (has a ``strand`` column) and the
    # atlas mapping are parsed correctly. ``strand`` is only used for keying when
    # the caller passes ``strands`` (parallel to mtx_paths); otherwise it is
    # collapsed to "" on both sides for backward-compatible last-wins behaviour.
    if strands is not None and len(strands) != len(mtx_paths):
        raise ValueError(
            f"strands length ({len(strands)}) must match mtx_paths ({len(mtx_paths)})"
        )
    strand_mode = strands is not None
    mapping: dict[tuple[str, str, int], str] = {}
    with open(mapping_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        i_ds = col.get("dataset_id", 0)
        i_strand = col.get("strand")  # None if absent (atlas legacy layout)
        i_old = col.get("old_pasnumber", 1 if i_strand is None else 2)
        i_new = col.get("new_pas_id", 2 if i_strand is None else 3)
        _need = max(x for x in (i_ds, i_old, i_new, i_strand) if x is not None)
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= _need:
                continue
            ds_id = parts[i_ds]
            strand_val = (
                parts[i_strand]
                if (strand_mode and i_strand is not None and i_strand < len(parts))
                else ""
            )
            old_pas_str, new_pas_id = parts[i_old], parts[i_new]
            try:
                old_pas_int = int(old_pas_str)
            except ValueError:
                # Non-integer pasnumber — store as-is with sentinel
                old_pas_int = hash(old_pas_str)
            mapping[(ds_id, strand_val, old_pas_int)] = new_pas_id

    # Step 2: Collect all unique new_pas_ids and sort them.
    # Primary sort: numeric suffix of "PAS_N"; fall back to lexicographic for
    # atlas-style IDs that don't follow the "PAS_<int>" pattern.
    all_new_pas_ids: list[str] = sorted(
        set(mapping.values()),
        key=_pas_sort_key,
    )
    new_pas_id_to_row: dict[str, int] = {pid: i for i, pid in enumerate(all_new_pas_ids)}

    # Step 3–5: Iterate datasets, parse MTX files, build COO data.
    #
    # Columns are keyed by (dataset_id, cb) so a cell that contributes a
    # separate mtx per strand (the pos and neg peak-call matrices share ONE
    # cb.tsv per BAM — see main.py "same cb.tsv path as pos") collapses into a
    # SINGLE column carrying both its +-strand and −-strand PAS (which B1's
    # strand-keyed mapping already routes to distinct rows). Without this the
    # same barcode landed in two columns — one all-pos, one all-neg — so every
    # two-strand cell was double-counted as two "cells" downstream (inflated
    # n_obs, split clustering). coo duplicate-summing then merges the two
    # sparse column vectors losslessly.
    # Vectorized per file (matrices are large — no Python per-triple loop):
    # read the whole MTX with pandas, map old-PAS->unified-row and
    # local-col->canonical-(dataset,cb)-col with numpy, accumulate arrays.
    row_chunks: list[np.ndarray] = []
    col_chunks: list[np.ndarray] = []
    val_chunks: list[np.ndarray] = []
    merged_cbs: list[str] = []
    col_key_to_idx: dict[tuple[str, str], int] = {}

    for _mtx_i, (mtx_path, dataset_id, cb_path) in enumerate(
        zip(mtx_paths, dataset_ids, cb_paths)
    ):
        mtx_path = Path(mtx_path)
        cb_path = Path(cb_path)
        strand_val = strands[_mtx_i] if strand_mode else ""

        # CB list: column N in MTX = line N in cb.tsv (0-based). Resolve each
        # local column to its canonical (dataset_id, cb) column index once,
        # into a lookup array (local 0-based col -> canonical col), so pos/neg
        # of the same cell reuse ONE column.
        cb_list = [line.strip() for line in open(cb_path) if line.strip()]
        local_canon = np.empty(len(cb_list), dtype=np.int64)
        for _local_c, _cb in enumerate(cb_list):
            _key = (dataset_id, _cb)
            _idx = col_key_to_idx.get(_key)
            if _idx is None:
                _idx = len(merged_cbs)
                col_key_to_idx[_key] = _idx
                merged_cbs.append(_cb)
            local_canon[_local_c] = _idx

        # Per-file old-pasnumber -> unified row map (dataset_id + strand fixed
        # for this file), as parallel numpy arrays for a vectorized lookup.
        pairs = [
            (old_int, new_pas_id_to_row[new_pas_id])
            for (ds, st, old_int), new_pas_id in mapping.items()
            if ds == dataset_id and st == strand_val
        ]
        if not pairs or len(cb_list) == 0:
            continue
        old_arr = np.fromiter((p[0] for p in pairs), dtype=np.int64, count=len(pairs))
        row_arr = np.fromiter((p[1] for p in pairs), dtype=np.int64, count=len(pairs))

        triples = _read_mtx_triples(mtx_path)  # (N,3) int64: pas, col, count
        if triples.size == 0:
            continue
        r = triples[:, 0]
        c = triples[:, 1] - 1  # 1-based -> 0-based local column
        v = triples[:, 2]

        # Map PAS rows via a Series (NaN where the PAS didn't survive the merge).
        new_row = pd.Series(r).map(pd.Series(row_arr, index=old_arr)).to_numpy()
        keep = ~np.isnan(new_row) & (c >= 0) & (c < len(cb_list))
        if not keep.any():
            continue
        row_chunks.append(new_row[keep].astype(np.int32))
        col_chunks.append(local_canon[c[keep]].astype(np.int32))
        val_chunks.append(v[keep].astype(np.int32))

    # Steps 6–7: Build sparse matrix and sum duplicates (same PAS merged from
    # multiple original pasnumbers hitting the same cell column, AND the pos/neg
    # column collapse above).
    n_rows = len(new_pas_id_to_row)
    n_cols = len(merged_cbs)
    all_rows = np.concatenate(row_chunks) if row_chunks else np.empty(0, dtype=np.int32)
    all_cols = np.concatenate(col_chunks) if col_chunks else np.empty(0, dtype=np.int32)
    all_vals = np.concatenate(val_chunks) if val_chunks else np.empty(0, dtype=np.int32)

    matrix = sp.coo_matrix(
        (all_vals, (all_rows, all_cols)),
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


def _read_mtx_triples(path: Path) -> np.ndarray:
    """Read a MatrixMarket-or-legacy COO ``.mtx`` into an ``(N, 3)`` int64
    array ``[pas_row, cell_col, count]`` using the pandas C parser (no
    per-triple Python loop — these matrices are large).

    Handles both formats concat_matrices accepts: MatrixMarket (``%`` comment
    lines + a ``n_rows n_cols nnz`` dimension header) and the legacy headerless
    triple format written by peakcalling. Returns an empty ``(0, 3)`` array for
    an empty/degenerate file rather than raising.
    """
    with open(path) as _peek:
        has_mm_header = _peek.read(1) == "%"
    try:
        df = pd.read_csv(
            path, sep=r"\s+", comment="%", header=None,
            skip_blank_lines=True, dtype=np.int64,
        )
    except (ValueError, pd.errors.EmptyDataError):
        return np.empty((0, 3), dtype=np.int64)
    if df.shape[0] == 0 or df.shape[1] < 3:
        return np.empty((0, 3), dtype=np.int64)
    if has_mm_header:
        df = df.iloc[1:]  # drop the "n_rows n_cols nnz" dimension line
    return df.iloc[:, :3].to_numpy(dtype=np.int64)


def _pas_sort_key(pas_id: str) -> tuple[int, str]:
    """Sort key for PAS IDs.

    For "PAS_<int>" style IDs (from merge_pas_beds), sorts numerically by the
    integer suffix.  For atlas-style or other non-numeric IDs, falls back to
    lexicographic order via a high sentinel integer so they sort after numeric
    entries consistently.
    """
    # Bare integer (current format from merge_pas_beds and atlas_snap)
    try:
        return (int(pas_id), "")
    except ValueError:
        pass
    # "PAS_<int>" style (legacy)
    parts = pas_id.rsplit("_", 1)
    if len(parts) == 2:
        try:
            return (int(parts[1]), "")
        except ValueError:
            pass
    # Non-numeric — place after all numeric entries.
    return (10**18, pas_id)
