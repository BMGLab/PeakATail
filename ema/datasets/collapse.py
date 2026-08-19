"""A4: collapse samtools-merge RG-suffix run tags in a unified count matrix.

``samtools merge`` renames colliding ``@RG`` IDs by appending an 8-hex suffix
(``GSM123-StageI-4B7F9BA8``, one per input run). Peak calling then reads those
per-run RG tags, so ``cb = f"{rg}_{barcode}"`` gets a PER-RUN prefix and the
runs of one library stop pooling into a single cell namespace. Peaks (matrix
rows) are unaffected — only the cell columns are mis-split.

This module rewrites the unified count matrix so cells pool correctly, WITHOUT
re-running peak calling:

  * strip the trailing ``-<8hex>`` run tag from each cb's library prefix
      ``GSM123-StageI-4B7F9BA8_AAAC…``  ->  ``GSM123-StageI_AAAC…``
  * SUM the columns that now share a canonical cb (same barcode across a
    library's runs = the same cell — exactly what ``before``-merge should pool)
  * DIFFERENT libraries keep distinct prefixes, so they never merge.

The pure functions :func:`canonical_cb` and :func:`collapse_columns` carry all
the logic and are unit-testable without touching the filesystem;
:func:`collapse_run` is the thin on-disk orchestration used by ``ema collapse``.
Promoted from ``scripts/collapse_run_suffixes.py``.
"""
from __future__ import annotations

import re
import shutil
from pathlib import Path

# samtools RG-collision suffix: hyphen + exactly 8 hex chars at the end of the id.
_RUN_SUFFIX = re.compile(r"-[0-9A-Fa-f]{8}$")


def canonical_cb(cb: str) -> str:
    """Map a per-run cb to its library-level canonical cb.

    ``GSM-Stage-<8hex>_<barcode>`` -> ``GSM-Stage_<barcode>``. A cb without a
    ``_`` separator (no barcode) is returned unchanged. Only the prefix (before
    the LAST ``_``) has the run suffix stripped, so barcodes that happen to end
    in ``-<8hex>`` are never touched, and library prefixes that themselves
    contain ``_`` (``lib_A-4B7F9BA8_<barcode>``) are stripped correctly instead
    of being left un-pooled.
    """
    i = cb.rfind("_")
    if i == -1:
        return cb
    prefix, barcode = cb[:i], cb[i + 1:]
    return f"{_RUN_SUFFIX.sub('', prefix)}_{barcode}"


def collapse_columns(matrix, cbs):
    """Collapse matrix columns that share a canonical cb by summing them.

    Args:
        matrix: sparse ``PAS x cells`` matrix (scipy.sparse), ``matrix.shape[1]``
            must equal ``len(cbs)``.
        cbs: per-column cell barcodes, aligned to ``matrix`` columns.

    Returns:
        ``(collapsed, uniq_cbs)`` where ``collapsed`` is a ``PAS x n_uniq``
        sparse matrix (columns summed across shared canonical cbs) and
        ``uniq_cbs`` is the first-seen-ordered list of canonical cbs.

    Raises:
        ValueError: if ``matrix.shape[1] != len(cbs)``.
    """
    import numpy as np
    import scipy.sparse as sp

    if matrix.shape[1] != len(cbs):
        raise ValueError(
            f"matrix has {matrix.shape[1]} columns but {len(cbs)} cbs"
        )
    canon = [canonical_cb(c) for c in cbs]
    uniq = list(dict.fromkeys(canon))          # preserve first-seen order
    col_of = {c: j for j, c in enumerate(uniq)}

    # Aggregation matrix A (n_cells x n_uniq) with a single 1 per row; then
    # matrix @ A sums the columns sharing a canonical cb — a clean sparse
    # collapse that never densifies.
    rows = np.arange(len(canon))
    cols = np.fromiter((col_of[c] for c in canon), dtype=np.int64, count=len(canon))
    A = sp.csc_matrix(
        (np.ones(len(canon), dtype=matrix.dtype), (rows, cols)),
        shape=(len(canon), len(uniq)),
    )
    collapsed = (matrix.tocsr() @ A).tocsc()
    return collapsed, uniq


def _library(cb: str) -> str:
    """Library prefix = everything before the LAST ``_``.

    The barcode half of a composite cb never contains ``_`` (see
    :func:`ema.countmatrix.indexing.split_cb`), so the last underscore is the
    separator and library ids containing ``_`` survive intact.
    """
    return cb.rsplit("_", 1)[0]


def collapse_run(in_run: str | Path, out_run: str | Path) -> dict:
    """Rewrite a completed before-merge run so per-run cbs pool per library.

    Reads ``<in_run>/unified/concatenated.mtx`` + ``concatenated_cbs.tsv``,
    collapses run-suffixed columns, and writes a corrected run dir at
    ``out_run`` (unified matrix + cbs, plus the RG-independent posbed/negbed
    copied verbatim so a re-annotate step can branch from it).

    Returns:
        dict summary: ``n_cells_before``, ``n_cells_after``, ``n_pooled``,
        ``n_libraries_before``, ``n_libraries_after``.

    Raises:
        FileNotFoundError: if a required input artifact is missing.
        ValueError: on a matrix/cb column-count mismatch.
    """
    import scipy.io as sio

    in_run = Path(in_run)
    out_run = Path(out_run)
    umtx = in_run / "unified" / "concatenated.mtx"
    ucbs = in_run / "unified" / "concatenated_cbs.tsv"
    required = [umtx, ucbs, in_run / "posbed.bed", in_run / "negbed.bed"]
    for p in required:
        if not p.exists():
            raise FileNotFoundError(
                f"missing {p} — is --in-run a completed multi-sample run?"
            )

    cbs = [ln.strip() for ln in ucbs.read_text().splitlines() if ln.strip()]
    matrix = sio.mmread(str(umtx)).tocsc()          # PAS (rows) x cells (cols)
    collapsed, uniq = collapse_columns(matrix, cbs)

    libs_before = {_library(c) for c in cbs}
    libs_after = {_library(c) for c in uniq}

    out_unified = out_run / "unified"
    out_unified.mkdir(parents=True, exist_ok=True)
    sio.mmwrite(str(out_unified / "concatenated.mtx"),
                collapsed.tocoo().astype(int), field="integer")
    (out_unified / "concatenated_cbs.tsv").write_text("\n".join(uniq) + "\n")
    for name in ("posbed.bed", "negbed.bed"):
        shutil.copy(in_run / name, out_run / name)

    return {
        "n_cells_before": len(cbs),
        "n_cells_after": len(uniq),
        "n_pooled": len(cbs) - len(uniq),
        "n_libraries_before": len(libs_before),
        "n_libraries_after": len(libs_after),
    }
