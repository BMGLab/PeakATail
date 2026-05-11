"""Functions for writing BED and MatrixMarket count-matrix output.

``matrix_write`` is the hot-path function called once per emitted peak.
It uses :meth:`~ema.countmatrix.indexing.BarcodeIndex.get_indices_batch`
to bulk-resolve all CB strings in the peak's ``cb_dict`` in a single
vectorised numpy call, rather than calling :func:`indexing` (the string-parse
+ encode + dict-lookup path) once per CB per peak.

The public API is unchanged: external callers pass the same ``cb_dict`` and
``pasnumber`` arguments as before.  New callers may additionally pass an
explicit ``index`` :class:`~ema.countmatrix.indexing.BarcodeIndex` instance
so the module-level singleton is not mutated -- important for parallel
workers that each maintain their own independent index.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ema.countmatrix.indexing import BarcodeIndex

strand_char: dict[bool, str] = {True: "-", False: "+"}
score: int = 0


def pas_write(
    chro: str,
    peak_start: int,
    l_end: int,
    strand: bool,
    pasnumber: int | str,
    output,
) -> None:
    """Write a single peak record to a BED-format file.

    Args:
        chro: Chromosome name.
        peak_start: Start coordinate of the peak (0-based, BED convention).
        l_end: Last read-end coordinate contributing to this peak.
        strand: ``True`` for reverse strand (``-``), ``False`` for forward (``+``).
        pasnumber: Peak identifier written to the BED name field.
        output: Writable file-like object.
    """
    bed_start = min(peak_start, l_end)
    bed_end = max(peak_start, l_end)
    peak_bed = (
        f"{chro}\t{bed_start}\t{bed_end}\t{pasnumber}\t{score}\t{strand_char[strand]}\n"
    )
    output.write(peak_bed)


def matrix_write(
    cb_dict: dict,
    pasnumber: int,
    output,
    index: "BarcodeIndex | None" = None,
) -> None:
    """Write count-matrix rows for one peak in MatrixMarket format.

    Uses :meth:`~ema.countmatrix.indexing.BarcodeIndex.get_indices_batch` to
    resolve all CB strings in *cb_dict* in a single batched call, avoiding
    repeated per-CB string parsing and encoding overhead.

    Args:
        cb_dict: Mapping of CB string -> read count for this peak.
        pasnumber: 1-based peak row index in the count matrix.
        output: Writable file-like object.
        index: Optional :class:`~ema.countmatrix.indexing.BarcodeIndex`
            instance to use for CB-to-column resolution.  When ``None``
            (default), the module-level singleton from
            :mod:`ema.countmatrix.indexing` is used -- preserving backward
            compatibility for all existing callers.  New code (e.g. parallel
            workers) should pass an explicit instance so the singleton is not
            touched.
    """
    if not cb_dict:
        return

    if index is None:
        from ema.countmatrix.indexing import _index as index  # type: ignore[assignment]

    cbs: list[str] = list(cb_dict.keys())
    counts: list[int] = [cb_dict[cb] for cb in cbs]

    # Bulk-resolve all CB strings to column indices in one numpy-accelerated call
    cols: list[int] = index.get_indices_batch(cbs)

    lines: list[str] = [
        f"{pasnumber} {col} {count}\n"
        for col, count in zip(cols, counts)
    ]
    output.write("".join(lines))
