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

#: Columns of the per-PAS poly(A) support sidecar (``pas_support.tsv``).
#: BED6 stays BED6 — the raw clip-read count and the ``-F 3844`` counts live
#: here instead of being packed into the name field (which downstream
#: parsers read as the PAS id) or the score column (one number only).
SUPPORT_COLUMNS: tuple[str, ...] = (
    "pas_id",             # BED column 4 — joins the sidecar to the BED
    "clip_reads",         # raw poly(A) clip reads at this PAS
    "clip_umis",          # distinct (CB, UMI) molecules == BED column 5
    "clip_reads_f3844",   # clip reads passing samtools -F 3844
    "clip_umis_f3844",    # distinct molecules among those reads
    "window_reads",       # reads counted into this PAS's matrix row
    "tier",               # 1 = clip-seeded cluster, 2 = coverage candidate
)

#: peakAtail-prime, ``--pas-features on``: two more columns written by the
#: caller itself, describing the clip cluster's SHAPE rather than its size.
#: They cost one ``len()`` and one subtraction per emitted PAS -- the cluster
#: member positions are already in hand -- and they are not derivable from
#: anything else the run writes, which is why they are emitted here and not
#: reconstructed later.  APPEND ONLY: they go after every v2 column, so a
#: reader indexing the v2 columns by position is unaffected.
#: See :mod:`ema.countmatrix.pas_features` for the rest of the feature set,
#: which is appended at the internal-priming seam.
CALL_FEATURE_COLUMNS: tuple[str, ...] = (
    "clip_positions",     # distinct poly(A) clip POSITIONS backing this PAS
    "clip_span",          # bp between the first and last of them (0 if <2)
    "clip_offset_mean",   # read-weighted mean(member - call), transcript
                          # orientation, + = downstream (TASK E); NA on tier 2
)


#: Accepted values of ``--pas-features``.  ``"off"`` is v2.
PAS_FEATURE_MODES: tuple[str, ...] = ("off", "on")
#: The v2-compatibility value of ``--pas-features``.
V2_PAS_FEATURES: str = "off"


def support_columns(features: bool = False) -> tuple[str, ...]:
    """Sidecar header for this run.

    *features* is ``--pas-features on``.  ``False`` is the v2 column set,
    byte-for-byte.
    """
    return SUPPORT_COLUMNS + CALL_FEATURE_COLUMNS if features else SUPPORT_COLUMNS


def support_path_for(bedfilepath) -> str:
    """Sidecar path for a caller BED (``x.bed`` -> ``x.support.tsv``)."""
    s = str(bedfilepath)
    if s.endswith(".bed"):
        s = s[:-4]
    return s + ".support.tsv"


def open_support(bedfilepath, features: bool = False):
    """Open the sidecar for *bedfilepath* and write its header row.

    *features* must be the same value every :func:`support_write` call for
    this handle passes, or the header and the rows disagree; it is threaded
    explicitly (rather than read from ``variable_config`` here) because the
    parallel callers run in spawned children whose globals are back at module
    defaults -- the same hazard that made ``--read-geometry`` silently differ
    between the tile path and the monolithic path.  ``tests/
    test_pas_features.py::test_every_calling_path_agrees_on_the_sidecar_header``
    pins the agreement.
    """
    fh = open(support_path_for(bedfilepath), "w")
    fh.write("\t".join(support_columns(features)) + "\n")
    return fh


def support_write(output, pasnumber, support: dict, features: bool = False) -> None:
    """Append one sidecar row.  *support* is a dict as produced by
    :meth:`~ema.countmatrix.polya.ClipSeeder.flush` (``support_out=``).

    With *features* false the row is v2's seven columns, byte-for-byte.
    """
    if output is None:
        return
    row = (
        f"{pasnumber}\t{support['clip_reads']}\t{support['clip_umis']}\t"
        f"{support['clip_reads_f3844']}\t{support['clip_umis_f3844']}\t"
        f"{support['window_reads']}\t{support['tier']}"
    )
    if features:
        row += "\t" + "\t".join(
            str(support.get(c, "NA")) for c in CALL_FEATURE_COLUMNS
        )
    output.write(row + "\n")


def pas_write(
    chro: str,
    peak_start: int,
    l_end: int,
    strand: bool,
    pasnumber: int | str,
    output,
    score: int = 0,
) -> None:
    """Write a single peak record to a BED-format file.

    Args:
        chro: Chromosome name.
        peak_start: Start coordinate of the peak (0-based, BED convention).
        l_end: Last read-end coordinate contributing to this peak.
        strand: ``True`` for reverse strand (``-``), ``False`` for forward (``+``).
        pasnumber: Peak identifier written to the BED name field.
        output: Writable file-like object.
        score: BED column 5. Historically hardcoded to 0; with
            ``--polya-evidence on`` it carries the PAS's poly(A) clip-read
            support (0 == coverage-only, >=1 == clip-supported), which every
            existing downstream reader already names ``score`` and ignores.
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
