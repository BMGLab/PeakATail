"""Barcode-to-column index mapping for the count matrix.

Internal representation uses ``(sample_id_str, cb_int64)`` tuples as dict
keys.  ``cb_int64`` is the 2-bit packed integer produced by
:func:`ema.countmatrix.cb_encode.encode_cb`.  Integer-keyed dicts have
roughly the same memory footprint as the old string-keyed dicts for short
keys, but the hash operation on an int64 is a single CPU instruction whereas
hashing a 16-byte string walks the whole string — this matters when
:func:`get_index` is called once per read in the hot path.

**Backwards-compatible API**:

  - :func:`indexing(cb_str)` — module-level convenience, still takes a
    ``"sample_id_barcode"`` string, same as before.
  - :func:`get_mapping()` — returns ``dict[str, int]`` (string keys) so
    :mod:`ema.matrixfilter` never needs changing.
  - :func:`reset_index()` — resets the singleton.
  - :func:`get_default_index()` — returns the module-level singleton instance
    directly so callers can pass it explicitly to downstream functions.

**New API** (used by :mod:`ema.countmatrix.paswrite`):

  - :meth:`BarcodeIndex.get_indices_batch(cb_list)` — bulk resolve a list of
    CB strings in one call; vectorises the encode step via
    :func:`encode_cb_batch`.

.. note::

    These module-level functions wrap a process-global ``BarcodeIndex``
    singleton.  They work correctly under spawn-based multiprocessing (each
    worker has its own copy of the module) but **ARE NOT THREAD-SAFE**.
    New code should create explicit ``BarcodeIndex`` instances and pass them
    through the call chain instead of relying on the module singleton.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from ema.countmatrix.cb_encode import encode_cb, encode_cb_batch, decode_cb


def split_cb(cb_str: str) -> tuple[str, str]:
    """Split a composite ``"<sample_id>_<barcode>"`` CB into its two parts.

    **The split is on the LAST underscore, and that is what makes the
    composite bijective.**  The barcode half is a fixed-length nucleotide
    string (``ACGTN``) emitted by
    :func:`ema.countmatrix.read.read_check`, so it can never contain an
    underscore; the sample half therefore owns every underscore in the
    string and is recovered verbatim, whatever it contains.

    Splitting on the *first* underscore (what this module did until the
    ``fix/cellranger-input-compat`` repair) silently truncated any sample id
    containing ``_`` — e.g. ``"pbmc_10k_v3_AAACCCAAGAAACCCA"`` parsed as
    sample ``"pbmc"`` + barcode ``"10k_v3_AAACCCAAGAAACCCA"``, which
    :func:`~ema.countmatrix.cb_encode.encode_cb` maps to ``-1``.  Every cell
    of such a run then collapsed onto the single column ``("pbmc", -1)``.
    ``ema merge`` stamps ``RG = dataset_id`` (``DatasetManager._tag_bam_with_rg``)
    and ``samtools merge`` derives RG ids from file names, so underscore-bearing
    sample ids are the norm, not an edge case.

    Sanitising the sample id instead (``"_"`` -> ``"-"``) is **not** a valid
    alternative: it is not injective (``"a_b"`` and ``"a-b"`` both become
    ``"a-b"``, merging two samples' cells into one column) and it breaks the
    per-dataset column selector ``cb.startswith(f"{ds_id}_")`` used by
    ``ema/main.py`` and ``ema/reannotate.py``.

    Args:
        cb_str: Composite CB string ``"<sample_id>_<barcode>"``.

    Returns:
        ``(sample_id, barcode_str)``.  When *cb_str* contains no underscore
        the sample id is ``""`` and the whole string is treated as the
        barcode (degenerate case, preserved for backward compatibility).
    """
    idx = cb_str.rfind("_")
    if idx == -1:
        return ("", cb_str)
    return cb_str[:idx], cb_str[idx + 1 :]


def _parse_cb_str(cb_str: str) -> tuple[str, int]:
    """Split a composite CB string and 2-bit encode the barcode half.

    Thin wrapper over :func:`split_cb`; see that function for why the split
    is on the *last* underscore.

    Returns:
        ``(sample_id, cb_int)`` where *cb_int* is the 2-bit encoding of the
        barcode portion, or -1 for a barcode with non-ACGT characters.
    """
    sample_id, barcode_str = split_cb(cb_str)
    return sample_id, encode_cb(barcode_str)


class BarcodeIndex:
    """Encapsulates barcode-to-column mapping state for one run.

    Internally keys on ``(sample_id: str, cb_int: int)`` tuples for fast
    dict hashing.  All public methods accept and return string CBs in the
    legacy ``"<sample_id>_<barcode>"`` format so callers need no changes.
    """

    def __init__(self) -> None:
        # Internal dict: (sample_id_str, cb_int64) → column_index (1-based)
        self._cb_total: dict[tuple[str, int], int] = {}
        self._cb_index: int = 0

    # ------------------------------------------------------------------
    # Core API
    # ------------------------------------------------------------------

    def get_index(self, cb: str) -> int:
        """Return the column index for *cb*, assigning a new index on first use.

        Args:
            cb: Prefixed CB string ``"<sample_id>_<barcode>"``.

        Returns:
            1-based column index (consistent with MatrixMarket convention).
        """
        sample_id, cb_int = _parse_cb_str(cb)
        key: tuple[str, int] = (sample_id, cb_int)
        try:
            return self._cb_total[key]
        except KeyError:
            self._cb_index += 1
            self._cb_total[key] = self._cb_index
            return self._cb_index

    def get_index_int(self, sample_id: str, cb_int: int) -> int:
        """Fast-path index lookup that skips string parsing and re-encoding.

        Use this when the caller already holds a ``(sample_id, cb_int)`` pair
        (e.g. inside :meth:`get_indices_batch`).

        Args:
            sample_id: Sample identifier string.
            cb_int: 2-bit-encoded barcode integer from
                :func:`~ema.countmatrix.cb_encode.encode_cb`.

        Returns:
            1-based column index.
        """
        key: tuple[str, int] = (sample_id, cb_int)
        try:
            return self._cb_total[key]
        except KeyError:
            self._cb_index += 1
            self._cb_total[key] = self._cb_index
            return self._cb_index

    def get_indices_batch(self, cbs: list[str]) -> list[int]:
        """Bulk-resolve a list of CB strings to column indices.

        Encodes the barcode portion of every CB in one vectorised numpy call
        (:func:`~ema.countmatrix.cb_encode.encode_cb_batch`), then does a
        plain Python dict lookup per entry.  The encoding step is O(n) in
        numpy rather than O(n) in pure Python, saving roughly 3–5 µs per CB
        on typical hardware.

        Args:
            cbs: List of prefixed CB strings ``"<sample_id>_<barcode>"``.

        Returns:
            List of 1-based column indices in the same order as *cbs*.
        """
        if not cbs:
            return []

        # Split each CB into (sample_id, barcode_str) without re-encoding.
        # rfind, not find: the barcode half is underscore-free, the sample
        # half is not (see split_cb).
        split_idx = [cb.rfind("_") for cb in cbs]
        sample_ids: list[str] = []
        barcode_strs: list[str] = []
        for cb, idx in zip(cbs, split_idx):
            if idx == -1:
                sample_ids.append("")
                barcode_strs.append(cb)
            else:
                sample_ids.append(cb[:idx])
                barcode_strs.append(cb[idx + 1 :])

        # Vectorised 2-bit encode — all barcodes in one numpy call
        cb_ints: np.ndarray = encode_cb_batch(barcode_strs)

        # Resolve indices using the fast-path (no re-parsing)
        return [
            self.get_index_int(sid, int(ci))
            for sid, ci in zip(sample_ids, cb_ints)
        ]

    def reset(self) -> None:
        """Reset state — call between independent peak-calling runs."""
        self._cb_total.clear()
        self._cb_index = 0

    @property
    def mapping(self) -> dict[str, int]:
        """Return a ``dict[str, int]`` copy with reconstructed string keys.

        Reconstructs ``"<sample_id>_<barcode>"`` strings from internal tuple
        keys so :mod:`ema.matrixfilter` and other callers that expect the
        legacy string-keyed mapping need no changes.

        The CB length used for decode defaults to 16 (10x Chromium v2/v3).
        If the pipeline ever supports different CB lengths the decode call
        should use a stored per-instance length; for now 16 is correct for
        all supported chemistries.
        """
        result: dict[str, int] = {}
        for (sample_id, cb_int), col_idx in self._cb_total.items():
            barcode_str = decode_cb(cb_int, length=16)
            if sample_id:
                key = f"{sample_id}_{barcode_str}"
            else:
                key = barcode_str
            result[key] = col_idx
        return result


# ---------------------------------------------------------------------------
# Module-level singleton — maintains backward compatibility
# ---------------------------------------------------------------------------

_index = BarcodeIndex()


def indexing(cb: str) -> int:
    """Return column index for barcode *cb*. Assigns a new index on first use.

    Args:
        cb: Prefixed CB string ``"<sample_id>_<barcode>"``.

    Returns:
        1-based column index.
    """
    return _index.get_index(cb)


def reset_index() -> None:
    """Reset the module-level barcode index. Call between independent runs."""
    _index.reset()


def get_mapping() -> dict[str, int]:
    """Return a copy of the current barcode-to-index mapping.

    Keys are reconstructed ``"<sample_id>_<barcode>"`` strings so callers
    that expect the legacy string-keyed dict work unchanged.
    """
    return _index.mapping


def get_default_index() -> BarcodeIndex:
    """Return the module-level singleton BarcodeIndex instance.

    Exposes direct access to the singleton so new call sites can pass it
    explicitly as an ``index`` parameter rather than using the module-level
    shim functions.  Under spawn-based multiprocessing each worker gets its
    own copy, so this is safe.  Under threading it is NOT thread-safe.

    New code should prefer creating a fresh :class:`BarcodeIndex` instance
    per worker/invocation rather than sharing this singleton.

    Returns:
        The process-global :class:`BarcodeIndex` singleton.
    """
    return _index
