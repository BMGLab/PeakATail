"""2-bit DNA encoding for cell barcodes.

A=0, C=1, G=2, T=3 — a 16 nt CB fits in 32 bits, well within an int64.

Provides:
  - Single-CB encode/decode (encode_cb / decode_cb)
  - Vectorised batch encoding via numpy (encode_cb_batch)

Non-ATCG characters (e.g. N, ambiguous IUPAC) return the sentinel value -1
so callers can filter invalid reads without raising exceptions.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Lookup table: ASCII ordinal → 2-bit value (0–3), -1 for non-ATCG
# ---------------------------------------------------------------------------

_BASE_TO_2BIT: np.ndarray = np.full(256, -1, dtype=np.int8)
_BASE_TO_2BIT[ord("A")] = 0
_BASE_TO_2BIT[ord("C")] = 1
_BASE_TO_2BIT[ord("G")] = 2
_BASE_TO_2BIT[ord("T")] = 3
_BASE_TO_2BIT[ord("a")] = 0
_BASE_TO_2BIT[ord("c")] = 1
_BASE_TO_2BIT[ord("g")] = 2
_BASE_TO_2BIT[ord("t")] = 3

_INT_TO_BASE: dict[int, str] = {0: "A", 1: "C", 2: "G", 3: "T"}

# Fast single-CB encoding via str.translate + int(s, 4).
# str.translate is C-level; int(s, 4) is also C. Together: ~250 ns/call vs
# ~2.5 μs for the numpy path. Used for hot-path single-CB lookups in indexing.
_TRANS_TABLE = str.maketrans({
    "A": "0", "C": "1", "G": "2", "T": "3",
    "a": "0", "c": "1", "g": "2", "t": "3",
})
_VALID_DIGITS = frozenset("0123")
_VALID_BASES = frozenset("ACGTacgt")


# ---------------------------------------------------------------------------
# Single-CB functions
# ---------------------------------------------------------------------------


def encode_cb(cb: str) -> int:
    """Encode a single CB string to a non-negative int via 2-bit packing.

    Uses an ASCII lookup table so the inner loop runs over a numpy int8 array
    rather than calling ``ord()`` character-by-character.

    Args:
        cb: Cell barcode string (expected ATCG only, any length).

    Returns:
        A non-negative integer with 2 bits per base, MSB first; or ``-1`` if
        *cb* is empty or contains any non-ATCG character.
    """
    if not cb:
        return -1
    # Pre-check via frozenset membership (C-level): rejects any char not in
    # ACGT/acgt — including digit chars like '0' which would otherwise pass
    # through translate unchanged and be misinterpreted as a base-4 digit.
    if not _VALID_BASES.issuperset(cb):
        return -1
    # Fast path: str.translate is C-level char remap; int(s, 4) parses base-4.
    return int(cb.translate(_TRANS_TABLE), 4)


def decode_cb(cb_int: int, length: int = 16) -> str:
    """Decode a 2-bit-packed integer back to a CB string.

    This is the inverse of :func:`encode_cb` and is used only when writing
    CB strings to output files (one decode per cell, never per read).

    Args:
        cb_int: Integer produced by :func:`encode_cb`.
        length: Number of bases in the original CB (default 16 for 10x).

    Returns:
        The original CB string, or ``""`` if *cb_int* is negative (sentinel).
    """
    if cb_int < 0:
        return ""
    chars: list[str] = []
    val = int(cb_int)
    for _ in range(length):
        chars.append(_INT_TO_BASE[val & 0x3])
        val >>= 2
    return "".join(reversed(chars))


# ---------------------------------------------------------------------------
# Vectorised batch encoding
# ---------------------------------------------------------------------------


def encode_cb_batch(cb_list: list[str] | np.ndarray) -> np.ndarray:
    """Vectorised encoding of a list of CB strings to an int64 array.

    All CBs **must** have the same length (validated against ``cb_list[0]``).
    The function builds a single contiguous byte buffer for the entire batch,
    applies the lookup table once, then accumulates bit-shifts over the base
    axis — avoiding per-string Python loops.

    Args:
        cb_list: Sequence of CB strings of equal length.

    Returns:
        ``np.ndarray`` of shape ``(len(cb_list),)`` with dtype ``int64``.
        Entries with non-ATCG characters receive the value ``-1``.
    """
    if len(cb_list) == 0:
        return np.array([], dtype=np.int64)

    n = len(cb_list)
    cb_len = len(cb_list[0])

    if cb_len == 0:
        return np.full(n, -1, dtype=np.int64)

    # Build (n × cb_len) byte matrix via a single frombuffer call.
    # Joining all strings and frombuffer-ing is faster than a Python loop.
    joined = "".join(cb_list)
    buf = np.frombuffer(joined.encode("ascii"), dtype=np.uint8).reshape(n, cb_len)

    # Apply lookup table: shape (n, cb_len), dtype int8, -1 for bad chars
    decoded: np.ndarray = _BASE_TO_2BIT[buf]

    # Track which rows have at least one invalid character
    invalid: np.ndarray = (decoded < 0).any(axis=1)

    # Bit-shift accumulate left-to-right across the base axis
    # Cast to int64 first so shifts don't overflow int8
    out = np.zeros(n, dtype=np.int64)
    for i in range(cb_len):
        out = (out << 2) | decoded[:, i].astype(np.int64)

    # Replace invalid rows with sentinel
    out[invalid] = -1
    return out
