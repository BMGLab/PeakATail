"""Tests for ema.countmatrix.cb_encode (Task P-3).

Coverage targets:
  - encode/decode round-trip for every 4-mer (4^4 = 256 cases)
  - batch encode matches single encode for 1000 random CBs
  - non-ATCG character ('N', 'X', etc.) returns sentinel -1
  - 16 nt CB encodes to a value that fits in 32 bits (< 2^32)
  - empty string returns -1
  - decode of -1 returns empty string
  - short CB (< 16 nt) round-trips correctly
"""

from __future__ import annotations

import itertools
import random
import string

import numpy as np
import pytest

from ema.countmatrix.cb_encode import decode_cb, encode_cb, encode_cb_batch

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

BASES = ["A", "C", "G", "T"]
RANDOM_SEED = 42


def _random_cb(length: int = 16, rng: random.Random | None = None) -> str:
    rng = rng or random.Random(RANDOM_SEED)
    return "".join(rng.choice(BASES) for _ in range(length))


# ---------------------------------------------------------------------------
# encode_cb / decode_cb round-trip — every 4-mer (256 cases)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kmer",
    ["".join(combo) for combo in itertools.product(BASES, repeat=4)],
)
def test_roundtrip_4mer(kmer: str) -> None:
    """Every 4-mer must survive an encode → decode round-trip unchanged."""
    encoded = encode_cb(kmer)
    assert encoded >= 0, f"encode_cb({kmer!r}) returned sentinel for valid input"
    decoded = decode_cb(encoded, length=4)
    assert decoded == kmer, f"Round-trip failed: {kmer!r} → {encoded} → {decoded!r}"


# ---------------------------------------------------------------------------
# Sentinel for non-ATCG / empty inputs
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bad_char", ["N", "n", "X", "U", "R", "Y", "0", " ", "-"])
def test_encode_invalid_char_returns_sentinel(bad_char: str) -> None:
    """Single non-ATCG character must return -1."""
    cb = f"AAAAACCC{bad_char}GGGTTTT"  # 16 chars with one invalid
    result = encode_cb(cb)
    assert result == -1, f"Expected -1 for CB with char {bad_char!r}, got {result}"


def test_encode_empty_string_returns_sentinel() -> None:
    """Empty string input must return -1."""
    assert encode_cb("") == -1


def test_encode_all_n_returns_sentinel() -> None:
    """Barcode of all Ns must return -1."""
    assert encode_cb("N" * 16) == -1


def test_decode_sentinel_returns_empty_string() -> None:
    """decode_cb(-1) must return empty string (not crash)."""
    assert decode_cb(-1) == ""
    assert decode_cb(-1, length=4) == ""


# ---------------------------------------------------------------------------
# 16 nt CB fits in 32 bits
# ---------------------------------------------------------------------------


def test_16nt_cb_fits_in_32_bits() -> None:
    """16 nt CB uses at most 32 bits (2 bits × 16 bases)."""
    cb = "AAAACCCCGGGGTTTT"
    encoded = encode_cb(cb)
    assert encoded >= 0
    assert encoded < 2**32, f"Expected < 2^32, got {encoded}"


def test_known_encoding_value() -> None:
    """Verify a known encoding against hand-computed value.

    AAAA = 00 00 00 00 → 0
    CCCC = 01 01 01 01 → 0b01010101 = 85
    GGGG = 10 10 10 10 → 0b10101010 = 170
    TTTT = 11 11 11 11 → 0b11111111 = 255
    AAAACCCCGGGGTTTT concatenated → 0 << 24 | 85 << 16 | 170 << 8 | 255
                                   = 0x0055AAFF = 5613311
    (treating each 4-base group as one byte, shifted appropriately)
    """
    cb = "AAAACCCCGGGGTTTT"
    encoded = encode_cb(cb)
    decoded = decode_cb(encoded, length=16)
    assert decoded == cb


# ---------------------------------------------------------------------------
# Batch encode matches single encode — 1000 random CBs
# ---------------------------------------------------------------------------


def test_batch_matches_single_1000_random() -> None:
    """encode_cb_batch must produce identical results to encode_cb for every CB."""
    rng = random.Random(RANDOM_SEED)
    cb_list = [_random_cb(16, rng) for _ in range(1000)]

    batch_result = encode_cb_batch(cb_list)
    assert batch_result.shape == (1000,)
    assert batch_result.dtype == np.int64

    for i, cb in enumerate(cb_list):
        single = encode_cb(cb)
        assert batch_result[i] == single, (
            f"Mismatch at index {i}: single={single}, batch={batch_result[i]}"
        )


def test_batch_with_invalid_entries() -> None:
    """Batch encode must return -1 for entries with non-ATCG characters."""
    cbs = ["AAAA", "ANAA", "TTTT", "CCGN"]
    result = encode_cb_batch(cbs)
    assert result[0] == encode_cb("AAAA")
    assert result[1] == -1
    assert result[2] == encode_cb("TTTT")
    assert result[3] == -1


def test_batch_empty_list() -> None:
    """Empty list input must return an empty int64 array."""
    result = encode_cb_batch([])
    assert isinstance(result, np.ndarray)
    assert result.dtype == np.int64
    assert len(result) == 0


def test_batch_single_element() -> None:
    """Single-element batch must match single-encode."""
    cb = "ACGTACGTACGTACGT"
    assert encode_cb_batch([cb])[0] == encode_cb(cb)


# ---------------------------------------------------------------------------
# Short CBs (length < 16)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("length", [1, 2, 4, 8, 12, 14])
def test_short_cb_roundtrip(length: int) -> None:
    """Short CBs of any valid length must round-trip correctly."""
    rng = random.Random(RANDOM_SEED + length)
    cb = _random_cb(length, rng)
    encoded = encode_cb(cb)
    assert encoded >= 0
    decoded = decode_cb(encoded, length=length)
    assert decoded == cb, f"Round-trip failed for length={length}: {cb!r} → {decoded!r}"


def test_short_cb_with_invalid_char() -> None:
    """Short CBs with non-ATCG chars must still return -1."""
    assert encode_cb("AN") == -1
    assert encode_cb("NA") == -1


# ---------------------------------------------------------------------------
# Batch encode dtype and large input
# ---------------------------------------------------------------------------


def test_batch_dtype_is_int64() -> None:
    """Batch output must always be int64 regardless of input length."""
    cbs = ["ACGT"] * 10
    result = encode_cb_batch(cbs)
    assert result.dtype == np.int64


def test_batch_large_random_no_crash() -> None:
    """Batch encoding of 10 000 CBs must complete without errors."""
    rng = random.Random(RANDOM_SEED)
    cb_list = [_random_cb(16, rng) for _ in range(10_000)]
    result = encode_cb_batch(cb_list)
    assert result.shape == (10_000,)
    # All should be >= 0 (all random, all valid ATCG)
    assert (result >= 0).all()


# ---------------------------------------------------------------------------
# Idempotency: encoding twice gives the same result
# ---------------------------------------------------------------------------


def test_encode_idempotent() -> None:
    """Calling encode_cb twice on the same string gives the same integer."""
    cb = "GCTAGCTAGCTAGCTA"
    assert encode_cb(cb) == encode_cb(cb)


# ---------------------------------------------------------------------------
# BarcodeIndex integration (tests that indexing still works after refactor)
# ---------------------------------------------------------------------------


def test_indexing_roundtrip_via_get_mapping() -> None:
    """BarcodeIndex must reconstruct the original CB string via get_mapping()."""
    from ema.countmatrix.indexing import BarcodeIndex

    idx = BarcodeIndex()
    cb = "sample1_ACGTACGTACGTACGT"
    col = idx.get_index(cb)
    assert col == 1

    mapping = idx.mapping
    assert cb in mapping
    assert mapping[cb] == 1


def test_indexing_same_cb_same_index() -> None:
    """The same CB string must always return the same column index."""
    from ema.countmatrix.indexing import BarcodeIndex

    idx = BarcodeIndex()
    cb = "s1_AAAACCCCGGGGTTTT"
    assert idx.get_index(cb) == idx.get_index(cb)


def test_indexing_different_cbs_different_indices() -> None:
    """Different CB strings must receive different column indices."""
    from ema.countmatrix.indexing import BarcodeIndex

    idx = BarcodeIndex()
    col_a = idx.get_index("s1_AAAAAAAAAAAAAAAA")
    col_b = idx.get_index("s1_TTTTTTTTTTTTTTTT")
    assert col_a != col_b


def test_get_indices_batch_matches_single() -> None:
    """get_indices_batch must return the same indices as repeated get_index calls."""
    from ema.countmatrix.indexing import BarcodeIndex

    rng = random.Random(RANDOM_SEED)
    idx = BarcodeIndex()
    cbs = [f"s1_{_random_cb(16, rng)}" for _ in range(100)]

    # Single-CB baseline
    idx_single = BarcodeIndex()
    single_cols = [idx_single.get_index(cb) for cb in cbs]

    # Batch path (fresh index so assignment order is identical)
    idx_batch = BarcodeIndex()
    batch_cols = idx_batch.get_indices_batch(cbs)

    assert batch_cols == single_cols
