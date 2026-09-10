"""E3 integrity invariant: pas_ledger survivors == n_vars(clusters.h5ad).

The provenance ledger records both survivors (dropped_at == "") and drops at
every drop site. The cross-repo contract invariant is that the count of
surviving PAS rows equals the number of variables (PAS) in the final
clusters.h5ad. This test exercises the reusable checker
(:func:`ema.provenance.check_survivor_invariant`) and the on-disk survivor
counter (:func:`ema.provenance.surviving_count_from_tsv`) on a synthetic
ledger. The *live-pipeline* value of ``n_vars`` needs a real run to validate
end-to-end (flagged); this test locks the checker logic itself.
"""
from __future__ import annotations

import pytest

from ema.provenance import (
    InvariantError,
    ProvenanceLedger,
    check_survivor_invariant,
    surviving_count_from_tsv,
)


def _seed_ledger(tmp_path):
    """3 surviving PAS, 2 dropped (atlas_snap + pas_gene). Survivors == 3."""
    led = ProvenanceLedger(tmp_path)
    # survivors
    for i in range(3):
        led.record_pas(orig_pas_key=f"ds::+::{i}", chrom="chr1", start=i * 10,
                       end=i * 10 + 5, strand="+", unified_pas_id=i,
                       last_stage="clusters")
    # drops at two different sites
    led.record_drop("pas", dropped_at="atlas_snap", drop_reason="no atlas PAS within 50bp",
                    orig_pas_key="ds::+::7", chrom="chr1")
    led.record_drop("pas", dropped_at="pas_gene", drop_reason="no gene within max_distance",
                    orig_pas_key="ds::-::9", chrom="chr2", strand="-")
    return led


def test_count_surviving_in_memory(tmp_path):
    led = _seed_ledger(tmp_path)
    assert led.count_surviving("pas") == 3
    assert led.surviving_pas_count == 3


def test_surviving_count_from_tsv_roundtrip(tmp_path):
    led = _seed_ledger(tmp_path)
    led.flush()
    path = tmp_path / "provenance" / "pas_ledger.tsv"
    assert surviving_count_from_tsv(path, kind="pas") == 3


def test_invariant_holds_when_counts_match(tmp_path):
    led = _seed_ledger(tmp_path)
    # clusters.h5ad has exactly 3 vars → invariant holds.
    assert check_survivor_invariant(led.surviving_pas_count, 3) is True


def test_invariant_raises_on_mismatch(tmp_path):
    led = _seed_ledger(tmp_path)
    with pytest.raises(InvariantError) as exc:
        check_survivor_invariant(led.surviving_pas_count, 5)
    # error message reports both counts and the delta.
    msg = str(exc.value)
    assert "surviving=3" in msg and "n_vars=5" in msg and "delta=-2" in msg


def test_invariant_soft_mode_returns_bool(tmp_path):
    led = _seed_ledger(tmp_path)
    assert check_survivor_invariant(led.surviving_pas_count, 4, raise_on_fail=False) is False


def test_empty_ledger_counts_zero(tmp_path):
    led = ProvenanceLedger(tmp_path)
    led.flush()
    path = tmp_path / "provenance" / "pas_ledger.tsv"
    assert surviving_count_from_tsv(path, kind="pas") == 0
    assert check_survivor_invariant(0, 0) is True


def test_tsv_without_dropped_at_raises(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("orig_pas_key\tchrom\nx\tchr1\n")
    with pytest.raises(ValueError):
        surviving_count_from_tsv(bad, kind="pas")
