"""E3 sidecar-then-reconcile on a SYNTHETIC 2-dataset mini-run.

Each per-dataset worker writes its own ledger sidecar; the parent reconciles
after the pool joins and self-checks the per-dataset invariant
``surviving PAS == n_vars(clusters.h5ad)``. This exercises the reconcile on a
constructed two-dataset fixture — NO real cohort run.
"""
from __future__ import annotations

from ema.provenance import (
    ProvenanceLedger,
    concat_ledger_tsvs,
    reconcile_dataset_ledgers,
    record_cell_drops,
    record_pas_drops,
    surviving_count_from_tsv,
)


def _write_dataset_sidecar(root, ds_id, input_pas, annotated_pas, final_pas,
                           input_cbs, kept_cbs, final_cbs):
    """Emulate one worker: build its per-dataset ledger + flush to a sidecar."""
    sidecar_parent = root / "07_clustering" / ds_id
    led = ProvenanceLedger(sidecar_parent)
    record_pas_drops(led, input_pas_ids=input_pas, annotated_pas_ids=annotated_pas,
                     final_pas_ids=final_pas, dataset_id=ds_id)
    record_cell_drops(led, input_cbs=input_cbs, kept_cbs=kept_cbs,
                      final_cbs=final_cbs, dataset_id=ds_id)
    led.flush()
    return sidecar_parent / "provenance", len(final_pas)


def _two_dataset_run(root):
    # dsA: 6 input PAS → 5 annotated → 3 final ; dsB: 4 → 4 → 2 final.
    a_dir, a_nvars = _write_dataset_sidecar(
        root, "dsA",
        input_pas=[f"a{i}" for i in range(6)],
        annotated_pas=[f"a{i}" for i in range(5)],
        final_pas=["a0", "a1", "a2"],
        input_cbs=[f"ca{i}" for i in range(5)], kept_cbs=[f"ca{i}" for i in range(4)],
        final_cbs=["ca0", "ca1", "ca2"],
    )
    b_dir, b_nvars = _write_dataset_sidecar(
        root, "dsB",
        input_pas=[f"b{i}" for i in range(4)],
        annotated_pas=[f"b{i}" for i in range(4)],
        final_pas=["b0", "b1"],
        input_cbs=[f"cb{i}" for i in range(3)], kept_cbs=[f"cb{i}" for i in range(3)],
        final_cbs=["cb0", "cb1"],
    )
    return [("dsA", a_dir, a_nvars), ("dsB", b_dir, b_nvars)]


def test_reconcile_invariant_holds_per_dataset(tmp_path):
    datasets = _two_dataset_run(tmp_path)
    cohort = tmp_path / "provenance" / "by_dataset"
    summary = reconcile_dataset_ledgers(datasets, cohort_dir=cohort)

    assert summary["n_datasets"] == 2
    assert summary["all_ok"] is True
    by_id = {d["ds_id"]: d for d in summary["datasets"]}
    assert by_id["dsA"]["surviving_pas"] == 3 == by_id["dsA"]["n_vars"]
    assert by_id["dsB"]["surviving_pas"] == 2 == by_id["dsB"]["n_vars"]


def test_reconcile_cohort_concatenation(tmp_path):
    datasets = _two_dataset_run(tmp_path)
    cohort = tmp_path / "provenance" / "by_dataset"
    reconcile_dataset_ledgers(datasets, cohort_dir=cohort)

    # cohort pas ledger = union of both datasets' survivors (3 + 2 = 5).
    cohort_pas = cohort / "pas_ledger.tsv"
    assert cohort_pas.exists()
    assert surviving_count_from_tsv(cohort_pas, kind="pas") == 5
    # cohort cell ledger survivors = 3 + 2 = 5.
    assert surviving_count_from_tsv(cohort / "cell_ledger.tsv", kind="cell") == 5


def test_reconcile_detects_a_bad_dataset(tmp_path):
    datasets = _two_dataset_run(tmp_path)
    # Corrupt dsB's declared n_vars so it no longer matches its 2 survivors.
    datasets[1] = ("dsB", datasets[1][1], 99)
    summary = reconcile_dataset_ledgers(datasets)
    assert summary["all_ok"] is False
    bad = [d for d in summary["datasets"] if not d["invariant_ok"]]
    assert len(bad) == 1 and bad[0]["ds_id"] == "dsB"
    assert bad[0]["surviving_pas"] == 2 and bad[0]["n_vars"] == 99


def test_reconcile_missing_sidecar_counts_zero(tmp_path):
    # A worker that never wrote a sidecar (crashed) → 0 survivors → flagged
    # unless its n_vars is also 0.
    summary = reconcile_dataset_ledgers([("dsGhost", tmp_path / "nope", 5)])
    assert summary["datasets"][0]["surviving_pas"] == 0
    assert summary["all_ok"] is False


def test_concat_skips_missing_and_empty(tmp_path):
    good = tmp_path / "a" / "provenance"
    good.mkdir(parents=True)
    led = ProvenanceLedger(tmp_path / "a")
    led.record_pas(orig_pas_key="x", dropped_at="")
    led.flush()
    out = tmp_path / "merged.tsv"
    n = concat_ledger_tsvs(
        [good / "pas_ledger.tsv", tmp_path / "missing.tsv"], out, kind="pas")
    assert n == 1 and out.exists()
