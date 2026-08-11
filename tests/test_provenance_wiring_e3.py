"""E3 drop-site wiring + integrity invariant on a SYNTHETIC mini-run.

Exercises the per-dataset drop accounting the downstream worker performs at the
pas_gene (annotate) + preprocess (min_cells) PAS drop sites and the cb_filter
(+ preprocess) cell drop sites, and asserts the headline provenance invariant:

    rows(pas_ledger where dropped_at == "") == n_vars(clusters.h5ad)

No real cohort run is needed — a tiny constructed PAS/cell funnel drives the
same helpers the worker calls, so the invariant is proven deterministically.
"""
from __future__ import annotations

import pytest

from ema.provenance import (
    ProvenanceLedger,
    check_survivor_invariant,
    record_cell_drops,
    record_pas_drops,
    surviving_count_from_tsv,
)


# A tiny synthetic dataset funnel:
#   input PAS      : p0..p5   (6)
#   annotate keeps : p0..p4   (p5 has no gene → pas_gene drop)
#   preprocess keeps (min_cells): p0..p2  (p3,p4 below min_cells → preprocess drop)
#   → clusters.h5ad n_vars = 3
INPUT_PAS = ["p0", "p1", "p2", "p3", "p4", "p5"]
ANNOTATED_PAS = ["p0", "p1", "p2", "p3", "p4"]
FINAL_PAS = ["p0", "p1", "p2"]
GENE_OF = {"p0": "gA", "p1": "gA", "p2": "gB", "p3": "gB", "p4": "gC"}


def test_pas_accounting_conserves_and_invariant_holds(tmp_path):
    led = ProvenanceLedger(tmp_path)
    surviving = record_pas_drops(
        led,
        input_pas_ids=INPUT_PAS,
        annotated_pas_ids=ANNOTATED_PAS,
        final_pas_ids=FINAL_PAS,
        dataset_id="dsX",
        gene_of=GENE_OF,
    )
    # survivors == final == would-be n_vars(clusters.h5ad)
    assert surviving == len(FINAL_PAS) == 3
    n_vars = 3
    assert check_survivor_invariant(surviving, n_vars) is True

    # Conservation: input == pas_gene_drops + preprocess_drops + survivors.
    pas_rows = led._rows["pas"]
    dropped_at_idx = 11  # PAS_LEDGER_COLUMNS.index("dropped_at")
    pas_gene = sum(1 for r in pas_rows if r[dropped_at_idx] == "pas_gene")
    preprocess = sum(1 for r in pas_rows if r[dropped_at_idx] == "preprocess")
    survivors = sum(1 for r in pas_rows if r[dropped_at_idx] == "")
    assert pas_gene == 1 and preprocess == 2 and survivors == 3
    assert pas_gene + preprocess + survivors == len(INPUT_PAS)


def test_pas_survivors_persist_and_reload(tmp_path):
    led = ProvenanceLedger(tmp_path)
    record_pas_drops(led, input_pas_ids=INPUT_PAS, annotated_pas_ids=ANNOTATED_PAS,
                     final_pas_ids=FINAL_PAS, gene_of=GENE_OF)
    led.flush()
    path = tmp_path / "provenance" / "pas_ledger.tsv"
    assert surviving_count_from_tsv(path, kind="pas") == 3


def test_cell_accounting_cb_filter_and_preprocess(tmp_path):
    led = ProvenanceLedger(tmp_path)
    #   input cbs : c0..c4 (5)
    #   cb_filter keeps : c0..c3 (c4 below min_read)
    #   preprocess keeps: c0..c2 (c3 below min_genes)
    input_cbs = ["c0", "c1", "c2", "c3", "c4"]
    kept = ["c0", "c1", "c2", "c3"]
    final = ["c0", "c1", "c2"]
    surviving = record_cell_drops(
        led, input_cbs=input_cbs, kept_cbs=kept, final_cbs=final, dataset_id="dsX",
    )
    assert surviving == 3
    cell_rows = led._rows["cell"]
    dropped_at_idx = 4  # CELL_LEDGER_COLUMNS.index("dropped_at")
    cb_drop = sum(1 for r in cell_rows if r[dropped_at_idx] == "cb_filter")
    pp_drop = sum(1 for r in cell_rows if r[dropped_at_idx] == "preprocess")
    surv = sum(1 for r in cell_rows if r[dropped_at_idx] == "")
    assert cb_drop == 1 and pp_drop == 1 and surv == 3
    assert cb_drop + pp_drop + surv == len(input_cbs)


def test_cell_accounting_without_preprocess_stage(tmp_path):
    led = ProvenanceLedger(tmp_path)
    surviving = record_cell_drops(
        led, input_cbs=["c0", "c1", "c2"], kept_cbs=["c0", "c1"], dataset_id="dsX",
    )
    assert surviving == 2  # no final_cbs → survivors = kept


def test_invariant_violation_is_detectable(tmp_path):
    # Simulate a BUG: a drop site not accounted for (final has a PAS never
    # recorded as surviving) → surviving != n_vars → invariant catches it.
    led = ProvenanceLedger(tmp_path)
    surviving = record_pas_drops(
        led, input_pas_ids=INPUT_PAS, annotated_pas_ids=ANNOTATED_PAS,
        final_pas_ids=FINAL_PAS, gene_of=GENE_OF,
    )
    # clusters.h5ad claims 4 vars but only 3 survived accounting → mismatch.
    assert check_survivor_invariant(surviving, 4, raise_on_fail=False) is False


def test_no_drops_all_survive(tmp_path):
    # Degenerate funnel: nothing dropped → survivors == input == n_vars.
    led = ProvenanceLedger(tmp_path)
    surviving = record_pas_drops(
        led, input_pas_ids=["p0", "p1"], annotated_pas_ids=["p0", "p1"],
        final_pas_ids=["p0", "p1"],
    )
    assert surviving == 2
    assert check_survivor_invariant(surviving, 2) is True
    # survivors carry gene_id when provided, else empty.
    assert sum(1 for r in led._rows["pas"] if r[11] == "") == 2
