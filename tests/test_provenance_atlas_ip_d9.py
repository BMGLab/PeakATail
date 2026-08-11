"""D9: atlas_match / atlas_distance_bp / internal_priming ledger columns.

Covers the provenance-layer half of "atlas-snap and internal-priming become
ANNOTATE (never drop by default)":

  * The three new ``PAS_LEDGER_COLUMNS`` are additive at the end and default
    to ``""`` -- a run with neither atlas nor ip enabled produces the SAME
    ledger values as before, plus the three new empty columns.
  * ``record_pas_drops(..., atlas_of=..., ip_of=...)`` populates them for
    every row it writes (drops AND survivors), keyed on the pas_id.
"""
from __future__ import annotations

from ema.provenance import PAS_LEDGER_COLUMNS, ProvenanceLedger, record_pas_drops

INPUT_PAS = ["p0", "p1", "p2", "p3", "p4", "p5"]
ANNOTATED_PAS = ["p0", "p1", "p2", "p3", "p4"]
FINAL_PAS = ["p0", "p1", "p2"]
GENE_OF = {"p0": "gA", "p1": "gA", "p2": "gB", "p3": "gB", "p4": "gC"}


def test_new_columns_are_last_and_additive():
    assert PAS_LEDGER_COLUMNS[-3:] == ["atlas_match", "atlas_distance_bp", "internal_priming"]
    # Everything before them is byte-identical to the pre-D9 contract.
    assert PAS_LEDGER_COLUMNS[:-3] == [
        "orig_pas_key", "chrom", "start", "end", "strand", "unified_pas_id",
        "snap_distance_bp", "gene_id", "gene_distance_bp", "tier", "last_stage",
        "dropped_at", "drop_reason",
    ]


def test_no_atlas_no_ip_produces_empty_new_columns(tmp_path):
    """Default-off: neither atlas_of nor ip_of passed -> every row's three
    new columns are "" -- identical, aside from those additive empty
    columns, to the pre-D9 ledger."""
    led = ProvenanceLedger(tmp_path)
    record_pas_drops(
        led, input_pas_ids=INPUT_PAS, annotated_pas_ids=ANNOTATED_PAS,
        final_pas_ids=FINAL_PAS, gene_of=GENE_OF,
    )
    led.flush()
    rows = tmp_path.joinpath("provenance", "pas_ledger.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    assert header == PAS_LEDGER_COLUMNS
    for line in rows[1:]:
        vals = dict(zip(header, line.split("\t")))
        assert vals["atlas_match"] == ""
        assert vals["atlas_distance_bp"] == ""
        assert vals["internal_priming"] == ""


def test_atlas_of_and_ip_of_populate_every_row(tmp_path):
    atlas_of = {
        "p0": (True, 5), "p1": (True, 12), "p2": (False, 200),
        "p3": (False, ""), "p4": (True, 3), "p5": (False, ""),
    }
    ip_of = {"p0": False, "p1": True, "p2": False, "p3": True, "p4": False, "p5": True}

    led = ProvenanceLedger(tmp_path)
    surviving = record_pas_drops(
        led, input_pas_ids=INPUT_PAS, annotated_pas_ids=ANNOTATED_PAS,
        final_pas_ids=FINAL_PAS, gene_of=GENE_OF, atlas_of=atlas_of, ip_of=ip_of,
    )
    assert surviving == 3

    led.flush()
    rows = tmp_path.joinpath("provenance", "pas_ledger.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    by_key = {}
    for line in rows[1:]:
        vals = dict(zip(header, line.split("\t")))
        by_key[vals["orig_pas_key"]] = vals

    # p5 dropped at pas_gene -- still gets its atlas/ip status.
    assert by_key["p5"]["dropped_at"] == "pas_gene"
    assert by_key["p5"]["atlas_match"] == "False"
    assert by_key["p5"]["internal_priming"] == "True"

    # p3, p4 dropped at preprocess -- still get their status.
    assert by_key["p3"]["dropped_at"] == "preprocess"
    assert by_key["p3"]["atlas_match"] == "False"
    assert by_key["p3"]["atlas_distance_bp"] == ""
    assert by_key["p4"]["atlas_match"] == "True"
    assert by_key["p4"]["atlas_distance_bp"] == "3"

    # p0, p1, p2 survive -- status carried through too.
    assert by_key["p0"]["dropped_at"] == ""
    assert by_key["p0"]["atlas_match"] == "True"
    assert by_key["p0"]["atlas_distance_bp"] == "5"
    assert by_key["p0"]["internal_priming"] == "False"
    assert by_key["p1"]["internal_priming"] == "True"
    assert by_key["p2"]["atlas_match"] == "False"


def test_atlas_of_ip_of_default_to_empty_for_unmapped_pas(tmp_path):
    """A pas_id absent from atlas_of/ip_of (e.g. atlas/ip didn't cover it)
    still gets "" defaults, never a KeyError."""
    led = ProvenanceLedger(tmp_path)
    record_pas_drops(
        led, input_pas_ids=["p0", "p1"], annotated_pas_ids=["p0", "p1"],
        final_pas_ids=["p0", "p1"], atlas_of={"p0": (True, 1)}, ip_of={"p0": True},
    )
    led.flush()
    rows = tmp_path.joinpath("provenance", "pas_ledger.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    by_key = {dict(zip(header, line.split("\t")))["orig_pas_key"]: dict(zip(header, line.split("\t")))
              for line in rows[1:]}
    assert by_key["p1"]["atlas_match"] == ""
    assert by_key["p1"]["atlas_distance_bp"] == ""
    assert by_key["p1"]["internal_priming"] == ""
