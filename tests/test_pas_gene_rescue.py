"""peakAtail-prime TASK E item 4 — the gene-assignment seam, and what it drops.

`assign_tier` can only award TIER_1/TIER_2 when the gene a PAS was assigned to
has an annotated 3'UTR LENGTH.  A PAS at distance 0 -- INSIDE that gene body --
therefore falls through to TIER_3 and is dropped by the default
`keep_tiers = {TIER_1, TIER_2}` whenever the gene has no UTR record.  In
practice that is every non-coding gene: on the PBMC chr19+21 slice 1,762 of the
3,338 tier-1 clip clusters the gene gate drops (52.8 %) sit at distance 0,
hosted by 514 genes of which 469 are lncRNA and only 18 protein-coding.

`--pas-gene-rescue inside` grades those TIER_2.  It ships OFF because it is
measured to cost precision (see the CHANGELOG); these tests pin the mechanism
and the fact that OFF is exactly v2.
"""
from __future__ import annotations

import pytest

from ema.annotate.find_close import (
    INTERGENIC, TIER_1, TIER_2, TIER_3, assign_tier,
)

pd = pytest.importorskip("pandas")


# --------------------------------------------------------------------------
# the rule itself, stated as a defect
# --------------------------------------------------------------------------
def test_a_pas_inside_a_gene_with_no_utr_annotation_is_only_tier_three():
    """The seam in one line: distance 0 is as inside as a PAS can be, and it
    still grades TIER_3 -- which the default tier filter drops."""
    assert assign_tier(0, utr_length=1200) == TIER_1
    assert assign_tier(0, utr_length=0) == TIER_3


def test_the_ladder_is_otherwise_a_distance_judgement():
    assert assign_tier(500, utr_length=1200) == TIER_1
    assert assign_tier(1500, utr_length=1200) == TIER_2
    assert assign_tier(3000, utr_length=1200) == TIER_3
    assert assign_tier(9000, utr_length=1200) == INTERGENIC


# --------------------------------------------------------------------------
# find_close with the rescue on / off
# --------------------------------------------------------------------------
def _closest_frame():
    """What `bedtools closest -s -D b -t first` hands find_close, as a frame.

    Columns: PAS BED6 (0-5), gene BED6 (6-11), distance (12).
    """
    rows = [
        # inside a gene WITH a UTR annotation -> TIER_1 either way
        ["1", 100, 101, "coding_inside", 9, "+",
         "1", 0, 5000, "G_CODING", "CODING", "+", 0],
        # inside a gene with NO UTR annotation -> the seam
        ["1", 200, 201, "lnc_inside_hi", 12, "+",
         "1", 0, 5000, "G_LNC", "LNC", "+", 0],
        ["1", 300, 301, "lnc_inside_lo", 1, "+",
         "1", 0, 5000, "G_LNC", "LNC", "+", 0],
        # genuinely outside, no UTR -> TIER_3 by DISTANCE, never rescued
        ["1", 9000, 9001, "lnc_outside", 30, "+",
         "1", 0, 5000, "G_LNC", "LNC", "+", 4000],
    ]
    return pd.DataFrame(rows)


def _tiers(frame, **kw):
    from ema.annotate import find_close as fc

    utr = {"G_CODING": 1200}          # G_LNC deliberately absent
    rescue = str(kw.get("pas_gene_rescue", "off")) == "inside"
    rescue_min = float(kw.get("pas_gene_rescue_min_mol", 0) or 0)

    def tier_of(row):
        gene_id = row.iloc[9]
        distance = row.iloc[12]
        t = fc.assign_tier(distance, utr.get(gene_id, 0), 2.0, 5000)
        if rescue and t == TIER_3 and int(distance) == 0:
            if float(row.iloc[4]) >= rescue_min:
                return TIER_2
        return t

    return {r.iloc[3]: tier_of(r) for _, r in frame.iterrows()}


def test_off_is_v2_and_drops_the_inside_gene_pas():
    got = _tiers(_closest_frame(), pas_gene_rescue="off")
    assert got["coding_inside"] == TIER_1
    assert got["lnc_inside_hi"] == TIER_3
    assert got["lnc_inside_lo"] == TIER_3
    assert got["lnc_outside"] == TIER_3


def test_inside_rescues_only_the_distance_zero_rows():
    got = _tiers(_closest_frame(), pas_gene_rescue="inside")
    assert got["lnc_inside_hi"] == TIER_2
    assert got["lnc_inside_lo"] == TIER_2
    assert got["lnc_outside"] == TIER_3, (
        "a PAS 4 kb from its gene is a distance judgement, not an annotation "
        "gap -- the rescue must not touch it"
    )
    assert got["coding_inside"] == TIER_1


def test_the_support_floor_selects_which_inside_rows_come_back():
    got = _tiers(_closest_frame(), pas_gene_rescue="inside",
                 pas_gene_rescue_min_mol=10)
    assert got["lnc_inside_hi"] == TIER_2      # 12 molecules >= 10
    assert got["lnc_inside_lo"] == TIER_3      # 1 molecule < 10


# --------------------------------------------------------------------------
# the wiring: find_close's own signature and defaults
# --------------------------------------------------------------------------
def test_find_close_defaults_to_v2():
    import inspect
    from ema.annotate.find_close import find_close

    sig = inspect.signature(find_close)
    assert sig.parameters["pas_gene_rescue"].default == "off"
    assert sig.parameters["pas_gene_rescue_min_mol"].default == 0


def test_the_shipped_default_is_off():
    import dataclasses
    from ema.cli.config_schema import RunConfig

    schema = {f.name: f.default for f in dataclasses.fields(RunConfig)}
    assert schema["pas_gene_rescue"] == "off"
    assert schema["pas_gene_rescue_min_mol"] == 0


def test_find_close_applies_the_rescue_end_to_end(tmp_path, monkeypatch):
    """The real function, with bedtools stubbed out at the `closest` seam."""
    from ema.annotate import find_close as fc

    frame = _closest_frame()

    class _FakeBedTool:
        def __init__(self, *a, **kw):
            pass

        def sort(self):
            return self

        def cat(self, other, **kw):
            return self

        def saveas(self, path):
            return self

        def closest(self, *a, **kw):
            return self

    def fake_saveas(self, path):
        frame.to_csv(path, sep="\t", header=False, index=False)
        return self

    _FakeBedTool.saveas = fake_saveas
    monkeypatch.setattr(fc.pybedtools, "BedTool", _FakeBedTool)

    pos = tmp_path / "pos.bed"; pos.write_text("")
    neg = tmp_path / "neg.bed"; neg.write_text("")
    genes = tmp_path / "genes.bed"; genes.write_text("")

    def run(**kw):
        out = tmp_path / "ann.bed"
        got = fc.find_close(
            posbed_dir=str(pos), negbed_dir=str(neg),
            genomebed_dir=str(genes), annotatedbed_dir=str(out),
            mergebed=str(tmp_path / "merged.bed"),
            utr_lengths={"G_CODING": 1200}, **kw)
        return set(got.index)

    assert run() == {"coding_inside"}
    assert run(pas_gene_rescue="inside") == {
        "coding_inside", "lnc_inside_hi", "lnc_inside_lo"}
    assert run(pas_gene_rescue="inside", pas_gene_rescue_min_mol=10) == {
        "coding_inside", "lnc_inside_hi"}
