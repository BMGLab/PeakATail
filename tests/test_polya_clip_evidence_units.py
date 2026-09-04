"""Stage 1c: the unit of poly(A) clip evidence (molecules, not reads).

The Stage-2 sweep verifier's finding, in one sentence: BED column 5 carried
RAW CLIP READS — PCR duplicates and secondary alignments included — while
``--polya-min-reads`` documented (and the gate actually used) distinct
molecules, so 23.7% of the "``>=2`` read" sites were one molecule counted
twice.  These tests pin the corrected contract:

* column 5 is distinct ``(CB, UMI)`` MOLECULES, the same unit
  ``--polya-min-umis`` gates on;
* the raw clip-read count and the ``samtools -F 3844`` counts survive in the
  ``pas_support.tsv`` sidecar, so BED6 stays BED6;
* ``clip_read_ok`` is exactly ``-F 3844``, and ``--polya-clip-filter f3844``
  makes those reads stop being evidence at all;
* ``--polya-min-reads`` still works as a deprecated alias;
* the tier-1 clip fallback no longer counts a read the neighbouring
  cluster's window already counted (measured at 13:43135248+, pas#23403).
"""
from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import pysam
import pytest

from ema.countmatrix.paswrite import SUPPORT_COLUMNS, support_path_for
from ema.countmatrix.polya import (
    CLIP_EXCLUDE_FLAGS,
    ClipAccumulator,
    ClipSeeder,
    clip_read_ok,
    cluster_clip_sites,
    site_molecules,
)


def _flush(seeder: ClipSeeder):
    """``(records, support_rows)`` — the sidecar rows stay aligned."""
    sup: list[dict] = []
    recs = seeder.flush(support_out=sup)
    assert len(recs) == len(sup)
    return recs, sup


# ---------------------------------------------------------------------------
# 1. UMI de-duplication
# ---------------------------------------------------------------------------

def test_three_reads_of_one_molecule_score_one_and_report_three_reads():
    s = ClipSeeder(direction=False)
    for _ in range(3):                     # one PCR stack: same CB, same UMI
        s.add_clip(1000, cb="CB1", umi="UMI1", end1=1091)
    (rec,), (sup,) = _flush(s)
    assert rec[2] == 1, "BED column 5 must be MOLECULES, not reads"
    assert sup["clip_umis"] == 1
    assert sup["clip_reads"] == 3, "the raw read count stays available"
    assert sup["tier"] == 1


def test_same_umi_in_two_cells_is_two_molecules():
    s = ClipSeeder(direction=False)
    s.add_clip(1000, cb="CB1", umi="SAME", end1=1091)
    s.add_clip(1000, cb="CB2", umi="SAME", end1=1091)
    (rec,), _ = _flush(s)
    assert rec[2] == 2, "molecules are keyed by (CB, UMI), not UMI alone"


def test_reads_without_a_umi_tag_count_as_one_molecule_each():
    s = ClipSeeder(direction=False)
    for _ in range(4):
        s.add_clip(1000, cb="CB1", umi=None, end1=1091)
    (rec,), (sup,) = _flush(s)
    assert rec[2] == 4 and sup["clip_reads"] == 4


def test_min_umis_gate_rejects_a_duplicate_stack():
    s = ClipSeeder(direction=False, min_umis=2)
    for _ in range(10):
        s.add_clip(1000, cb="CB1", umi="SAME", end1=1091)
    assert s.flush() == [], "10 PCR duplicates are not 2 molecules"


def test_molecules_are_keyed_by_the_barcode_not_the_read_group():
    """``read_check`` prefixes every barcode with its RG, so a 10x BAM with
    one read group per LANE would report one molecule once per lane
    (+4.7% on the PBMC chr21 (+) slice).  Molecule identity drops the
    prefix; matrix cell identity does not."""
    from ema.countmatrix.polya import molecule_cb

    assert molecule_cb("pbmc:0:1:HFWFVDMXX:1_ACGTACGTACGTACGT") == \
        "ACGTACGTACGTACGT"
    assert molecule_cb("ACGTACGTACGTACGT") == "ACGTACGTACGTACGT"

    s = ClipSeeder(direction=False)
    s.add_clip(1000, cb="lane1_ACGTACGTACGTACGT", umi="u", end1=1091)
    s.add_clip(1000, cb="lane2_ACGTACGTACGTACGT", umi="u", end1=1091)
    (rec,), (sup,) = _flush(s)
    assert rec[2] == 1, "one molecule seen in two lanes is one molecule"
    assert sup["clip_reads"] == 2
    assert rec[3] == {"lane1_ACGTACGTACGTACGT": 1,
                      "lane2_ACGTACGTACGTACGT": 1}, (
        "matrix cell identity must keep the read-group prefix"
    )


# ---------------------------------------------------------------------------
# 2. Alignment filter (samtools -F 3844)
# ---------------------------------------------------------------------------

class _FlagRead:
    def __init__(self, flag: int) -> None:
        self.flag = flag


@pytest.mark.parametrize("flag,ok", [
    (0, True),          # primary forward
    (16, True),         # reverse strand is not an exclusion
    (99, True),         # paired/proper-pair bits are not exclusions
    (4, False),         # unmapped
    (256, False),       # secondary
    (512, False),       # qcfail
    (1024, False),      # duplicate
    (2048, False),      # supplementary
    (1040, False),      # reverse + duplicate
])
def test_clip_read_ok_is_exactly_minus_F_3844(flag, ok):
    assert CLIP_EXCLUDE_FLAGS == 3844 == (4 | 256 | 512 | 1024 | 2048)
    assert clip_read_ok(_FlagRead(flag)) is ok
    # equivalent to what `samtools view -F 3844` keeps
    assert clip_read_ok(_FlagRead(flag)) == (not flag & 3844)


def test_flagged_reads_are_excluded_from_the_f3844_columns_only():
    s = ClipSeeder(direction=False)
    s.add_clip(1000, cb="CB1", umi="u1", end1=1091, primary=True)
    s.add_clip(1000, cb="CB2", umi="u2", end1=1091, primary=False)  # secondary
    s.add_clip(1000, cb="CB3", umi="u3", end1=1091, primary=False)  # duplicate
    (rec,), (sup,) = _flush(s)
    assert rec[2] == 3, "the default unit counts every accepted clip read"
    assert (sup["clip_reads"], sup["clip_umis"]) == (3, 3)
    assert (sup["clip_reads_f3844"], sup["clip_umis_f3844"]) == (1, 1)


def test_clip_filter_f3844_scores_and_gates_on_the_filtered_molecules():
    s = ClipSeeder(direction=False, clip_filter="f3844")
    s.add_clip(1000, cb="CB1", umi="u1", end1=1091, primary=True)
    s.add_clip(1000, cb="CB2", umi="u2", end1=1091, primary=False)
    (rec,), (sup,) = _flush(s)
    assert rec[2] == 1, "under f3844 the score is the filtered molecule count"
    assert sup["clip_reads"] == 2, "the sidecar still reports both units"


def test_clip_filter_f3844_drops_a_cluster_with_only_flagged_evidence():
    """The measured cost of the strict filter: on PBMC 8.3% (chr19+) to
    27.0% (chr21+) of tier-1 clusters have no -F 3844 clip read at all, so
    the strict unit is a CALL-SET change and cannot be the default without
    re-running the benchmark."""
    kw = dict(direction=False, window=100)
    default = ClipSeeder(**kw)
    strict = ClipSeeder(**kw, clip_filter="f3844")
    for s in (default, strict):
        s.add_clip(1000, cb="CB1", umi="u1", end1=1091, primary=False)
        s.add_coverage_pas(5000, 5200, {"CB1": 4})
    assert [r[:2] for r in default.flush()] == [(1000, 1001), (5000, 5200)]
    assert [r[:2] for r in strict.flush()] == [(5000, 5200)], (
        "a cluster whose only evidence is secondary/duplicate alignments "
        "does not survive --polya-clip-filter f3844"
    )


def test_accumulator_tracks_both_units_per_site():
    acc = ClipAccumulator()
    acc.add(500, cb="A", umi="u1", end1=591, primary=True)
    acc.add(500, cb="A", umi="u1", end1=591, primary=False)   # dup of the same
    acc.add(500, cb="B", umi="u2", end1=592, primary=False)
    rec = acc.sites[500]
    assert rec[0] == 3
    assert site_molecules(rec) == 2
    assert site_molecules(rec, strict=True) == 1
    assert rec[4] == {591: {"A": 2}, 592: {"B": 1}}
    (clu,) = cluster_clip_sites(acc.sites, 25)
    assert (clu["nreads"], clu["numis"]) == (3, 2)
    assert (clu["nreads_f3844"], clu["numis_f3844"]) == (1, 1)


# ---------------------------------------------------------------------------
# 3. Flag rename / deprecated alias
# ---------------------------------------------------------------------------

def test_seeder_accepts_the_deprecated_min_reads_kwarg():
    s = ClipSeeder(direction=False, min_reads=3)
    assert s.min_umis == 3 and s.min_reads == 3


def test_cli_accepts_both_spellings_and_yaml_alias():
    import click
    from click.testing import CliRunner

    from ema.cli.config_schema import RunConfig, click_options_from_schema

    seen: dict = {}

    @click.command()
    @click_options_from_schema(RunConfig)
    def _cmd(**kwargs):
        seen.update(kwargs)

    runner = CliRunner()
    assert runner.invoke(_cmd, ["--polya-min-umis", "4"]).exit_code == 0
    assert seen["polya_min_umis"] == 4
    assert runner.invoke(_cmd, ["--polya-min-reads", "7"]).exit_code == 0, (
        "--polya-min-reads must keep working as a deprecated alias"
    )
    assert seen["polya_min_umis"] == 7
    assert RunConfig.from_yaml_dict({"polya_min_reads": 5}).polya_min_umis == 5
    assert RunConfig.from_yaml_dict({"polya_min_umis": 6}).polya_min_umis == 6
    assert RunConfig().polya_min_umis == 1, "default must stay 1"
    assert RunConfig().polya_clip_filter == "none", "default must stay 'none'"


def test_peak_calling_maps_the_deprecated_kwarg_to_molecules(monkeypatch):
    """``peak_calling(polya_min_reads=N)`` must land on the molecule gate."""
    import ema.countmatrix.peackcalling as pc

    captured: dict = {}
    real = pc.ClipSeeder

    def _spy(*a, **kw):
        captured.update(kw)
        return real(*a, **kw)

    monkeypatch.setattr(pc, "ClipSeeder", _spy)
    with pytest.raises(Exception):
        # no BAM: we only need the setup block to run
        pc.peak_calling(False, "/dev/null", "/dev/null",
                        bamfile_dir="/nonexistent.bam",
                        strategy=__import__(
                            "ema.strategies", fromlist=["get_strategy"]
                        ).get_strategy("clip_seeded"),
                        polya_min_reads=9)
    assert captured.get("min_umis") == 9


# ---------------------------------------------------------------------------
# 4. Clip fallback must not count a neighbour's reads (regression)
# ---------------------------------------------------------------------------

def test_clip_fallback_excludes_reads_already_counted_by_a_neighbour():
    """Regression for 13:43135248+ (pas#23403) on testis mouse1.

    That cluster's partition slice and midpoint-clipped window were both
    empty, so it fell back to its 37 clip reads — whose ``end1`` (43135267
    -43135294) lie past the midpoint to the 43135276 cluster and were
    therefore ALSO counted in that cluster's window.  The row summed to 371
    where the shipped candidate + remainder was 334.
    """
    SEQ_LEN = 98
    A, B = 43_135_248, 43_135_276          # cluster anchors ('+': mode)
    s = ClipSeeder(direction=False, window=100, count_window=(-1, 25))
    ends = list(range(43_135_267, 43_135_295))   # 28 reads, all past midpoint
    for i, end1 in enumerate(ends):
        s.add_read(end1 - SEQ_LEN, end1, f"CB{i}")
        s.add_clip(A, cb=f"CB{i}", umi=f"u{i}", end1=end1)   # cluster A's clips
    s.add_clip(B, cb="CBB", umi="ub", end1=43_135_290)       # cluster B exists
    s.add_read(43_135_290 - SEQ_LEN, 43_135_290, "CBB")

    recs, sup = _flush(s)
    by_start = {r[0]: r for r in recs}
    assert set(by_start) == {A, B}
    midpoint = (A + B) // 2
    assert all(e > midpoint for e in ends), "fixture must reproduce the case"

    total = sum(sum(r[3].values()) for r in recs)
    assert total == len(ends) + 1, (
        f"every read end must be counted exactly once, got {total}"
    )
    assert by_start[A][3] == {}, (
        "cluster A's clip reads belong to B's window — counting them again "
        "is the double count this test pins"
    )
    assert sum(by_start[B][3].values()) == len(ends) + 1
    assert by_start[A][2] == len(ends), "the BED score still reports A's molecules"
    assert s.stats["empty_tier1_rows"] == 1
    assert s.stats["clip_fallback_rows"] == 0


def test_clip_fallback_keeps_reads_inside_its_own_territory():
    """The fallback is not disabled — only midpoint-clipped: a clip read in
    the cluster's own territory (here, past its count window on the 3' side
    but well before the midpoint) is still its own."""
    SEQ_LEN = 98
    A, B = 10_000, 20_000
    s = ClipSeeder(direction=False, window=100, count_window=(SEQ_LEN, 5))
    # '+' count window is [site - SEQ_LEN, site + 5]; end1 10_200 is past it
    # but far short of the midpoint to B (15_000), so it is A's own read.
    s.add_clip(A, cb="CBA", umi="ua", end1=10_200)
    s.add_read(10_200 - SEQ_LEN, 10_200, "CBA")
    # B's read end sits inside B's window, so B never reaches the fallback.
    s.add_clip(B, cb="CBB", umi="ub", end1=19_950)
    s.add_read(19_950 - SEQ_LEN, 19_950, "CBB")
    recs, _ = _flush(s)
    by_start = {r[0]: r for r in recs}
    assert by_start[A][3] == {"CBA": 1}, "A's own clip read is not discarded"
    assert by_start[B][3] == {"CBB": 1}
    assert s.stats["clip_fallback_rows"] == 1
    assert s.stats["reads_from_window"] == 1


def test_clip_fallback_ignores_reads_inside_a_coverage_candidate():
    """A clip read whose end1 sits in a coverage candidate is already in
    that candidate's partition."""
    SEQ_LEN = 98
    s = ClipSeeder(direction=False, window=10, count_window=(0, 0))
    s.add_clip(5_000, cb="CBA", umi="ua", end1=7_000)
    s.add_read(7_000 - SEQ_LEN, 7_000, "CBA")
    s.add_coverage_pas(6_900, 7_100, {"CBA": 1})   # far from the cluster
    recs, _ = _flush(s)
    by_start = {r[0]: r for r in recs}
    assert by_start[5_000][3] == {}, "already counted in the candidate"
    assert by_start[6_900][3] == {"CBA": 1}


# ---------------------------------------------------------------------------
# 5. The sidecar, end to end on a BAM
# ---------------------------------------------------------------------------

CHROM = "1"
SEQ_LEN = 91
SITE = 2_400


def _mk_bam(path: Path) -> Path:
    """One clipped locus: 40 summit reads, of which 10 are one PCR stack
    (same CB+UMI) and 6 are secondary alignments."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": 50_000}],
              "RG": [{"ID": "s", "SM": "s"}]}
    recs = []
    n = 0
    specs = [(SITE - 350 + i * 4, 80, 0, 0, i) for i in range(60)]
    for i in range(40):
        start = SITE - 80 + (i % 10)
        flag = 256 if i >= 34 else 0                 # 6 secondary alignments
        dup = i < 10                                  # 10 reads, one molecule
        specs.append((start, SITE - start + 1, 12, flag, -1 if dup else i))
    for start, aligned, clip, flag, cb_i in specs:
        a = pysam.AlignedSegment()
        a.query_name = f"r{n}"
        n += 1
        a.query_sequence = "C" * aligned + "A" * clip
        a.flag = flag
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = 60
        a.cigartuples = [(0, aligned)] + ([(4, clip)] if clip else [])
        a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
        a.set_tag("CB", "ACGTACGTACGTACGT" if cb_i < 0
                  else f"{'ACGT' * 4}"[:12] + f"{cb_i:04d}".replace("0", "A")
                       .replace("1", "C").replace("2", "G").replace("3", "T"))
        a.set_tag("UB", "UMIDUP" if cb_i < 0 else f"UMI{n:06d}")
        a.set_tag("RG", "s")
        recs.append(a)
    recs.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for r in recs:
            out.write(r)
    pysam.index(str(path))
    return path


@pytest.fixture(scope="module")
def clip_bam(tmp_path_factory) -> Path:
    return _mk_bam(tmp_path_factory.mktemp("polya_units") / "units.bam")


def _run(bam: Path, out: Path, tag: str, **kwargs):
    from ema.config import variable_config
    from ema.countmatrix.indexing import reset_index
    from ema.countmatrix.peackcalling import peak_calling
    from ema.countmatrix.peak import Peak
    from ema.strategies import get_strategy

    saved = {k: getattr(variable_config, k) for k in
             ("seqlen", "cb_len", "barcode_tag", "ignore_chro",
              "default_threshold", "merge_len")}
    variable_config.seqlen = SEQ_LEN
    variable_config.cb_len = 16
    variable_config.barcode_tag = "CB"
    variable_config.ignore_chro = []
    variable_config.default_threshold = 5
    variable_config.merge_len = 100
    reset_index()
    Peak.reset_pasnumber()
    bed = out / f"{tag}.bed"
    try:
        with patch.object(sys, "argv", ["ema", "--sequenceLen", str(SEQ_LEN),
                                        "--CellBarcodeLen", "16",
                                        "--BarcodeTag", "CB"]):
            peak_calling(False, bedfilepath=str(bed),
                         matrixpath=str(out / f"{tag}.mtx"),
                         bamfile_dir=str(bam),
                         strategy=get_strategy("clip_seeded"),
                         default_threshold=5, merge_len=100,
                         min_pas_spacing=SEQ_LEN, min_pas_prominence=5.0,
                         **kwargs)
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)
    return bed


def _read_support(bed: Path) -> dict[str, dict]:
    """The v2 columns of the sidecar, whatever else follows them.

    peakAtail-prime's ``--pas-features on`` APPENDS columns (paswrite's
    ``CALL_FEATURE_COLUMNS``, and more at the internal-priming seam).  This
    module is about the v2 columns' VALUES, so it pins their names and
    positions -- which is the guarantee that matters -- and ignores the tail.
    """
    lines = Path(support_path_for(bed)).read_text().splitlines()
    header = lines[0].split("\t")
    assert header[:len(SUPPORT_COLUMNS)] == list(SUPPORT_COLUMNS), (
        "the sidecar's v2 columns must keep their names AND their positions; "
        "new columns are APPENDED, never inserted or reordered"
    )
    out = {}
    for line in lines[1:]:
        f = line.split("\t")[:len(SUPPORT_COLUMNS)]
        out[f[0]] = dict(zip(SUPPORT_COLUMNS[1:], [int(x) for x in f[1:]]))
    return out


def test_sidecar_is_written_and_joins_the_bed_by_pas_id(clip_bam, tmp_path):
    bed = _run(clip_bam, tmp_path, "sidecar")
    rows = [line.split("\t") for line in bed.read_text().splitlines() if line.strip()]
    sup = _read_support(bed)
    assert len(rows) == len(sup) and rows, "one sidecar row per BED row"
    for r in rows:
        assert len(r) == 6, "pasbed.bed must stay plain BED6"
        s = sup[r[3]]
        assert int(r[4]) == s["clip_umis"], "BED column 5 == distinct molecules"
        assert s["tier"] == (1 if int(r[4]) > 0 else 2)
    tier1 = [s for s in sup.values() if s["tier"] == 1]
    assert tier1, "the clipped locus must produce a tier-1 PAS"
    top = max(tier1, key=lambda s: s["clip_reads"])
    assert top["clip_reads"] > top["clip_umis"], (
        "the PCR stack must collapse: 10 reads of one molecule"
    )
    assert top["clip_reads_f3844"] < top["clip_reads"], (
        "the 6 secondary alignments must be excluded from the -F 3844 count"
    )
    assert top["clip_umis_f3844"] <= top["clip_umis"]
    assert top["window_reads"] > 0


def test_clip_filter_flag_end_to_end_keeps_coordinates(clip_bam, tmp_path):
    """Switching the support unit may drop a PAS but must never move one."""
    a = _run(clip_bam, tmp_path, "unit_none")
    b = _run(clip_bam, tmp_path, "unit_f3844", polya_clip_filter="f3844")
    coords_a = {tuple(l.split("\t")[:3]) for l in a.read_text().splitlines()}
    coords_b = {tuple(l.split("\t")[:3]) for l in b.read_text().splitlines()}
    assert coords_b <= coords_a, "f3844 may only remove PAS, never move them"
