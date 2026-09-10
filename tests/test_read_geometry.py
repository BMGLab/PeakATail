"""``--read-geometry`` (peakAtail-prime TASK A): the read acceptance geometry.

v2 (``fixed``) normalises every accepted read to exactly ``--seq-len`` bp of
REFERENCE span: a read whose span is longer is DISCARDED, a shorter one has
its end rewritten to ``start + seq_len``.  Measured on the PBMC chr19+21 dev
slice (``results/prime/read_geometry_census_pbmc_slice.tsv``):

* the discard removes 11,768,752 / 48,647,964 = **24.19 %** of valid-CB reads
  before BOTH the poly(A) clip detector and the count matrix, 98.08 % of them
  spliced (the genome-wide A2 census says 13.74 % / 96 %; the slice is
  spliced-richer, so slice numbers here are an over-estimate of the genome);
* the rewrite fabricates the 3' end of a further 4,752,308 reads (9.77 %),
  pushing it a mean 13.53 bp DOWNSTREAM of the read's real last aligned base.

This module pins the three modes, and above all the property that makes the
change safe to reason about: **``keep`` and ``true`` accept every read
``fixed`` accepts.**  A change justified by "we were throwing reads away"
must not throw a different set away, and the obvious spelling of the fix (a
query-length test) does exactly that -- it drops short alignments carrying a
long terminal soft clip, which is the shape of a poly(A) clip read.
"""
from __future__ import annotations

import itertools
from pathlib import Path

import pysam
import pytest

from ema.config import variable_config
from ema.countmatrix.polya import ClipSeeder, ClipStream
from ema.countmatrix.read import (
    READ_GEOMETRIES,
    V2_READ_GEOMETRY,
    read_check,
    reset_composite_cache,
)

BARCODE = "AAAACCCCGGGGTTTT"
SEQ_LEN = 91
CHROM = "1"
CHROM_LEN = 200_000

# CIGAR op codes
M, I, D, N, S = 0, 1, 2, 3, 4

_HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": CHROM, "LN": CHROM_LEN}],
    "RG": [{"ID": "s1", "SM": "s1"}],
})


def _seg(cigar, start=1_000, reverse=False, flag=None, cb=BARCODE, rg="s1"):
    """A pysam ``AlignedSegment`` with the given CIGAR (not a mock).

    Using the real class matters: ``read_check`` reads ``reference_end`` and
    ``cigartuples``, both of which pysam derives from the CIGAR in C.
    """
    a = pysam.AlignedSegment(header=_HEADER)
    a.query_name = "r"
    qlen = sum(n for op, n in cigar if op in (M, I, S))
    a.query_sequence = "C" * qlen
    a.flag = (16 if reverse else 0) if flag is None else flag
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.cigartuples = list(cigar)
    a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
    if cb is not None:
        a.set_tag("CB", cb)
    a.set_tag("RG", rg)
    return a


def _check(seg, geometry, *, direction=False, exclude_flags=0):
    return read_check(seg, direction=direction, barcode="CB", barcode_len=16,
                      seq_len=SEQ_LEN, ignore_chro=(), sample_id="ds",
                      geometry=geometry, exclude_flags=exclude_flags)


@pytest.fixture(autouse=True)
def _clean():
    reset_composite_cache()
    yield
    reset_composite_cache()


# ---------------------------------------------------------------------------
# the flag itself
# ---------------------------------------------------------------------------

def test_the_branch_default_is_single_valued():
    """One default, in one place — the bug this test exists for.

    An earlier revision left ``variable_config.read_geometry`` at v2's
    ``"fixed"`` and put the branch default only in ``RunConfig``.  ``peakatail run``
    then used one geometry and a direct library call the other, and
    ``RunConfig.apply_to_legacy_globals()`` leaked the branch value into the
    process globals part-way through a pytest session, making 15 unrelated
    tests order-dependent.
    """
    from ema.cli.config_schema import RunConfig

    assert variable_config.read_geometry in READ_GEOMETRIES
    assert RunConfig().read_geometry == variable_config.read_geometry
    assert V2_READ_GEOMETRY == "fixed"


def test_peak_calling_rejects_an_unknown_geometry(tmp_path):
    from ema.countmatrix.peackcalling import peak_calling
    from ema.strategies import get_strategy

    saved = variable_config.read_geometry
    variable_config.read_geometry = "nonsense"
    try:
        with pytest.raises(ValueError, match="read_geometry"):
            peak_calling(False, str(tmp_path / "b.bed"), str(tmp_path / "m.mtx"),
                         bamfile_dir=str(tmp_path / "missing.bam"),
                         strategy=get_strategy("clip_seeded"))
    finally:
        variable_config.read_geometry = saved


# ---------------------------------------------------------------------------
# the safety property
# ---------------------------------------------------------------------------

_SHAPES = [
    ("full length",                [(M, 91)]),
    ("short alignment",            [(M, 60)]),
    ("soft clip 3' (poly-A read)", [(M, 71), (S, 20)]),
    ("soft clip 5'",               [(S, 20), (M, 71)]),
    ("huge 3' soft clip",          [(M, 20), (S, 71)]),
    ("deletion inside",            [(M, 40), (D, 5), (M, 51)]),
    ("insertion inside",           [(M, 40), (I, 5), (M, 46)]),
    ("spliced, small intron",      [(M, 40), (N, 200), (M, 51)]),
    ("spliced, big intron",        [(M, 40), (N, 50_000), (M, 51)]),
    ("spliced + poly-A clip",      [(M, 30), (N, 4_000), (M, 41), (S, 20)]),
    ("two introns",                [(M, 30), (N, 900), (M, 30), (N, 900), (M, 31)]),
    ("long deletion (not intron)", [(M, 40), (D, 400), (M, 51)]),
]


@pytest.mark.parametrize("label,cigar", _SHAPES, ids=[s[0] for s in _SHAPES])
@pytest.mark.parametrize("reverse", [False, True])
def test_keep_and_true_accept_every_read_fixed_accepts(label, cigar, reverse):
    """The strict-superset property.  This is the whole safety argument."""
    seg = _seg(cigar, reverse=reverse)
    fixed = _check(seg, "fixed", direction=reverse)
    if fixed[0] == 0:
        return  # nothing to prove for a read v2 already rejected
    for mode in ("keep", "true"):
        got = _check(seg, mode, direction=reverse)
        assert got[0] != 0, (
            f"{mode} rejected a read that v2 ACCEPTS ({label}); the change is "
            "supposed to add reads, never remove any"
        )
        assert got[1] == fixed[1], "the 5'-most aligned base must not move"


@pytest.mark.parametrize("label,cigar", _SHAPES, ids=[s[0] for s in _SHAPES])
def test_a_query_length_rule_would_not_have_that_property(label, cigar):
    """Documents WHY the acceptance test is on the de-introned footprint.

    ``query_length > seq_len`` looks like the natural reading of "a
    span-vs-query-length test", but a short alignment with a long terminal
    soft clip has span <= seq_len (v2 keeps it) and query length > seq_len
    (a query-length rule drops it).  The `huge 3' soft clip` shape is exactly
    that, and it is the shape of a poly(A) clip read.
    """
    seg = _seg(cigar)
    span = seg.reference_end - seg.reference_start
    qlen = seg.query_length
    if span <= SEQ_LEN and qlen > SEQ_LEN:
        assert _check(seg, "fixed")[0] != 0
        assert _check(seg, "true")[0] != 0, "footprint rule must keep it"
        assert label == "huge 3' soft clip"


# ---------------------------------------------------------------------------
# what each mode actually returns
# ---------------------------------------------------------------------------

def test_fixed_is_v2_discard_and_pad():
    assert _check(_seg([(M, 60)]), "fixed")[1:3] == (1_000, 1_091)   # padded
    assert _check(_seg([(M, 91)]), "fixed")[1:3] == (1_000, 1_091)
    assert _check(_seg([(M, 40), (N, 200), (M, 51)]), "fixed")[0] == 0  # discarded


def test_keep_rescues_the_spliced_read_but_keeps_the_fixed_length_interval():
    spliced = _seg([(M, 40), (N, 200), (M, 51)])
    assert _check(spliced, "fixed")[0] == 0
    chro, start, end, _, _ = _check(spliced, "keep")
    assert (chro, start, end) == (CHROM, 1_000, 1_091), (
        "keep is the ablation arm: it stops discarding and changes nothing else"
    )
    # ... and it still pads a short read, exactly like v2
    assert _check(_seg([(M, 60)]), "keep")[1:3] == (1_000, 1_091)


def test_true_uses_the_real_aligned_footprint():
    # unspliced short read: the true last aligned base, not a fabricated one
    assert _check(_seg([(M, 60)]), "true")[1:3] == (1_000, 1_060)
    # the poly(A) clip is NOT part of the footprint (it is not aligned)
    assert _check(_seg([(M, 71), (S, 20)]), "true")[1:3] == (1_000, 1_071)
    # a 5' soft clip is not either, and reference_start is unaffected by it
    assert _check(_seg([(S, 20), (M, 71)]), "true")[1:3] == (1_000, 1_071)
    # D advances the reference and stays INSIDE the footprint
    assert _check(_seg([(M, 40), (D, 5), (M, 46)]), "true")[1:3] == (1_000, 1_091)
    # I consumes query only, so it does not
    assert _check(_seg([(M, 40), (I, 5), (M, 46)]), "true")[1:3] == (1_000, 1_086)


@pytest.mark.parametrize("intron", [200, 4_000, 50_000])
def test_true_removes_the_intron_from_the_span(intron):
    """The measurement that forced this: on 7.82 M valid-CB reads of the PBMC
    slice the introns carry 10.02 Gb, ~14x the real read mass over the same
    105 Mb.  A coverage state machine fed genomic spans becomes a gene-body
    detector, so N is removed from the interval."""
    seg = _seg([(M, 40), (N, intron), (M, 51)])
    chro, start, end, _, _ = _check(seg, "true")
    assert (start, end) == (1_000, 1_091), "40 + 51 aligned bases, no intron"
    assert seg.reference_end == 1_000 + 91 + intron, "fixture sanity"


def test_true_keeps_two_introns_out_as_well():
    seg = _seg([(M, 30), (N, 900), (M, 30), (N, 900), (M, 31)])
    assert _check(seg, "true")[1:3] == (1_000, 1_091)


def test_a_long_deletion_is_still_rejected_in_every_mode():
    """Only introns are rescued.  A read whose ALIGNED footprint really does
    exceed --seq-len is as wrong under prime as under v2."""
    seg = _seg([(M, 40), (D, 400), (M, 51)])
    for mode in READ_GEOMETRIES:
        assert _check(seg, mode)[0] == 0, mode


# ---------------------------------------------------------------------------
# --read-exclude-flags
# ---------------------------------------------------------------------------

def test_read_exclude_flags_defaults_to_counting_every_alignment():
    """v2 applies no ``-F``: a multimapper contributes one read of coverage
    per ALIGNMENT.  10.42 % of valid-CB reads on the PBMC chr19+21 slice are
    secondary alignments."""
    secondary = _seg([(M, 91)], flag=256)
    assert _check(secondary, "true")[0] != 0
    assert _check(secondary, "true", exclude_flags=256)[0] == 0
    assert _check(secondary, "true", exclude_flags=3844)[0] == 0
    assert _check(secondary, "fixed", exclude_flags=256)[0] == 0
    primary = _seg([(M, 91)])
    assert _check(primary, "true", exclude_flags=3844)[0] != 0


def test_read_exclude_flags_is_off_by_default_in_the_config():
    assert variable_config.read_exclude_flags == 0
    from ema.cli.config_schema import RunConfig
    assert RunConfig().read_exclude_flags == 0


# ---------------------------------------------------------------------------
# ClipStream: the invariant read_check used to guarantee
# ---------------------------------------------------------------------------

def test_clipstream_refuses_to_guess_seq_len_outside_v2_geometry():
    """v2 learns seq_len from ``end1 - start1`` of the first read.  Outside
    ``fixed`` that quantity varies per read, so guessing it would silently
    shift every minus-strand count window."""
    ClipStream()                                  # v2: fine, learns it
    for geometry in ("keep", "true"):
        with pytest.raises(ValueError, match="seq_len"):
            ClipStream(geometry=geometry, direction=True)
    ClipStream(seq_len=SEQ_LEN, geometry="true", direction=True)


def test_three_key_is_the_transcript_three_prime_end_in_end1_space():
    # + strand: the key IS end1 (which under `true` is the real last base)
    plus = ClipStream(seq_len=SEQ_LEN, geometry="true", direction=False)
    assert plus.three_key(1_000, 1_060) == 1_060
    # - strand: the transcript 3' end is start1, stored as start1 + seq_len so
    # it lands exactly on the cluster anchor `mode + seq_len`
    minus = ClipStream(seq_len=SEQ_LEN, geometry="true", direction=True)
    assert minus.three_key(1_000, 1_060) == 1_000 + SEQ_LEN
    # v2: both reduce to end1, because end1 == start1 + seq_len anyway
    v2 = ClipStream()
    assert v2.three_key(1_000, 1_091) == 1_091


def test_finalize_restores_the_sorted_invariant_count_ends_bisects_on():
    """Under ``true`` geometry on ``+`` a soft-clipped read's true end can
    fall behind its predecessor's.  ``count_ends`` bisects, so the arrays are
    repaired once, at flush."""
    st = ClipStream(seq_len=SEQ_LEN, geometry="true", direction=False)
    for end, cb in ((1_091, "cb1"), (1_055, "cb2"), (1_200, "cb3"), (1_060, "cb4")):
        st.add_read(end - 60, end, cb)
    assert list(st.ends) == [1_091, 1_055, 1_200, 1_060], "fixture: unsorted"
    st.finalize()
    assert list(st.ends) == sorted(st.ends)
    counts: dict = {}
    n = st.count_ends(1_050, 1_100, counts)
    assert n == 3 and counts == {"cb1": 1, "cb2": 1, "cb4": 1}


def test_finalize_is_a_no_op_when_the_stream_is_already_sorted():
    """It must be, or it would be free to move a v2 byte."""
    st = ClipStream()
    for i in range(5):
        st.add_read(1_000 + i, 1_091 + i, f"cb{i}")
    before = list(st.ends), list(st.cbids)
    st.finalize()
    assert (list(st.ends), list(st.cbids)) == before


def test_clipstream_survives_the_pipeline_pickle_with_its_geometry():
    """The 3-stage pipeline hands the finder's stream to the writer."""
    import pickle

    st = ClipStream(seq_len=SEQ_LEN, geometry="true", direction=True)
    st.add_read(1_000, 1_040, "cb1")
    back = pickle.loads(pickle.dumps(st))
    assert back.geometry == "true" and back.direction is True
    assert back.seq_len == SEQ_LEN
    assert list(back.ends) == [1_000 + SEQ_LEN]


def test_seeder_maps_a_clip_reads_end_into_the_same_space_as_add_read():
    """Both feed the same midpoint arithmetic in the tier-1 clip fallback."""
    seeder = ClipSeeder(True, geometry="true", seq_len=SEQ_LEN)
    seeder.add_read(1_000, 1_040, "cb1")
    seeder.add_clip(1_000, "cb1", "UMI1", 1_040, True, start1=1_000)
    (site_rec,) = seeder.stream.sites.values()
    assert list(site_rec[4]) == [1_000 + SEQ_LEN] == list(seeder.stream.ends)


# ---------------------------------------------------------------------------
# end to end on the committed fixture BAM
# ---------------------------------------------------------------------------

FIXTURE = Path(__file__).parent / "fixtures" / "cellranger_pbmc_tiny.bam"


def _run(tmp_path: Path, geometry: str) -> tuple[int, int]:
    """Return (BED rows, total matrix mass) for one geometry."""
    from ema.countmatrix.indexing import BarcodeIndex
    from ema.countmatrix.peackcalling import peak_calling
    from ema.countmatrix.peak_state import PeakCallingState
    from ema.strategies import get_strategy

    saved = {k: getattr(variable_config, k)
             for k in ("seqlen", "cb_len", "barcode_tag", "read_geometry")}
    variable_config.seqlen = SEQ_LEN
    variable_config.cb_len = 16
    variable_config.barcode_tag = "CB"
    variable_config.read_geometry = geometry
    try:
        index, state = BarcodeIndex(), PeakCallingState(pasnumber=0)
        rows = mass = 0
        for direction in (False, True):
            bed = tmp_path / f"{geometry}_{int(direction)}.bed"
            mtx = tmp_path / f"{geometry}_{int(direction)}.mtx"
            peak_calling(direction=direction, bedfilepath=str(bed),
                         matrixpath=str(mtx), bamfile_dir=str(FIXTURE),
                         index=index, state=state, sample_id="fixture",
                         strategy=get_strategy("clip_seeded"))
            rows += sum(1 for ln in bed.read_text().splitlines() if ln.strip())
            mass += sum(int(ln.split()[2])
                        for ln in mtx.read_text().splitlines() if ln.strip())
        return rows, mass
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)


@pytest.mark.skipif(not FIXTURE.exists(), reason="fixture BAM not present")
def test_the_recovered_reads_reach_the_count_matrix(tmp_path):
    """The point of TASK A: the discarded reads are lost to the MATRIX, not
    just to the caller.  Genome-wide that is 88.8 M reads (A2 census)."""
    fixed_rows, fixed_mass = _run(tmp_path, "fixed")
    keep_rows, keep_mass = _run(tmp_path, "keep")
    assert keep_mass > fixed_mass, (
        f"--read-geometry keep must add matrix mass ({keep_mass} vs "
        f"{fixed_mass}); it is v2 with the discard removed and nothing else"
    )
    assert fixed_rows and keep_rows


@pytest.mark.skipif(not FIXTURE.exists(), reason="fixture BAM not present")
def test_every_geometry_runs_end_to_end_and_they_differ(tmp_path):
    out = {g: _run(tmp_path, g) for g in READ_GEOMETRIES}
    assert all(rows > 0 for rows, _ in out.values())
    assert len({v for v in out.values()}) > 1, (
        "the three geometries produced identical output — the flag is inert"
    )


# ---------------------------------------------------------------------------
# the run record
# ---------------------------------------------------------------------------

def test_run_config_records_the_geometry_that_actually_ran():
    """``run_config.json`` must say which geometry produced the output.

    Without this the shipped/v2/prime comparison cannot be audited from the
    run tree, and the ``args`` block is no substitute: it is the argparse
    namespace and carries schema DEFAULTS for anything Click resolved, so a
    run invoked with ``--read-geometry fixed`` still shows
    ``args.read_geometry == "true"`` there (the same class of defect bug B0
    fixed for atlas/gtf).
    """
    from ema.outputs import build_resolved_run_config

    saved = variable_config.read_geometry
    try:
        for geometry in READ_GEOMETRIES:
            variable_config.read_geometry = geometry
            v = build_resolved_run_config()["variables"]
            assert v["read_geometry"] == geometry
            assert "read_exclude_flags" in v
    finally:
        variable_config.read_geometry = saved
