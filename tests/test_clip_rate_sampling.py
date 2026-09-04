"""peakAtail-prime TASK E item 2 — the poly(A) clip-rate QC must sample the file.

v2 read the first 200,000 CB reads of a coordinate-sorted BAM, i.e. the head of
the first contig, and on PBMC 10k v3 it returned **2.2565 %** where the
whole-file rate on the same denominator is **0.5364 %**.  A 4x over-estimate is
the dangerous direction: the estimator exists to shout when the poly(A)
evidence channel has been destroyed, and an over-estimate MASKS exactly that.

The fixture here is built to have the same property in miniature -- a head
that is not representative of the file -- so the test states the defect rather
than the implementation.
"""
from __future__ import annotations

import pytest

pysam = pytest.importorskip("pysam")

from ema.config import variable_config
from ema.countmatrix import polya
from ema.countmatrix.polya import (
    CLIP_RATE_SAMPLINGS,
    V2_CLIP_RATE_SAMPLING,
    check_clip_rate,
)

CLIPPED = "A" * 30          # a qualifying poly(A) soft clip on a '+' read
PLAIN = "C" * 30


def _make_bam(path, per_contig):
    """One BAM, two contigs; *per_contig* is {name: [(n_reads, clipped)]}."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": name, "LN": 100_000} for name in per_contig]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for tid, (name, blocks) in enumerate(per_contig.items()):
            pos = 100
            for n_reads, clipped in blocks:
                for i in range(n_reads):
                    a = pysam.AlignedSegment()
                    a.query_name = "%s_%d_%d" % (name, pos, i)
                    a.query_sequence = ("G" * 60) + (CLIPPED if clipped else PLAIN)
                    a.flag = 0
                    a.reference_id = tid
                    a.reference_start = pos
                    a.mapping_quality = 60
                    a.cigartuples = [(0, 60), (4, 30)]   # 60M30S
                    a.query_qualities = pysam.qualitystring_to_array("I" * 90)
                    a.set_tag("CB", "AAACCCAAGAAACACT-1")
                    a.set_tag("UB", "ACGT%06d" % i)
                    out.write(a)
                    pos += 1
    pysam.index(str(path))
    return path


@pytest.fixture
def unrepresentative_bam(tmp_path):
    """chrA is 100 % clipped, chrB is 0 %; the file's true rate is 10 %."""
    path = tmp_path / "head_lies.bam"
    _make_bam(path, {"chrA": [(200, True)], "chrB": [(1800, False)]})
    return path


@pytest.fixture(autouse=True)
def _clear_cache():
    polya._clip_rate_warned.clear()
    saved = getattr(variable_config, "clip_rate_sampling", "strided")
    yield
    polya._clip_rate_warned.clear()
    variable_config.clip_rate_sampling = saved


def test_the_v2_sampler_is_wrong_on_this_file_and_that_is_the_point(
        unrepresentative_bam):
    rate = check_clip_rate(str(unrepresentative_bam), max_reads=200,
                           sampling="head")
    assert rate == pytest.approx(1.0), (
        "v2 reads the head of the file, so it sees only the clipped contig"
    )


def test_the_strided_sampler_recovers_the_files_own_rate(unrepresentative_bam):
    rate = check_clip_rate(str(unrepresentative_bam), max_reads=2000,
                           sampling="strided")
    assert rate == pytest.approx(0.10, abs=0.02), (
        "strided sampling must land on the whole-file rate, not the head's"
    )


def test_the_strided_sampler_is_closer_at_the_same_budget(unrepresentative_bam):
    truth = 0.10
    polya._clip_rate_warned.clear()
    head = check_clip_rate(str(unrepresentative_bam), max_reads=400,
                           sampling="head")
    polya._clip_rate_warned.clear()
    strided = check_clip_rate(str(unrepresentative_bam), max_reads=400,
                              sampling="strided")
    assert abs(strided - truth) < abs(head - truth)


def test_the_default_comes_from_variable_config(unrepresentative_bam):
    variable_config.clip_rate_sampling = "head"
    assert check_clip_rate(str(unrepresentative_bam), max_reads=200) == \
        pytest.approx(1.0)
    polya._clip_rate_warned.clear()
    variable_config.clip_rate_sampling = "strided"
    assert check_clip_rate(str(unrepresentative_bam), max_reads=2000) == \
        pytest.approx(0.10, abs=0.02)


def test_an_unknown_sampling_is_an_error_not_a_silent_head_scan(
        unrepresentative_bam):
    with pytest.raises(ValueError):
        check_clip_rate(str(unrepresentative_bam), sampling="middle")


def test_the_v2_value_is_named_and_is_head():
    assert V2_CLIP_RATE_SAMPLING == "head"
    assert set(CLIP_RATE_SAMPLINGS) == {"head", "strided", "pass"}


def test_pass_does_not_estimate_at_all(unrepresentative_bam):
    """'pass' is not a sampler: the peak-calling loop counts the real thing, so
    check_clip_rate must do nothing rather than return a number nobody asked
    for (and must not put an estimate in the log)."""
    assert check_clip_rate(str(unrepresentative_bam), sampling="pass") is None


def test_the_exact_counter_is_exact():
    from ema.countmatrix.polya import ClipRateCounter

    c = ClipRateCounter()
    assert c.rate == 0.0 and c.n_reads == 0
    for has_clip in (True, False, False, False, False,
                     False, False, False, False, False):
        c.add(has_clip)
    assert c.n_reads == 10 and c.n_clip == 1
    assert c.rate == pytest.approx(0.10)


def test_a_single_contig_strand_job_does_not_raise_the_alarm(caplog):
    """chr21's '+' strand alone measures 0.2878 % on a healthy PBMC library --
    under the 0.3 % bar.  A per-job warning would fire on a good BAM, so the
    per-job line is information and the alarm belongs to the run-level total."""
    import logging

    from ema.countmatrix.polya import LOW_CLIP_RATE, ClipRateCounter

    c = ClipRateCounter()
    for i in range(100_000):
        c.add(i % 1000 == 0)             # 0.1 %, well under the bar
    assert c.rate < LOW_CLIP_RATE
    with caplog.at_level(logging.INFO):
        c.report("21+ strand")           # warn defaults to False
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert any("poly(A) clip rate (exact" in r.getMessage()
               for r in caplog.records)


def test_the_run_level_total_shouts_below_the_bar(caplog):
    import logging

    from ema.countmatrix.polya import report_clip_rate_total

    with caplog.at_level(logging.WARNING):
        rate = report_clip_rate_total(100_000, 100, "some.bam")
    assert rate == pytest.approx(0.001)
    msgs = [r.getMessage() for r in caplog.records]
    assert any("LOW POLY(A) CLIP RATE (exact" in m for m in msgs), msgs
    assert any("0.5730" in m for m in msgs), (
        "the warning must carry the reference rate on the SAME denominator"
    )


def test_the_run_level_total_is_quiet_above_the_bar(caplog):
    import logging

    from ema.countmatrix.polya import report_clip_rate_total

    with caplog.at_level(logging.WARNING):
        report_clip_rate_total(100_000, 600, "some.bam")
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


def test_the_counts_travel_to_the_dispatcher(caplog):
    """The two integers reach chrom_parallel without widening peak_calling's
    return type, and are consumed exactly once."""
    from ema.countmatrix.polya import ClipRateCounter, take_last_clip_counts

    take_last_clip_counts()                       # clear
    c = ClipRateCounter()
    for i in range(1_000):
        c.add(i % 100 == 0)
    c.report("19+ strand")
    assert take_last_clip_counts() == (1_000, 10)
    assert take_last_clip_counts() == (0, 0), "a second read must not double-count"


def test_the_exact_counter_is_quiet_above_the_bar(caplog):
    import logging

    from ema.countmatrix.polya import ClipRateCounter

    c = ClipRateCounter()
    for i in range(100_000):
        c.add(i % 100 == 0)              # 1 %
    with caplog.at_level(logging.WARNING):
        c.report("19+ strand", warn=True)
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


def test_the_branch_default_is_single_valued():
    import dataclasses

    from ema.cli.config_schema import RunConfig
    from ema.config import variable_config as vc

    schema = {f.name: f.default for f in dataclasses.fields(RunConfig)}
    assert schema["clip_rate_sampling"] == vc.clip_rate_sampling == "pass"


def test_an_unindexed_bam_falls_back_and_says_so(tmp_path, caplog):
    path = tmp_path / "noindex.bam"
    _make_bam(path, {"chrA": [(50, True)], "chrB": [(450, False)]})
    (tmp_path / "noindex.bam.bai").unlink()
    import logging
    with caplog.at_level(logging.WARNING):
        rate = check_clip_rate(str(path), max_reads=500, sampling="strided")
    assert rate is not None
    assert any("no usable index" in r.getMessage() for r in caplog.records), \
        "an unindexed BAM must say the sample could not be spread"


def test_the_result_is_cached_per_bam(unrepresentative_bam):
    first = check_clip_rate(str(unrepresentative_bam), max_reads=2000,
                            sampling="strided")
    # a second call with the v2 sampler must NOT rescan -- the cache wins
    assert check_clip_rate(str(unrepresentative_bam), max_reads=200,
                           sampling="head") == first


def test_a_broken_bam_never_fails_the_run(tmp_path):
    bad = tmp_path / "not.bam"
    bad.write_text("this is not a BAM\n")
    assert check_clip_rate(str(bad), sampling="strided") is None


# --------------------------------------------------------------------------
# the design property the FIRST implementation got wrong
# --------------------------------------------------------------------------
def test_the_estimate_does_not_depend_on_the_window_count(tmp_path):
    """A first implementation gave each stratum an equal READ quota starting at
    an evenly spaced coordinate.  That weights read-sparse coordinate regions
    equally with dense ones, and it showed: on the PBMC chr19+21 slice
    (truth 0.6169 %) it returned 0.9111 % at 25 strata and 0.1853 % at 2,500 --
    a 3x swing driven by nothing but the stratum count.  Taking EVERY read in a
    coordinate-uniform window instead estimates the same quantity whatever the
    count, which is what this test pins.
    """
    import pysam as _pysam

    from ema.countmatrix.polya import _sample_strided

    # dense and sparse blocks alternating, same clip rate in both, so the only
    # thing a window count can change is which blocks are sampled.
    path = tmp_path / "mixed.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chrA", "LN": 20_000}]}
    with _pysam.AlignmentFile(str(path), "wb", header=header) as out:
        n = 0
        for pos in range(0, 20_000, 10):
            dense = (pos // 1000) % 2 == 0
            for i in range(20 if dense else 2):
                a = _pysam.AlignedSegment()
                a.query_name = "r%d" % n
                clipped = (n % 10 == 0)          # exactly 10 % everywhere
                a.query_sequence = ("G" * 60) + (CLIPPED if clipped else PLAIN)
                a.flag = 0
                a.reference_id = 0
                a.reference_start = pos
                a.mapping_quality = 60
                a.cigartuples = [(0, 60), (4, 30)]
                a.query_qualities = _pysam.qualitystring_to_array("I" * 90)
                a.set_tag("CB", "AAACCCAAGAAACACT-1")
                a.set_tag("UB", "ACGT%06d" % i)
                out.write(a)
                n += 1
    _pysam.index(str(path))

    rates = []
    for strata in (5, 20, 80):
        with _pysam.AlignmentFile(str(path), "rb") as bam:
            n_cb, n_clip, _d = _sample_strided(bam, 6, 0.8, "CB", 20_000,
                                               strata=strata)
        assert n_cb > 0
        rates.append(n_clip / n_cb)
    assert max(rates) - min(rates) < 0.02, rates
    for r in rates:
        assert r == pytest.approx(0.10, abs=0.02), rates


def test_a_read_is_counted_in_exactly_one_window(tmp_path):
    """Windows are half-open on the read START, so a read overlapping two of
    them is counted once -- otherwise the denominator drifts with the window
    width and the rate with it."""
    import pysam as _pysam

    from ema.countmatrix.polya import _sample_strided

    path = tmp_path / "wide.bam"
    _make_bam(path, {"chrA": [(500, False)]})
    with _pysam.AlignmentFile(str(path), "rb") as bam:
        n_cb, _n_clip, _d = _sample_strided(bam, 6, 0.8, "CB", 10_000_000,
                                            strata=7)
    # budget far exceeds the file, so every read must be sampled exactly once
    assert n_cb == 500
