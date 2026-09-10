"""Per-(contig, strand) parallel peak calling must be BYTE-IDENTICAL to the
legacy single-process two-pass caller (ema/countmatrix/chrom_parallel.py).

A synthetic coordinate-sorted BAM with four contigs (one of them the ignored
``MT``), peaks on both strands, poly(A)-clipped loci, cell barcodes that are
first seen in a different order on every contig, and a few reads without an
RG tag is called both ways; the five files the pipeline consumes
(``pos.bed``, ``neg.bed``, ``pos.mtx``, ``neg.mtx``, ``cb.tsv``) and both
support sidecars are compared as bytes.  The merge helpers are also pinned
on hand-made parts.
"""
from __future__ import annotations

import filecmp
import sys
from pathlib import Path
from unittest.mock import patch

import pysam
import pytest

SEQ_LEN = 91
N_RAMP = 60
N_SUMMIT = 40
CONTIGS = (("1", 30_000), ("2", 30_000), ("3", 22_000), ("MT", 16_000))
#: (contig, strand_is_reverse, site, clipped)
LOCI = (
    ("1", False, 2_400, True), ("1", False, 12_400, False), ("1", True, 20_400, True),
    ("2", False, 5_400, False), ("2", True, 9_400, True), ("2", True, 21_400, False),
    ("3", False, 3_400, True), ("3", True, 15_400, False),
    ("MT", False, 4_400, True), ("MT", True, 9_400, False),
)


def _barcode(i: int) -> str:
    alphabet = "ACGT"
    return ("".join(alphabet[(i >> (2 * k)) & 3] for k in range(8)) + "TTGACCAG")[:16]


def _locus_specs(site: int, reverse: bool, clipped: bool):
    """``(start, aligned_len, clip_len)``; summit reads end (forward) or
    start (reverse) exactly at *site* so it is the cleavage coordinate.  The
    ramp lies on the transcript-body side of the site (left of it on ``+``,
    right of it on ``-``) so the pile is continuous in coordinate order."""
    clip = 12 if clipped else 0
    if reverse:
        specs = [(site + 30 + i * 4, 80, 0) for i in range(N_RAMP)]
        specs += [(site, 80 - (i % 10), clip) for i in range(N_SUMMIT)]
        return specs
    specs = [(site - 350 + i * 4, 80, 0) for i in range(N_RAMP)]
    for i in range(N_SUMMIT):
        start = site - 80 + (i % 10)
        specs.append((start, site - start + 1, clip))
    return specs


def write_multi_contig_bam(path: Path) -> Path:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": ln} for c, ln in CONTIGS],
        "RG": [{"ID": "testsample", "SM": "testsample"}],
    }
    tid = {c: i for i, (c, _) in enumerate(CONTIGS)}
    records = []
    n = 0
    for li, (contig, reverse, site, clipped) in enumerate(LOCI):
        for i, (start, aligned, clip) in enumerate(_locus_specs(site, reverse, clipped)):
            a = pysam.AlignedSegment()
            a.query_name = f"r{n}"
            a.reference_id = tid[contig]
            a.reference_start = start
            a.mapping_quality = 60
            if reverse:
                a.flag = 16
                a.query_sequence = "T" * clip + "C" * aligned
                a.cigartuples = ([(4, clip)] if clip else []) + [(0, aligned)]
            else:
                a.flag = 0
                a.query_sequence = "C" * aligned + "A" * clip
                a.cigartuples = [(0, aligned)] + ([(4, clip)] if clip else [])
            a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
            # a different first-seen barcode order on every locus
            a.set_tag("CB", _barcode((i * 7 + li * 11) % 37))
            a.set_tag("UB", f"UMI{n:06d}")
            if n % 13:  # a few reads without RG -> dataset-id fallback
                a.set_tag("RG", "testsample")
            records.append(a)
            n += 1
    records.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for rec in records:
            out.write(rec)
    pysam.index(str(path))
    return path


@pytest.fixture(scope="module")
def bam(tmp_path_factory) -> Path:
    return write_multi_contig_bam(tmp_path_factory.mktemp("chrom_parallel") / "multi.bam")


# spawn propagates sys.argv; the legacy ema.config shim re-parses it in the
# child (see tests/test_polya_three_path_agreement.py).
_ARGV = ["ema", "--sequenceLen", str(SEQ_LEN), "--CellBarcodeLen", "16", "--BarcodeTag", "CB"]

PEAK_KWARGS = dict(
    dynamic_threshold=False, floor_threshold=3, lambda_fold_change=2.0,
    lambda_window=5000, bam_threads=1, min_pas_spacing=-1,
    min_pas_prominence=5.0, polya_enabled=True, polya_min_clip=6,
    polya_min_purity=0.8, polya_window=100, polya_seed_window=25,
    polya_min_umis=1, polya_clip_filter="none", polya_count_window=(-1, 25),
)


@pytest.fixture(autouse=True)
def _config():
    from ema.config import variable_config

    # peakAtail-prime: read_geometry / pas_features / clip_rate_sampling are
    # process-global too, and the geometry test below deliberately moves one.
    keys = ("seqlen", "cb_len", "barcode_tag", "ignore_chro", "default_threshold",
            "merge_len", "read_geometry", "pas_features", "clip_rate_sampling")
    saved = {k: getattr(variable_config, k) for k in keys}
    variable_config.seqlen = SEQ_LEN
    variable_config.cb_len = 16
    variable_config.barcode_tag = "CB"
    variable_config.ignore_chro = ["MT"]
    variable_config.default_threshold = 5
    variable_config.merge_len = 100
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)


def _sequential(bam: Path, out: Path, strategy_name: str, *, reset_pasnumber: bool = True) -> dict[str, Path]:
    """What ema/main.py's legacy loop does: two full passes on one index.
    ``reset_pasnumber=False`` is what the loop does for the SECOND BAM of a
    multi-dataset run: the index is reset, ``Peak.pasnumber`` is not."""
    from ema.countmatrix.indexing import get_mapping, reset_index
    from ema.countmatrix.peackcalling import peak_calling
    from ema.countmatrix.peak import Peak
    from ema.countmatrix.read import set_default_sample_id
    from ema.strategies import get_strategy

    out.mkdir()
    reset_index()
    if reset_pasnumber:
        Peak.reset_pasnumber()
    set_default_sample_id("default")
    files = {k: out / f"default_0.{k}" for k in ("pos.bed", "neg.bed", "pos.mtx", "neg.mtx", "cb.tsv")}
    for direction, b, m in ((False, "pos.bed", "pos.mtx"), (True, "neg.bed", "neg.mtx")):
        peak_calling(
            direction, bedfilepath=str(files[b]), matrixpath=str(files[m]),
            bamfile_dir=str(bam), default_threshold=5, merge_len=100,
            strategy=get_strategy(strategy_name), **PEAK_KWARGS,
        )
    with open(files["cb.tsv"], "w") as fh:
        for cb, _ in sorted(get_mapping().items(), key=lambda kv: kv[1]):
            fh.write(cb + "\n")
    return files


def _parallel(bam: Path, out: Path, strategy_name: str, n_workers: int, *, reset_pasnumber: bool = True) -> tuple[dict[str, Path], dict]:
    """What ema/main.py's default path does for one BAM.  The dispatcher
    seeds the merged ids from ``Peak.pasnumber`` exactly like the legacy
    ``peak_calling()`` (0 in a fresh ``peakatail run`` process; ``_sequential``
    leaves it at its last id, so reset it here unless the test wants the
    multi-BAM continuation)."""
    from ema.countmatrix.chrom_parallel import run_chrom_parallel
    from ema.countmatrix.peak import Peak

    out.mkdir()
    if reset_pasnumber:
        Peak.reset_pasnumber()
    files = {k: out / f"default_0.{k}" for k in ("pos.bed", "neg.bed", "pos.mtx", "neg.mtx", "cb.tsv")}
    with patch.object(sys, "argv", _ARGV):
        summary = run_chrom_parallel(
            dataset_id="default", bam_path=str(bam),
            pos_bed=str(files["pos.bed"]), neg_bed=str(files["neg.bed"]),
            pos_mtx=str(files["pos.mtx"]), neg_mtx=str(files["neg.mtx"]),
            cb_tsv=str(files["cb.tsv"]), n_workers=n_workers,
            strategy_name=strategy_name, strategy_kwargs={},
            peak_kwargs=dict(PEAK_KWARGS), workdir=str(out / "_jobs"),
            thread_budget=4,
        )
    return files, summary


def _assert_same_bytes(a: Path, b: Path):
    assert a.exists() and b.exists(), (a, b)
    assert filecmp.cmp(a, b, shallow=False), f"{a.name} differs between sequential and parallel"


@pytest.mark.parametrize("strategy_name", ["lambda_gradient", "clip_seeded"])
def test_parallel_is_byte_identical_to_sequential(bam, tmp_path, strategy_name):
    from ema.countmatrix.paswrite import support_path_for

    seq = _sequential(bam, tmp_path / "seq", strategy_name)
    par, summary = _parallel(bam, tmp_path / "par", strategy_name, n_workers=3)

    # The synthetic data must actually exercise the merge: PAS on both
    # strands of several contigs, clip-supported rows, shared barcodes.
    pos_rows = seq["pos.bed"].read_text().splitlines()
    neg_rows = seq["neg.bed"].read_text().splitlines()
    assert len({r.split("\t")[0] for r in pos_rows}) >= 3
    assert len({r.split("\t")[0] for r in neg_rows}) >= 3
    assert not any(r.startswith("MT\t") for r in pos_rows + neg_rows)
    assert any(int(r.split("\t")[4]) > 0 for r in pos_rows + neg_rows)
    assert summary["n_pas"] == len(pos_rows) + len(neg_rows)
    assert len(summary["jobs"]) == 6  # 3 contigs x 2 strands, MT skipped
    assert not (tmp_path / "par" / "_jobs").exists()  # scratch removed

    for k in ("pos.bed", "neg.bed", "pos.mtx", "neg.mtx", "cb.tsv"):
        _assert_same_bytes(seq[k], par[k])
    for k in ("pos.bed", "neg.bed"):
        _assert_same_bytes(Path(support_path_for(seq[k])), Path(support_path_for(par[k])))


def test_pas_ids_continue_across_datasets_like_the_legacy_loop(bam, tmp_path):
    """main.py resets the barcode index per BAM but never Peak.pasnumber, so
    the second BAM of a multi-dataset run is numbered from the first one's
    last id.  The parallel path must seed its merge from the same counter and
    write the final value back, or the per-dataset raw files of BAM 2+ would
    differ from the legacy loop's."""
    from ema.countmatrix.peak import Peak
    from ema.countmatrix.paswrite import support_path_for

    seq1 = _sequential(bam, tmp_path / "seq1", "clip_seeded")
    n1 = Peak.pasnumber
    assert n1 == len(seq1["pos.bed"].read_text().splitlines()) + len(seq1["neg.bed"].read_text().splitlines())
    seq2 = _sequential(bam, tmp_path / "seq2", "clip_seeded", reset_pasnumber=False)
    assert Peak.pasnumber == 2 * n1
    assert seq2["pos.bed"].read_text().splitlines()[0].split("\t")[3] == str(n1 + 1)

    # Parallel "dataset 1" (fresh counter) then "dataset 2" (counter carried).
    par1, s1 = _parallel(bam, tmp_path / "par1", "clip_seeded", n_workers=2)
    assert Peak.pasnumber == n1 and s1["n_pas"] == n1
    par2, s2 = _parallel(bam, tmp_path / "par2", "clip_seeded", n_workers=2, reset_pasnumber=False)
    assert Peak.pasnumber == 2 * n1 and s2["n_pas"] == n1

    for seq, par in ((seq1, par1), (seq2, par2)):
        for k in ("pos.bed", "neg.bed", "pos.mtx", "neg.mtx", "cb.tsv"):
            _assert_same_bytes(seq[k], par[k])
        for k in ("pos.bed", "neg.bed"):
            _assert_same_bytes(Path(support_path_for(seq[k])), Path(support_path_for(par[k])))
    Peak.reset_pasnumber()


def test_strategy_kwargs_reach_the_workers(bam, tmp_path):
    """A tunable set on the CLI (e.g. --max-pas) must survive the spawn: the
    workers rebuild the strategy by name, so the kwargs travel in the job."""
    from ema.countmatrix.chrom_parallel import run_chrom_parallel
    from ema.countmatrix.indexing import reset_index
    from ema.countmatrix.peackcalling import peak_calling
    from ema.countmatrix.peak import Peak
    from ema.countmatrix.read import set_default_sample_id
    from ema.strategies import get_strategy

    kwargs = {"max_pas": 1, "min_prominence": 1.0}
    seq = tmp_path / "seq"
    seq.mkdir()
    reset_index()
    Peak.reset_pasnumber()
    set_default_sample_id("default")
    seq_bed = seq / "default_0.pos.bed"
    peak_calling(
        False, bedfilepath=str(seq_bed), matrixpath=str(seq / "default_0.pos.mtx"),
        bamfile_dir=str(bam), default_threshold=5, merge_len=100,
        strategy=get_strategy("lambda_gradient", **kwargs), **PEAK_KWARGS,
    )
    par = tmp_path / "par"
    par.mkdir()
    files = {k: par / f"default_0.{k}" for k in ("pos.bed", "neg.bed", "pos.mtx", "neg.mtx", "cb.tsv")}
    Peak.reset_pasnumber()  # the dispatcher seeds ids from the counter, like peak_calling()
    with patch.object(sys, "argv", _ARGV):
        run_chrom_parallel(
            dataset_id="default", bam_path=str(bam),
            pos_bed=str(files["pos.bed"]), neg_bed=str(files["neg.bed"]),
            pos_mtx=str(files["pos.mtx"]), neg_mtx=str(files["neg.mtx"]),
            cb_tsv=str(files["cb.tsv"]), n_workers=2,
            strategy_name="lambda_gradient", strategy_kwargs=kwargs,
            peak_kwargs=dict(PEAK_KWARGS), workdir=str(par / "_jobs"),
            thread_budget=2,
        )
    _assert_same_bytes(seq_bed, files["pos.bed"])


def test_singleton_index_matches_merged_column_order(bam, tmp_path):
    """The CB filter reads get_mapping() -- it must see the merged order."""
    from ema.countmatrix.chrom_parallel import load_index_from_cb_list
    from ema.countmatrix.indexing import get_mapping

    seq = _sequential(bam, tmp_path / "seq", "clip_seeded")
    expected = get_mapping()
    par, summary = _parallel(bam, tmp_path / "par", "clip_seeded", n_workers=2)
    load_index_from_cb_list(par["cb.tsv"].read_text().splitlines())
    assert get_mapping() == expected
    assert list(get_mapping()) == list(expected)


def test_parallel_without_polya_evidence_writes_no_support(bam, tmp_path):
    from ema.countmatrix.paswrite import support_path_for

    kwargs = dict(PEAK_KWARGS, polya_enabled=False)
    with patch.dict(PEAK_KWARGS, kwargs):
        seq = _sequential(bam, tmp_path / "seq", "lambda_gradient")
        par, _ = _parallel(bam, tmp_path / "par", "lambda_gradient", n_workers=2)
    for k in ("pos.bed", "neg.bed", "pos.mtx", "neg.mtx", "cb.tsv"):
        _assert_same_bytes(seq[k], par[k])
    assert not Path(support_path_for(par["pos.bed"])).exists()
    assert not Path(support_path_for(seq["pos.bed"])).exists()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def test_contigs_with_reads_skips_ignored_and_empty(bam):
    from ema.countmatrix.chrom_parallel import can_run_parallel, contigs_with_reads

    assert can_run_parallel(str(bam))
    got = contigs_with_reads(str(bam), ["MT"])
    assert [c for c, _, _ in got] == ["1", "2", "3"]  # header order
    assert all(n > 0 for _, _, n in got)
    assert [c for c, _, _ in contigs_with_reads(str(bam), [])] == ["1", "2", "3", "MT"]


def test_can_run_parallel_requires_index(tmp_path):
    from ema.countmatrix.chrom_parallel import can_run_parallel

    assert not can_run_parallel(str(tmp_path / "missing.bam"))


def test_worker_bam_threads_respects_budget():
    from ema.countmatrix.chrom_parallel import worker_bam_threads

    assert worker_bam_threads(4, 16, 16) == 1
    assert worker_bam_threads(4, 4, 16) == 4
    assert worker_bam_threads(4, 3, 16) == 4
    assert worker_bam_threads(4, 6, 16) == 2
    assert worker_bam_threads(4, 2, None) == 4
    assert worker_bam_threads(0, 2, None) == 1


def test_merge_renumbers_in_emission_order_and_keeps_first_write_cb_order(tmp_path):
    """Hand-made parts: '+' jobs first in header order, then '-'; pas ids
    continue across strands; the shared CB list is the first-write union."""
    from ema.countmatrix.chrom_parallel import merge_chrom_results
    from ema.countmatrix.paswrite import SUPPORT_COLUMNS

    def part(name, contig, direction, beds, mtx, cbs, sup):
        d = tmp_path / name
        d.mkdir()
        (d / "j.bed").write_text("".join(
            f"{contig}\t{s}\t{e}\t{i}\t{sc}\t{'-' if direction else '+'}\n"
            for i, (s, e, sc) in enumerate(beds, 1)))
        (d / "j.mtx").write_text("".join(f"{p} {c} {n}\n" for p, c, n in mtx))
        (d / "j.cb.tsv").write_text("".join(cb + "\n" for cb in cbs))
        (d / "j.support.tsv").write_text(
            "\t".join(SUPPORT_COLUMNS) + "\n" + "".join(f"{i}\t{r}\n" for i, r in enumerate(sup, 1)))
        return {"contig": contig, "direction": direction, "bed": str(d / "j.bed"),
                "mtx": str(d / "j.mtx"), "support": str(d / "j.support.tsv"), "cb": str(d / "j.cb.tsv")}

    results = [
        # arrival order deliberately scrambled
        part("b_neg", "2", True, [(10, 11, 3)], [(1, 1, 4)], ["s_B"], ["0\t3\t0\t3\t4\t1"]),
        part("a_pos", "1", False, [(100, 150, 0), (300, 301, 2)], [(1, 1, 5), (1, 2, 1), (2, 2, 7)],
             ["s_A", "s_B"], ["0\t0\t0\t0\t6\t2", "2\t2\t2\t2\t7\t1"]),
        part("b_pos", "2", False, [(20, 21, 1)], [(1, 1, 2), (1, 2, 3)], ["s_C", "s_A"], ["1\t1\t1\t1\t5\t1"]),
        part("a_neg", "1", True, [], [], [], []),
    ]
    out = {k: tmp_path / k for k in ("pos.bed", "neg.bed", "pos.mtx", "neg.mtx", "cb.tsv")}
    n_pas, cbs = merge_chrom_results(
        results, ["1", "2", "3"], str(out["pos.bed"]), str(out["neg.bed"]),
        str(out["pos.mtx"]), str(out["neg.mtx"]), str(out["cb.tsv"]),
    )
    assert n_pas == 4
    assert cbs == ["s_A", "s_B", "s_C"]
    assert out["cb.tsv"].read_text() == "s_A\ns_B\ns_C\n"
    assert out["pos.bed"].read_text() == (
        "1\t100\t150\t1\t0\t+\n1\t300\t301\t2\t2\t+\n2\t20\t21\t3\t1\t+\n")
    assert out["neg.bed"].read_text() == "2\t10\t11\t4\t3\t-\n"
    # columns: job b_pos's local 1 = s_C -> global 3, local 2 = s_A -> global 1
    assert out["pos.mtx"].read_text() == "1 1 5\n1 2 1\n2 2 7\n3 3 2\n3 1 3\n"
    assert out["neg.mtx"].read_text() == "4 2 4\n"
    pos_sup = (tmp_path / "pos.support.tsv").read_text().splitlines()
    assert pos_sup[0] == "\t".join(SUPPORT_COLUMNS)
    assert pos_sup[1:] == ["1\t0\t0\t0\t0\t6\t2", "2\t2\t2\t2\t2\t7\t1", "3\t1\t1\t1\t1\t5\t1"]
    assert (tmp_path / "neg.support.tsv").read_text().splitlines()[1:] == ["4\t0\t3\t0\t3\t4\t1"]


# ---------------------------------------------------------------------------
# peakAtail-prime: the spawned workers must inherit the prime knobs
# ---------------------------------------------------------------------------
# Only "true" is an arm here: on this fixture no read's reference span exceeds
# --seq-len, so nothing is discarded and "keep" is byte-equal to "fixed" -- a
# "keep" arm could not tell a broken ChromJob from a working one.
@pytest.mark.parametrize("geometry", ["true"])
def test_parallel_matches_sequential_under_a_non_default_read_geometry(
        bam, tmp_path, geometry):
    """``--read-geometry`` has to reach the SPAWNED per-contig workers.

    ``chrom_parallel`` is the path a real ``peakatail run --threads N`` takes, and
    its workers are spawned: ``ema.config``'s legacy globals come back at
    their MODULE defaults in the child.  ``read_geometry``'s module default is
    ``"fixed"``, and every other integration test in this suite runs at that
    same value -- so a ``ChromJob`` that forgot to carry the geometry would be
    invisible to all of them while the parent used one geometry and the
    children another.  The only symptom would be wrong coordinates.

    (This coverage existed implicitly while ``"true"`` was the branch default;
    it stopped existing when the measurement moved the default back to
    ``"fixed"``.  Hence an explicit arm.)
    """
    from ema.config import variable_config
    from ema.countmatrix.paswrite import support_path_for

    variable_config.read_geometry = geometry
    seq = _sequential(bam, tmp_path / "seq", "clip_seeded")
    par, summary = _parallel(bam, tmp_path / "par", "clip_seeded", n_workers=3)
    assert len(summary["jobs"]) == 6, "the merge was not exercised"
    for key in ("pos.bed", "neg.bed", "pos.mtx", "neg.mtx", "cb.tsv"):
        _assert_same_bytes(seq[key], par[key])
    for key in ("pos.bed", "neg.bed"):
        _assert_same_bytes(Path(support_path_for(str(seq[key]))),
                           Path(support_path_for(str(par[key]))))

    # ...and the arm is not vacuous: v2 geometry gives DIFFERENT bytes on this
    # fixture, so byte-identity above is a statement about the plumbing.
    variable_config.read_geometry = "fixed"
    v2 = _sequential(bam, tmp_path / "v2", "clip_seeded")
    assert v2["pos.bed"].read_bytes() != seq["pos.bed"].read_bytes(), (
        f"--read-geometry {geometry} produced v2's bytes on this fixture, so "
        "the test cannot tell 'the geometry reached the child' from 'nothing "
        "happened'"
    )


def test_parallel_matches_sequential_with_the_feature_columns_on(bam, tmp_path):
    """``--pas-features on`` decides the sidecar's COLUMN SET, and the merged
    sidecar's header is taken from the children.  Both halves have to agree
    across the spawn."""
    from ema.config import variable_config
    from ema.countmatrix.paswrite import (
        CALL_FEATURE_COLUMNS, SUPPORT_COLUMNS, support_path_for,
    )

    variable_config.pas_features = "on"
    seq = _sequential(bam, tmp_path / "seq", "clip_seeded")
    par, _ = _parallel(bam, tmp_path / "par", "clip_seeded", n_workers=3)
    for key in ("pos.bed", "neg.bed"):
        a = Path(support_path_for(str(seq[key])))
        b = Path(support_path_for(str(par[key])))
        _assert_same_bytes(a, b)
        header = a.read_text().splitlines()[0].split("\t")
        assert header == list(SUPPORT_COLUMNS) + list(CALL_FEATURE_COLUMNS), header
        rows = a.read_text().splitlines()[1:]
        assert rows, f"{key}: empty sidecar makes this check vacuous"
        for row in rows:
            assert len(row.split("\t")) == len(header)


def test_every_prime_knob_that_changes_output_travels_with_the_job():
    """A structural guard: a new prime option that a child would otherwise
    re-read from ``variable_config`` must be a ``ChromJob`` field AND be
    assigned in ``chrom_worker``.

    The failure mode this pins is silent and has already happened twice on
    this branch (``run_tiled``'s legacy dicts, and the sidecar header): the
    parent honours the flag, the spawned child does not, and the two halves of
    one run disagree.
    """
    import dataclasses
    import inspect

    from ema.countmatrix import chrom_parallel as CP

    fields = {f.name for f in dataclasses.fields(CP.ChromJob)}
    src = inspect.getsource(CP.chrom_worker)
    for knob in ("read_geometry", "read_exclude_flags", "pas_features",
                 "clip_rate_sampling"):
        assert knob in fields, f"ChromJob has no {knob} field"
        assert f"vc.{knob} = job.{knob}" in src, (
            f"chrom_worker never pushes {knob} into the child's variable_config"
        )
