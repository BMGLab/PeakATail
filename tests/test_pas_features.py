"""TASK C -- the per-site scoring features written into ``pas_support.tsv``.

What these tests protect, in order of what it would cost to get wrong:

1. **v2 byte-identity.**  ``--pas-features off`` must leave the sidecar exactly
   as v2 wrote it.  ``test_prime_v2_compat_golden.py`` hashes the whole file
   against goldens from the frozen v2 worktree; here we pin the narrower and
   more diagnostic statement: turning the flag ON changes the SIDECAR and
   NOTHING ELSE -- the BED and the count matrix stay byte-identical.
2. **The features are the offline recipe's, not a lookalike.**  A model is
   about to be trained on ``results/algo_headroom/``'s feature tables and then
   evaluated with the tool's own columns.  If the two disagree by a base the
   score is silently wrong at run time and nothing fails.  So the sequence
   features are checked against an INDEPENDENT re-implementation of
   ``A3_scoring_model/code/seq_features.py`` -- written from the symmetric
   ``bedtools getfasta -s`` +/-120 window that analysis actually used, i.e. a
   different code path to the same numbers -- over randomised sequence.
3. **Strand.**  #96 fixed the internal-priming window on '-'.  Every window
   here is transcript-relative and every '-' feature is computed on the
   reverse complement; a plus-strand site and the mirrored minus-strand site
   on the reverse-complemented genome must produce identical features.
4. **The internal-priming veto is untouched.**  Collecting features must not
   change which rows are flagged, dropped or written -- ``ip_tool_flag`` is a
   covariate emitted IN ADDITION to the veto (VERIFY §1.4, §7.1).
"""
from __future__ import annotations

import hashlib
import random
from pathlib import Path

import pytest

from ema.config import variable_config
from ema.countmatrix import pas_features as pf
from ema.countmatrix.indexing import BarcodeIndex
from ema.countmatrix.pas_features import FeatureCollector
from ema.countmatrix.paswrite import (
    CALL_FEATURE_COLUMNS,
    SUPPORT_COLUMNS,
    support_columns,
    support_path_for,
    support_write,
)
from ema.countmatrix.peackcalling import peak_calling
from ema.countmatrix.peak_state import PeakCallingState
from ema.experimental.internal_priming import filter_internal_priming
from ema.strategies import get_strategy

FIXTURE = Path(__file__).parent / "fixtures" / "cellranger_pbmc_tiny.bam"


@pytest.fixture(autouse=True)
def _restore_variable_config():
    keys = ("seqlen", "cb_len", "barcode_tag", "ignore_chro",
            "read_geometry", "pas_features")
    saved = {k: getattr(variable_config, k) for k in keys}
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)


# ---------------------------------------------------------------------------
# The independent reference: A3's seq_features.py, re-derived
# ---------------------------------------------------------------------------
# A3 fetched a SYMMETRIC +/-120 nt window with `bedtools getfasta -s` (so the
# sequence arrives already reverse-complemented on '-') and indexed it with
# `rel(seq, lo, hi) == seq[lo+W : hi+W+1]`, r=0 being the cleavage base.  The
# tool fetches r in [-40, +30] and orients it itself.  Two different windows,
# two different indexings, one set of numbers -- that is the point.
_W = 120


def _ref_window(genome: str, cleavage: int, strand: str) -> str:
    """A3's `bedtools getfasta -s` window: [c-120, c+120] in transcript order."""
    lo, hi = cleavage - _W, cleavage + _W + 1
    assert lo >= 0 and hi <= len(genome), "widen the synthetic contig"
    seq = genome[lo:hi].upper()
    return seq if strand == "+" else pf.reverse_complement(seq)


def _ref_rel(seq: str, lo: int, hi: int) -> str:
    return seq[lo + _W: hi + _W + 1]


def _ref_longest_run(s: str, ch: str) -> int:
    best = cur = 0
    for c in s:
        if c == ch:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def _ref_features(seq: str) -> dict:
    """A3 ``seq_features.py``'s columns, restricted to the ones the tool emits."""
    hexwin = _ref_rel(seq, -40, -5)
    hits = [h for h in pf.HEX12 if h in hexwin]
    best_off = 0
    for h in hits:
        off = -40 + hexwin.rfind(h) + 5
        if best_off == 0 or off > best_off:
            best_off = off
    strong_off = 0
    for h in pf.HEX_STRONG:
        p = hexwin.rfind(h)
        if p >= 0:
            off = -40 + p + 5
            if strong_off == 0 or off > strong_off:
                strong_off = off
    k18 = _ref_rel(seq, 1, 18)
    d30 = _ref_rel(seq, 1, 30)
    a18 = k18.count("A")
    return {
        "seq_ok": 1,
        "a_count_d18": a18,
        "a_frac_d18": round(a18 / 18, 4),
        "a_run_d18": _ref_longest_run(k18, "A"),
        "a_frac_d30": round(d30.count("A") / 30, 4),
        "a_run_d30": _ref_longest_run(d30, "A"),
        "kin_ip_flag": int(a18 >= 12),
        "hex_strong": int(any(h in pf.HEX_STRONG for h in hits)),
        "hex_any12": int(bool(hits)),
        "hex_n_types": len(hits),
        "hex_best_off": best_off,
        "hex_strong_off": strong_off,
    }


def _tool_features(genome: str, cleavage: int, strand: str) -> dict:
    """The tool's path: feature_window -> orient_window -> sequence_features."""
    pas_pos = cleavage + 1 if strand == "+" else cleavage
    g0, g1 = pf.feature_window(pas_pos, strand)
    tseq, r0 = pf.orient_window(genome[g0:g1], pas_pos, strand, g0)
    return pf.sequence_features(tseq, r0)


def _random_genome(n: int, rng: random.Random, alpha: str = "ACGT") -> str:
    return "".join(rng.choice(alpha) for _ in range(n))


# ---------------------------------------------------------------------------
# 2. the features are the offline recipe's
# ---------------------------------------------------------------------------
def test_sequence_features_match_the_offline_recipe_on_random_sequence():
    """500 random sites x both strands, tool vs an independent re-derivation."""
    rng = random.Random(20260822)
    genome = _random_genome(4000, rng)
    checked = 0
    for _ in range(500):
        c = rng.randrange(_W + 1, len(genome) - _W - 1)
        for strand in ("+", "-"):
            got = _tool_features(genome, c, strand)
            want = _ref_features(_ref_window(genome, c, strand))
            assert got == want, (
                f"cleavage {c} strand {strand}: tool {got} != A3 recipe {want}"
            )
            checked += 1
    assert checked == 1000


def test_sequence_features_match_the_offline_recipe_on_an_A_rich_sequence():
    """Random ACGT rarely contains a canonical hexamer or a long A-run.  Bias
    the alphabet so the hexamer and downstream-A columns are actually
    exercised rather than agreeing at zero."""
    rng = random.Random(7)
    genome = _random_genome(6000, rng, alpha="AAAATTCG")
    seen_hex = seen_kin = 0
    for _ in range(400):
        c = rng.randrange(_W + 1, len(genome) - _W - 1)
        for strand in ("+", "-"):
            got = _tool_features(genome, c, strand)
            assert got == _ref_features(_ref_window(genome, c, strand))
            seen_hex += got["hex_any12"]
            seen_kin += got["kin_ip_flag"]
    assert seen_hex > 20, f"only {seen_hex} hexamer hits -- the test is vacuous"
    assert seen_kin > 5, f"only {seen_kin} Kinnex-IP hits -- the test is vacuous"


# ---------------------------------------------------------------------------
# 3. strand
# ---------------------------------------------------------------------------
def test_a_minus_site_sees_the_reverse_complement_of_what_a_plus_site_sees():
    rng = random.Random(11)
    genome = _random_genome(2000, rng, alpha="AAATCG")
    rc = pf.reverse_complement(genome)
    n = len(genome)
    for c in range(200, 1800, 37):
        plus = _tool_features(genome, c, "+")
        # base c of `genome` is base n-1-c of `rc`, and transcript orientation
        # on '-' of rc runs the same way as '+' on genome.
        minus = _tool_features(rc, n - 1 - c, "-")
        assert plus == minus, f"cleavage {c}: {plus} != {minus}"


def test_downstream_on_minus_means_lower_genomic_coordinates():
    """The A-stretch that primes oligo-dT on '-' sits genomically BELOW the
    cleavage base and reads as a T-run on the forward strand."""
    g = list("CG" * 400)
    g[300 - 12:300 - 6] = list("T" * 6)      # r +7..+12 on '-' at cleavage 300
    genome = "".join(g)
    minus = _tool_features(genome, 300, "-")
    assert minus["a_run_d18"] == 6 and minus["a_count_d18"] == 6
    # the same T-run is UPSTREAM for a '+' site at the same base: invisible
    plus = _tool_features(genome, 300, "+")
    assert plus["a_run_d18"] == 0 and plus["a_count_d18"] == 0


# ---------------------------------------------------------------------------
# window edges
# ---------------------------------------------------------------------------
def test_a_truncated_window_reports_seq_ok_zero_and_NA_everywhere():
    genome = "ACGT" * 20  # 80 nt: no site can carry a full [-40, +30] window
    feats = _tool_features(genome, 5, "+")
    assert feats["seq_ok"] == 0
    assert all(feats[c] == pf.NA for c in pf.SEQ_FEATURE_COLUMNS if c != "seq_ok")
    # ...and the last base of a long contig on '-', where the truncation hits
    # the transcript-oriented HEAD rather than the tail
    long_genome = "ACGT" * 100
    assert _tool_features(long_genome, len(long_genome) - 3, "-")["seq_ok"] == 0


def test_a_site_exactly_at_the_window_boundary_is_complete():
    genome = "ACGT" * 100
    # '+': needs bases c-40 .. c+30
    assert _tool_features(genome, 40, "+")["seq_ok"] == 1
    assert _tool_features(genome, 39, "+")["seq_ok"] == 0
    # '-': needs bases c-30 .. c+40
    assert _tool_features(genome, 30, "-")["seq_ok"] == 1
    assert _tool_features(genome, 29, "-")["seq_ok"] == 0


# ---------------------------------------------------------------------------
# ip_tool_* -- the covariate the verifier singled out
# ---------------------------------------------------------------------------
def test_ip_covariates_describe_the_string_the_veto_tested():
    """VERIFY §6(c): the tool's own inverted ``ip_tool_afrac`` alone scores
    AUC 0.7790 on truth-vs-decoy, beating A3's whole model (0.7655).  It has
    to be the VETO's window, not a lookalike."""
    from ema.experimental.internal_priming import check_internal_priming

    rng = random.Random(3)
    genome = _random_genome(2000, rng, alpha="AAATCG")
    for strand in ("+", "-"):
        for c in range(300, 1700, 53):
            pas_pos = c + 1 if strand == "+" else c
            for left, right in ((10, 30), (5, 50), (25, 25)):
                flag, tested = check_internal_priming(
                    genome, pas_pos, strand, left, right)
                afrac, arun = pf.ip_covariates(tested)
                assert abs(afrac - tested.count("A") / len(tested)) < 1e-12
                assert arun == _ref_longest_run(tested, "A")
                # ...and the flag is exactly the veto's own answer
                assert flag == (("A" * 6 in tested) or afrac >= 0.7)


# ---------------------------------------------------------------------------
# local candidate context
# ---------------------------------------------------------------------------
def test_local_context_counts_the_right_neighbours():
    rows = [
        # (chrom, strand, cleavage, molecules, pas_id)
        ("1", "+", 1000, 5, "a"),
        ("1", "+", 1050, 2, "b"),
        ("1", "+", 1400, 9, "c"),
        ("1", "+", 9000, 1, "d"),      # isolated
        ("1", "-", 1010, 7, "e"),      # other strand: never a neighbour of a-d
        ("2", "+", 1005, 3, "f"),      # other contig: likewise
    ]
    ctx = pf.local_context(rows)
    assert ctx["a"]["d_prev_cand"] == pf.NO_NEIGHBOUR
    assert ctx["a"]["d_next_cand"] == 50
    assert ctx["a"]["n_cand_100"] == 1        # b only (c is 400 away)
    assert ctx["a"]["n_cand_500"] == 2        # b and c
    assert ctx["a"]["mol_500_sum"] == 5 + 2 + 9
    assert ctx["a"]["is_local_mol_max"] == 0  # c has 9
    assert ctx["a"]["mol_frac_local"] == round(5 / 16, 4)
    assert ctx["c"]["is_local_mol_max"] == 1
    assert ctx["d"]["n_cand_500"] == 0
    assert ctx["d"]["mol_500_sum"] == 1
    assert ctx["d"]["is_local_mol_max"] == 1
    assert ctx["d"]["mol_frac_local"] == 1.0
    assert ctx["d"]["d_prev_cand"] == 9000 - 1400
    assert ctx["d"]["d_next_cand"] == pf.NO_NEIGHBOUR
    # strand and contig isolation
    assert ctx["e"]["n_cand_500"] == 0 and ctx["f"]["n_cand_500"] == 0


def test_local_context_never_divides_by_zero():
    ctx = pf.local_context([("1", "+", 10, 0, "z"), ("1", "+", 20, 0, "y")])
    assert ctx["z"]["mol_500_sum"] == 0
    assert ctx["z"]["mol_frac_local"] == 0.0
    assert ctx["z"]["is_local_mol_max"] == 1   # 0 >= max(0, 0)


# ---------------------------------------------------------------------------
# 1. v2 byte-identity, and the flag being load-bearing
# ---------------------------------------------------------------------------
def _run_caller(features: str, out_dir: Path) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    variable_config.seqlen = 91
    variable_config.cb_len = 16
    variable_config.barcode_tag = "CB"
    variable_config.ignore_chro = ["MT", "mt"]
    variable_config.read_geometry = "fixed"
    variable_config.pas_features = features
    index = BarcodeIndex()
    state = PeakCallingState(pasnumber=0)
    digests: dict[str, str] = {}
    for direction in (False, True):
        bed = out_dir / f"pas_{int(direction)}.bed"
        matrix = out_dir / f"matrix_{int(direction)}.txt"
        peak_calling(
            direction=direction, bedfilepath=str(bed), matrixpath=str(matrix),
            bamfile_dir=str(FIXTURE), index=index, state=state,
            sample_id="fixture", strategy=get_strategy("clip_seeded"),
        )
        for path in (bed, matrix, Path(support_path_for(bed))):
            digests[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return digests


@pytest.mark.skipif(not FIXTURE.exists(), reason="fixture BAM not present")
def test_the_features_change_the_sidecar_and_nothing_else(tmp_path):
    off = _run_caller("off", tmp_path / "off")
    on = _run_caller("on", tmp_path / "on")
    assert set(off) == set(on)
    for name in off:
        if name.endswith(".support.tsv"):
            assert off[name] != on[name], (
                f"{name}: --pas-features on wrote the same bytes as off -- "
                "either the flag is dead or the fixture has no PAS"
            )
        else:
            assert off[name] == on[name], (
                f"{name} changed: --pas-features must touch the SIDECAR ONLY. "
                "The BED and the count matrix are what every downstream "
                "benchmark reads."
            )


@pytest.mark.skipif(not FIXTURE.exists(), reason="fixture BAM not present")
def test_the_sidecar_is_extended_by_appending(tmp_path):
    """Every v2 column keeps its position AND its value."""
    _run_caller("off", tmp_path / "off")
    _run_caller("on", tmp_path / "on")
    for strand in (0, 1):
        v2 = (tmp_path / "off" / f"pas_{strand}.support.tsv").read_text().splitlines()
        pr = (tmp_path / "on" / f"pas_{strand}.support.tsv").read_text().splitlines()
        assert len(v2) == len(pr) > 1
        assert pr[0].split("\t")[:len(SUPPORT_COLUMNS)] == list(SUPPORT_COLUMNS)
        assert pr[0].split("\t")[len(SUPPORT_COLUMNS):] == list(CALL_FEATURE_COLUMNS)
        for a, b in zip(v2, pr):
            assert b.split("\t")[:len(SUPPORT_COLUMNS)] == a.split("\t")


@pytest.mark.skipif(not FIXTURE.exists(), reason="fixture BAM not present")
def test_clip_positions_and_span_describe_the_cluster(tmp_path):
    _run_caller("on", tmp_path / "on")
    rows = []
    for strand in (0, 1):
        lines = (tmp_path / "on" / f"pas_{strand}.support.tsv").read_text().splitlines()
        cols = lines[0].split("\t")
        rows.extend(dict(zip(cols, ln.split("\t"))) for ln in lines[1:])
    assert rows, "fixture produced no PAS"
    for r in rows:
        n, span, reads = (int(r["clip_positions"]), int(r["clip_span"]),
                          int(r["clip_reads"]))
        assert n >= 0 and span >= 0
        # a cluster cannot have more distinct positions than clip reads, and a
        # single position has zero span
        assert n <= reads
        assert (span == 0) == (n <= 1)
        # single-linkage at --polya-seed-window 25 bounds the span
        assert span <= 25 * max(0, n - 1)
    assert any(int(r["clip_positions"]) > 1 for r in rows), (
        "no multi-position cluster in the fixture -- the columns are vacuous"
    )


def test_the_writer_and_the_header_cannot_disagree(tmp_path):
    """``open_support`` and ``support_write`` take the same flag; a row must
    always have exactly as many fields as the header."""
    row = {c: 1 for c in SUPPORT_COLUMNS[1:]}
    row.update({c: 2 for c in CALL_FEATURE_COLUMNS})
    for features in (False, True):
        path = tmp_path / f"f{int(features)}.tsv"
        with open(path, "w") as fh:
            fh.write("\t".join(support_columns(features)) + "\n")
            support_write(fh, 1, row, features)
        header, data = path.read_text().splitlines()
        assert len(header.split("\t")) == len(data.split("\t"))
        assert len(header.split("\t")) == 7 + (2 if features else 0)


def test_a_missing_call_feature_is_NA_not_a_zero(tmp_path):
    """A support dict without the feature keys must not silently write 0 --
    a zero clip span is a real measurement and an absent one is not."""
    path = tmp_path / "s.tsv"
    with open(path, "w") as fh:
        support_write(fh, 7, {c: 3 for c in SUPPORT_COLUMNS[1:]}, True)
    assert path.read_text().rstrip("\n").split("\t")[-2:] == ["NA", "NA"]


# ---------------------------------------------------------------------------
# the seam writer
# ---------------------------------------------------------------------------
def test_append_columns_is_append_only(tmp_path):
    src = tmp_path / "pas_support.tsv"
    src.write_text(
        "pas_id\tclip_reads\ttier\n"
        "1\t10\t1\n"
        "2\t0\t2\n"
        "3\t4\t1\n"
    )
    original = src.read_text().splitlines()
    feats = {"1": {"seq_ok": 1, "hex_strong": 1}, "3": {"seq_ok": 0}}
    n = pf.append_columns(src, feats, columns=("seq_ok", "hex_strong"))
    assert n == 3
    out = src.read_text().splitlines()
    assert out[0] == "pas_id\tclip_reads\ttier\tseq_ok\thex_strong"
    for before, after in zip(original, out):
        assert after.split("\t")[:len(before.split("\t"))] == before.split("\t")
    assert out[1].split("\t")[3:] == ["1", "1"]
    # a pas_id the seam never saw is NA, not a fabricated 0
    assert out[2].split("\t")[3:] == ["NA", "NA"]
    assert out[3].split("\t")[3:] == ["0", "NA"]


def test_write_features_tsv_is_sorted_by_numeric_pas_id(tmp_path):
    out = tmp_path / "pas_features.tsv"
    pf.write_features_tsv(out, {"10": {"seq_ok": 1}, "2": {"seq_ok": 0}},
                          columns=("seq_ok",))
    assert out.read_text().splitlines() == ["pas_id\tseq_ok", "2\t0", "10\t1"]


def test_collector_without_a_genome_says_NA_not_zero():
    coll = FeatureCollector()
    coll.add("1", "1", "+", 100, 5, None)
    feats = coll.finish()["1"]
    # every SEQUENCE column is NA, including seq_ok: "no genome was supplied"
    # is a different fact from "this site is at a contig edge" (seq_ok == 0)
    for c in pf.SEQ_FEATURE_COLUMNS:
        assert feats[c] == pf.NA
    # ...while the context columns are real -- they only need the BED
    assert feats["n_cand_100"] == 0
    assert feats["mol_500_sum"] == 5


def test_collect_from_bed_reads_the_cleavage_base_per_strand(tmp_path):
    bed = tmp_path / "pas.bed"
    bed.write_text(
        "1\t100\t101\tp1\t5\t+\n"
        "1\t200\t201\tm1\t3\t-\n"
        "1\t300\t400\tp2\t0\t+\n"     # a tier-2 interval, not a 1-bp point
    )
    coll = FeatureCollector()
    assert pf.collect_from_bed(bed, coll) == 3
    by_id = {r[4]: r for r in coll.rows}
    assert by_id["p1"][2] == 100      # end - 1
    assert by_id["m1"][2] == 200      # start
    assert by_id["p2"][2] == 399      # end - 1, wide interval
    assert by_id["p1"][3] == 5 and by_id["p2"][3] == 0


# ---------------------------------------------------------------------------
# the branch default lives in exactly one place
# ---------------------------------------------------------------------------
def test_the_branch_default_is_single_valued():
    """Same lesson as ``--read-geometry``: `ema run` reads RunConfig and a
    direct library call reads variable_config.  If they disagree the two entry
    points write different sidecars and neither looks wrong."""
    from ema.cli.config_schema import RunConfig

    assert RunConfig().pas_features == variable_config.pas_features


def test_the_flag_is_validated_not_guessed(tmp_path):
    from ema.countmatrix.paswrite import PAS_FEATURE_MODES, V2_PAS_FEATURES

    assert V2_PAS_FEATURES in PAS_FEATURE_MODES
    variable_config.seqlen = 91
    variable_config.cb_len = 16
    variable_config.barcode_tag = "CB"
    variable_config.pas_features = "yes-please"
    with pytest.raises(ValueError, match="pas_features must be one of"):
        peak_calling(
            direction=False, bedfilepath=str(tmp_path / "x.bed"),
            matrixpath=str(tmp_path / "x.mtx"), bamfile_dir=str(FIXTURE),
            index=BarcodeIndex(), state=PeakCallingState(pasnumber=0),
            sample_id="fixture", strategy=get_strategy("clip_seeded"),
        )


# ---------------------------------------------------------------------------
# 4. the veto, and the seam end to end
# ---------------------------------------------------------------------------
def _tiny_genome_and_bed(tmp_path: Path) -> tuple[Path, Path]:
    """A C/G background with implanted A/T runs and hexamers, plus a BED that
    exercises both strands, a contig edge and a missing contig."""
    g = list(("CG" * 400)[:800])
    for pos, motif in ((410, "A" * 8), (275, "T" * 8), (555, "AATAAA"),
                       (600, "A" * 14)):
        g[pos:pos + len(motif)] = list(motif)
    seq = "".join(g)
    fa = tmp_path / "genome.fa"
    with open(fa, "w") as fh:
        for name, s in (("chrT", seq), ("chrE", "ACGT" * 5)):
            fh.write(f">{name}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    bed = tmp_path / "pas.bed"
    bed.write_text(
        "chrT\t399\t400\tP1\t5\t+\n"
        "chrT\t299\t300\tM1\t3\t-\n"
        "chrT\t549\t550\tP2\t0\t+\n"
        "chrT\t594\t595\tP3\t9\t+\n"
        "chrE\t5\t6\tEDGE\t1\t+\n"
        "chrMissing\t100\t101\tGONE\t1\t+\n"
    )
    return fa, bed


def test_collecting_features_does_not_change_the_veto(tmp_path):
    pytest.importorskip("pyfaidx")
    genome, bed = _tiny_genome_and_bed(tmp_path)
    out_a, out_b = tmp_path / "a.bed", tmp_path / "b.bed"
    plain = filter_internal_priming(str(bed), str(genome), str(out_a), mode="filter")
    coll = FeatureCollector()
    withf = filter_internal_priming(str(bed), str(genome), str(out_b),
                                    mode="filter", features=coll)
    assert plain["flags"] == withf["flags"]
    assert out_a.read_bytes() == out_b.read_bytes()
    assert withf["total"] == plain["total"] and withf["filtered"] == plain["filtered"]
    # every scanned row got a feature dict, flagged or not
    assert set(coll.features) == set(plain["flags"])
    for pas_id, flagged in plain["flags"].items():
        emitted = coll.features[pas_id]["ip_tool_flag"]
        if emitted != pf.NA:
            assert emitted == int(flagged)


def test_a_features_only_scan_writes_no_bed(tmp_path):
    """``output_path=None``: the pass --pas-features makes when --ip-filter is
    off.  It must read the genome once and touch nothing else."""
    pytest.importorskip("pyfaidx")
    genome, bed = _tiny_genome_and_bed(tmp_path)
    coll = FeatureCollector()
    stats = filter_internal_priming(str(bed), str(genome), None,
                                    mode="annotate", features=coll)
    assert stats["filtered"] == 0
    # no BED was written -- pyfaidx's own .fai index is the only new file
    assert sorted(p.name for p in tmp_path.iterdir()
                  if p.suffix == ".bed") == ["pas.bed"]
    assert len(coll.features) == stats["total"] > 0
    with pytest.raises(ValueError, match="output_path=None"):
        filter_internal_priming(str(bed), str(genome), None, mode="filter")


def test_end_to_end_seam_features_are_real(tmp_path):
    pytest.importorskip("pyfaidx")
    genome, bed = _tiny_genome_and_bed(tmp_path)
    coll = FeatureCollector()
    filter_internal_priming(str(bed), str(genome), None, mode="annotate",
                            features=coll)
    feats = coll.finish()
    # P1: cleavage 399 on '+', an 8-A run implanted at 410 => r +11..+18
    p1 = feats["P1"]
    assert p1["seq_ok"] == 1
    assert p1["a_run_d18"] == 8 and p1["a_count_d18"] == 8
    assert p1["ip_tool_flag"] == 1          # 8 >= --ip-a-stretch 6
    assert p1["ip_tool_arun"] == 8
    # M1: cleavage 299 on '-', an 8-T run at genomic [275, 283) => transcript
    # r = 299 - x, i.e. r +17..+24 -- straddling the +1..+18 boundary, which
    # is exactly the kind of off-by-one an untested strand convention hides.
    m1 = feats["M1"]
    assert m1["a_run_d30"] == 8
    assert m1["a_run_d18"] == 2 and m1["a_count_d18"] == 2
    # P2: AATAAA implanted at 555, cleavage 549 => last base at r +11, NOT in
    # the -40..-5 hexamer window, so it must NOT be counted
    assert feats["P2"]["hex_strong"] == 0
    # P3: cleavage 594, the same AATAAA at 555..560 => last base r -34
    assert feats["P3"]["hex_strong"] == 1
    assert feats["P3"]["hex_strong_off"] == -34
    assert feats["P3"]["hex_best_off"] == -34
    # a contig edge and a missing contig degrade rather than crash
    assert feats["EDGE"]["seq_ok"] == 0
    assert feats["GONE"]["seq_ok"] == 0
    assert feats["GONE"]["ip_tool_afrac"] == pf.NA
    # ...and context is real for every row, even the ones with no sequence
    for row in feats.values():
        assert row["n_cand_100"] != pf.NA
