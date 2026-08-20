"""End-to-end poly(A) evidence tests on a synthetic BAM (Stage 1).

Builds a small BAM in-process (no fixture file, so `.gitignore:15 *.bam`
cannot silently skip these — see plan Stage 0f) containing two constructed
loci on the forward strand:

  * ``CLIPPED`` locus at ~2000: a coverage peak whose reads carry genuine
    terminal poly(A) soft clips at one cleavage site.
  * ``BARE`` locus at ~12000: an equally deep coverage peak with no clip
    evidence anywhere near it.

The tests then prove the three Phase-1/Phase-2 contracts:

1. BED column 5 carries per-PAS clip-read support (it was hardcoded 0).
2. ``--peak-strategy clip_seeded`` emits BOTH tiers — a clip-cluster PAS at
   the cleavage site (score >= 1) and a coverage-only PAS (score == 0) —
   and stays BED6-parseable.
3. ``--polya-mode filter`` drops ONLY the coverage-only tier.
4. Two identical runs are byte-identical (determinism).
"""
from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import pysam
import pytest

CHROM = "1"
CHROM_LEN = 200_000
SEQ_LEN = 91          # matches the PBMC benchmark's --seq-len
N_RAMP = 60           # background reads building the coverage ramp
N_SUMMIT = 40         # reads whose 3' end defines the cleavage site

CLIP_SITE = 2_400     # cleavage coordinate the clipped locus should report
BARE_SITE = 12_400    # coverage-only locus, no clip evidence


def _barcode(i: int) -> str:
    """Distinct valid 16-mer barcodes."""
    alphabet = "ACGT"
    out = [alphabet[(i >> (2 * k)) & 3] for k in range(8)]
    return ("".join(out) + "ACGTACGT")[:16]


def _locus_specs(site: int, clipped: bool, n_summit: int = N_SUMMIT):
    """``(start, aligned_len, clip_len)`` triples for one synthetic locus.

    The ramp reads spread their starts over ~240 bp so the caller's
    ``start1 > i_end`` slice fires inside the locus and ``peak_list`` is
    actually built; every aligned span stays <= --seq-len so ``read_check``
    does not drop them.  The summit reads all terminate at *site*, which is
    what makes *site* the true cleavage coordinate.
    """
    specs = [(site - 350 + i * 4, 80, 0) for i in range(N_RAMP)]
    for i in range(n_summit):
        start = site - 80 + (i % 10)
        specs.append((start, site - start + 1, 12 if clipped else 0))
    return specs


def _write_bam(path: Path, *, clipped_summit: int = N_SUMMIT,
               bare_summit: int = N_SUMMIT) -> Path:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": CHROM_LEN}],
        "RG": [{"ID": "testsample", "SM": "testsample"}],
    }
    records = []
    n = 0
    loci = (
        (CLIP_SITE, True, clipped_summit),
        (BARE_SITE, False, bare_summit),
    )
    for site, clipped, n_summit in loci:
        for i, (start, aligned, clip) in enumerate(
            _locus_specs(site, clipped, n_summit)
        ):
            a = pysam.AlignedSegment()
            a.query_name = f"r{n}"
            n += 1
            a.query_sequence = "C" * aligned + "A" * clip
            a.flag = 0                     # forward, primary, mapped
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = [(0, aligned)] + ([(4, clip)] if clip else [])
            a.query_qualities = pysam.qualitystring_to_array(
                "I" * len(a.query_sequence)
            )
            a.set_tag("CB", _barcode(i))
            a.set_tag("UB", f"UMI{n:06d}")
            a.set_tag("RG", "testsample")
            records.append(a)

    records.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for rec in records:
            out.write(rec)
    pysam.index(str(path))
    return path


@pytest.fixture(scope="module")
def synthetic_bam(tmp_path_factory) -> Path:
    d = tmp_path_factory.mktemp("polya_bam")
    return _write_bam(d / "synthetic.bam")


@pytest.fixture(scope="module")
def peak_calling():
    """Point ``variable_config`` at the synthetic BAM's parameters.

    ``read_check`` looks these up at CALL time (not at import time — see its
    docstring), so no module reload is needed.  The values are restored on
    teardown: reloading or leaking config here silently changes every test
    that runs afterwards.
    """
    from ema.config import variable_config
    from ema.countmatrix.peackcalling import peak_calling as _pc

    saved = {
        k: getattr(variable_config, k)
        for k in ("seqlen", "cb_len", "barcode_tag", "ignore_chro",
                  "default_threshold", "merge_len")
    }
    variable_config.seqlen = SEQ_LEN
    variable_config.cb_len = 16
    variable_config.barcode_tag = "CB"
    variable_config.ignore_chro = []
    variable_config.default_threshold = 5
    variable_config.merge_len = 100
    try:
        yield _pc
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)


def _run(peak_calling, bam: Path, out_dir: Path, tag: str, **kwargs):
    bed = out_dir / f"{tag}.bed"
    mtx = out_dir / f"{tag}.mtx"
    from ema.countmatrix.peak import Peak
    from ema.strategies import get_strategy
    strategy_name = kwargs.pop("strategy_name", "lambda_gradient")
    # Each _run is an independent "run": reset the class-level pasnumber
    # counter the same way ema.main does between datasets, so PAS ids start
    # at 1 and two runs are comparable byte-for-byte.
    Peak.reset_pasnumber()
    with patch.object(sys, "argv", ["ema"]):
        peak_calling(
            False,
            bedfilepath=str(bed),
            matrixpath=str(mtx),
            bamfile_dir=str(bam),
            strategy=get_strategy(strategy_name),
            **kwargs,
        )
    return bed, mtx


def _rows(bed: Path):
    rows = []
    for line in bed.read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        assert len(f) == 6, f"pasbed must stay BED6-parseable, got {len(f)}: {f}"
        rows.append({
            "chrom": f[0], "start": int(f[1]), "end": int(f[2]),
            "name": f[3], "score": int(f[4]), "strand": f[5],
        })
    return rows


def _near(rows, site, window=150):
    return [r for r in rows if abs(r["end"] - 1 - site) <= window
            or abs(r["start"] - site) <= window]


# ---------------------------------------------------------------------------
# Phase 1: BED column 5 carries clip support
# ---------------------------------------------------------------------------

def test_annotate_mode_writes_clip_support_into_bed_score(synthetic_bam, peak_calling, tmp_path):
    bed, _ = _run(peak_calling, synthetic_bam, tmp_path, "annot")
    rows = _rows(bed)
    assert rows, "expected at least one PAS from the synthetic BAM"

    clipped = _near(rows, CLIP_SITE)
    bare = _near(rows, BARE_SITE)
    assert clipped, "no PAS called near the clipped locus"
    assert bare, "no PAS called near the bare locus"

    assert max(r["score"] for r in clipped) >= 1, (
        "PAS at the clipped locus must carry non-zero clip support in BED col 5"
    )
    assert all(r["score"] == 0 for r in bare), (
        "PAS at the clip-free locus must have score 0"
    )


def test_polya_evidence_off_keeps_score_zero_everywhere(synthetic_bam, peak_calling, tmp_path):
    """--polya-evidence off must restore the pre-Stage-1 output."""
    bed, _ = _run(peak_calling, synthetic_bam, tmp_path, "off", polya_enabled=False)
    rows = _rows(bed)
    assert rows
    assert all(r["score"] == 0 for r in rows)


def test_polya_off_and_on_agree_on_coordinates(synthetic_bam, peak_calling, tmp_path):
    """Annotate mode must not move a single coordinate — only column 5."""
    bed_on, _ = _run(peak_calling, synthetic_bam, tmp_path, "cmp_on")
    bed_off, _ = _run(peak_calling, synthetic_bam, tmp_path, "cmp_off",
                      polya_enabled=False)
    on = [(r["chrom"], r["start"], r["end"], r["strand"]) for r in _rows(bed_on)]
    off = [(r["chrom"], r["start"], r["end"], r["strand"]) for r in _rows(bed_off)]
    assert on == off


# ---------------------------------------------------------------------------
# Phase 2: two-tier clip_seeded output
# ---------------------------------------------------------------------------

def test_clip_seeded_emits_two_tiers(synthetic_bam, peak_calling, tmp_path):
    bed, _ = _run(peak_calling, synthetic_bam, tmp_path, "seeded",
                  strategy_name="clip_seeded")
    rows = _rows(bed)
    assert rows, "clip_seeded produced no PAS"

    tier1 = [r for r in rows if r["score"] >= 1]
    tier2 = [r for r in rows if r["score"] == 0]
    assert tier1, "expected a clip-supported (tier 1) PAS"
    assert tier2, "expected a coverage-only (tier 2) PAS"

    # tier 1 must sit at the cleavage site, as a 1-bp interval
    seeds = [r for r in tier1 if r["start"] == CLIP_SITE]
    assert seeds, (
        f"expected a tier-1 PAS at the cleavage site {CLIP_SITE}, "
        f"got {[r['start'] for r in tier1]}"
    )
    assert seeds[0]["end"] == CLIP_SITE + 1, "clip cluster must be a 1-bp BED interval"
    assert seeds[0]["score"] == N_SUMMIT, (
        f"tier-1 score must be the cluster's clip-read count ({N_SUMMIT})"
    )

    # tier 2 must be the bare locus
    assert _near(tier2, BARE_SITE), "coverage-only tier must contain the bare locus"

    # output must be coordinate-sorted within the chromosome
    starts = [r["start"] for r in rows]
    assert starts == sorted(starts), "clip_seeded output must be coordinate-sorted"


def test_clip_seeded_beats_coverage_only_on_cleavage_accuracy(synthetic_bam, peak_calling, tmp_path):
    """The point of seeding: the emitted 3' base should BE the cleavage site,
    not the coverage summit's outer edge."""
    seeded, _ = _run(peak_calling, synthetic_bam, tmp_path, "acc_seeded",
                     strategy_name="clip_seeded")
    cov, _ = _run(peak_calling, synthetic_bam, tmp_path, "acc_cov")

    def best_offset(rows):
        cands = [abs(r["end"] - 1 - CLIP_SITE) for r in rows]
        return min(cands) if cands else 10**9

    assert best_offset(_rows(seeded)) <= best_offset(_rows(cov))
    assert best_offset(_rows(seeded)) == 0


def test_clip_seeded_requires_polya_evidence(synthetic_bam, peak_calling, tmp_path):
    """Contradictory flags must fail loudly, not silently drop to coverage."""
    with pytest.raises(ValueError, match="clip_seeded"):
        _run(peak_calling, synthetic_bam, tmp_path, "bad",
             strategy_name="clip_seeded", polya_enabled=False)


def test_clip_seeded_matrix_rows_match_bed_rows(synthetic_bam, peak_calling, tmp_path):
    """Every emitted PAS id must be a matrix row id and vice versa — the
    pasnumber/matrix pairing is what Stage 0a/0g showed is easy to break."""
    bed, mtx = _run(peak_calling, synthetic_bam, tmp_path, "keys",
                    strategy_name="clip_seeded")
    bed_ids = {r["name"] for r in _rows(bed)}
    mtx_ids = {
        line.split()[0] for line in mtx.read_text().splitlines() if line.strip()
    }
    assert mtx_ids, "clip_seeded wrote no matrix rows"
    assert mtx_ids <= bed_ids, f"matrix rows not in BED: {mtx_ids - bed_ids}"


def test_polya_min_reads_gates_tier1(synthetic_bam, peak_calling, tmp_path):
    """A threshold above the locus's molecule count must remove tier 1 and
    leave the (now unsuppressed) coverage candidate in its place."""
    bed, _ = _run(peak_calling, synthetic_bam, tmp_path, "gate",
                  strategy_name="clip_seeded", polya_min_reads=N_SUMMIT + 10)
    rows = _rows(bed)
    assert not [r for r in rows if r["start"] == CLIP_SITE and r["end"] == CLIP_SITE + 1], (
        "cluster below --polya-min-reads must not be emitted as tier 1"
    )
    assert _near(rows, CLIP_SITE), (
        "with tier 1 gated out, the coverage candidate must survive as tier 2"
    )


# ---------------------------------------------------------------------------
# --polya-mode filter drops only the coverage-only tier
# ---------------------------------------------------------------------------

def test_polya_mode_filter_drops_only_coverage_only_peaks(synthetic_bam, peak_calling, tmp_path):
    """The gate is a pure BED-score filter at the _apply_pas_filters seam, so
    exercise it on the real BED the caller produced."""
    bed, _ = _run(peak_calling, synthetic_bam, tmp_path, "filt",
                  strategy_name="clip_seeded")
    before = _rows(bed)
    kept_expected = [r for r in before if r["score"] >= 1]
    dropped_expected = [r for r in before if r["score"] == 0]
    assert kept_expected and dropped_expected, "need both tiers to test the gate"

    import types
    from unittest.mock import patch as _patch

    import ema.main as main_mod

    posbed = tmp_path / "posbed.bed"
    negbed = tmp_path / "negbed.bed"
    posbed.write_text(bed.read_text())
    negbed.write_text("")

    fake_args = types.SimpleNamespace(
        polya_evidence="on", polya_mode="filter",
        ip_filter=False, annot_filter=False,
    )

    mgr = types.SimpleNamespace(
        path=lambda *a: str(tmp_path / "polya_gate_stats.json")
    )
    dircfg = types.SimpleNamespace(posbed=posbed, negbed=negbed)

    with _patch.object(main_mod, "args", fake_args), \
         _patch.object(main_mod, "directory_config", dircfg):
        stats = main_mod._apply_polya_gate(mgr)

    assert stats is not None
    after = _rows(posbed)
    assert all(r["score"] >= 1 for r in after), "coverage-only PAS must be dropped"
    assert {r["name"] for r in after} == {r["name"] for r in kept_expected}
    assert stats["pos"]["dropped"] == len(dropped_expected)


def test_polya_mode_require_raises_when_no_support(tmp_path):
    """require must fail the run rather than ship an unsupported call set."""
    import types
    from unittest.mock import patch as _patch

    import ema.main as main_mod

    posbed = tmp_path / "posbed.bed"
    negbed = tmp_path / "negbed.bed"
    posbed.write_text("1\t100\t200\t1\t0\t+\n1\t300\t400\t2\t0\t+\n")
    negbed.write_text("")

    fake_args = types.SimpleNamespace(
        polya_evidence="on", polya_mode="require",
        ip_filter=False, annot_filter=False,
    )

    mgr = types.SimpleNamespace(path=lambda *a: str(tmp_path / "s.json"))
    dircfg = types.SimpleNamespace(posbed=posbed, negbed=negbed)

    with _patch.object(main_mod, "args", fake_args), \
         _patch.object(main_mod, "directory_config", dircfg):
        with pytest.raises(RuntimeError, match="require"):
            main_mod._apply_polya_gate(mgr)


def test_polya_gate_is_noop_in_annotate_mode(tmp_path):
    import types
    from unittest.mock import patch as _patch

    import ema.main as main_mod

    posbed = tmp_path / "posbed.bed"
    negbed = tmp_path / "negbed.bed"
    original = "1\t100\t200\t1\t0\t+\n"
    posbed.write_text(original)
    negbed.write_text("")

    fake_args = types.SimpleNamespace(
        polya_evidence="on", polya_mode="annotate",
        ip_filter=False, annot_filter=False,
    )

    mgr = types.SimpleNamespace(path=lambda *a: str(tmp_path / "s.json"))
    dircfg = types.SimpleNamespace(posbed=posbed, negbed=negbed)

    with _patch.object(main_mod, "args", fake_args), \
         _patch.object(main_mod, "directory_config", dircfg):
        assert main_mod._apply_polya_gate(mgr) is None
    assert posbed.read_text() == original


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("strategy_name", ["lambda_gradient", "clip_seeded"])
def test_two_runs_are_byte_identical(synthetic_bam, peak_calling, tmp_path, strategy_name):
    from ema.countmatrix.indexing import reset_index

    reset_index()
    bed_a, mtx_a = _run(peak_calling, synthetic_bam, tmp_path,
                        f"det_a_{strategy_name}", strategy_name=strategy_name)
    reset_index()
    bed_b, mtx_b = _run(peak_calling, synthetic_bam, tmp_path,
                        f"det_b_{strategy_name}", strategy_name=strategy_name)

    assert bed_a.read_bytes() == bed_b.read_bytes(), "BED output is not deterministic"
    assert mtx_a.read_bytes() == mtx_b.read_bytes(), "MTX output is not deterministic"


def test_clip_rate_warning_fires_on_clip_free_bam(tmp_path, caplog):
    """R2 mitigation: a BAM whose poly(A) evidence was destroyed upstream must
    produce a loud startup warning, not a silent unsupported call set."""
    import logging

    from ema.countmatrix.polya import check_clip_rate, _clip_rate_warned

    bam = _write_bam(tmp_path / "noclip.bam", clipped_summit=0)
    _clip_rate_warned.pop(str(bam), None)
    with caplog.at_level(logging.WARNING, logger="ema.countmatrix.polya"):
        rate = check_clip_rate(str(bam))
    assert rate == 0.0
    assert any("LOW POLY(A) CLIP RATE" in r.message for r in caplog.records)


def test_clip_rate_no_warning_on_clip_rich_bam(tmp_path, caplog):
    import logging

    from ema.countmatrix.polya import check_clip_rate, _clip_rate_warned

    bam = _write_bam(tmp_path / "clippy.bam", clipped_summit=N_SUMMIT)
    _clip_rate_warned.pop(str(bam), None)
    with caplog.at_level(logging.WARNING, logger="ema.countmatrix.polya"):
        rate = check_clip_rate(str(bam))
    assert rate > 0.003
    assert not any("LOW POLY(A) CLIP RATE" in r.message for r in caplog.records)
