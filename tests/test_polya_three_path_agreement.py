"""The three peak-calling paths must agree on poly(A) evidence (plan risk R5).

``read_check``'s 5-tuple is shared by the monolithic loop, the 3-stage
pipeline and the tile workers, and the benchmark exercises only the
monolithic path — so a divergence in the other two would never be caught by
re-running the benchmark (see the Stage-0 stop signal in
``manuscript/10_caller_fix_plan.md``).  Stage 1 adds a per-read measurement
to all three; these tests pin that they produce the same BED.

The pipeline path spawns three subprocesses and the tile path a process
pool, so both are slow — but they are the only way to prove the mirrors are
real rather than aspirational.
"""
from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from tests.test_polya_two_tier_integration import (  # noqa: E402
    BARE_SITE,
    CLIP_SITE,
    SEQ_LEN,
    _rows,
    _write_bam,
)


@pytest.fixture(scope="module")
def bam(tmp_path_factory) -> Path:
    d = tmp_path_factory.mktemp("polya_paths")
    return _write_bam(d / "synthetic.bam")


# The pipeline and tile paths use multiprocessing 'spawn', and spawn
# propagates the parent's sys.argv into every child — where ema.config's
# lazy args shim re-parses it (ema/cli/__init__.py::cli honours the legacy
# --sequenceLen / --CellBarcodeLen / --BarcodeTag flags).  So these flags
# must stay on sys.argv for the whole call, not only across the import.
_ARGV = [
    "ema",
    "--sequenceLen", str(SEQ_LEN),
    "--CellBarcodeLen", "16",
    "--BarcodeTag", "CB",
]


@pytest.fixture(autouse=True)
def _slice_config():
    """Set the synthetic BAM's read parameters, and restore them afterwards.

    ``read_check`` reads ``variable_config`` at call time, so setting the
    values is enough — and restoring them keeps this module from silently
    reconfiguring every test that runs after it.
    """
    from ema.config import variable_config

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
        yield
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)


def _call(bam: Path, out: Path, tag: str, strategy_name: str, **kwargs):
    from ema.countmatrix.indexing import reset_index
    from ema.countmatrix.peackcalling import peak_calling
    from ema.countmatrix.peak import Peak
    from ema.strategies import get_strategy

    reset_index()
    Peak.reset_pasnumber()
    bed = out / f"{tag}.bed"
    mtx = out / f"{tag}.mtx"
    with patch.object(sys, "argv", _ARGV):
        peak_calling(
            False,
            bedfilepath=str(bed),
            matrixpath=str(mtx),
            bamfile_dir=str(bam),
            strategy=get_strategy(strategy_name),
            **kwargs,
        )
    return bed


def _coords_and_scores(bed: Path):
    return [(r["chrom"], r["start"], r["end"], r["score"], r["strand"])
            for r in _rows(bed)]


@pytest.mark.parametrize("strategy_name", ["lambda_gradient", "clip_seeded"])
def test_pipeline_path_matches_monolithic(bam, tmp_path, strategy_name):
    mono = _call(bam, tmp_path, f"mono_{strategy_name}", strategy_name)
    pipe = _call(bam, tmp_path, f"pipe_{strategy_name}", strategy_name,
                 use_pipeline=True)
    assert _coords_and_scores(pipe) == _coords_and_scores(mono), (
        "--pipeline diverged from the monolithic path on poly(A) evidence"
    )


@pytest.mark.parametrize("strategy_name", ["lambda_gradient", "clip_seeded"])
def test_tile_path_matches_monolithic(bam, tmp_path, strategy_name):
    mono = _call(bam, tmp_path, f"mono_t_{strategy_name}", strategy_name)
    # tile_size below the chromosome length forces the real pool path
    # (2 tiles) rather than the single-tile short-circuit.
    tiled = _call(bam, tmp_path, f"tile_{strategy_name}", strategy_name,
                  use_tiles=True, tile_size=8_000, tile_overlap=2_000,
                  n_workers=2)
    assert _coords_and_scores(tiled) == _coords_and_scores(mono), (
        "--tiles diverged from the monolithic path on poly(A) evidence"
    )


def test_all_three_paths_find_the_cleavage_site_when_seeded(bam, tmp_path):
    """The measurement that matters: every path must place a clip-supported
    PAS exactly on the cleavage site, and keep the bare locus as tier 2."""
    beds = {
        "monolithic": _call(bam, tmp_path, "all_mono", "clip_seeded"),
        "pipeline": _call(bam, tmp_path, "all_pipe", "clip_seeded",
                          use_pipeline=True),
        "tiles": _call(bam, tmp_path, "all_tile", "clip_seeded",
                       use_tiles=True, tile_size=8_000, tile_overlap=2_000,
                       n_workers=2),
    }
    for label, bed in beds.items():
        rows = _rows(bed)
        seeds = [r for r in rows
                 if r["start"] == CLIP_SITE and r["end"] == CLIP_SITE + 1]
        assert seeds, f"{label}: no tier-1 PAS at the cleavage site"
        assert seeds[0]["score"] >= 1, f"{label}: tier-1 PAS has no clip support"
        bare = [r for r in rows
                if r["score"] == 0 and r["start"] <= BARE_SITE <= r["end"] + 150]
        assert bare, f"{label}: coverage-only tier missing at the bare locus"
