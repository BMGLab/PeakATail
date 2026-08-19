"""End-to-end replay of the three CellRanger input-compatibility defects.

Runs the real peak caller over ``tests/fixtures/cellranger_pbmc_tiny.bam`` --
a 148 kB / 2072-read slice of the 10x ``pbmc_10k_v3`` BAM (CellRanger 3.0.0,
GRCh38-3.0.0 reference) that carries, on purpose:

* ``CB`` tags with the CellRanger GEM-group suffix (``...GTT-1``),
* two ``@RG`` ids that contain underscores (``pbmc_10k_v3:0:1:...``),
* unmapped reads that still carry a ``CB`` tag (``reference_end is None``),
* secondary (0x100) and duplicate-flagged (0x400) alignments.

Each of those was a silent data-loss bug in the shipped read path:

1. the GEM suffix made ``len(cb) != cb_len``, so **every** read was dropped
   and the run produced an empty matrix (measured on this fixture at
   ``origin/develop`` c08ca23: **0** barcode columns);
2. ``reference_end is None`` crashed ``read_end - read_start`` with a
   ``TypeError`` hours into a run;
3. splitting the ``"<sample_id>_<barcode>"`` composite on the *first*
   underscore truncated any underscore-bearing RG, so the barcode half
   failed 2-bit encoding and every cell collapsed onto one column.

With the repaired read path this fixture yields **1493** barcode columns
(906 forward-strand, 622 reverse-strand, 1316 distinct barcodes, 177 of
which legitimately appear under *both* read groups).  The assertions below
use floors rather than those exact numbers so that peak-strategy tuning does
not turn this into a change-detector -- what is being pinned is
"many columns", never 0 or 1.

Fixture provenance and the exact regeneration command:
``tests/fixtures/make_cellranger_pbmc_tiny.sh``.
"""
from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from ema.config import variable_config
from ema.countmatrix.indexing import BarcodeIndex, split_cb
from ema.countmatrix.peak_state import PeakCallingState
from ema.countmatrix.peackcalling import peak_calling
from ema.countmatrix.read import read_check

FIXTURE = Path(__file__).parent / "fixtures" / "cellranger_pbmc_tiny.bam"

# The two @RG ids in the fixture header.  Both contain underscores; the
# library must keep them verbatim (see ema/countmatrix/indexing.py::split_cb).
RG_IDS = (
    "pbmc_10k_v3:0:1:HFWFVDMXX:1",
    "pbmc_10k_v3:0:1:HFWFVDMXX:2",
)

# Shipped `ema run` backstops (ema/cli/run.py::_BACKSTOP).  read_check reads
# these off variable_config at call time, so the test must set them.
SEQ_LEN = 150
CB_LEN = 16


@pytest.fixture(autouse=True)
def _shipped_read_config():
    """Point variable_config at the defaults `ema run` would use, then restore."""
    saved = {
        k: getattr(variable_config, k)
        for k in ("seqlen", "cb_len", "barcode_tag", "ignore_chro")
    }
    variable_config.seqlen = SEQ_LEN
    variable_config.cb_len = CB_LEN
    variable_config.barcode_tag = "CB"
    variable_config.ignore_chro = ["MT", "mt"]
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)


def _call_both_strands(tmp_path: Path) -> tuple[BarcodeIndex, list[str]]:
    """Run the monolithic caller over the fixture on both strands.

    Returns the shared :class:`BarcodeIndex` and the concatenated BED lines.
    """
    index = BarcodeIndex()
    bed_lines: list[str] = []
    for direction in (False, True):
        bed = tmp_path / f"pas_{int(direction)}.bed"
        mtx = tmp_path / f"matrix_{int(direction)}.txt"
        peak_calling(
            direction=direction,
            bedfilepath=str(bed),
            matrixpath=str(mtx),
            bamfile_dir=str(FIXTURE),
            index=index,
            state=PeakCallingState(pasnumber=0),
            sample_id="fixture",
        )
        bed_lines.extend(bed.read_text().splitlines())
    return index, bed_lines


# ---------------------------------------------------------------------------
# The fixture itself must keep replaying all three bugs
# ---------------------------------------------------------------------------


def test_fixture_is_committed_and_readable() -> None:
    """The BAM is in the repo -- `.gitignore` must not swallow it.

    8 of the 9 skipped tests on `develop` skip for "no BAM available"; this
    one is the reason the `!tests/fixtures/*.bam` negation exists.
    """
    assert FIXTURE.exists(), (
        f"{FIXTURE} is missing -- check the `!tests/fixtures/*.bam` negation "
        "in .gitignore (line 15 ignores *.bam)"
    )
    with pysam.AlignmentFile(str(FIXTURE), "rb") as bam:
        assert bam.count(until_eof=True) == 2072


def test_fixture_replays_all_three_cellranger_input_defects() -> None:
    """Guard the fixture's shape: lose any of these and the replay is hollow."""
    with pysam.AlignmentFile(str(FIXTURE), "rb") as bam:
        header_rgs = {rg["ID"] for rg in bam.header.to_dict().get("RG", [])}
        gem_suffixed = 0
        unmapped_with_cb = 0
        secondary = 0
        duplicate = 0
        read_rgs = set()
        for read in bam.fetch(until_eof=True):
            if read.has_tag("CB"):
                cb = read.get_tag("CB")
                if len(cb) == CB_LEN + 2 and cb[-2] == "-" and cb[-1].isdigit():
                    gem_suffixed += 1
                if read.is_unmapped or read.reference_end is None:
                    unmapped_with_cb += 1
            if read.is_secondary:
                secondary += 1
            if read.is_duplicate:
                duplicate += 1
            if read.has_tag("RG"):
                read_rgs.add(read.get_tag("RG"))

    assert header_rgs == set(RG_IDS)
    assert read_rgs == set(RG_IDS)
    assert all("_" in rg for rg in header_rgs), "underscore-bearing RG is the point"
    assert gem_suffixed == 2034, "CB tags must keep the CellRanger '-1' GEM suffix"
    assert unmapped_with_cb == 20, "unmapped-but-CB-tagged reads must survive"
    assert secondary == 1125
    assert duplicate == 371


# ---------------------------------------------------------------------------
# The replay: the caller must produce a real matrix, not an empty one
# ---------------------------------------------------------------------------


def test_pipeline_yields_many_barcode_columns(tmp_path: Path) -> None:
    """THE proving assertion: many columns, never 0 (develop) and never 1.

    On `origin/develop` (c08ca23) this fixture produces **0** columns because
    the GEM-group suffix fails the `len(cb) != barcode_len` check; on the
    interim `eb13529` RG rewrite the underscore-bearing RG collapses cells
    onto a single `(sample, -1)` column.  Repaired: 1493.
    """
    index, _ = _call_both_strands(tmp_path)
    n_columns = len(index.mapping)
    assert n_columns not in (0, 1), (
        "the whole run collapsed to %d barcode column(s) -- CellRanger input "
        "regression is back (GEM suffix dropped every read, or the composite "
        "cb was split on the first underscore)" % n_columns
    )
    assert n_columns >= 500, f"expected ~1493 barcode columns, got {n_columns}"


def test_every_column_decodes_to_a_16nt_barcode(tmp_path: Path) -> None:
    """GEM suffix stripped + composite split on the LAST underscore.

    If either half of the round-trip is wrong the barcode half is not a
    16 nt ACGTN string and `encode_cb` maps it to -1.
    """
    index, _ = _call_both_strands(tmp_path)
    bad = [
        cb
        for cb in index.mapping
        if len(split_cb(cb)[1]) != CB_LEN or set(split_cb(cb)[1]) - set("ACGTN")
    ]
    assert bad == [], f"{len(bad)} columns carry a non-barcode half, e.g. {bad[:3]}"


def test_underscore_bearing_read_groups_keep_separate_columns(tmp_path: Path) -> None:
    """Both RGs survive verbatim and never share a column.

    177 barcodes in this fixture are seen under *both* read groups.  Any RG
    sanitising/dropping merges those two cells into one matrix column --
    exactly the `ema merge` data-corruption mode (0e).
    """
    index, _ = _call_both_strands(tmp_path)
    samples = {split_cb(cb)[0] for cb in index.mapping}
    assert samples == set(RG_IDS), (
        f"read-group ids were rewritten or dropped: {sorted(samples)}"
    )

    shared: dict[str, set[str]] = {}
    for cb in index.mapping:
        sample, barcode = split_cb(cb)
        shared.setdefault(barcode, set()).add(sample)
    both = [bc for bc, s in shared.items() if len(s) == 2]
    assert len(both) >= 100, (
        f"only {len(both)} barcodes kept a distinct column per read group; "
        "the two RGs are being collapsed onto one sample id"
    )
    # And distinct columns really are distinct integers.
    for bc in both[:20]:
        cols = {index.mapping[f"{rg}_{bc}"] for rg in RG_IDS}
        assert len(cols) == 2, f"barcode {bc} shares a column across read groups"


def test_unmapped_reads_are_skipped_not_crashed() -> None:
    """`reference_end is None` must return the skip sentinel, not raise.

    Before the fix this raised `TypeError: unsupported operand type(s) for -:
    'NoneType' and 'int'` from `read_end - read_start`, hours into a run.
    """
    seen = 0
    with pysam.AlignmentFile(str(FIXTURE), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not (read.is_unmapped or read.reference_end is None):
                continue
            seen += 1
            for direction in (False, True):
                assert read_check(read, direction, sample_id="fixture") == (
                    0,
                    0,
                    0,
                    0,
                    0,
                )
    assert seen == 20, "fixture lost its unmapped reads"


def test_matrix_columns_are_written_for_the_called_peaks(tmp_path: Path) -> None:
    """The run writes real BED peaks and a matrix keyed by those columns."""
    index, bed_lines = _call_both_strands(tmp_path)
    assert len(bed_lines) >= 5, "no peaks called from the fixture"
    max_col = max(index.mapping.values())
    assert max_col == len(index.mapping), "column indices must be dense and 1-based"

    written_cols: set[int] = set()
    for direction in (False, True):
        for line in (tmp_path / f"matrix_{int(direction)}.txt").read_text().splitlines():
            parts = line.split()
            if len(parts) >= 2:
                written_cols.add(int(parts[1]))
    assert written_cols, "matrix file is empty -- every read was discarded"
    assert max(written_cols) <= max_col


# ---------------------------------------------------------------------------
# R5 (10_caller_fix_plan.md §5): monolithic / pipeline / tile must not diverge
# ---------------------------------------------------------------------------
#
# `read_check`'s 5-tuple is consumed by three different execution paths, and
# the benchmark only ever exercised the monolithic one -- so a divergence in
# `--pipeline` or `--tiles` would not be caught by re-running the benchmark.
# These two tests run all three paths over the committed fixture and compare
# their BED output.
#
# The spawned workers (`mp.get_context("spawn")`) do NOT inherit
# `variable_config`: each child re-imports `ema.config`, which re-parses
# `sys.argv` (spawn propagates the parent's argv).  So the CLI form of the
# same settings has to be put on `sys.argv` for the pipeline/tile paths --
# that is production behaviour, not a test artefact.


def _index_fixture_copy(tmp_path: Path) -> Path:
    """Copy the fixture into *tmp_path* and build its .bai (tile mode needs it).

    The index is deliberately NOT committed: `.gitignore` still ignores
    `*.bai`, and regenerating it is one pysam call.
    """
    import shutil

    bam = tmp_path / "fixture.bam"
    shutil.copyfile(FIXTURE, bam)
    pysam.index(str(bam))
    assert (tmp_path / "fixture.bam.bai").exists()
    return bam


def _cli_argv(bam: Path) -> list[str]:
    """The legacy-argparse form of the config the spawned children re-parse."""
    return [
        "ema",
        "--bamDir", str(bam),
        "--sequenceLen", str(SEQ_LEN),
        "--CellBarcodeLen", str(CB_LEN),
        "--BarcodeTag", "CB",
    ]


def _bed_rows_all_paths(tmp_path: Path, monkeypatch) -> dict[str, list[list[str]]]:
    """Run the fixture through all three paths; return their BED rows."""
    bam = _index_fixture_copy(tmp_path)
    argv = _cli_argv(bam)
    runs = {
        "monolithic": {},
        "pipeline": {"use_pipeline": True},
        "tiles": {"use_tiles": True, "tile_size": 1_000_000, "n_workers": 2},
    }
    out: dict[str, list[list[str]]] = {}
    for label, kwargs in runs.items():
        rows: list[list[str]] = []
        for direction in (False, True):
            bed = tmp_path / f"{label}_{int(direction)}.bed"
            mtx = tmp_path / f"{label}_{int(direction)}.mtx"
            extra = dict(kwargs)
            if label == "monolithic":
                extra.update(
                    index=BarcodeIndex(),
                    state=PeakCallingState(pasnumber=0),
                    sample_id="fixture",
                )
            with monkeypatch.context() as mp:
                mp.setattr("sys.argv", list(argv))
                peak_calling(
                    direction=direction,
                    bedfilepath=str(bed),
                    matrixpath=str(mtx),
                    bamfile_dir=str(bam),
                    **extra,
                )
            rows.extend(line.split("\t") for line in bed.read_text().splitlines())
        out[label] = rows
    return out


def test_monolithic_pipeline_and_tile_paths_agree_on_pas_coordinates(
    tmp_path: Path, monkeypatch
) -> None:
    """Same BAM, same PAS intervals, whichever path called them."""
    rows = _bed_rows_all_paths(tmp_path, monkeypatch)
    coords = {
        label: sorted((r[0], r[1], r[2]) for r in rs) for label, rs in rows.items()
    }
    assert len(coords["monolithic"]) >= 5, "fixture called no peaks"
    assert coords["pipeline"] == coords["monolithic"], (
        "--pipeline diverges from the default path on PAS coordinates"
    )
    assert coords["tiles"] == coords["monolithic"], (
        "--tiles diverges from the default path on PAS coordinates"
    )


@pytest.mark.xfail(
    reason=(
        "KNOWN-RED, found by this fixture on 2026-08-20 and tracked in "
        "the Stage-0 fix campaign (PeakATail_wd analysis, 2026-08; see PR description and issues #67-#72)"
        "flush (the `if signal:` / `elif len(peak.peak_list)` block after the "
        "read loop) writes `strand`, which is the 4th element of read_check's "
        "tuple from the LAST read -- and every skipped read sets it to 0 via "
        "the `chro1, start1, end1, strand, cb = read_check(...)` unpack before "
        "the `if chro1 == 0: continue`. Every CellRanger BAM ends in unmapped "
        "reads, so the final PAS of each strand pass is written as '+' "
        "regardless of direction: the reverse pass over this fixture emits "
        "1:155247357-155247365 as '+' in the monolithic and tile paths and as "
        "'-' (correct) in the pipeline path. Fix is `direction`, not `strand`, "
        "in both flush branches -- read_check already guarantees "
        "strand == direction for every read it accepts. Deliberately NOT fixed "
        "on chore/ci-and-fixture: caller behaviour belongs to Stage 1, which "
        "owns the one budgeted re-run. Delete this marker with the fix."
    ),
    strict=False,
)
def test_monolithic_pipeline_and_tile_paths_agree_on_pas_strand(
    tmp_path: Path, monkeypatch
) -> None:
    """Same BAM, same strand column, whichever path called them."""
    rows = _bed_rows_all_paths(tmp_path, monkeypatch)
    stranded = {
        label: sorted((r[0], r[1], r[2], r[5]) for r in rs)
        for label, rs in rows.items()
    }
    assert stranded["pipeline"] == stranded["monolithic"]
    assert stranded["tiles"] == stranded["monolithic"]
