"""The `peakAtail-prime` compatibility guarantee, asserted on a committed fixture.

**The cardinal rule of this branch** (``PRIME_PLAN.md``, ``manuscript/24_prime_preregistration.md``
§1): the branch may change the caller's defaults, but v2 -- commit ``9dfdefb``, the code that
produced the numbers currently in the manuscript -- must stay reproducible **byte-for-byte**
through explicit flags.  The manuscript's shipped / v2 / prime comparison is meaningless the
moment that stops being true.

This test is the cheap, always-on half of that guarantee.  It runs the real caller
(``ema.countmatrix.peackcalling.peak_calling``, i.e. ``read.py`` -> ``polya.py`` ->
the peak state machine -> the strategy) over the 2,072-read fixture
``tests/fixtures/cellranger_pbmc_tiny.bam`` on **v2 settings**, and asserts the SHA-256 of every
output byte -- the BED, the count matrix **and** the ``pas_support.tsv`` sidecar, which TASK C
(``--pas-features``) appends columns to.  The expensive half is ``scripts/prime/identity_check.py`` on the real PBMC
chr19+21 slice, which compares whole run trees and is run by hand after every behavioural change.

Two arms, chosen because between them they cover every module Changes 1-5 touch:

* ``clip_seeded_91`` -- the manuscript arm's geometry (``--seq-len 91``, ``--peak-strategy
  clip_seeded``): the poly(A) clip detector, the clip-cluster seeding, the (CB,UMI) molecule
  counts in BED column 5, and the ``pas_support.tsv`` sidecar.
* ``original_150`` -- the coverage path at the shipped ``peakatail run`` backstop ``--seq-len 150``:
  the read-acceptance geometry of ``read.py`` and the coverage state machine, with poly(A)
  evidence in annotate-only mode.

The arms are sensitive to the read geometry by construction: the same fixture at
``seqlen`` 91 vs 150 gives different matrices, which is exactly what Change 1
(``--read-span-mode``) will move.

WHEN THIS TEST FAILS AFTER A PRIME CHANGE
-----------------------------------------
**Do not regenerate the goldens.**  A failure means either

1. the change is not flag-gated (fix the change: every behavioural change is behind a flag), or
2. it is flag-gated but its v2 value is not what ``_v2_settings()`` pins (fix ``_v2_settings()``
   -- that function is the single place where a new compat knob is pinned to its v2 value, and
   it is what ``--compat v2`` must be equivalent to).

The goldens themselves were produced by the **frozen v2 worktree**
``tools/pa-polya-run-9dfdefb3`` (``git rev-parse HEAD`` = ``9dfdefb3eb35...``), not by this
branch, on fixture md5 ``331574d31cf861e1c7a8bf84364d09ad``.  Regenerate them only if the
fixture itself is regenerated, and then from that frozen worktree, never from prime.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from ema.config import args, variable_config
from ema.countmatrix.indexing import BarcodeIndex
from ema.countmatrix.peak_state import PeakCallingState
from ema.countmatrix.paswrite import support_path_for
from ema.countmatrix.peackcalling import peak_calling
from ema.strategies import get_strategy

FIXTURE = Path(__file__).parent / "fixtures" / "cellranger_pbmc_tiny.bam"
FIXTURE_MD5 = "331574d31cf861e1c7a8bf84364d09ad"
V2_COMMIT = "9dfdefb3eb353b0817ef79c4eb9ace6d6c8aab53"

# SHA-256 of every byte the v2 caller writes for each arm.  Produced on
# tools/pa-polya-run-9dfdefb3 (see the module docstring).
V2_GOLDEN: dict[str, dict[str, str]] = {
    "clip_seeded_91": {
        "matrix_0.txt": "722e3001cef8250446fd01ff4eafebac0486c2a2f8ffc8792a7a05e6826621a6",
        "matrix_1.txt": "ddcadbffedc66bb51c23739c385d9511ab7f31e394cc974d4238ffd0ab03306c",
        "pas_0.bed": "2518f62b96263cd3416a98410199112ba0b03f3fda2a4bff620f6659e7314825",
        "pas_1.bed": "5cd77ab689e9379642200f89a75d24a02c0a541a3d2031e909cacfc10e2e4c6a",
        "pas_0.support.tsv": "1cb1abe1eba0a8fd25a18172a41f99fb1f2d64c81394a6305eb3329faa2d40c8",
        "pas_1.support.tsv": "9cf8efcd9d99b03a6502e417314f148521861c1b7608a1d8329bcf4c5de9eaf8",
    },
    "original_150": {
        "matrix_0.txt": "b9b063eb881bf68232fa1bb6410f743ce475bb7dd0bd22a188bcc80dc63ebb1e",
        "matrix_1.txt": "b7d3cdc34760ec90b3322fb07ba0a8ba478e3d90f1f8aee9fa3867cd814842db",
        "pas_0.bed": "1efb9d0f9a3c44c38aa0805955cd03a32f266d74a0b4573caeabc64472d34b69",
        "pas_1.bed": "e79b2fa134e12610003e46b21b0f43f04847777138a10b8637faf47e71b26cf4",
        "pas_0.support.tsv": "32aa26512a43a6040e4589f7724ba0f9e33078f1d2b72fd83b8ae85a17d13ef6",
        "pas_1.support.tsv": "27eb45b7684fb183d051bab13c5685926b471cc1779bfd5a69f688b343560218",
    },
}

# (arm name, peak strategy, --seq-len)
ARMS = (
    ("clip_seeded_91", "clip_seeded", 91),
    ("original_150", "original", 150),
)


def _v2_settings(seq_len: int) -> None:
    """Pin every knob this branch may move to its **v2** value.

    THE COMPAT SURFACE LIVES HERE.  Each of Changes 1-5 in ``PRIME_PLAN.md`` adds one line to
    this function setting its new option to the v2 value, and ``--compat v2`` must be exactly
    equivalent to what this function does.  Nothing else in this file changes.
    """
    variable_config.seqlen = seq_len
    variable_config.cb_len = 16
    variable_config.barcode_tag = "CB"
    variable_config.ignore_chro = ["MT", "mt"]
    # --- v2 values of prime options, added as each change lands -------------
    # TASK A (--read-geometry): v2 discards a read whose REFERENCE span exceeds
    # --seq-len and rewrites a shorter read's end to start + seq_len.  The
    # branch default is "true"; this line is what makes these two arms a
    # compatibility test rather than a change detector.
    variable_config.read_geometry = "fixed"
    # TASK A (--read-exclude-flags): v2 applies no SAM-flag filter on the
    # coverage channel.  This is already the branch default; pinned anyway so
    # the compat surface is complete and readable in one place.
    variable_config.read_exclude_flags = 0
    # TASK C (--pas-features): v2 writes seven sidecar columns.  The branch
    # default is "on", which APPENDS two more at call time (clip_positions,
    # clip_span) and more still at the internal-priming seam.  Appending
    # columns still changes the sidecar's BYTES, which is why it is a flag
    # and why it is pinned here.
    variable_config.pas_features = "off"
    # TASK D (--pas-score): v2 has no per-site score and no pas_score column.
    # "none" is ALSO the branch default (manuscript/24 3.3 ships a score
    # flag-off until 3.1 is met), so this line is not what makes these arms
    # pass today -- it is here so the compat surface stays complete in one
    # place if that default ever moves, and because --compat v2 must be
    # exactly equivalent to this function.
    variable_config.pas_score = "none"
    variable_config.pas_score_model = "prime1"
    variable_config.pas_score_min = -1.0
    # TASK E (--cleavage-offset): v2 spelled "no shift" as the integer 0; the
    # branch spells it "none" and both parse to the same no-op.  Pinned as the
    # int so this function keeps stating the v2 SPELLING, not the branch's.
    variable_config.cleavage_offset = 0
    variable_config.auto_cleavage_offset = False
    # TASK E (--emit-inferred-cleavage): v2 has no inferred_cleavage column.
    # The branch default is "on" and it changes the sidecar's bytes, so this
    # line is load-bearing for the two sidecar goldens below.
    variable_config.emit_inferred_cleavage = "off"
    # TASK E (--clip-rate-sampling): v2 sampled the head of the file.  The QC
    # is log-only so this cannot move a golden, but the compat surface must be
    # complete: --compat v2 has to reproduce v2's LOG line too.
    variable_config.clip_rate_sampling = "head"
    # TASK E (--ip-filter-default): v2 ran the internal-priming veto only when
    # --ip-filter was passed.  The veto lives above peak_calling(), so it
    # cannot move these goldens either; pinned for the same reason.
    variable_config.ip_filter_default = "off"
    # TASK E (--ip-filter-mode): v2's literal default was "annotate" -- keep
    # every flagged PAS, drop none.  The branch default is now the sentinel
    # "auto", which resolves to "filter", because leaving it at "annotate"
    # made --ip-filter-default auto a no-op (the veto ran and dropped
    # nothing).  This is the FIRST compat knob that lives on the legacy
    # `args` namespace rather than variable_config, which is why
    # test_prime_compat_flags.py now reads pins from both.
    args.ip_filter_mode = "annotate"
    # Change 4 (--polya-genomic-a-gate):  variable_config.polya_genomic_a_gate = False


@pytest.fixture(autouse=True)
def _restore_variable_config():
    """Every knob `_v2_settings` touches is process-global; put it back afterwards."""
    keys = ("seqlen", "cb_len", "barcode_tag", "ignore_chro",
            "read_geometry", "read_exclude_flags", "pas_features",
            "pas_score", "pas_score_model", "pas_score_min",
            "cleavage_offset", "auto_cleavage_offset",
            "emit_inferred_cleavage", "clip_rate_sampling",
            "ip_filter_default")
    saved = {k: getattr(variable_config, k) for k in keys}
    _ns = args._get()
    _MISSING = object()
    saved_args = {k: getattr(_ns, k, _MISSING) for k in ("ip_filter_mode",)}
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)
        for k, v in saved_args.items():
            if v is _MISSING:
                try:
                    delattr(_ns, k)
                except AttributeError:
                    pass
            else:
                setattr(_ns, k, v)


#: Every knob `_v2_settings` pins that is a peakAtail-prime option, i.e. every
#: one whose SHIPPED default may differ from the v2 value pinned above.  Used
#: by :func:`_branch_default_settings` to state the branch defaults without
#: hand-typing a single one of them.
_PRIME_KNOBS = (
    "read_geometry", "read_exclude_flags", "pas_features", "pas_score",
    "pas_score_model", "pas_score_min", "cleavage_offset",
    "auto_cleavage_offset", "emit_inferred_cleavage", "clip_rate_sampling",
    "ip_filter_default",
)


def _branch_default_settings(seq_len: int) -> None:
    """The opposite of :func:`_v2_settings`: every prime knob at its SHIPPED default.

    Only the four knobs that describe the FIXTURE BAM are set by hand; every
    prime option is read out of ``RunConfig``'s own field defaults, so this
    really is "what a user gets by typing nothing" and not a second hand-kept
    list that can drift from the schema.
    """
    from ema.cli.config_schema import RunConfig

    variable_config.seqlen = seq_len
    variable_config.cb_len = 16
    variable_config.barcode_tag = "CB"
    variable_config.ignore_chro = ["MT", "mt"]
    cfg = RunConfig()
    for name in _PRIME_KNOBS:
        setattr(variable_config, name, getattr(cfg, name))
    args.ip_filter_mode = cfg.ip_filter_mode


def _run_arm(strategy: str, seq_len: int, tmp_path: Path,
             settings=_v2_settings) -> dict[str, str]:
    """Run both strands of the caller over the fixture; return {filename: sha256}."""
    settings(seq_len)
    index = BarcodeIndex()
    state = PeakCallingState(pasnumber=0)
    digests: dict[str, str] = {}
    for direction in (False, True):
        bed = tmp_path / f"pas_{int(direction)}.bed"
        matrix = tmp_path / f"matrix_{int(direction)}.txt"
        peak_calling(
            direction=direction,
            bedfilepath=str(bed),
            matrixpath=str(matrix),
            bamfile_dir=str(FIXTURE),
            index=index,
            state=state,
            sample_id="fixture",
            strategy=get_strategy(strategy),
        )
        support = Path(support_path_for(bed))
        for path in (bed, matrix, support):
            if path.exists():
                digests[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return digests


def test_fixture_is_the_one_the_goldens_were_made_from() -> None:
    """A golden is only meaningful against a known input."""
    assert FIXTURE.exists(), f"{FIXTURE} missing (see .gitignore's !tests/fixtures/*.bam)"
    assert hashlib.md5(FIXTURE.read_bytes()).hexdigest() == FIXTURE_MD5, (
        "the fixture BAM changed; the v2 goldens in this file were computed from "
        f"md5 {FIXTURE_MD5} and must be regenerated from tools/pa-polya-run-9dfdefb3"
    )


@pytest.mark.parametrize(("arm", "strategy", "seq_len"), ARMS, ids=[a[0] for a in ARMS])
def test_v2_output_is_reproduced_byte_for_byte(
    arm: str, strategy: str, seq_len: int, tmp_path: Path
) -> None:
    """v2 (9dfdefb) bytes, on v2 settings, from this branch.

    See the module docstring before touching the golden values: a failure here is a
    finding about the change, not about the goldens.
    """
    got = _run_arm(strategy, seq_len, tmp_path)
    want = V2_GOLDEN[arm]
    assert set(got) == set(want), f"arm {arm} wrote {sorted(got)}, expected {sorted(want)}"
    bad = {k: (want[k], got[k]) for k in want if want[k] != got[k]}
    assert not bad, (
        "v2 byte-identity LOST on arm %s for %s.\n"
        "Every behavioural change on peakAtail-prime must be behind a flag whose v2 value is "
        "pinned in _v2_settings(); do NOT regenerate these goldens.\n%s"
        % (arm, ", ".join(sorted(bad)), bad)
    )


@pytest.mark.parametrize(("arm", "strategy", "seq_len"), ARMS, ids=[a[0] for a in ARMS])
def test_branch_defaults_still_write_v2s_beds_and_matrix(
    arm: str, strategy: str, seq_len: int, tmp_path: Path
) -> None:
    """The PR's headline property, which nothing else pinned.

    The compat test above proves v2 is REACHABLE (with 14 flags).  This proves
    the thing the PR body actually leads with: at the caller, **typing nothing**
    already gives v2's BED and count-matrix bytes -- the internal-priming veto
    that does move the call set lives above ``peak_calling()``, so every prime
    default at this seam is call-set-neutral by construction, and this is the
    test that says so instead of asserting it in prose.

    ``pas_support.tsv`` is deliberately EXCLUDED and asserted to differ:
    ``--pas-features on`` and ``--emit-inferred-cleavage on`` append sidecar
    columns at the branch defaults.  Claiming the flagless output is
    byte-identical to v2's *without* that exclusion is the CHANGELOG error this
    test exists to keep out of the tree.
    """
    got = _run_arm(strategy, seq_len, tmp_path, settings=_branch_default_settings)
    want = V2_GOLDEN[arm]
    call_set = {k: v for k, v in got.items() if not k.endswith(".support.tsv")}
    assert call_set == {k: v for k, v in want.items() if not k.endswith(".support.tsv")}, (
        "a peakAtail-prime DEFAULT moved the call set at the caller on arm %s. "
        "Every default this branch flips below the veto is supposed to be "
        "call-set-neutral; one of them is not." % arm
    )
    sidecar = {k: v for k, v in got.items() if k.endswith(".support.tsv")}
    assert sidecar and sidecar != {k: v for k, v in want.items()
                                   if k.endswith(".support.tsv")}, (
        "pas_support.tsv is byte-identical to v2's at the branch defaults, so "
        "the CHANGELOG may say the flagless output is byte-identical without "
        "scoping it -- but --pas-features on is supposed to append columns"
    )


def test_the_compat_pin_is_load_bearing(tmp_path: Path) -> None:
    """`_v2_settings` must be doing work, not agreeing with whatever is set.

    If the pinned values reproduced the v2 goldens no matter what, the two
    tests above would pass while proving nothing.  ``--read-geometry true``
    restores 24.19 % of the reads v2 discards, so it MUST move bytes on a
    fixture containing spliced or soft-clipped reads.  If this test starts
    failing, either the branch has no behavioural change left or the fixture
    stopped covering it -- both worth knowing.

    Deliberately independent of which value the branch DEFAULTS to: the
    default is an adoption decision that measurement can flip either way,
    while the compatibility guarantee must hold regardless.
    """
    from ema.countmatrix.read import READ_GEOMETRIES, V2_READ_GEOMETRY

    non_v2 = [g for g in READ_GEOMETRIES if g != V2_READ_GEOMETRY]
    assert non_v2, "there is no non-v2 geometry left to gate"
    branch_dir = tmp_path / "branch"
    branch_dir.mkdir()
    index = BarcodeIndex()
    state = PeakCallingState(pasnumber=0)
    digests: dict[str, str] = {}
    # NB: deliberately does NOT call _v2_settings().  The three BAM-shape
    # knobs are set by hand and the geometry is forced to a NON-v2 value, so
    # this test says the same thing whichever value the branch defaults to.
    variable_config.seqlen = 91
    variable_config.cb_len = 16
    variable_config.barcode_tag = "CB"
    variable_config.read_geometry = non_v2[-1]
    variable_config.pas_features = "on"
    for direction in (False, True):
        bed = branch_dir / f"pas_{int(direction)}.bed"
        matrix = branch_dir / f"matrix_{int(direction)}.txt"
        peak_calling(
            direction=direction, bedfilepath=str(bed), matrixpath=str(matrix),
            bamfile_dir=str(FIXTURE), index=index, state=state,
            sample_id="fixture", strategy=get_strategy("clip_seeded"),
        )
        support = Path(support_path_for(bed))
        for path in (bed, matrix, support):
            if path.exists():
                digests[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digests != V2_GOLDEN["clip_seeded_91"], (
        "the branch default reproduced the v2 goldens byte-for-byte, so the "
        "compat test above cannot distinguish 'flag-gated' from 'no change'"
    )


def test_the_arms_actually_carry_signal(tmp_path: Path) -> None:
    """Guard against a hollow golden: empty files hash consistently too."""
    clip_dir, plain_dir = tmp_path / "clip", tmp_path / "plain"
    clip_dir.mkdir()
    plain_dir.mkdir()
    clip = _run_arm("clip_seeded", 91, clip_dir)
    plain = _run_arm("original", 150, plain_dir)
    for name, digests, root in (("clip_seeded_91", clip, clip_dir),
                                ("original_150", plain, plain_dir)):
        for fname in digests:
            size = (root / fname).stat().st_size
            assert size > 0, f"{name}/{fname} is empty -- the golden would be hollow"
        assert (root / "matrix_0.txt").stat().st_size > 1000, f"{name}: matrix_0 suspiciously small"
    assert clip != plain, (
        "the two arms produced identical bytes -- they are not covering different code paths"
    )
