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
output byte.  The expensive half is ``scripts/prime/identity_check.py`` on the real PBMC
chr19+21 slice, which compares whole run trees and is run by hand after every behavioural change.

Two arms, chosen because between them they cover every module Changes 1-5 touch:

* ``clip_seeded_91`` -- the manuscript arm's geometry (``--seq-len 91``, ``--peak-strategy
  clip_seeded``): the poly(A) clip detector, the clip-cluster seeding, the (CB,UMI) molecule
  counts in BED column 5.
* ``original_150`` -- the coverage path at the shipped ``ema run`` backstop ``--seq-len 150``:
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

from ema.config import variable_config
from ema.countmatrix.indexing import BarcodeIndex
from ema.countmatrix.peak_state import PeakCallingState
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
    },
    "original_150": {
        "matrix_0.txt": "b9b063eb881bf68232fa1bb6410f743ce475bb7dd0bd22a188bcc80dc63ebb1e",
        "matrix_1.txt": "b7d3cdc34760ec90b3322fb07ba0a8ba478e3d90f1f8aee9fa3867cd814842db",
        "pas_0.bed": "1efb9d0f9a3c44c38aa0805955cd03a32f266d74a0b4573caeabc64472d34b69",
        "pas_1.bed": "e79b2fa134e12610003e46b21b0f43f04847777138a10b8637faf47e71b26cf4",
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
    # Change 1 (--read-span-mode):        variable_config.read_span_mode = "truncate"
    # Change 2 (--polya-genomic-a-gate):  variable_config.polya_genomic_a_gate = False
    # Change 3 (--pas-score):             variable_config.pas_score = "none"
    # Change 5 (--cleavage-edge-refine):  variable_config.cleavage_edge_refine = False


@pytest.fixture(autouse=True)
def _restore_variable_config():
    """Every knob `_v2_settings` touches is process-global; put it back afterwards."""
    keys = ("seqlen", "cb_len", "barcode_tag", "ignore_chro")
    saved = {k: getattr(variable_config, k) for k in keys}
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)


def _run_arm(strategy: str, seq_len: int, tmp_path: Path) -> dict[str, str]:
    """Run both strands of the caller over the fixture; return {filename: sha256}."""
    _v2_settings(seq_len)
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
        for path in (bed, matrix):
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
