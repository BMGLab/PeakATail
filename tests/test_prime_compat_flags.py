"""The v2-compat FLAG LIST and the v2-compat PIN must not drift apart.

The branch's cardinal rule is that v2 (commit ``9dfdefb``) stays reproducible
byte-for-byte.  There are two expressions of that promise and they were not
connected to each other:

* ``tests/test_prime_v2_compat_golden.py::_v2_settings()`` -- the LIBRARY-level
  pin, checked against goldens on the committed fixture; and
* the CLI incantation a reviewer actually types, which existed only as prose in
  ``CHANGELOG.md`` and was checked by nothing.

A new prime option with a non-v2 default, added to ``_v2_settings()`` and
forgotten in the command line, would leave the unit test green while every
*published* compat run silently stopped being v2 -- which is exactly the
failure the cardinal rule exists to prevent, and the hardest kind to notice
after the fact because the run tree still looks fine.

``ema.cli.config_schema.V2_COMPAT_FLAGS`` is now the canonical copy.  This
module ties it to ``_v2_settings()`` in BOTH directions and checks the prose
against it.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from ema.cli.config_schema import RunConfig, V2_COMPAT_FLAGS, field_specs
from ema.config import args, variable_config

_MISSING = object()

import importlib.util as _ilu

CHANGELOG = Path(__file__).resolve().parents[1] / "CHANGELOG.md"

#: ``_v2_settings()`` also sets the four knobs that describe the FIXTURE BAM
#: (its read length, barcode length, barcode tag and ignored contigs).  Those
#: are not compat knobs -- they are the same in v2 and on this branch, and a
#: real run states them on its own command line -- so they are excluded here.
#: Listed by name rather than inferred, so that a genuinely new compat knob
#: cannot hide behind a heuristic.
_FIXTURE_SHAPE = {"seq_len", "cb_len", "barcode_tag", "ignore_chro"}

#: peakAtail-prime options bridged to the legacy ``args`` namespace instead of
#: ``variable_config``, and which ``_v2_settings()`` does NOT pin.  The fixture
#: golden is blind to those, so the test below requires their branch default to
#: BE their v2 value -- the only thing that makes the blindness safe.
#:
#: ``--ip-filter-mode`` is deliberately NOT in here.  It is args-bridged AND
#: its branch default is not v2's, which is exactly the combination that let
#: the no-op default ship: the veto turned itself on while the mode literal
#: stayed at v2's ``annotate``, so nothing was dropped.  ``_v2_settings()``
#: now pins it on ``args`` and :func:`_pinned_value` reads pins from both
#: holders, so the command line and the library pin are checked against each
#: other for args-bridged options too.
#: ``--dynamic-threshold-clamp`` qualifies: bridged to ``args``, default
#: ``off``, and the fixture golden cannot see it because the fixture never
#: turns ``--dynamic-threshold`` on.  Since the merge with ``develop`` the
#: flag is an accepted NO-OP -- the look-back bound is unconditional (issue
#: #101) -- so ``off`` is trivially still safe to leave on the documented
#: compat command line.  ``tests/test_prime_dynamic_threshold_clamp.py`` pins
#: the no-op contract.
_ARGS_BRIDGED = {"--pas-gene-rescue", "--pas-gene-rescue-min-mol",
                 "--dynamic-threshold-clamp"}


def _pinned_value(name, spec):
    """The value ``_v2_settings()`` pins for a field -- wherever it lives.

    Returns ``_MISSING`` when the pin does not touch the field at all.
    """
    lda = spec.legacy_dataclass_attr or ""
    if lda.startswith("variable_config."):
        return getattr(variable_config, lda.split(".", 1)[1], _MISSING)
    if not lda:
        return getattr(args._get(), spec.legacy_args_attr or name, _MISSING)
    return _MISSING


def _load_v2_settings():
    """Import ``_v2_settings`` from its own module without assuming that
    ``tests/`` is on ``sys.path`` (it is not under the default rootdir)."""
    src = Path(__file__).with_name("test_prime_v2_compat_golden.py")
    spec = _ilu.spec_from_file_location("_prime_compat_golden", src)
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod._v2_settings


_v2_settings = _load_v2_settings()


def _flag_pairs() -> dict[str, str]:
    assert len(V2_COMPAT_FLAGS) % 2 == 0, "V2_COMPAT_FLAGS is flag/value pairs"
    return dict(zip(V2_COMPAT_FLAGS[::2], V2_COMPAT_FLAGS[1::2]))


def _by_cli_flag() -> dict[str, tuple[str, object]]:
    """``{cli_flag: (field_name, spec)}`` for every RunConfig field."""
    return {spec.cli_flag: (name, spec)
            for name, spec in field_specs(RunConfig).items() if spec.cli_flag}


def _coerce(spec, text: str, default):
    if isinstance(default, bool):
        return text.lower() in ("1", "true", "yes", "on")
    if isinstance(default, int) and not isinstance(default, bool):
        return int(text)
    if isinstance(default, float):
        return float(text)
    return text


@pytest.fixture(autouse=True)
def _restore():
    keys = [n for n, s in field_specs(RunConfig).items()
            if s.legacy_dataclass_attr
            and s.legacy_dataclass_attr.startswith("variable_config.")]
    attrs = [field_specs(RunConfig)[k].legacy_dataclass_attr.split(".", 1)[1]
             for k in keys]
    saved = {a: getattr(variable_config, a) for a in attrs
             if hasattr(variable_config, a)}
    # ...and the args-bridged pins (`_v2_settings()` sets ip_filter_mode).
    ns = args._get()
    args_attrs = [s_.legacy_args_attr or n for n, s_ in field_specs(RunConfig).items()
                  if not (s_.legacy_dataclass_attr or "")]
    saved_args = {a: getattr(ns, a, _MISSING) for a in args_attrs}
    try:
        yield
    finally:
        for a, v in saved.items():
            setattr(variable_config, a, v)
        for a, v in saved_args.items():
            if v is _MISSING:
                try:
                    delattr(ns, a)
                except AttributeError:
                    pass
            else:
                setattr(ns, a, v)


def test_every_compat_flag_is_a_real_option():
    known = _by_cli_flag()
    for flag in _flag_pairs():
        assert flag in known, (
            f"{flag} is in V2_COMPAT_FLAGS but is not a RunConfig cli_flag; a "
            "typo here makes the documented compat command silently wrong"
        )


def test_every_knob_the_pin_moves_is_on_the_command_line():
    """`_v2_settings()` -> `V2_COMPAT_FLAGS`.

    Anything the library-level pin has to CHANGE from the branch default is a
    knob whose v2 value a reviewer must type; if it is missing from the
    documented command, a published compat run is not v2.
    """
    specs = field_specs(RunConfig)
    defaults = {f.name: f.default for f in RunConfig.__dataclass_fields__.values()}
    # branch defaults into the legacy globals, then the v2 pin over the top.
    # BOTH holders: --ip-filter-mode is bridged to `args`, and it is the knob
    # whose v2/branch difference (annotate vs the "auto" sentinel that drops)
    # a variable_config-only sweep could not see.
    ns = args._get()
    for name, spec in specs.items():
        lda = spec.legacy_dataclass_attr or ""
        if lda.startswith("variable_config."):
            setattr(variable_config, lda.split(".", 1)[1], defaults[name])
        elif not lda:
            setattr(ns, spec.legacy_args_attr or name, defaults[name])
    _v2_settings(91)

    flags = _flag_pairs()
    by_flag = _by_cli_flag()
    field_of = {name: flag for flag, (name, _s) in by_flag.items()}
    missed = []
    for name, spec in specs.items():
        lda = spec.legacy_dataclass_attr or ""
        if lda and not lda.startswith("variable_config."):
            continue          # directory_config / filter_config: not compat knobs
        pinned = _pinned_value(name, spec)
        if pinned is _MISSING:
            continue
        branch = defaults[name]
        if pinned == branch:
            continue                      # branch default is already v2's
        if name in _FIXTURE_SHAPE:
            continue                      # describes the fixture, not v2
        if name == "cleavage_offset":
            # v2 spelled "no shift" as the int 0, the branch spells it "none";
            # both parse to the same no-op, which is what has to match.
            from ema.countmatrix.cleavage_offset import parse_cleavage_offset
            typed = flags.get(field_of.get(name, ""))
            assert typed is not None, "cleavage_offset missing from the command"
            assert parse_cleavage_offset(typed) == parse_cleavage_offset(pinned)
            continue
        flag = field_of.get(name)
        if flag is None or flag not in flags:
            missed.append((name, branch, pinned))
            continue
        got = _coerce(spec, flags[flag], branch)
        assert got == pinned, (
            f"{flag} is documented as {flags[flag]!r} but _v2_settings pins "
            f"{name} to {pinned!r}"
        )
    assert not missed, (
        "these prime options have a non-v2 default and are pinned by "
        "_v2_settings(), but the documented compat command line does not set "
        "them -- a published compat run would NOT be v2: %r" % (missed,)
    )


def test_every_flag_on_the_command_line_is_pinned_by_the_library():
    """`V2_COMPAT_FLAGS` -> `_v2_settings()`, the other direction.

    A flag typed on the command line whose value the library-level pin does not
    also set means the fixture golden is not testing what the slice run does.
    """
    specs = field_specs(RunConfig)
    by_flag = _by_cli_flag()
    defaults = {f.name: f.default for f in RunConfig.__dataclass_fields__.values()}
    _v2_settings(91)
    for flag, text in _flag_pairs().items():
        name, spec = by_flag[flag]
        lda = spec.legacy_dataclass_attr or ""
        if not lda.startswith("variable_config.") and \
                _pinned_value(name, spec) is not _MISSING and \
                _pinned_value(name, spec) != defaults[name]:
            # args-bridged AND actually pinned by _v2_settings() to something
            # other than the branch default -- check the command line against
            # the pin exactly as for a variable_config knob.  This is the
            # branch of --ip-filter-mode.
            assert _coerce(spec, text, defaults[name]) == _pinned_value(name, spec), (
                f"{flag} {text!r} on the command line vs "
                f"{_pinned_value(name, spec)!r} in _v2_settings"
            )
            continue
        if not lda.startswith("variable_config."):
            # STRUCTURAL GAP, declared rather than hidden: `_v2_settings()`
            # only touches variable_config, so an option bridged to the legacy
            # `args` namespace cannot be pinned by the fixture golden at all.
            # That is safe TODAY only because such an option's branch default
            # IS its v2 value -- i.e. there is nothing to undo.  Assert exactly
            # that, so the day someone gives one a non-v2 default this fails
            # and forces either a variable_config bridge or a wider pin.
            assert flag in _ARGS_BRIDGED, (
                f"{flag} is bridged to `args`, not variable_config, so "
                "_v2_settings() cannot pin it and the fixture golden cannot "
                "see it. Add it to _ARGS_BRIDGED after checking its default."
            )
            assert _coerce(spec, text, defaults[name]) == defaults[name], (
                f"{flag}'s branch default {defaults[name]!r} is NOT its v2 "
                f"value {text!r}, and it is bridged to `args` where "
                "_v2_settings() cannot reach it -- the fixture golden is blind "
                "to this option. Bridge it to variable_config, or extend the "
                "pin mechanism, before shipping a non-v2 default."
            )
            continue
        attr = lda.split(".", 1)[1]
        pinned = getattr(variable_config, attr)
        if name == "cleavage_offset":
            from ema.countmatrix.cleavage_offset import parse_cleavage_offset
            assert parse_cleavage_offset(text) == parse_cleavage_offset(pinned)
            continue
        assert _coerce(spec, text, defaults[name]) == pinned, (
            f"{flag} {text!r} on the command line vs {pinned!r} in _v2_settings"
        )


def test_the_changelog_prose_matches_the_tuple():
    """The paragraph a reviewer copies from must be the tuple, verbatim."""
    text = CHANGELOG.read_text()
    flags = _flag_pairs()
    # the compat paragraph, normalised over its line breaks
    m = re.search(r"\(`(--read-geometry[^`]*)`", text, re.S)
    assert m, "CHANGELOG no longer documents a v2-compat command line"
    typed = dict(zip(*[iter(m.group(1).split())] * 2))
    assert typed == flags, (
        "the CHANGELOG's compat command and V2_COMPAT_FLAGS disagree:\n"
        f"  prose : {typed}\n  tuple : {flags}"
    )


# ---------------------------------------------------------------------------
# completeness, checked against v2's RunConfig rather than a hand-kept list
#
# THE DEFECT THIS SECTION USED TO HAVE.  v2's defaults were read only from a
# worktree at one absolute path on one machine, behind a ``skipif``, so the two
# completeness checks below -- the branch's strongest guarantee -- SKIPPED
# everywhere else, CI included: a green checkmark that had never run them.
#
# So there are now two sources, in this order:
#
# * a live v2 worktree, when there is one (``PEAKATAIL_V2_WORKTREE``, or the
#   author's historical path); and
# * :data:`V2_DEFAULTS_SNAPSHOT`, a committed JSON dump of the same thing,
#   which makes the checks run unconditionally.
#
# The snapshot is a CACHE, never the authority: whenever a live worktree IS
# present, ``test_the_frozen_v2_snapshot_still_matches_the_v2_worktree`` asserts
# the two agree, so a stale snapshot is a failure rather than a quiet fiction.
# ---------------------------------------------------------------------------
#: Environment variable naming a checkout of v2 (commit ``9dfdefb``).  Set it
#: to re-derive v2's defaults live instead of trusting the snapshot.
V2_WORKTREE_ENV = "PEAKATAIL_V2_WORKTREE"

#: Where the frozen v2 worktree lives on the author's machine.  Used only when
#: it happens to exist; nothing depends on it any more.
V2_WORKTREE_FALLBACK = Path("/mnt/ssd1/Projects/PeakATail_wd/tools/pa-polya-run-9dfdefb3")

#: The committed dump of v2's ``RunConfig`` defaults.  Regeneration command is
#: inside the file, under ``_how_to_regenerate``.
V2_DEFAULTS_SNAPSHOT = Path(__file__).parent / "fixtures" / "v2_runconfig_defaults.json"


def _v2_worktree() -> Path | None:
    """A checkout of v2 to read defaults out of, or ``None`` to use the snapshot.

    Raises:
        AssertionError: if ``PEAKATAIL_V2_WORKTREE`` is set but does not point
            at a tree containing ``ema/cli/config_schema.py``.  Someone who
            asked for the live check must not silently get the cached one.
    """
    import os

    raw = os.environ.get(V2_WORKTREE_ENV)
    root = Path(raw) if raw else V2_WORKTREE_FALLBACK
    ok = (root / "ema" / "cli" / "config_schema.py").exists()
    assert ok or not raw, (
        f"{V2_WORKTREE_ENV}={raw!r} does not contain ema/cli/config_schema.py; "
        "unset it to fall back to tests/fixtures/v2_runconfig_defaults.json"
    )
    return root if ok else None


def _snapshot() -> dict:
    import json

    assert V2_DEFAULTS_SNAPSHOT.exists(), (
        f"{V2_DEFAULTS_SNAPSHOT} is missing -- without it the completeness "
        "checks below silently stop running anywhere without a v2 worktree, "
        "which is the exact failure they were added to end"
    )
    return json.loads(V2_DEFAULTS_SNAPSHOT.read_text())


def _v2_field_defaults_from_worktree(root: Path) -> dict:
    """``{field_name: default}`` for RunConfig as it stands in a v2 checkout."""
    import json
    import subprocess
    import sys

    out = subprocess.run(
        [sys.executable, "-c",
         f"import sys; sys.path.insert(0, {str(root)!r}); "
         "import dataclasses, json; from ema.cli.config_schema import RunConfig; "
         "print(json.dumps({f.name: f.default for f in dataclasses.fields(RunConfig)}, "
         "default=str))"],
        capture_output=True, text=True,
        env={"PYTHONDONTWRITEBYTECODE": "1", "PATH": "/usr/bin:/bin",
             "HOME": str(Path.home())},
    )
    assert out.returncode == 0, out.stderr[-2000:]
    return json.loads(out.stdout.strip().splitlines()[-1])


def _v2_field_defaults() -> dict:
    """``{field_name: default}`` for v2's RunConfig -- live tree if any, else snapshot."""
    root = _v2_worktree()
    if root is None:
        return _snapshot()["fields"]
    return _v2_field_defaults_from_worktree(root)


def test_the_frozen_v2_snapshot_is_the_v2_this_branch_claims_to_reproduce():
    """The snapshot must name the commit the goldens were made from, and be real.

    A cache nobody can trace back to a commit is worse than no cache: it would
    let the completeness checks below pass against whatever happened to be
    dumped into it.
    """
    from tests.test_prime_v2_compat_golden import V2_COMMIT

    snap = _snapshot()
    assert snap["v2_commit"] == V2_COMMIT, (
        "the frozen defaults were dumped from %r but the compat goldens are "
        "v2 = %r" % (snap["v2_commit"], V2_COMMIT)
    )
    fields = snap["fields"]
    assert len(fields) > 50, f"only {len(fields)} fields -- that is not v2's RunConfig"
    # The one field whose moved default IS the branch's behavioural change.
    # If this ever reads "auto", the snapshot was dumped from prime, not v2.
    assert fields["ip_filter_mode"] == "annotate", fields["ip_filter_mode"]
    # ...and the prime-only fields must be ABSENT, or the dump came from prime.
    assert "ip_filter_default" not in fields
    assert "pas_features" not in fields


def test_the_frozen_v2_snapshot_still_matches_the_v2_worktree():
    """Where a real v2 checkout exists, the cache must agree with it.

    Skipped (and only skipped) on a machine that has no v2 tree -- which is
    precisely why the two checks below no longer are.
    """
    root = _v2_worktree()
    if root is None:
        pytest.skip(
            f"no v2 worktree (set {V2_WORKTREE_ENV}); the frozen snapshot is "
            "used instead and the completeness checks still run"
        )
    live = _v2_field_defaults_from_worktree(root)
    frozen = _snapshot()["fields"]
    assert {k: str(v) for k, v in frozen.items()} == {k: str(v) for k, v in live.items()}, (
        f"{V2_DEFAULTS_SNAPSHOT} has drifted from {root}; regenerate it with "
        "the command in the file's _how_to_regenerate key"
    )


#: Fields that exist in v2's ``RunConfig`` but CANNOT appear on the documented
#: v2-compat command line, because that command line is an ``ema run``
#: invocation and these options belong to a different subcommand.  Declared by
#: name, with the reason, rather than inferred silently -- the whole point of
#: the check below is that a moved default must never go unnoticed.
#:
#: ``marker_top_n`` is a ``switch diff`` option.  ``develop`` moved its default
#: 200 -> 0 (issue #94: pre-selecting the tested PAS with the SAME cluster
#: labels the test then contrasts is a double-dip, and 0 was the only
#: FDR-controlled setting measured).  ``ema run`` does not accept
#: ``--marker-top-n``, so restoring v2's 200 on the compat RUN command is
#: impossible AND meaningless: no run-level output depends on it.  The default
#: itself is guarded by ``tests/test_marker_top_n_double_dip_i94.py``.
_NOT_ON_THE_RUN_COMMAND_LINE = {"marker_top_n"}


def _is_run_option(spec) -> bool:
    """True when a field can be typed on an ``ema run`` command line."""
    return not spec.applies_to or "run" in spec.applies_to


def test_every_v2_option_whose_default_this_branch_moved_is_on_the_command_line():
    """The check that would have caught the no-op default.

    The completeness test below asks "which fields are NEW on this branch?"
    and requires each to state its v2 value.  A field that already existed in
    v2 and whose DEFAULT this branch quietly moved is invisible to it -- and
    that is exactly what happened in reverse with ``--ip-filter-mode``: the
    branch added ``--ip-filter-default auto`` (a new field, duly listed) which
    turns the internal-priming veto ON, while ``--ip-filter-mode`` -- an
    *existing* v2 field -- kept v2's ``annotate``, so the veto dropped nothing
    and the branch's one behavioural default emitted v2's exact call set.
    Nothing checked the pair.

    So: diff the DEFAULTS, not just the field names.  Any v2 field whose
    branch default differs must appear on the documented compat command line
    with its v2 value, or a reviewer typing that command does not get v2.
    """
    import json

    v2_defaults = _v2_field_defaults()
    specs = field_specs(RunConfig)
    branch = {f.name: f.default for f in RunConfig.__dataclass_fields__.values()}
    flags = _flag_pairs()

    moved, unstated = [], []
    for name, spec in specs.items():
        if name not in v2_defaults or not spec.cli_flag:
            continue
        if str(branch[name]) == str(v2_defaults[name]):
            continue
        moved.append((name, v2_defaults[name], branch[name]))
        if not _is_run_option(spec):
            # Not typeable on an `ema run` line at all; must be declared.
            assert name in _NOT_ON_THE_RUN_COMMAND_LINE, (
                f"{spec.cli_flag} belongs to {sorted(spec.applies_to)} rather "
                "than to `ema run`, so the v2-compat RUN command cannot "
                "restore it. Add it to _NOT_ON_THE_RUN_COMMAND_LINE with the "
                "reason, and make sure something else guards its default."
            )
            continue
        text = flags.get(spec.cli_flag)
        if text is None:
            unstated.append((name, spec.cli_flag, v2_defaults[name], branch[name]))
            continue
        if name == "cleavage_offset":
            from ema.countmatrix.cleavage_offset import parse_cleavage_offset
            assert parse_cleavage_offset(text) == parse_cleavage_offset(
                v2_defaults[name])
            continue
        assert str(_coerce(spec, text, branch[name])) == str(v2_defaults[name]), (
            f"{spec.cli_flag} is documented as {text!r} but v2's default is "
            f"{v2_defaults[name]!r}"
        )

    assert not unstated, (
        "these options exist in v2 and this branch MOVED their default, but "
        "the documented compat command line does not restore them -- a "
        "published compat run would not be v2, and (as with --ip-filter-mode) "
        "nothing else would notice: %r" % (unstated,)
    )
    # A guard on the guard: this test is only meaningful while it has
    # something to check.  --ip-filter-mode is that something today.
    assert any(n == "ip_filter_mode" for n, _v2, _b in moved), (
        "--ip-filter-mode's branch default is v2's again; if that is "
        "deliberate, the branch has no behavioural default left"
    )


def test_every_option_this_branch_adds_is_on_the_compat_command_line():
    """Self-maintaining completeness: diff RunConfig against v2's.

    A hand-kept list of "the prime options" is exactly the thing that goes
    stale.  This asks the frozen v2 worktree which fields are new and requires
    each one to be either a flag/value pair in ``V2_COMPAT_FLAGS`` or a
    declared boolean flag in ``V2_COMPAT_OMITTED_FLAGS`` (whose v2 behaviour is
    "do not pass it").  A thirteenth prime option cannot then be added without
    someone stating its v2 value.
    """
    from ema.cli.config_schema import V2_COMPAT_OMITTED_FLAGS

    v2_fields = set(_v2_field_defaults())

    specs = field_specs(RunConfig)
    listed = set(_flag_pairs()) | set(V2_COMPAT_OMITTED_FLAGS)
    missing = []
    for name, spec in specs.items():
        if name in v2_fields or not spec.cli_flag:
            continue
        if spec.cli_flag not in listed:
            missing.append((name, spec.cli_flag))
    assert not missing, (
        "these options exist on peakAtail-prime and not in v2, and neither "
        "V2_COMPAT_FLAGS nor V2_COMPAT_OMITTED_FLAGS says what their v2 value "
        "is -- a reviewer typing the documented compat command would not "
        "restore them: %r" % (missing,)
    )
