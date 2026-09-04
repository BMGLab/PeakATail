"""The branch's one behavioural default must not be a no-op.

THE DEFECT THIS FILE PINS
-------------------------
peakAtail-prime's TASK E item 3 makes the internal-priming veto run by itself:
``--ip-filter-default auto`` turns it on whenever a readable ``--genome-fasta``
is available, and the branch advertises that as "the largest measured accuracy
lift in this caller" (+7.7 % to +12.8 % relative recall at matched atlas
precision; the verifier's independent sweep puts the same quantity at +6.8 % to
+9.1 %).  That lift is the lift of **dropping** flagged sites.

Whether they are dropped was decided by a *different* flag, ``--ip-filter-mode``,
whose literal default was v2's ``"annotate"`` -- keep every PAS, record a flag.
So the branch shipped with the veto ON and the drop OFF.  Measured on the PBMC
chr19+21 slice with the branch defaults, a genome FASTA and no other flag::

    peak_filters_stats.json:  mode "annotate", flagged 5015, filtered 0
    pas 18865 | tier1 10886 | tier2 7979 | tier1>=2mol 3643

which is v2's call set to the row.  The branch's only behavioural default
changed nothing, and every published claim about it was measured on runs that
also passed ``--ip-filter-mode filter``.

THE FIX THESE TESTS PIN
-----------------------
``RunConfig.ip_filter_mode``'s default is now the SENTINEL
``"auto"`` (:data:`ema.cli.config_schema.IP_FILTER_MODE_UNSET`), meaning "the
user did not name a mode", and
:func:`ema.cli.config_schema.resolve_ip_filter_mode` turns it into ``"filter"``.
An explicit ``--ip-filter-mode annotate`` is honoured verbatim and is what
``--compat v2`` pins (``V2_COMPAT_FLAGS`` now carries it).

Every test below fails on the pre-fix literal.  ``test_..._drops_ip_flagged_sites``
is the one that matters: it drives the real seam, ``ema.main._apply_pas_filters``,
with ``args`` seeded from the schema exactly the way ``ema run`` seeds it for a
run with no flags at all, and asserts a flagged site is GONE from the BED on
disk.  With ``default="annotate"`` restored it fails with the flagged row still
present and ``filtered == 0`` -- i.e. it reproduces the defect.
"""
from __future__ import annotations

import json
import logging
from dataclasses import fields as _dc_fields
from pathlib import Path

import pytest

from ema.cli.config_schema import (
    IP_FILTER_MODES_EFFECTIVE,
    IP_FILTER_MODE_UNSET,
    IP_FILTER_MODE_V2,
    IP_FILTER_MODE_WHEN_UNSET,
    RunConfig,
    V2_COMPAT_FLAGS,
    field_specs,
    resolve_ip_filter_mode,
)
from ema.config import args, directory_config, set_directory_config, variable_config
from ema.main import _ip_filter_decision, _ip_filter_resolved, _resolve_ip_filter
from ema.outputs import OutputManager


# ---------------------------------------------------------------------------
# "the branch's DEFAULTS" -- seeded from the schema through the same bridge
# metadata ema/cli/__init__.py::cli() and RunConfig.apply_to_legacy_globals()
# use, so no value here is hand-typed and the test really does exercise the
# literal in ema/cli/config_schema.py.
# ---------------------------------------------------------------------------
def _schema_default_targets() -> list[tuple[object, str, object]]:
    """``[(holder, attr, default), ...]`` for every args/variable_config field."""
    specs = field_specs(RunConfig)
    ns = args._get()
    out: list[tuple[object, str, object]] = []
    for f in _dc_fields(RunConfig):
        spec = specs[f.name]
        if spec is None or getattr(spec, "skip_legacy_bridge", False):
            continue
        lda = spec.legacy_dataclass_attr or ""
        if lda.startswith("variable_config."):
            out.append((variable_config, lda.split(".", 1)[1], f.default))
        elif not lda:
            out.append((ns, spec.legacy_args_attr or f.name, f.default))
    return out


_MISSING = object()


@pytest.fixture
def branch_defaults():
    """Put ``args`` / ``variable_config`` in the state a no-flag run has."""
    targets = _schema_default_targets()
    saved = [(holder, attr, getattr(holder, attr, _MISSING))
             for holder, attr, _d in targets]
    for holder, attr, default in targets:
        setattr(holder, attr, default)
    _ip_filter_resolved.clear()
    try:
        yield args._get()
    finally:
        for holder, attr, old in saved:
            if old is _MISSING:
                try:
                    delattr(holder, attr)
                except AttributeError:
                    pass
            else:
                setattr(holder, attr, old)
        _ip_filter_resolved.clear()


def _genome(tmp_path: Path) -> Path:
    """200 bp of C with a 10 bp poly-A stretch at [100, 110)."""
    seq = list("C" * 200)
    seq[100:110] = list("A" * 10)
    fa = tmp_path / "genome.fa"
    fa.write_text(">chr1\n" + "".join(seq) + "\n")
    return fa


def _run_dir(tmp_path: Path) -> OutputManager:
    run_dir = tmp_path / "run"
    set_directory_config(output_dir=run_dir, gtf_dir=None)
    mgr = OutputManager(base_dir=str(run_dir))
    mgr.setup()
    return mgr


def _write_two_pas() -> None:
    """PAS '1' is clean (window all C); PAS '2' sits on the poly-A stretch."""
    with open(directory_config.posbed, "w") as fh:
        fh.write("chr1\t40\t50\t1\t0\t+\n")
        fh.write("chr1\t90\t100\t2\t0\t+\n")
    with open(directory_config.negbed, "w") as fh:
        fh.write("")


def _stats(mgr: OutputManager) -> dict:
    return json.loads(Path(mgr.path("pas_gene", "peak_filters_stats.json")).read_text())


def _surviving() -> list[str]:
    return [ln.split("\t")[3]
            for ln in directory_config.posbed.read_text().splitlines() if ln.strip()]


# ---------------------------------------------------------------------------
# 1. the literal itself
# ---------------------------------------------------------------------------
def test_the_schema_default_is_the_sentinel_not_a_mode():
    """FAILS on the pre-fix tree, where the default was the string 'annotate'."""
    assert RunConfig().ip_filter_mode == IP_FILTER_MODE_UNSET
    spec = field_specs(RunConfig)["ip_filter_mode"]
    assert IP_FILTER_MODE_UNSET in spec.choice
    assert set(IP_FILTER_MODES_EFFECTIVE) <= set(spec.choice)
    assert IP_FILTER_MODE_UNSET not in IP_FILTER_MODES_EFFECTIVE, (
        "the sentinel must never be a mode the filter itself accepts"
    )


def test_the_sentinel_resolves_to_a_mode_that_drops():
    assert resolve_ip_filter_mode(IP_FILTER_MODE_UNSET)[0] == IP_FILTER_MODE_WHEN_UNSET
    assert IP_FILTER_MODE_WHEN_UNSET == "filter"
    # ...and both real modes survive round-tripping unchanged.
    for mode in IP_FILTER_MODES_EFFECTIVE:
        assert resolve_ip_filter_mode(mode)[0] == mode
    with pytest.raises(ValueError):
        resolve_ip_filter_mode("keep-everything-please")


# ---------------------------------------------------------------------------
# 2. THE REGRESSION TEST: a run at the branch defaults must actually drop.
# ---------------------------------------------------------------------------
def test_a_branch_default_run_with_a_fasta_drops_ip_flagged_sites(
        branch_defaults, tmp_path):
    """The whole finding, at the seam that implements it.

    Nothing is set here but ``genome_fasta`` -- the one thing a real command
    line supplies beyond the input paths.  ``ip_filter``, ``no_ip_filter``,
    ``ip_filter_mode`` and ``ip_filter_default`` all come from the schema.

    On the pre-fix literal this fails twice over: ``surviving == ['1', '2']``
    and ``stats['ip_filter_mode'] == 'annotate'``.
    """
    pytest.importorskip("pyfaidx")
    from ema.main import _apply_pas_filters, _validate_pas_filter_config

    mgr = _run_dir(tmp_path)
    _write_two_pas()
    branch_defaults.genome_fasta = str(_genome(tmp_path))

    _validate_pas_filter_config()
    result = _apply_pas_filters(mgr)

    stats = _stats(mgr)
    assert stats["ip_filter"] is True, "the branch default must run the veto"
    assert stats["ip_filter_mode"] == "filter", (
        "the branch default resolved to %r: the veto ran and dropped nothing"
        % (stats["ip_filter_mode"],)
    )
    assert stats["ip_filter_mode_requested"] == IP_FILTER_MODE_UNSET
    assert _surviving() == ["1"], (
        "the internally-primed PAS survived a DEFAULT run: the branch's only "
        "behavioural default is a no-op again"
    )
    assert result is not None and result["ip_of"] == {"1": False, "2": True}


def test_an_explicit_annotate_is_still_honoured(branch_defaults, tmp_path):
    """The D9 behaviour stays reachable, and stays one flag away."""
    pytest.importorskip("pyfaidx")
    from ema.main import _apply_pas_filters

    mgr = _run_dir(tmp_path)
    _write_two_pas()
    branch_defaults.genome_fasta = str(_genome(tmp_path))
    branch_defaults.ip_filter_mode = "annotate"

    _apply_pas_filters(mgr)

    assert _stats(mgr)["ip_filter_mode"] == "annotate"
    assert _surviving() == ["1", "2"], (
        "an explicit --ip-filter-mode annotate must keep every PAS"
    )


def test_the_v2_pin_still_keeps_every_pas(branch_defaults, tmp_path):
    """``V2_COMPAT_FLAGS`` must still reproduce v2 at this seam.

    v2 ran the veto only when asked (``--ip-filter-default off``) and kept
    every flagged PAS when it did (``--ip-filter-mode annotate``).  Both are
    now on the documented compat command line, so applying it here must leave
    the BEDs byte-identical.
    """
    pytest.importorskip("pyfaidx")
    from ema.main import _apply_pas_filters

    pairs = dict(zip(V2_COMPAT_FLAGS[::2], V2_COMPAT_FLAGS[1::2]))
    assert pairs["--ip-filter-mode"] == IP_FILTER_MODE_V2, (
        "the v2 compat command no longer pins the mode; a published compat "
        "run would drop sites v2 kept"
    )
    variable_config.ip_filter_default = pairs["--ip-filter-default"]
    branch_defaults.ip_filter_mode = pairs["--ip-filter-mode"]
    branch_defaults.genome_fasta = str(_genome(tmp_path))

    mgr = _run_dir(tmp_path)
    _write_two_pas()
    before = directory_config.posbed.read_bytes()

    assert _resolve_ip_filter() is False, "v2 ran the veto only on request"
    _apply_pas_filters(mgr)
    assert directory_config.posbed.read_bytes() == before


# ---------------------------------------------------------------------------
# 3. the resolved mode must be VISIBLE -- log and run record.
# ---------------------------------------------------------------------------
def test_the_run_log_names_the_resolved_mode(branch_defaults, tmp_path, caplog):
    branch_defaults.genome_fasta = str(_genome(tmp_path))
    with caplog.at_level(logging.INFO, logger="ema.main"):
        assert _resolve_ip_filter() is True
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "mode=filter" in msgs, msgs
    assert "DROPPED" in msgs, msgs


def test_the_run_log_says_so_when_nothing_will_be_dropped(
        branch_defaults, tmp_path, caplog):
    branch_defaults.genome_fasta = str(_genome(tmp_path))
    branch_defaults.ip_filter_mode = "annotate"
    with caplog.at_level(logging.INFO, logger="ema.main"):
        assert _resolve_ip_filter() is True
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "NOTHING IS DROPPED" in msgs, msgs


def test_the_mode_is_part_of_the_one_time_notice_key(
        branch_defaults, tmp_path, caplog):
    """The notice is memoised; the memo must not swallow a mode change."""
    branch_defaults.genome_fasta = str(_genome(tmp_path))
    with caplog.at_level(logging.INFO, logger="ema.main"):
        _resolve_ip_filter()
        branch_defaults.ip_filter_mode = "annotate"
        _resolve_ip_filter()
    msgs = [r.getMessage() for r in caplog.records]
    assert any("mode=filter" in m for m in msgs), msgs
    assert any("NOTHING IS DROPPED" in m for m in msgs), msgs


def test_the_run_record_carries_the_resolved_mode(branch_defaults, tmp_path):
    """run_config.json / run_manifest.json must show what actually happened.

    ``args`` alone cannot: it carries the unresolved sentinel.  This is the
    surface on which the defect hid for a whole verification pass.
    """
    from ema.outputs import build_resolved_run_config

    _run_dir(tmp_path)
    branch_defaults.genome_fasta = str(_genome(tmp_path))
    rec = build_resolved_run_config()["internal_priming"]
    assert rec["ip_filter"] is True
    assert rec["ip_filter_mode_requested"] == IP_FILTER_MODE_UNSET
    assert rec["ip_filter_mode_resolved"] == "filter"
    assert "filter" in rec["ip_filter_mode_why"]


def test_the_run_record_reports_no_mode_when_the_veto_does_not_run(
        branch_defaults, tmp_path):
    from ema.outputs import build_resolved_run_config

    _run_dir(tmp_path)
    variable_config.ip_filter_default = "off"
    rec = build_resolved_run_config()["internal_priming"]
    assert rec["ip_filter"] is False
    assert rec["ip_filter_mode_resolved"] is None, (
        "naming a mode for a veto that never ran is the same class of untruth"
    )


# ---------------------------------------------------------------------------
# 4. the sentinel must never reach the code that acts on it.
# ---------------------------------------------------------------------------
def test_the_sentinel_never_reaches_the_filter(tmp_path):
    from ema.experimental.internal_priming import filter_internal_priming

    fa = tmp_path / "g.fa"
    fa.write_text(">1\n" + ("C" * 100 + "A" * 40 + "C" * 100) + "\n")
    bed = tmp_path / "in.bed"
    bed.write_text("1\t99\t100\tp1\t5\t+\n")
    with pytest.raises(ValueError):
        filter_internal_priming(str(bed), str(fa), str(tmp_path / "o.bed"),
                                mode=IP_FILTER_MODE_UNSET)


def test_annotate_mode_really_keeps_every_internally_primed_pas(tmp_path):
    """The other half of the claim, at the function that implements it."""
    from ema.experimental.internal_priming import filter_internal_priming

    fa = tmp_path / "g.fa"
    fa.write_text(">1\n" + ("C" * 100 + "A" * 40 + "C" * 100) + "\n")
    bed = tmp_path / "in.bed"
    bed.write_text("1\t99\t100\tp1\t5\t+\n1\t20\t21\tp2\t5\t+\n")
    out = tmp_path / "out.bed"
    stats = filter_internal_priming(str(bed), str(fa), str(out), mode="annotate")
    assert stats["flagged"] >= 1, stats
    assert stats["filtered"] == 0, stats
    assert len(out.read_text().splitlines()) == 2, (
        "annotate mode dropped a PAS; the D9 contract is that it keeps every one"
    )
    out2 = tmp_path / "out2.bed"
    stats2 = filter_internal_priming(str(bed), str(fa), str(out2), mode="filter")
    assert stats2["filtered"] == stats["flagged"] >= 1, stats2
    assert len(out2.read_text().splitlines()) == 2 - stats2["filtered"]
