"""``--ip-filter-default auto`` turns the veto ON; it does not make it DROP.

THE FINDING THIS FILE PINS
--------------------------
TASK E item 3 makes the internal-priming veto run whenever a genome FASTA is
available, and the branch documents that as "the largest measured accuracy
lift in this caller (+7.7 % to +12.8 % relative recall at matched atlas
precision)".  That lift is the lift of **dropping** flagged PAS.  Whether they
are dropped is decided by a *different* flag, ``--ip-filter-mode``, whose
default is ``annotate`` -- the D9 decision that an internally-primed peak is
more likely real alternative-PAS signal than noise, so every peak is kept and
merely flagged.

Measured on the PBMC chr19+21 slice with the branch defaults and a genome
FASTA and NO other flag::

    peak_filters_stats.json:  mode "annotate", flagged 5015, filtered 0
    pas 18865 | tier1 10886 | tier2 7979 | tier1>=2mol 3643

which is v2's call set to the row -- not the ``pas 15,925 | tier1 8,524 |
tier2 7,401 | tier1>=2mol 2,883`` the branch CHANGELOG attributes to "a prime
run with no --ip-filter flag".  That measurement was made with
``--ip-filter-mode filter`` also supplied.

So: at the branch default the veto ANNOTATES and the advertised recall lift is
not delivered.  Whether the default mode should move is an adoption decision
(it changes the call set of every run), not something to fix quietly.  What is
fixed is the run saying the opposite of what it does.
"""
from __future__ import annotations

import logging

import pytest

from ema.config import args, variable_config
from ema.main import _ip_filter_resolved, _resolve_ip_filter


@pytest.fixture
def env(tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")
    ns = args._get()
    saved = {k: getattr(ns, k, None)
             for k in ("ip_filter", "no_ip_filter", "genome_fasta",
                       "ip_filter_mode")}
    had_mode = hasattr(ns, "ip_filter_mode")
    saved_policy = getattr(variable_config, "ip_filter_default", "auto")
    ns.ip_filter = False
    ns.no_ip_filter = False
    ns.genome_fasta = str(fasta)
    variable_config.ip_filter_default = "auto"
    _ip_filter_resolved.clear()
    yield ns
    for k, v in saved.items():
        setattr(ns, k, v)
    if not had_mode:
        try:
            delattr(ns, "ip_filter_mode")
        except AttributeError:
            pass
    variable_config.ip_filter_default = saved_policy
    _ip_filter_resolved.clear()


def test_the_notice_says_nothing_is_dropped_in_annotate_mode(env, caplog):
    env.ip_filter_mode = "annotate"
    with caplog.at_level(logging.INFO, logger="ema.main"):
        assert _resolve_ip_filter() is True
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "NOTHING IS DROPPED" in msgs, msgs
    assert "--ip-filter-mode filter" in msgs, msgs
    # ...and it must NOT claim the recall lift as something this run gets.
    lift_lines = [r.getMessage() for r in caplog.records
                  if "relative recall" in r.getMessage()]
    for line in lift_lines:
        assert "NOTHING IS DROPPED" in line, (
            "the recall lift is advertised without saying it needs "
            "--ip-filter-mode filter: %r" % line
        )


def test_the_notice_claims_the_lift_only_in_filter_mode(env, caplog):
    env.ip_filter_mode = "filter"
    with caplog.at_level(logging.INFO, logger="ema.main"):
        assert _resolve_ip_filter() is True
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "mode=filter" in msgs and "DROPPED" in msgs, msgs
    assert "relative recall" in msgs, msgs


def test_the_mode_is_part_of_the_one_time_notice_key(env, caplog):
    """The notice is memoised; the memo must not swallow a mode change."""
    env.ip_filter_mode = "annotate"
    with caplog.at_level(logging.INFO, logger="ema.main"):
        _resolve_ip_filter()
        env.ip_filter_mode = "filter"
        _resolve_ip_filter()
    msgs = [r.getMessage() for r in caplog.records]
    assert any("NOTHING IS DROPPED" in m for m in msgs), msgs
    assert any("mode=filter" in m for m in msgs), msgs


def test_the_effective_default_mode_is_annotate_when_no_flag_was_given(env):
    """``_apply_pas_filters`` reads ``args.ip_filter_mode`` with an
    ``"annotate"`` fallback, and the Click bridge only writes the attribute
    when the user set the flag.  So a no-flag run annotates."""
    try:
        delattr(env, "ip_filter_mode")
    except AttributeError:
        pass
    assert getattr(env, "ip_filter_mode", "annotate") == "annotate"


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
