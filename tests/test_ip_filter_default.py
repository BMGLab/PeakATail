"""peakAtail-prime TASK E item 3 — the internal-priming veto turns itself on.

The veto is the single largest measured accuracy lift in the caller: +7.7 % to
+12.8 % relative recall at matched atlas precision and +17.7 % to +22.9 % at
matched long-read precision (results/algo_headroom/VERIFY/tables/
v8_ipveto_value.tsv).  In v2 it was an opt-in flag, so every benchmark arm in
the manuscript ran with it and every user who did not read the flag list did
not.  On this branch it runs whenever a genome FASTA is available -- and says
LOUDLY when there is none, because that is the case where the run silently
loses the lift.
"""
from __future__ import annotations

import logging

import pytest

from ema.config import args, variable_config
from ema.main import _resolve_ip_filter, _ip_filter_resolved


@pytest.fixture
def env(tmp_path):
    """A namespace + a real FASTA path, restored afterwards."""
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")
    ns = args._get()
    saved = {k: getattr(ns, k, None)
             for k in ("ip_filter", "no_ip_filter", "genome_fasta")}
    saved_policy = getattr(variable_config, "ip_filter_default", "auto")
    ns.ip_filter = False
    ns.no_ip_filter = False
    ns.genome_fasta = None
    _ip_filter_resolved.clear()
    yield {"ns": ns, "fasta": str(fasta), "missing": str(tmp_path / "nope.fa")}
    for k, v in saved.items():
        setattr(ns, k, v)
    variable_config.ip_filter_default = saved_policy
    _ip_filter_resolved.clear()


def test_auto_turns_it_on_when_a_fasta_is_available(env):
    variable_config.ip_filter_default = "auto"
    env["ns"].genome_fasta = env["fasta"]
    assert _resolve_ip_filter() is True


def test_auto_cannot_turn_it_on_without_a_fasta(env):
    variable_config.ip_filter_default = "auto"
    env["ns"].genome_fasta = None
    assert _resolve_ip_filter() is False


def test_a_fasta_path_that_does_not_exist_is_not_a_fasta(env):
    variable_config.ip_filter_default = "auto"
    env["ns"].genome_fasta = env["missing"]
    assert _resolve_ip_filter() is False


def test_the_v2_policy_needs_the_flag_even_with_a_fasta(env):
    variable_config.ip_filter_default = "off"
    env["ns"].genome_fasta = env["fasta"]
    assert _resolve_ip_filter() is False
    env["ns"].ip_filter = True
    _ip_filter_resolved.clear()
    assert _resolve_ip_filter() is True


def test_no_ip_filter_wins_over_everything(env):
    variable_config.ip_filter_default = "auto"
    env["ns"].genome_fasta = env["fasta"]
    env["ns"].ip_filter = True
    env["ns"].no_ip_filter = True
    assert _resolve_ip_filter() is False


def test_it_shouts_when_it_cannot_run_for_lack_of_a_fasta(env, caplog):
    variable_config.ip_filter_default = "auto"
    env["ns"].genome_fasta = None
    with caplog.at_level(logging.WARNING):
        assert _resolve_ip_filter() is False
    msgs = [r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING]
    assert any("INTERNAL-PRIMING FILTER CANNOT RUN" in m for m in msgs), msgs
    # the warning has to name what it costs, or it is not a loud warning
    assert any("relative recall" in m for m in msgs), msgs


def test_it_does_not_shout_when_the_user_said_no(env, caplog):
    variable_config.ip_filter_default = "auto"
    env["ns"].no_ip_filter = True
    with caplog.at_level(logging.WARNING):
        assert _resolve_ip_filter() is False
    assert not [r for r in caplog.records
                if "INTERNAL-PRIMING FILTER CANNOT RUN" in r.getMessage()]


def test_it_does_not_shout_under_the_v2_policy(env, caplog):
    variable_config.ip_filter_default = "off"
    env["ns"].genome_fasta = None
    with caplog.at_level(logging.WARNING):
        assert _resolve_ip_filter() is False
    assert not [r for r in caplog.records
                if "INTERNAL-PRIMING FILTER CANNOT RUN" in r.getMessage()]


def test_the_notice_is_emitted_once_per_resolution_not_per_call(env, caplog):
    variable_config.ip_filter_default = "auto"
    env["ns"].genome_fasta = None
    with caplog.at_level(logging.WARNING):
        for _ in range(5):
            _resolve_ip_filter()
    hits = [r for r in caplog.records
            if "INTERNAL-PRIMING FILTER CANNOT RUN" in r.getMessage()]
    assert len(hits) == 1, "a per-call warning would drown the log"


def test_explicit_ip_filter_still_validates_the_fasta(env):
    """--ip-filter must remain an ERROR without a FASTA; only the DEFAULT is
    allowed to degrade quietly to v2 behaviour."""
    from ema.main import _validate_pas_filter_config

    variable_config.ip_filter_default = "auto"
    env["ns"].ip_filter = True
    env["ns"].genome_fasta = None
    env["ns"].annot_filter = False
    with pytest.raises(ValueError, match="--ip-filter is enabled"):
        _validate_pas_filter_config()


def test_the_default_that_ships_is_auto():
    from ema.cli.config_schema import RunConfig
    import dataclasses

    schema = {f.name: f.default for f in dataclasses.fields(RunConfig)}
    assert schema["ip_filter_default"] == "auto"
    assert schema["ip_filter_default"] == variable_config.ip_filter_default, (
        "one default, one place: RunConfig and variable_config must agree"
    )
    assert schema["no_ip_filter"] is False and schema["ip_filter"] is False
