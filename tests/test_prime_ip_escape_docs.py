"""The documented escape from the breaking default must be an escape that works.

THE DEFECT THIS FILE PINS
-------------------------
``peakAtail-prime`` flips one behavioural default: the internal-priming veto
runs, and drops, when the user types nothing but ``--genome-fasta``.  The
CHANGELOG and ``docs/cli/run.md`` named ONE affected population (the flagless
user) and offered ``--ip-filter-default off`` as the escape.

There are TWO populations, and that escape does not work for the second one.
``ema.main._ip_filter_decision()`` resolves ``forced_on`` (``--ip-filter``)
BEFORE ``policy`` (``ip_filter_default``), so a v2 user who typed
``--ip-filter`` and never named a mode -- who got ``annotate`` in v2, i.e. the
veto ran, flagged, and dropped nothing -- now gets ``filter`` and loses the
same ~17 % of calls, and adding ``--ip-filter-default off`` changes nothing for
them because their own ``--ip-filter`` wins.  Their escape is
``--ip-filter-mode annotate`` (or ``--no-ip-filter``, or the full
``V2_COMPAT_FLAGS`` line).

The precedence itself is deliberate and is NOT changed here: ``--ip-filter``
means "run the veto", ``--ip-filter-default`` only decides what happens when
nobody said, and ``--ip-filter-mode`` is the single flag that decides whether
the veto DROPS.  Making an explicit ``--ip-filter`` imply ``annotate`` would
rebuild the exact trap this branch removed -- two flags silently disagreeing
about whether anything is dropped -- and would break the v2 compat golden.
What changes is the documentation, so the tests below are half behaviour (the
precedence, pinned so the prose stays true) and half prose (the two documents
must name both populations and the right escape for each).

The second finding pinned here: the CHANGELOG claimed the no-FASTA flagless
output "stays byte-identical to v2's".  The BEDs and the count matrix are;
``pas_support.tsv`` is NOT, because ``--pas-features on`` and
``--emit-inferred-cleavage on`` are branch defaults that append columns to it
whatever the FASTA situation is.  (The executable half of that statement is
``tests/test_prime_v2_compat_golden.py::test_branch_defaults_still_write_v2s_beds_and_matrix``.)

Third, non-blocking: the ``mode=filter`` notice was ``log.info`` while the two
LESS consequential notices next to it were warnings, so ``--quiet`` hid the
only one that removes calls.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import pytest

from ema.config import args, variable_config
from ema.main import _ip_filter_decision, _ip_filter_resolved, _resolve_ip_filter

REPO = Path(__file__).resolve().parents[1]
CHANGELOG = REPO / "CHANGELOG.md"
RUN_DOC = REPO / "docs" / "cli" / "run.md"

_MISSING = object()


def _flat(text: str) -> str:
    """Markdown emphasis and line wrapping removed, so prose checks are not
    hostage to where a paragraph happens to wrap."""
    return re.sub(r"\s+", " ", text.replace("*", "").replace("\n", " "))


@pytest.fixture
def ip_args(tmp_path):
    """``args`` / ``variable_config`` seeded as the branch defaults leave them."""
    ns = args._get()
    keys = ("ip_filter", "no_ip_filter", "ip_filter_mode", "genome_fasta")
    saved = {k: getattr(ns, k, _MISSING) for k in keys}
    saved_default = getattr(variable_config, "ip_filter_default", _MISSING)
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "C" * 200 + "\n")
    ns.ip_filter = False
    ns.no_ip_filter = False
    ns.ip_filter_mode = "auto"
    ns.genome_fasta = str(fasta)
    variable_config.ip_filter_default = "auto"
    _ip_filter_resolved.clear()
    try:
        yield ns
    finally:
        for k, v in saved.items():
            if v is _MISSING:
                try:
                    delattr(ns, k)
                except AttributeError:
                    pass
            else:
                setattr(ns, k, v)
        if saved_default is not _MISSING:
            variable_config.ip_filter_default = saved_default
        _ip_filter_resolved.clear()


# ---------------------------------------------------------------------------
# 1. the behaviour the documents have to describe
# ---------------------------------------------------------------------------
def test_the_flagless_user_with_a_fasta_gets_the_veto_and_the_drop(ip_args):
    d = _ip_filter_decision()
    assert (d["ip_filter"], d["mode"]) == (True, "filter")


def test_ip_filter_default_off_is_the_escape_for_the_FLAGLESS_user(ip_args):
    variable_config.ip_filter_default = "off"
    assert _ip_filter_decision()["ip_filter"] is False


def test_ip_filter_default_off_is_NOT_an_escape_for_the_explicit_ip_filter_user(ip_args):
    """The precedence the documented escape was wrong about.

    ``--ip-filter`` beats ``ip_filter_default``, so this user still runs the
    veto -- and, with no ``--ip-filter-mode``, still DROPS, where v2 annotated.
    """
    ip_args.ip_filter = True
    variable_config.ip_filter_default = "off"
    d = _ip_filter_decision()
    assert d["ip_filter"] is True, "--ip-filter no longer beats the policy"
    assert d["why"] == "--ip-filter"
    assert d["mode"] == "filter", (
        "if this is now 'annotate', the precedence changed and the CHANGELOG / "
        "docs / V2_COMPAT_FLAGS prose in this branch describe the old rule"
    )


def test_ip_filter_mode_annotate_IS_the_escape_for_that_user(ip_args):
    """v2's behaviour for the explicit-``--ip-filter`` user, exactly: veto runs,
    flags recorded, nothing dropped."""
    ip_args.ip_filter = True
    ip_args.ip_filter_mode = "annotate"
    d = _ip_filter_decision()
    assert (d["ip_filter"], d["mode"]) == (True, "annotate")


def test_no_ip_filter_beats_an_explicit_ip_filter(ip_args):
    ip_args.ip_filter = True
    ip_args.no_ip_filter = True
    d = _ip_filter_decision()
    assert d["ip_filter"] is False and d["mode"] is None


def test_the_notice_that_removes_calls_is_a_warning_not_info(ip_args, caplog):
    """``--quiet`` must not be able to hide the only notice that drops 17 % of
    the call set, while the two less consequential ones are already warnings."""
    with caplog.at_level(logging.INFO, logger="ema.main"):
        assert _resolve_ip_filter() is True
    dropping = [r for r in caplog.records if "mode=filter" in r.getMessage()]
    assert dropping, [r.getMessage() for r in caplog.records]
    assert all(r.levelno >= logging.WARNING for r in dropping), (
        "the mode=filter notice is emitted at %r; it removes 17.09 %% of the "
        "call set and must survive --quiet"
        % [logging.getLevelName(r.levelno) for r in dropping]
    )


# ---------------------------------------------------------------------------
# 2. the documents
# ---------------------------------------------------------------------------
def _breaking_section() -> str:
    text = CHANGELOG.read_text()
    m = re.search(r"^### BREAKING.*?(?=^#{2,3} )", text, re.S | re.M)
    assert m, "the CHANGELOG no longer has a BREAKING section for the default flip"
    return m.group(0)


def test_the_changelog_breaking_section_names_both_populations():
    flat = _flat(_breaking_section())
    assert "--genome-fasta" in flat, flat
    assert "`--ip-filter`" in flat, (
        "the BREAKING section names only the flagless population; the user who "
        "typed --ip-filter without --ip-filter-mode loses the same ~17 % of "
        "calls and is not mentioned"
    )
    assert "`--ip-filter-mode annotate`" in flat, (
        "the escape that actually works for the explicit --ip-filter user is "
        "not offered"
    )
    assert "`--no-ip-filter`" in flat and "V2_COMPAT_FLAGS" in flat, flat
    assert re.search(r"`--ip-filter-default off` does NOT", flat), (
        "the BREAKING section still offers --ip-filter-default off without "
        "saying that it does NOT restore v2 for a user who typed --ip-filter"
    )


def test_the_changelog_no_fasta_claim_is_scoped_to_the_beds_and_the_matrix():
    """`pas_support.tsv` is not byte-identical to v2's at the branch defaults."""
    section = _breaking_section()
    para = [p for p in section.split("\n\n") if "no readable `--genome-fasta`" in p]
    assert para, "the BREAKING section no longer states the no-FASTA case"
    flat = _flat(para[-1])
    assert "pas_support.tsv" in flat, (
        "the no-FASTA paragraph claims byte-identity without excluding "
        "pas_support.tsv, which --pas-features on / --emit-inferred-cleavage "
        "on append columns to at the branch defaults: %s" % flat
    )
    assert "count matrix" in flat or "matrix" in flat, flat


def _d9_section() -> str:
    text = RUN_DOC.read_text()
    m = re.search(r"^### Internal-priming annotation \(D9\).*?(?=^#{2,3} )",
                  text, re.S | re.M)
    assert m, "docs/cli/run.md no longer has the D9 internal-priming section"
    return m.group(0)


def test_the_run_docs_name_both_populations_and_the_right_escape():
    flat = _flat(_d9_section())
    assert "`--ip-filter-mode annotate`" in flat
    assert re.search(r"`--ip-filter-default off` does not help", flat, re.I), (
        "docs/cli/run.md still presents --ip-filter-default off as the escape "
        "without saying it is ignored by a run that typed --ip-filter"
    )
    assert "pas_support.tsv" in flat, (
        "docs/cli/run.md claims the no-FASTA output is byte-identical without "
        "excluding the sidecar"
    )
    assert "the internal-priming filter is off unless enabled" not in flat, (
        "docs/cli/run.md still describes v2's defaults as this branch's"
    )
