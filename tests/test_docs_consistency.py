"""Regression guards for documentation/code consistency (issue #71).

These tests pin the doc claims that were previously contradictory to the
actual behaviour of the code:

1. UMI handling -- ``ema/countmatrix/read.py`` reads only the ``CB`` barcode
   tag and the ``RG`` read-group tag; it never reads a ``UB``/UMI tag and
   counts raw read 3' ends. The user-facing docs must therefore NOT state
   that a ``UB``/UMI tag is *required* or *used*.
2. ``docs/cli/run.md`` must not carry the stale "Currently a no-op"
   ``--ip-filter`` / ``--genome-fasta`` rows -- those flags are wired into
   ``ema run`` (D9, ``ema/main.py::_apply_pas_filters``).
3. The README command table must list every ``ema`` subcommand.
"""
from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
README = REPO_ROOT / "README.md"
RUN_MD = REPO_ROOT / "docs" / "cli" / "run.md"
DATA_FLOW_MD = REPO_ROOT / "docs" / "concepts" / "data-flow.md"
READ_PY = REPO_ROOT / "ema" / "countmatrix" / "read.py"


def test_read_parser_does_not_read_umi_tag():
    """Ground truth: the hot-path read parser never touches a UB/UMI tag."""
    src = READ_PY.read_text()
    tags = re.findall(r"get_tag\(\s*['\"]([A-Za-z0-9]+)['\"]\s*\)", src)
    # ``barcode`` is a variable (defaults to CB) resolved at call time, so the
    # only string-literal tag read directly is the RG read-group tag.
    assert "UB" not in tags, f"read.py unexpectedly reads a UMI tag: {tags!r}"


def test_docs_do_not_claim_umi_is_required():
    """No user-facing doc may state a UB/UMI tag is required or used."""
    forbidden = [
        "must carry `CB:Z` (corrected barcode) and `UB:Z`",
        "must have `CB:Z` (corrected barcode) and `UB:Z`",
        "carries `CB:Z` (corrected cell barcode) and `UB:Z`",
        "cell-barcode (`CB`) and UMI (`UB`) tags",
    ]
    for md in (README, REPO_ROOT / "docs" / "index.md"):
        text = md.read_text()
        for phrase in forbidden:
            assert phrase not in text, (
                f"{md.name} still claims UMI tags are required: {phrase!r}"
            )


def test_data_flow_states_umis_unused():
    """data-flow.md is the correct doc and must keep the truthful claim."""
    text = DATA_FLOW_MD.read_text()
    assert "UMI tags are not used" in text


def test_run_md_has_no_stale_no_op_ip_filter_rows():
    text = RUN_MD.read_text()
    assert "Currently a no-op" not in text, (
        "run.md still contains stale 'Currently a no-op' ip-filter rows"
    )
    # The correct D9 annotate-mode description must remain.
    assert "Enable the internal-priming check" in text


def test_readme_lists_all_cli_commands():
    text = README.read_text()
    for cmd in (
        "`ema run`",
        "`ema reannotate`",
        "`ema switch diff`",
        "`ema switch length`",
        "`ema switch trend`",
        "`ema switch combine`",
        "`ema switch match`",
        "`ema switch geneview`",
        "`ema collapse`",
        "`ema merge`",
        "`ema parse-gtf`",
        "`ema wizard`",
    ):
        assert cmd in text, f"README command table is missing {cmd}"
