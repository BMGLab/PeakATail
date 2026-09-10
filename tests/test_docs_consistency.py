"""Regression guards for documentation/code consistency (issue #71).

These tests pin the doc claims that were previously contradictory to the
actual behaviour of the code:

1. UMI handling -- ``ema/countmatrix/read.py`` reads only the ``CB`` barcode
   tag and the ``RG`` read-group tag; it never reads a ``UB``/UMI tag and
   counts raw read 3' ends. The user-facing docs must therefore NOT state
   that a ``UB``/UMI tag is *required* or *used*.
2. ``docs/cli/run.md`` must not carry the stale "Currently a no-op"
   ``--ip-filter`` / ``--genome-fasta`` rows -- those flags are wired into
   ``peakatail run`` (D9, ``ema/main.py::_apply_pas_filters``).
3. The README command table must list every ``ema`` subcommand.
4. The two user-visible behaviour changes of issue #65 -- ``peakatail reannotate``
   refusing a ``--out`` that another live run holds (hard error, exit code 1)
   and the atomic artifact write -- must be documented in
   ``docs/cli/reannotate.md`` and in ``CHANGELOG.md``'s ``Unreleased``
   section, since a user who suddenly hits the refusal has nowhere else to
   look.
"""
from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
README = REPO_ROOT / "README.md"
RUN_MD = REPO_ROOT / "docs" / "cli" / "run.md"
DATA_FLOW_MD = REPO_ROOT / "docs" / "concepts" / "data-flow.md"
READ_PY = REPO_ROOT / "ema" / "countmatrix" / "read.py"
REANNOTATE_MD = REPO_ROOT / "docs" / "cli" / "reannotate.md"
CHANGELOG = REPO_ROOT / "CHANGELOG.md"
REANNOTATE_PY = REPO_ROOT / "ema" / "reannotate.py"
SWITCH_DIFF_MD = REPO_ROOT / "docs" / "cli" / "switch-diff.md"


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
        "`peakatail run`",
        "`peakatail reannotate`",
        "`peakatail switch diff`",
        "`peakatail switch length`",
        "`peakatail switch trend`",
        "`peakatail switch combine`",
        "`peakatail switch match`",
        "`peakatail switch geneview`",
        "`peakatail collapse`",
        "`peakatail merge`",
        "`peakatail parse-gtf`",
        "`peakatail wizard`",
    ):
        assert cmd in text, f"README command table is missing {cmd}"


# --------------------------------------------------------------------- #
# Issue #65: the duplicate-``--out`` refusal and the atomic write are both
# user-visible behaviour changes, so both must be documented.
# --------------------------------------------------------------------- #

def _unreleased_section() -> str:
    """Every ``## Unreleased`` section of CHANGELOG.md, concatenated.

    The changelog convention in this repo is one themed section per change
    (``## Unreleased -- <topic>``), so several sit stacked at the top at any
    time.  Reading only the FIRST one made this helper order-dependent: the
    next merge to prepend a section silently hid every earlier entry from the
    coverage tests below.  Collect them all instead.
    """
    text = CHANGELOG.read_text()
    parts, idx = [], text.find("## Unreleased")
    while idx != -1:
        rest = text[idx + len("## Unreleased"):]
        nxt = rest.find("\n## ")
        parts.append(rest if nxt == -1 else rest[:nxt])
        nxt_abs = text.find("\n## ", idx + 1)
        if nxt_abs == -1:
            break
        idx = text.find("## Unreleased", nxt_abs)
    return "\n".join(parts)


def test_reannotate_md_documents_the_duplicate_out_refusal():
    """A user hitting the hard error must find it, and its exit code, here."""
    text = REANNOTATE_MD.read_text()
    # The lock file name is the ground truth in the code.
    lock_name = re.search(
        r'OUT_DIR_LOCK_NAME\s*=\s*["\']([^"\']+)["\']', REANNOTATE_PY.read_text()
    ).group(1)
    assert lock_name in text, (
        f"reannotate.md never names the lock file {lock_name!r} that the "
        "refusal points users at"
    )
    assert "exit code 1" in text, (
        "reannotate.md does not state the exit code of the duplicate---out refusal"
    )
    assert "refuses to start" in text, (
        "reannotate.md does not say that a duplicate --out is refused up front"
    )


def test_reannotate_md_documents_atomic_writes():
    text = REANNOTATE_MD.read_text()
    assert "os.replace" in text and "atomic" in text.lower(), (
        "reannotate.md does not document that artifacts are written atomically"
    )


def test_changelog_unreleased_covers_issue_65():
    section = _unreleased_section()
    for phrase in ("issue #65", "--out", "exits 1", "atomic"):
        assert phrase in section, (
            f"CHANGELOG.md's Unreleased section is missing {phrase!r} "
            "(CONTRIBUTING.md requires a changelog entry per change)"
        )


# --------------------------------------------------------------------- #
# Issue #94: the ``switch diff`` per-pair TSV schema on the CLI page drifted
# from the code -- it documented a ``statistic`` column no strategy produces
# and omitted most columns that are actually written.  The doc table must
# match, exactly and in order, the columns ``ema switch diff`` writes for the
# default ``fisher`` strategy.
# --------------------------------------------------------------------- #

# Lead columns prepended by ``_augment_diff_df`` in ema/switch_test/runner.py
# (per_gene scope; ``between_utr`` drops the four coordinate columns).
DIFF_LEAD_COLS = [
    "pas_id", "gene_id", "chrom", "start", "end", "strand",
    "cluster1", "cluster2",
]


def _fisher_stat_columns() -> list[str]:
    """Ground truth: the stat columns the fisher strategy actually returns."""
    import numpy as np
    import pandas as pd

    from ema.switch_test.strategies import get_diff_strategy

    rng = np.random.default_rng(0)
    cells = [f"c{i}" for i in range(40)]
    pas = ["P1", "P2", "P3", "P4"]
    count_matrix = pd.DataFrame(
        rng.integers(0, 6, size=(len(cells), len(pas))), index=cells, columns=pas
    )
    labels = pd.Series(["A"] * 20 + ["B"] * 20, index=cells)
    result = get_diff_strategy("fisher").test(
        count_matrix=count_matrix,
        cluster_labels=labels,
        cluster1="A",
        cluster2="B",
        min_cells_per_group=1,
        pas_gene_map={"P1": "G1", "P2": "G1", "P3": "G2", "P4": "G2"},
    )
    assert not result.empty, "fisher strategy returned no rows for the fixture"
    return list(result.columns)


def _documented_diff_columns() -> list[str]:
    """Column names of the per-pair TSV table in docs/cli/switch-diff.md."""
    text = SWITCH_DIFF_MD.read_text()
    anchor = text.index("**`differential/<strategy>_<c1>_vs_<c2>.tsv`**")
    header = text.index("| Column | Type | Description |", anchor)
    cols = []
    for line in text[header:].splitlines()[2:]:  # skip header + separator row
        if not line.startswith("|"):
            break
        cell = line.split("|")[1].strip()
        match = re.fullmatch(r"`([A-Za-z0-9_]+)`", cell)
        assert match, f"unparsable column cell in switch-diff.md table: {cell!r}"
        cols.append(match.group(1))
    return cols


def test_switch_diff_md_documents_the_real_tsv_columns():
    """The documented per-pair schema must equal what the code writes."""
    expected = DIFF_LEAD_COLS + _fisher_stat_columns()
    documented = _documented_diff_columns()
    assert documented == expected, (
        "docs/cli/switch-diff.md's per-pair TSV table is out of sync with "
        f"ema/switch_test/.\n  documented: {documented}\n  actual:     {expected}\n"
        f"  invented:   {sorted(set(documented) - set(expected))}\n"
        f"  missing:    {sorted(set(expected) - set(documented))}"
    )


def test_switch_diff_md_documents_the_opposite_sign_conventions():
    """``delta_proportion`` and ``log2fc`` point opposite ways -- say so."""
    text = SWITCH_DIFF_MD.read_text()
    for phrase in (
        "prop(cluster1) - prop(cluster2)",
        "log2(prop(cluster2) / prop(cluster1))",
        "opposite sign",
    ):
        assert phrase in text, (
            f"switch-diff.md never states {phrase!r}; a reader cannot tell "
            "that delta_proportion and log2fc use opposite sign conventions"
        )


def test_changelog_unreleased_covers_issue_94_doc_fix():
    section = _unreleased_section()
    for phrase in ("issue #94", "per-pair TSV", "`statistic` column",
                   "opposite sign"):
        assert phrase in section, (
            f"CHANGELOG.md's Unreleased section is missing {phrase!r} "
            "(CONTRIBUTING.md requires a changelog entry per change)"
        )
