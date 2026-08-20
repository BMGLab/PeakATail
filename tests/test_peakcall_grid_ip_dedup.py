"""Issue #69 regression guard: the Laughney peak-call grid must not contain
two arms that produce a byte-identical ``pasbed.bed``.

Background
----------
``experiments/laughney/main.nf::GRID_RUN`` turns each row of
``peakcall_grid.tsv`` into one ``ema run``. The internal-priming CLI it emits
is (see main.nf)::

    ip_mode == 'off'   -> (no --ip-filter flag)                # keep-all
    ip_mode == 'annotate' -> --ip-filter --ip-filter-mode annotate  # keep-all
    ip_mode == 'filter'   -> --ip-filter --ip-filter-mode filter     # DROPS PAS

``annotate`` mode never drops or edits a ``pasbed.bed`` row -- the
``internal_priming`` flag is written only to ``annotatedpas.bed`` (see
``ema.experimental.internal_priming.filter_internal_priming`` and
``tests/test_pasbed_status_columns_d9.py``). So at the ``pasbed.bed`` /
``benchmark_vs_polyasite_v3.json`` level, ``annotate`` and ``off`` are the SAME
run. Grid arm ``lg_ip_off`` (ip_mode=off) was therefore byte-identical to
``lg_annotate`` (ip_mode=annotate) and got removed; this test stops it (or any
equivalent duplicate) from silently coming back.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

GRID = (
    Path(__file__).resolve().parents[1]
    / "experiments" / "laughney" / "peakcall_grid.tsv"
)

# Only ip_mode=='filter' changes pasbed.bed; every other value is keep-all.
_IP_DROPS = "filter"


def _active_rows() -> list[dict]:
    if not GRID.exists():  # pragma: no cover - grid is an experiment artifact
        pytest.skip(f"peakcall_grid.tsv not present at {GRID}")
    rows: list[dict] = []
    with open(GRID) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            name = (row.get("run_name") or "").strip()
            # main.nf skips '#'-prefixed and blank run_names in both parsers.
            if not name or name.startswith("#"):
                continue
            rows.append(row)
    return rows


def _pasbed_identity(row: dict) -> tuple[str, str, bool]:
    """Everything that determines the bytes of pasbed.bed for a grid arm."""
    return (
        row["peak_strategy"].strip(),
        row["atlas_mode"].strip(),
        row["ip_mode"].strip() == _IP_DROPS,  # ip_drops: only 'filter' differs
    )


def test_grid_has_rows():
    assert _active_rows(), "peakcall_grid.tsv has no active (non-comment) rows"


def test_no_two_arms_produce_identical_pasbed():
    """The exact #69 defect: two arms whose (strategy, atlas, ip_drops) match
    write byte-identical pasbed.bed and thus an identical benchmark JSON."""
    seen: dict[tuple[str, str, bool], str] = {}
    for row in _active_rows():
        ident = _pasbed_identity(row)
        assert ident not in seen, (
            f"grid arms {seen[ident]!r} and {row['run_name']!r} resolve to the "
            f"same pasbed.bed identity {ident} (peak_strategy, atlas_mode, "
            f"ip_drops) -- they will produce byte-identical output (issue #69). "
            f"'annotate' and 'off' ip_mode are both keep-all; only 'filter' "
            f"changes pasbed.bed."
        )
        seen[ident] = row["run_name"]


def test_no_ip_off_arm_shadows_a_keep_all_arm():
    """An explicit ip_mode='off' arm is only meaningful if no keep-all
    (annotate/off) arm already exists for the same strategy+atlas."""
    keep_all: dict[tuple[str, str], list[str]] = {}
    for row in _active_rows():
        if row["ip_mode"].strip() != _IP_DROPS:  # keep-all arm
            key = (row["peak_strategy"].strip(), row["atlas_mode"].strip())
            keep_all.setdefault(key, []).append(row["run_name"])
    dupes = {k: v for k, v in keep_all.items() if len(v) > 1}
    assert not dupes, (
        f"multiple keep-all IP arms for the same (peak_strategy, atlas_mode): "
        f"{dupes} -- these differ only in ip_mode annotate-vs-off, which does "
        f"not change pasbed.bed (issue #69)."
    )
