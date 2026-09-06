"""Regression tests for `switch diff --isoform-agg within_utr`.

Both defects covered here shipped in a released version and neither was caught,
because no test exercised `within_utr` or `between_utr` at all:

1. `_build_diff_isoform_groups` read `gene_fallback_map`, a name bound only in
   `run_length`, so every `within_utr` invocation raised NameError while building
   groups -- before any statistical test ran.

2. Once (1) was fixed, a spliced 3'UTR (several `three_prime_utr` records for one
   transcript) appended the same PAS to one UTR group repeatedly. The duplicate
   columns reach the fisher strategy as `int(agg1[p])` on a Series, raising
   TypeError.

Both tests fail on the unfixed source and pass on the fixed source. Defect 1 is
checked by static analysis rather than a functional run because reaching that line
needs a real matrix, a GTF and a pasbed, and a NameError there is invisible to any
test that stops short of them.
"""

from __future__ import annotations

import ast
import builtins
import inspect
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from ema.switch_test import runner


def _free_names(fn) -> set[str]:
    """Names the function LOADS that it never binds and that are not globals."""
    tree = ast.parse(textwrap.dedent(inspect.getsource(fn)))
    fdef = tree.body[0]

    bound: set[str] = set()
    for node in ast.walk(fdef):
        # parameters of this function AND of any nested def/lambda
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
            a = node.args
            for group in (a.args, a.posonlyargs, a.kwonlyargs):
                bound |= {x.arg for x in group}
            for extra in (a.vararg, a.kwarg):
                if extra is not None:
                    bound.add(extra.arg)
    for node in ast.walk(fdef):
        if isinstance(node, ast.Name) and isinstance(node.ctx, (ast.Store, ast.Del)):
            bound.add(node.id)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            bound.add(node.name)
        elif isinstance(node, ast.comprehension):
            for t in ast.walk(node.target):
                if isinstance(t, ast.Name):
                    bound.add(t.id)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            for al in node.names:
                bound.add((al.asname or al.name).split(".")[0])
        elif isinstance(node, ast.ExceptHandler) and node.name:
            bound.add(node.name)

    module_globals = set(vars(runner)) | set(dir(builtins))
    loads = {n.id for n in ast.walk(fdef)
             if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
    return loads - bound - module_globals


def test_within_utr_group_builder_has_no_unbound_name():
    """Defect 1. `gene_fallback_map` was read here but bound only in run_length."""
    free = _free_names(runner._build_diff_isoform_groups)
    assert not free, (
        "_build_diff_isoform_groups reads names it never binds: "
        f"{sorted(free)}. Every --isoform-agg within_utr run dies with NameError "
        "at group construction."
    )


def test_run_diff_entrypoint_has_no_unbound_name():
    """The same class of defect on the function that calls it."""
    free = _free_names(runner.run_diff)
    assert not free, f"run_diff reads names it never binds: {sorted(free)}"


def test_spliced_utr_does_not_duplicate_a_pas_in_its_group():
    """Defect 2, on the shipped source rather than a copy of its logic.

    A transcript whose 3'UTR spans two exons yields two entries naming the same
    (gene, transcript) for one PAS; the group must still hold that PAS once.
    """
    src = inspect.getsource(runner._build_diff_isoform_groups)
    start = src.index("for pas_id, entries in pas_isoform_map_raw.items():")
    end = src.index("orphans_by_gene", start)
    loop = textwrap.dedent(src[start:end])

    utr_pas_members: dict[str, list] = {}
    gene_of_utr: dict[str, str] = {}
    ns = {
        "pas_isoform_map_raw": {
            "PAS1": [("G1", "T1", 0, 1, 2), ("G1", "T1", 1, 1, 2)],  # spliced 3'UTR
            "PAS2": [("G1", "T1", 0, 2, 2)],
        },
        "diff_cols": {"PAS1", "PAS2"},
        "utr_pas_members": utr_pas_members,
        "gene_of_utr": gene_of_utr,
    }
    exec(loop, ns)

    members = utr_pas_members["G1::T1"]
    assert members == ["PAS1", "PAS2"], (
        f"UTR group holds {members}; a spliced 3'UTR must not append a PAS twice. "
        "Duplicate columns raise TypeError in the fisher strategy."
    )
