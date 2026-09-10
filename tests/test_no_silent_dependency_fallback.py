"""A missing dependency must never silently change the science.

`ema/main.py` states the contract: "Never silently skips a requested filter."
Two paths broke it by catching ImportError and degrading:

* `filter_internal_priming` logged, copied the input through **unfiltered**,
  and returned ``{"error": ...}`` -- a key no caller reads. The run exited 0
  with a call set roughly 17 % larger than intended, and the internal-priming
  veto is on by default, so this was the default path.
* `estimate_cleavage_offset` substituted a constant for the data-driven
  estimate, shifting every reported PAS 3' end.

Both now raise. A broken install is not a data condition.
"""
from __future__ import annotations

import builtins

import pytest


def _hide(monkeypatch, name: str) -> None:
    """Make `import <name>` raise ImportError, as a broken install would."""
    real_import = builtins.__import__

    def fake(mod, globals=None, locals=None, fromlist=(), level=0):
        if mod == name or mod.startswith(name + "."):
            raise ImportError(f"No module named {name!r} (simulated)")
        return real_import(mod, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake)


def test_internal_priming_refuses_to_run_without_pyfaidx(tmp_path, monkeypatch):
    from ema.experimental.internal_priming import filter_internal_priming

    bed = tmp_path / "pas.bed"
    bed.write_text("chr1\t100\t102\tP1\t0\t+\n", encoding="utf-8")
    out = tmp_path / "out.bed"

    _hide(monkeypatch, "pyfaidx")
    with pytest.raises(RuntimeError, match="internal-priming"):
        filter_internal_priming(
            bed_path=str(bed), genome_fasta=str(tmp_path / "g.fa"),
            output_path=str(out), mode="filter",
        )
    assert not out.exists(), (
        "the filter wrote an output file despite being unable to filter; the "
        "old behaviour copied the input through unfiltered"
    )


def test_internal_priming_error_mentions_how_to_proceed(tmp_path, monkeypatch):
    """An error a user cannot act on just moves the problem."""
    from ema.experimental.internal_priming import filter_internal_priming

    bed = tmp_path / "pas.bed"
    bed.write_text("chr1\t100\t102\tP1\t0\t+\n", encoding="utf-8")

    _hide(monkeypatch, "pyfaidx")
    with pytest.raises(RuntimeError) as excinfo:
        filter_internal_priming(
            bed_path=str(bed), genome_fasta=str(tmp_path / "g.fa"),
            output_path=str(tmp_path / "o.bed"), mode="filter",
        )
    message = str(excinfo.value)
    assert "pyfaidx" in message
    assert "--no-ip-filter" in message, "the escape hatch must be named"


def test_auto_cleavage_offset_refuses_to_substitute_a_constant(tmp_path, monkeypatch):
    from ema.countmatrix import cleavage_offset as mod

    fasta = tmp_path / "g.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    bed = tmp_path / "p.bed"
    bed.write_text("chr1\t1\t2\tP1\t0\t+\n", encoding="utf-8")

    _hide(monkeypatch, "pyfaidx")
    with pytest.raises(RuntimeError, match="auto-cleavage-offset"):
        mod.resolve_cleavage_offset(
            0, auto=True, bed_paths=[str(bed)], genome_fasta=str(fasta),
        )


def test_auto_cleavage_offset_refuses_a_missing_genome_fasta(tmp_path):
    """Opt-in feature + absent input = misconfiguration, not a fallback."""
    from ema.countmatrix import cleavage_offset as mod

    with pytest.raises(FileNotFoundError, match="auto-cleavage-offset"):
        mod.resolve_cleavage_offset(
            0, auto=True, bed_paths=[], genome_fasta=str(tmp_path / "absent.fa"),
        )


def test_manual_cleavage_offset_is_unaffected(tmp_path):
    """auto=False is the default path and must keep working untouched."""
    from ema.countmatrix import cleavage_offset as mod

    offset, diag = mod.resolve_cleavage_offset(73, auto=False)
    assert offset == 73 and diag is None
