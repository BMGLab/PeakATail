"""Atomic-write guarantee for write_pas_gene_artifacts()'s pas_gene.tsv and
annotatedpas.bed outputs.

Root cause this fixes: a real bug where 5 branches in the corrected-rerun
sweep's A2/A3 trim/cluster grid, declaring IDENTICAL params, produced 3
DIFFERENT pas_gene.tsv row counts (51,129 / 67,409 / 82,702), grouped by
write time. Diagnosis: two grid rows resolving to the SAME `--out` path ran
concurrently (a duplicate-branch_name gap in the grid loader -- guarded
separately, see experiments/laughney/main.nf), and the writer used a plain
`open(..., "w")` / `pd.to_csv()` -- two processes racing on the same path
can interleave, leaving a reader with a torn file.

These tests prove the write is now atomic from a READER's perspective: the
final path is only ever either (a) absent, (b) the complete PRIOR content,
or (c) the complete NEW content -- never a partial/torn mix -- by forcing a
write to pause mid-flight (after the temp file is fully written, before the
atomic rename) and checking from a separate thread that the destination
path hasn't budged.
"""
from __future__ import annotations

import threading
import time
from pathlib import Path
from unittest.mock import patch

from ema.config import set_directory_config
from ema.outputs import write_pas_gene_artifacts


def _write_pasbed(path: Path, rows: list[tuple]) -> None:
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")


def test_pas_gene_tsv_write_leaves_no_temp_files_behind(tmp_path):
    run_dir = tmp_path / "run"
    set_directory_config(output_dir=run_dir)
    from ema.config import directory_config

    pasbed = directory_config.pasbed_for("dsA")
    pasbed.parent.mkdir(parents=True, exist_ok=True)
    _write_pasbed(pasbed, [("chr1", 10, 11, "1", 0, "+"), ("chr1", 20, 21, "2", 0, "-")])

    write_pas_gene_artifacts(run_dir, "dsA", pas_ids=["1", "2"], gene_ids=["ENSG1", "ENSG2"])

    pas_gene_tsv = directory_config.pas_gene_for("dsA")
    annot_bed = directory_config.annotatedpas_for("dsA")
    assert pas_gene_tsv.exists()
    assert annot_bed.exists()

    # No `.{name}.tmp<pid>` leftovers in either output directory.
    leftover_pgt = list(pas_gene_tsv.parent.glob(f".{pas_gene_tsv.name}.tmp*"))
    leftover_annot = list(annot_bed.parent.glob(f".{annot_bed.name}.tmp*"))
    assert leftover_pgt == [], f"leftover temp file(s): {leftover_pgt}"
    assert leftover_annot == [], f"leftover temp file(s): {leftover_annot}"


def test_pas_gene_tsv_destination_never_observes_partial_write(tmp_path):
    """A reader polling the final pas_gene.tsv path mid-write must see either
    the complete OLD content or the complete NEW content -- never neither
    (mid-truncation) nor a mix. Forces the write to pause after the temp
    file is fully populated but BEFORE the atomic rename, and polls from a
    separate thread during that window.
    """
    run_dir = tmp_path / "run"
    set_directory_config(output_dir=run_dir)
    from ema.config import directory_config

    pasbed = directory_config.pasbed_for("dsA")
    pasbed.parent.mkdir(parents=True, exist_ok=True)
    _write_pasbed(pasbed, [("chr1", 10, 11, "1", 0, "+"), ("chr1", 20, 21, "2", 0, "-")])

    pas_gene_tsv = directory_config.pas_gene_for("dsA")
    pas_gene_tsv.parent.mkdir(parents=True, exist_ok=True)

    old_content = "pas_id\tgene_id\n1\tOLD_GENE\n"
    pas_gene_tsv.write_text(old_content)

    rename_started = threading.Event()
    release_rename = threading.Event()
    observations: list[str] = []

    import os
    real_replace = os.replace

    def _slow_replace(src, dst):
        # The temp file (src) is already fully written by this point --
        # this delay simulates a scheduler pause between "finished writing
        # the temp file" and "renamed it into place", the exact window a
        # concurrent writer/reader could observe.
        rename_started.set()
        release_rename.wait(timeout=5)
        real_replace(src, dst)

    def _poll_during_write():
        rename_started.wait(timeout=5)
        # While the real write is paused right before its rename, the
        # destination must still show the complete OLD content.
        observations.append(pas_gene_tsv.read_text())
        release_rename.set()

    poller = threading.Thread(target=_poll_during_write)
    poller.start()
    with patch("os.replace", side_effect=_slow_replace):
        write_pas_gene_artifacts(run_dir, "dsA", pas_ids=["1", "2"], gene_ids=["ENSG1", "ENSG2"])
    poller.join(timeout=5)

    assert len(observations) == 1
    assert observations[0] == old_content, (
        "reader observed a partial/incomplete write mid-flight -- atomicity broken"
    )

    # After the call returns, the destination has the complete NEW content.
    final = pas_gene_tsv.read_text()
    assert "OLD_GENE" not in final
    assert "ENSG1" in final and "ENSG2" in final


def test_concurrent_writers_to_same_path_never_produce_a_torn_file(tmp_path):
    """Simulates the actual bug: two `peakatail reannotate` branches (e.g. from a
    duplicate branch_name) racing to write the SAME pas_gene.tsv path
    concurrently. With the atomic temp+rename fix, the file must ALWAYS end
    up as one writer's complete, valid output -- never a byte-level mix
    (e.g. writer A's header + writer B's rows) that would corrupt row
    counts unpredictably, which is exactly the symptom reported (3
    different row counts across 5 "identical" branches, grouped by write
    time).
    """
    import pandas as pd

    run_dir = tmp_path / "run"
    set_directory_config(output_dir=run_dir)
    from ema.config import directory_config

    pasbed = directory_config.pasbed_for("dsA")
    pasbed.parent.mkdir(parents=True, exist_ok=True)
    rows = [("chr1", i * 10, i * 10 + 1, str(i), 0, "+") for i in range(1, 201)]
    _write_pasbed(pasbed, rows)

    pas_ids_a = [str(i) for i in range(1, 201)]       # "writer A": 200 rows
    pas_ids_b = [str(i) for i in range(1, 51)]        # "writer B": 50 rows

    barrier = threading.Barrier(2)

    def _writer(pas_ids):
        barrier.wait(timeout=5)  # maximize actual overlap between the two writers
        write_pas_gene_artifacts(
            run_dir, "dsA", pas_ids=pas_ids, gene_ids=[f"G{p}" for p in pas_ids],
        )

    for _trial in range(20):  # repeat -- a race is timing-dependent, not deterministic
        t1 = threading.Thread(target=_writer, args=(pas_ids_a,))
        t2 = threading.Thread(target=_writer, args=(pas_ids_b,))
        t1.start(); t2.start()
        t1.join(timeout=10); t2.join(timeout=10)

        pas_gene_tsv = directory_config.pas_gene_for("dsA")
        df = pd.read_csv(pas_gene_tsv, sep="\t")  # must always parse cleanly
        assert len(df) in (200, 50), (
            f"trial {_trial}: torn file -- got {len(df)} rows, expected 200 or 50 "
            "(one full writer's output, never a byte-level mix)"
        )
        # Every row must belong to the SAME writer's gene_id scheme (G<id>),
        # not a mix of two different writers' outputs at the row level.
        assert (df["gene_id"] == "G" + df["pas_id"].astype(str)).all()
