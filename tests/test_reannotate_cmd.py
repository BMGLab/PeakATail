"""Tests for `peakatail reannotate` (ema/cli/reannotate.py + ema/reannotate.py).

Builds a tiny synthetic BASE ``peakatail run`` output dir with the raw artifacts
`peakatail reannotate` needs:

    <base>/posbed.bed                         (unified +strand PAS)
    <base>/negbed.bed                         (unified -strand PAS)
    <base>/unified/concatenated.mtx           (PAS x cell counts)
    <base>/unified/concatenated_cbs.tsv       (namespaced barcodes)
    <base>/01_peak_calling/ds1/pasbed.bed     (per-dataset raw PAS, for
                                                annotatedpas.bed)

plus a minimal 2-gene GTF, then drives the real ``find_close`` ->
``annotate`` -> ``preprocessing`` chain end to end (real bedtools, real
pandas/scipy) through the ``peakatail reannotate`` CLI.  Only the clustering step
is stubbed -- exactly like ``tests/test_downstream_parallel.py`` does --
since Leiden clustering correctness is exercised by the dedicated clustering
test suite, not here; this file is about reannotate-specific plumbing: trim
wiring, the chainable run-dir contract, and NOT re-running peak calling.
"""
from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

# Several downstream modules import ema.config at module load time (see
# tests/test_downstream_parallel.py for the same precaution); pre-load them
# under a patched argv before any test runs.
with patch.object(sys, "argv", ["ema"]):
    import ema.matrixfilter as _mf  # noqa: F401
    import ema.annotate.annotate as _ann  # noqa: F401
    import ema.clustering.clustering as _clust_mod

import pandas as pd
import pytest
import scipy.io as sci
import scipy.sparse as sp
from click.testing import CliRunner

pytest.importorskip("anndata")
pytest.importorskip("pybedtools")

from ema.cli import main  # noqa: E402

N_CELLS = 10

# Two genes:
#   ENSG...1 has a 500bp 3'UTR -> PAS "1" (inside it) is TIER_1, always kept.
#   ENSG...2 has no UTR entry  -> PAS "2" (~1950bp past its extended end) is
#   only TIER_3, kept ONLY when --include-extended AND --max-gene-distance
#   covers that distance.  This is the knob the differential test flips.
GTF_TEXT = (
    'chr1\tHAVANA\tgene\t1000\t1050\t.\t+\t.\tgene_id "ENSG00000000001"; '
    'gene_type "protein_coding"; gene_name "GENE1";\n'
    'chr1\tHAVANA\tthree_prime_utr\t1000\t1500\t.\t+\t.\tgene_id "ENSG00000000001"; '
    'gene_type "protein_coding"; gene_name "GENE1";\n'
    'chr1\tHAVANA\tgene\t20000\t20050\t.\t+\t.\tgene_id "ENSG00000000002"; '
    'gene_type "protein_coding"; gene_name "GENE2";\n'
)

PASBED_TEXT = (
    "chr1\t1100\t1101\t1\t0\t+\n"
    "chr1\t27000\t27001\t2\t0\t+\n"
)


def _fake_clustering(adata, output_h5ad=None, **_kwargs):
    """Stand-in for ema.clustering.clustering.clustering.

    Mutates ``adata`` in place (adds a fake 'leiden' column, matching what
    the real strategy does) and writes it to ``output_h5ad`` -- exactly the
    side effects run_one_dataset_downstream relies on. Real Leiden/TF-IDF is
    exercised elsewhere; this test is about the reannotate-specific wiring
    around it (trim -> annotate -> preprocess -> provenance -> manifest).
    """
    adata.obs["leiden"] = pd.Categorical([str(i % 2) for i in range(adata.n_obs)])
    if output_h5ad is not None:
        Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
        adata.write(output_h5ad)
    return adata


def _build_base_run(tmp_path: Path) -> Path:
    """Synthetic BASE `peakatail run` output: 1 dataset ('ds1'), 2 PAS, 10 cells."""
    base = tmp_path / "base_run"
    (base / "unified").mkdir(parents=True)
    (base / "01_peak_calling" / "ds1").mkdir(parents=True)

    (base / "posbed.bed").write_text(PASBED_TEXT)
    (base / "negbed.bed").write_text("")
    # Per-dataset raw pasbed (peak-calling artifact) -- reannotate copies this
    # into --out so annotatedpas.bed can be written without re-peak-calling.
    (base / "01_peak_calling" / "ds1" / "pasbed.bed").write_text(PASBED_TEXT)

    cbs = [f"ds1_CB{i}" for i in range(N_CELLS)]
    (base / "unified" / "concatenated_cbs.tsv").write_text("\n".join(cbs) + "\n")

    # PAS "1" (row 0) present in all 10 cells; PAS "2" (row 1) in 5 of them.
    rows = [0] * N_CELLS + [1] * 5
    cols = list(range(N_CELLS)) + list(range(5))
    data = [5] * N_CELLS + [3] * 5
    mat = sp.coo_matrix((data, (rows, cols)), shape=(2, N_CELLS))
    sci.mmwrite(str(base / "unified" / "concatenated.mtx"), mat, field="integer")

    return base


def _invoke_reannotate(base: Path, out: Path, gtf: Path, *, extra_args=()):
    args = [
        "reannotate",
        "--base-run", str(base),
        "--out", str(out),
        "--gtf", str(gtf),
        "--min-read", "1",
        "--min-cells", "1",
        "--min-pas-per-cell", "1",
        "--no-log-file",
        "--no-progress",
        *extra_args,
    ]
    with patch.object(_clust_mod, "clustering", side_effect=_fake_clustering):
        return CliRunner().invoke(main, args)


@pytest.fixture(autouse=True)
def _isolated_gtf_cache(tmp_path, monkeypatch):
    # Keep the global GTF cache (~/.cache/peakatail/gtf/) out of this test's
    # way -- hermetic, and avoids collisions with a developer's real cache.
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg_cache"))


@pytest.fixture()
def gtf(tmp_path) -> Path:
    p = tmp_path / "genes.gtf"
    p.write_text(GTF_TEXT)
    return p


# --------------------------------------------------------------------------- #
# CLI wiring
# --------------------------------------------------------------------------- #


class TestReannotateCliWiring:
    def test_help_lists_required_options(self):
        result = CliRunner().invoke(main, ["reannotate", "--help"])
        assert result.exit_code == 0
        # Every trim/filter/cluster knob run_one_dataset_downstream/find_close
        # accept must be a CLI flag -- nothing hardcoded, so an OFAT sweep can
        # vary any of them purely via CLI flags.
        for flag in ("--base-run", "--out", "--gtf",
                     "--max-gene-distance", "--utr-multiplier", "--include-extended",
                     "--cluster-method", "--resolution", "--n-neighbors", "--n-pcs",
                     "--n-svd-components", "--n-top-hvg", "--random-seed",
                     "--tfidf-scale-factor", "--depth-corr-threshold",
                     "--external-clusters",
                     "--min-read", "--min-cells", "--min-pas-per-cell", "--threads"):
            assert flag in result.output, f"{flag} missing from --help output"

    def test_missing_required_options_errors(self):
        result = CliRunner().invoke(main, ["reannotate"])
        assert result.exit_code != 0

    def test_out_equal_base_run_errors_clearly(self, tmp_path, gtf):
        base = _build_base_run(tmp_path)
        result = _invoke_reannotate(base, base, gtf)
        assert result.exit_code != 0

    def test_missing_base_artifact_errors_before_heavy_work(self, tmp_path, gtf):
        empty_base = tmp_path / "empty_base"
        empty_base.mkdir()
        result = _invoke_reannotate(empty_base, tmp_path / "out", gtf)
        assert result.exit_code != 0
        assert "posbed.bed" in str(result.output) + str(result.exception)


# --------------------------------------------------------------------------- #
# End-to-end: chainable run dir + no peak calling
# --------------------------------------------------------------------------- #


class TestReannotateProducesChainableRunDir:
    def test_default_invocation_produces_full_run_dir_without_peak_calling(
        self, tmp_path, gtf,
    ):
        base = _build_base_run(tmp_path)
        out = tmp_path / "branch_default"

        import ema.countmatrix.peackcalling as pc_mod
        with patch.object(
            pc_mod, "peak_calling",
            side_effect=AssertionError("peak calling must NOT run under peakatail reannotate"),
        ):
            result = _invoke_reannotate(base, out, gtf)

        assert result.exit_code == 0, result.output

        # 07_clustering/<ds>/clusters.h5ad — the new matrix + clustering.
        h5ad = out / "07_clustering" / "ds1" / "clusters.h5ad"
        assert h5ad.exists()

        # E3 provenance: reconciled per-dataset ledgers + summary.
        assert (out / "provenance" / "by_dataset" / "pas_ledger.tsv").exists()
        assert (out / "provenance" / "by_dataset" / "cell_ledger.tsv").exists()
        assert (out / "provenance" / "reconcile_summary.json").exists()

        # E2 run manifest.
        assert (out / "run_manifest.json").exists()

        # Per-dataset annotatedpas.bed + the branch's own root pasbed.bed.
        assert (out / "03_gtf_annotation" / "ds1" / "annotatedpas.bed").exists()
        assert (out / "pasbed.bed").exists()

        # The branch's own human-readable manifest.
        assert (out / "branch_manifest.json").exists()

    def test_run_dir_loads_via_ema_data_run(self, tmp_path, gtf):
        base = _build_base_run(tmp_path)
        out = tmp_path / "branch_loadable"
        result = _invoke_reannotate(base, out, gtf)
        assert result.exit_code == 0, result.output

        from ema.data import Run
        run = Run.from_dir(out)
        h5ad_paths = list((out / "07_clustering").glob("*/clusters.h5ad"))
        assert h5ad_paths
        # Manifest-driven artifact resolution works, same as a base run.
        assert run.manifest["root"] == str(out.resolve())

    def test_run_dir_is_valid_switch_diff_input(self, tmp_path, gtf):
        """`peakatail switch diff` can consume the reannotated matrix directly —
        proves the branch is a chainable input to the next pipeline step."""
        base = _build_base_run(tmp_path)
        out = tmp_path / "branch_switch"
        result = _invoke_reannotate(
            base, out, gtf,
            extra_args=["--include-extended", "--max-gene-distance", "5000"],
        )
        assert result.exit_code == 0, result.output

        h5ad = out / "07_clustering" / "ds1" / "clusters.h5ad"
        assert h5ad.exists()

        diff_result = CliRunner().invoke(main, [
            "switch", "diff",
            "-i", str(h5ad),
            "--cluster-key", "leiden",
            "--strategy", "fisher",
            "-o", str(tmp_path / "switch_diff_out"),
            "--no-log-file",
            "--no-progress",
        ])
        assert diff_result.exit_code == 0, diff_result.output


# --------------------------------------------------------------------------- #
# The trim actually takes effect
# --------------------------------------------------------------------------- #


class TestMaxGeneDistanceChangesSurvivingPas:
    def test_wider_max_gene_distance_keeps_more_pas(self, tmp_path, gtf):
        import anndata as ad

        base = _build_base_run(tmp_path)

        out_wide = tmp_path / "branch_wide"
        result_wide = _invoke_reannotate(
            base, out_wide, gtf,
            extra_args=["--include-extended", "--max-gene-distance", "5000"],
        )
        assert result_wide.exit_code == 0, result_wide.output

        out_narrow = tmp_path / "branch_narrow"
        result_narrow = _invoke_reannotate(
            base, out_narrow, gtf,
            extra_args=["--include-extended", "--max-gene-distance", "1000"],
        )
        assert result_narrow.exit_code == 0, result_narrow.output

        n_vars_wide = ad.read_h5ad(
            out_wide / "07_clustering" / "ds1" / "clusters.h5ad"
        ).n_vars
        n_vars_narrow = ad.read_h5ad(
            out_narrow / "07_clustering" / "ds1" / "clusters.h5ad"
        ).n_vars

        # PAS "2" (~1950bp past gene2's extended boundary) survives at
        # max_gene_distance=5000 but is dropped at max_gene_distance=1000 —
        # proves the trim parameter actually changes the surviving PAS set,
        # not just that the command runs.
        assert n_vars_wide == 2
        assert n_vars_narrow == 1
        assert n_vars_wide > n_vars_narrow

    def test_default_reuses_base_peak_calls_single_pas(self, tmp_path, gtf):
        """Default (--include-extended off) keeps only the always-TIER_1 PAS,
        proving the base run's raw peak calls (not a re-call) are what feeds
        the branch."""
        import anndata as ad

        base = _build_base_run(tmp_path)
        out = tmp_path / "branch_default2"
        result = _invoke_reannotate(base, out, gtf)
        assert result.exit_code == 0, result.output

        n_vars = ad.read_h5ad(out / "07_clustering" / "ds1" / "clusters.h5ad").n_vars
        assert n_vars == 1


# --------------------------------------------------------------------------- #
# Every clustering knob is a real CLI flag, not a hidden default
# --------------------------------------------------------------------------- #


class TestClusteringKnobsReachClustering:
    def test_all_cluster_kwargs_flow_through_to_clustering(self, tmp_path, gtf):
        """CLI values for every clustering knob must reach the real
        ``clustering()`` call unchanged -- proves none of them are silently
        pinned to a hidden default inside reannotate_run, so an OFAT sweep
        can vary any of them purely via CLI flags."""
        base = _build_base_run(tmp_path)
        out = tmp_path / "branch_knobs"

        captured: dict = {}

        def _capturing_clustering(adata, output_h5ad=None, **kwargs):
            captured.update(kwargs)
            return _fake_clustering(adata, output_h5ad=output_h5ad, **kwargs)

        args = [
            "reannotate",
            "--base-run", str(base),
            "--out", str(out),
            "--gtf", str(gtf),
            "--min-read", "1", "--min-cells", "1", "--min-pas-per-cell", "1",
            "--no-log-file", "--no-progress",
            "--cluster-method", "leiden_tfidf",
            "--resolution", "0.42",
            "--n-neighbors", "3",
            "--n-pcs", "5",
            "--n-svd-components", "6",
            "--n-top-hvg", "7",
            "--random-seed", "99",
            "--tfidf-scale-factor", "12345",
            "--depth-corr-threshold", "0.33",
        ]
        with patch.object(_clust_mod, "clustering", side_effect=_capturing_clustering):
            result = CliRunner().invoke(main, args)
        assert result.exit_code == 0, result.output

        assert captured["method"] == "leiden_tfidf"
        assert captured["resolution"] == 0.42
        assert captured["n_neighbors"] == 3
        assert captured["n_pcs"] == 5
        assert captured["n_svd_components"] == 6
        assert captured["n_top_hvg"] == 7
        assert captured["random_seed"] == 99
        assert captured["tfidf_scale_factor"] == 12345.0
        assert captured["depth_corr_threshold"] == 0.33


# --------------------------------------------------------------------------- #
# PAS labels + clustering-only exclusion mask — the FILTER-EFFECT experiment's
# mechanism: LABEL every PAS (atlas_match / internal_priming / in_3utr) in
# annotatedpas.bed, always, regardless of exclusion; independently, the
# exclude_* flags narrow ONLY the matrix preprocessing()/clustering() sees.
# --------------------------------------------------------------------------- #

import json  # noqa: E402


def _atlas_status_tsv(matches: dict) -> str:
    lines = ["unified_pas_id\tatlas_match\tatlas_distance_bp"]
    for pas_id, match in matches.items():
        lines.append(f"{pas_id}\t{match}\t{0 if match else ''}")
    return "\n".join(lines) + "\n"


def _annotatedpas_rows(out: Path, ds: str = "ds1") -> dict:
    """{pas_id: [trailing columns after strand]} from a dataset's annotatedpas.bed."""
    path = out / "03_gtf_annotation" / ds / "annotatedpas.bed"
    rows = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        rows[fields[3]] = fields[6:]
    return rows


# Both synthetic PAS need --include-extended + a wide enough --max-gene-distance
# to survive find_close()'s gene assignment (PAS "2" is TIER_3-only by default,
# per the GTF fixture comment above) -- these tests need BOTH PAS present so
# the clustering-exclusion mask has something to narrow.
_BOTH_PAS_ARGS = ["--include-extended", "--max-gene-distance", "5000"]


class TestPasLabelsAndClusteringMask:
    def test_no_labels_no_mask_reproduces_default_exactly(self, tmp_path, gtf):
        base = _build_base_run(tmp_path)
        out = tmp_path / "branch_nolabel"
        result = _invoke_reannotate(base, out, gtf, extra_args=_BOTH_PAS_ARGS)
        assert result.exit_code == 0, result.output

        manifest = json.loads((out / "branch_manifest.json").read_text())
        pm = manifest["pas_labels_and_mask"]
        assert pm["exclude_atlas_nonmatch"] is False
        assert pm["exclude_internal_priming"] is False
        assert pm["exclude_not_in_3utr"] is False
        assert pm["n_excluded_for_clustering"] == 0
        assert pm["n_kept_in_results"] == pm["n_used_for_clustering"] == 2

        # Base run's own PAS results are untouched.
        assert (base / "posbed.bed").read_text() == PASBED_TEXT
        # No --out posbed/negbed copies are written when nothing is labeled.
        assert not (out / "posbed.bed").exists()

        import anndata as ad
        a = ad.read_h5ad(out / "07_clustering" / "ds1" / "clusters.h5ad")
        assert sorted(a.var_names) == ["1", "2"]

    def test_atlas_exclude_masks_clustering_but_keeps_full_annotatedpas(self, tmp_path, gtf):
        base = _build_base_run(tmp_path)
        (base / "unified" / "atlas_status.tsv").write_text(
            _atlas_status_tsv({"1": True, "2": False})
        )
        out = tmp_path / "branch_atlas_excl"
        result = _invoke_reannotate(
            base, out, gtf, extra_args=[*_BOTH_PAS_ARGS, "--exclude-atlas-nonmatch"],
        )
        assert result.exit_code == 0, result.output

        manifest = json.loads((out / "branch_manifest.json").read_text())
        pm = manifest["pas_labels_and_mask"]
        assert pm["exclude_atlas_nonmatch"] is True
        assert pm["n_kept_in_results"] == 2          # BOTH PAS still in the results
        assert pm["n_used_for_clustering"] == 1       # only the matched one clustered
        assert pm["n_excluded_for_clustering"] == 1
        assert pm["labels"]["atlas"]["n_matched"] == 1

        # annotatedpas.bed carries BOTH PAS, both labeled (never dropped there).
        rows = _annotatedpas_rows(out)
        assert set(rows) == {"1", "2"}
        assert rows["1"][1] == "True"   # atlas_match
        assert rows["2"][1] == "False"

        # Base run untouched.
        assert (base / "posbed.bed").read_text() == PASBED_TEXT

        # Only PAS "1" (atlas-matched) survives into the clustering matrix.
        import anndata as ad
        a = ad.read_h5ad(out / "07_clustering" / "ds1" / "clusters.h5ad")
        assert list(a.var_names) == ["1"]

    def test_atlas_label_without_exclude_flag_labels_but_does_not_mask(self, tmp_path, gtf):
        """--atlas alone (no --exclude-atlas-nonmatch) labels every PAS but
        clusters on the full set -- labeling and exclusion are independent."""
        base = _build_base_run(tmp_path)
        (base / "unified" / "atlas_status.tsv").write_text(
            _atlas_status_tsv({"1": True, "2": False})
        )
        out = tmp_path / "branch_atlas_label_only"
        # cache is auto-detected, no exclude flag needed
        result = _invoke_reannotate(base, out, gtf, extra_args=_BOTH_PAS_ARGS)
        assert result.exit_code == 0, result.output

        rows = _annotatedpas_rows(out)
        assert rows["1"][1] == "True"
        assert rows["2"][1] == "False"

        manifest = json.loads((out / "branch_manifest.json").read_text())
        assert manifest["pas_labels_and_mask"]["n_used_for_clustering"] == 2

    def test_exclude_atlas_nonmatch_without_cache_or_atlas_errors_clearly(self, tmp_path, gtf):
        base = _build_base_run(tmp_path)  # no unified/atlas_status.tsv
        out = tmp_path / "branch_atlas_excl_err"
        result = _invoke_reannotate(base, out, gtf, extra_args=["--exclude-atlas-nonmatch"])
        assert result.exit_code != 0
        assert "atlas" in str(result.output).lower() + str(result.exception).lower()

    def test_exclude_not_in_3utr_masks_clustering_but_keeps_full_annotatedpas(
        self, tmp_path, gtf,
    ):
        pytest.importorskip("pybedtools")
        base = _build_base_run(tmp_path)
        # Region BED covers only PAS "1"'s locus (chr1:1090-1110), not PAS "2"'s.
        region_bed = tmp_path / "three_prime_utr.bed"
        region_bed.write_text("chr1\t1090\t1110\tregion\t0\t+\n")

        out = tmp_path / "branch_3utr_excl"
        result = _invoke_reannotate(
            base, out, gtf,
            extra_args=[*_BOTH_PAS_ARGS, "--annotation-bed", str(region_bed),
                        "--exclude-not-in-3utr"],
        )
        assert result.exit_code == 0, result.output

        manifest = json.loads((out / "branch_manifest.json").read_text())
        pm = manifest["pas_labels_and_mask"]
        assert pm["exclude_not_in_3utr"] is True
        assert pm["n_kept_in_results"] == 2
        assert pm["n_used_for_clustering"] == 1

        rows = _annotatedpas_rows(out)
        assert set(rows) == {"1", "2"}
        assert rows["1"][4] == "True"    # in_3utr
        assert rows["2"][4] == "False"

        import anndata as ad
        a = ad.read_h5ad(out / "07_clustering" / "ds1" / "clusters.h5ad")
        assert list(a.var_names) == ["1"]

    def test_exclude_not_in_3utr_without_annotation_bed_errors_clearly(self, tmp_path, gtf):
        base = _build_base_run(tmp_path)
        out = tmp_path / "branch_3utr_excl_err"
        result = _invoke_reannotate(base, out, gtf, extra_args=["--exclude-not-in-3utr"])
        assert result.exit_code != 0
        assert "3utr" in str(result.output).lower() + str(result.exception).lower() or \
               "annotation-bed" in str(result.output).lower() + str(result.exception).lower()

    def test_exclude_internal_priming_masks_clustering_but_keeps_full_annotatedpas(
        self, tmp_path, gtf,
    ):
        pytest.importorskip("pyfaidx")
        base = _build_base_run(tmp_path)

        # chr1 filler with no long A-run, except an A-rich stretch right
        # after PAS "2" (chr1:27000-27001, + strand -> checked window is
        # [26991:27031)) so only PAS "2" gets flagged internal-priming.
        seq = list(("ACGT" * 10000)[:30000])
        for i in range(27000, 27030):
            seq[i] = "A"
        fasta_path = tmp_path / "genome.fa"
        fasta_path.write_text(">chr1\n" + "".join(seq) + "\n")

        out = tmp_path / "branch_ip_excl"
        result = _invoke_reannotate(
            base, out, gtf,
            extra_args=[*_BOTH_PAS_ARGS, "--genome-fasta", str(fasta_path),
                        "--exclude-internal-priming"],
        )
        assert result.exit_code == 0, result.output

        manifest = json.loads((out / "branch_manifest.json").read_text())
        pm = manifest["pas_labels_and_mask"]
        assert pm["exclude_internal_priming"] is True
        assert pm["n_kept_in_results"] == 2
        assert pm["n_used_for_clustering"] == 1

        rows = _annotatedpas_rows(out)
        assert set(rows) == {"1", "2"}
        assert rows["1"][3] == "False"   # internal_priming
        assert rows["2"][3] == "True"

        import anndata as ad
        a = ad.read_h5ad(out / "07_clustering" / "ds1" / "clusters.h5ad")
        assert list(a.var_names) == ["1"]

    def test_exclude_internal_priming_without_genome_fasta_errors_clearly(self, tmp_path, gtf):
        base = _build_base_run(tmp_path)
        out = tmp_path / "branch_ip_excl_err"
        result = _invoke_reannotate(base, out, gtf, extra_args=["--exclude-internal-priming"])
        assert result.exit_code != 0
        assert "genome-fasta" in str(result.output).lower() + str(result.exception).lower() or \
               "internal_priming" in str(result.output).lower() + str(result.exception).lower()

    def test_help_lists_pas_label_and_mask_options(self):
        result = CliRunner().invoke(main, ["reannotate", "--help"])
        assert result.exit_code == 0
        for flag in ("--atlas", "--atlas-distance", "--genome-fasta", "--annotation-bed",
                     "--exclude-atlas-nonmatch", "--exclude-internal-priming",
                     "--exclude-not-in-3utr"):
            assert flag in result.output, f"{flag} missing from --help output"
