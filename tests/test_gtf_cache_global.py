"""Tests for the global GTF cache tier in ema.annotate.gtf_cache.

Covers:
- global_cache_dir() respects XDG_CACHE_HOME
- gtf_global_fingerprint() is stable across calls and changes on content/mtime
- lookup_global_cache() returns None on cold cache
- populate_global_cache() writes a valid entry and lookup then hits
- process_gtf_cached() tries global cache first, then per-output-dir, then parses
- Warm global cache skips parse (speed / call count assertion)
- Atomic write: partial tmp dir is cleaned on error
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

import ema.annotate.gtf_cache as gtf_cache_mod
from ema.annotate.gtf_cache import (
    global_cache_dir,
    gtf_global_fingerprint,
    lookup_global_cache,
    populate_global_cache,
    process_gtf_cached,
)


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

MINIMAL_GTF = """\
##format: GTF
chr1\tENSEMBL\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG001"; gene_name "G1"; gene_biotype "protein_coding";
chr1\tENSEMBL\tthree_prime_utr\t1800\t2000\t.\t+\t.\tgene_id "ENSG001"; transcript_id "T001"; gene_name "G1";
"""


@pytest.fixture()
def gtf_file(tmp_path: Path) -> Path:
    """Write a minimal GTF and return its path."""
    p = tmp_path / "test.gtf"
    p.write_text(MINIMAL_GTF)
    return p


@pytest.fixture()
def cache_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Redirect global_cache_dir() to a temp directory for isolation."""
    root = tmp_path / "global_cache"
    root.mkdir()
    monkeypatch.setattr(gtf_cache_mod, "global_cache_dir", lambda: root)
    return root


@pytest.fixture()
def output_dir(tmp_path: Path) -> Path:
    """Per-output-dir for process_gtf_cached tests."""
    d = tmp_path / "out"
    d.mkdir()
    return d


# ---------------------------------------------------------------------------
# global_cache_dir()
# ---------------------------------------------------------------------------

class TestGlobalCacheDir:
    def test_default_uses_home_cache(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
        result = global_cache_dir()
        assert str(result).endswith(os.path.join(".cache", "peakatail", "gtf"))

    def test_respects_xdg_cache_home(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        custom = tmp_path / "xdg_cache"
        monkeypatch.setenv("XDG_CACHE_HOME", str(custom))
        result = global_cache_dir()
        assert result == custom / "peakatail" / "gtf"

    def test_empty_xdg_falls_back_to_home(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("XDG_CACHE_HOME", "   ")
        result = global_cache_dir()
        assert "peakatail" in str(result)
        assert ".cache" in str(result)


# ---------------------------------------------------------------------------
# gtf_global_fingerprint()
# ---------------------------------------------------------------------------

class TestGtfGlobalFingerprint:
    def test_stable_across_calls(self, gtf_file: Path) -> None:
        fp1 = gtf_global_fingerprint(gtf_file)
        fp2 = gtf_global_fingerprint(gtf_file)
        assert fp1 == fp2

    def test_contains_sha_mtime_size(self, gtf_file: Path) -> None:
        fp = gtf_global_fingerprint(gtf_file)
        parts = fp.split("_")
        assert len(parts) == 3
        # SHA prefix: 16 hex chars
        assert len(parts[0]) == 16
        assert all(c in "0123456789abcdef" for c in parts[0])
        # mtime_ns and size are numeric
        assert parts[1].isdigit()
        assert parts[2].isdigit()

    def test_changes_on_content_change(self, gtf_file: Path) -> None:
        fp_before = gtf_global_fingerprint(gtf_file)
        gtf_file.write_text(MINIMAL_GTF + "# extra line\n")
        fp_after = gtf_global_fingerprint(gtf_file)
        assert fp_before != fp_after

    def test_changes_on_mtime_change(self, gtf_file: Path) -> None:
        fp_before = gtf_global_fingerprint(gtf_file)
        # Set future mtime without changing content
        future = os.stat(gtf_file).st_mtime + 10
        os.utime(gtf_file, (future, future))
        fp_after = gtf_global_fingerprint(gtf_file)
        assert fp_before != fp_after


# ---------------------------------------------------------------------------
# lookup_global_cache() — cold cache
# ---------------------------------------------------------------------------

class TestLookupGlobalCacheCold:
    def test_returns_none_when_no_entry(
        self, gtf_file: Path, cache_root: Path
    ) -> None:
        result = lookup_global_cache(gtf_file)
        assert result is None

    def test_returns_none_on_corrupt_manifest(
        self, gtf_file: Path, cache_root: Path
    ) -> None:
        fp = gtf_global_fingerprint(gtf_file)
        entry = cache_root / fp
        entry.mkdir()
        (entry / "manifest.json").write_text("NOT JSON")
        assert lookup_global_cache(gtf_file) is None

    def test_returns_none_on_wrong_fingerprint_in_manifest(
        self, gtf_file: Path, cache_root: Path
    ) -> None:
        fp = gtf_global_fingerprint(gtf_file)
        entry = cache_root / fp
        entry.mkdir()
        (entry / "manifest.json").write_text(json.dumps({"fingerprint": "wrong"}))
        assert lookup_global_cache(gtf_file) is None

    def test_returns_none_when_utr_json_missing(
        self, gtf_file: Path, cache_root: Path
    ) -> None:
        fp = gtf_global_fingerprint(gtf_file)
        entry = cache_root / fp
        entry.mkdir()
        manifest = {
            "fingerprint": fp,
            "gtf_path": str(gtf_file),
            "endbed": "gene_ends.bed",
            "features": "raw_features.tsv",
            "utr_lengths_tsv": "utr_lengths.tsv",
        }
        (entry / "manifest.json").write_text(json.dumps(manifest))
        # No utr_lengths.json
        assert lookup_global_cache(gtf_file) is None


# ---------------------------------------------------------------------------
# populate_global_cache() + lookup_global_cache() — round trip
# ---------------------------------------------------------------------------

class TestPopulateAndLookup:
    def _make_artifacts(self, tmp_path: Path) -> tuple[str, str, str]:
        endbed = tmp_path / "gene_ends.bed"
        features = tmp_path / "raw_features.tsv"
        utr_tsv = tmp_path / "utr_lengths.tsv"
        endbed.write_text("chr1\t999\t2000\tENSG001\t0\t+\n")
        features.write_text("ENSG001\tG1\n")
        utr_tsv.write_text("ENSG001\t200\n")
        return str(endbed), str(features), str(utr_tsv)

    def test_populate_then_lookup_returns_data(
        self, gtf_file: Path, cache_root: Path, tmp_path: Path
    ) -> None:
        endbed, features, utr_tsv = self._make_artifacts(tmp_path)
        utr_lengths = {"ENSG001": 200}

        populate_global_cache(
            gtf_path=gtf_file,
            utr_lengths=utr_lengths,
            endbed_path=endbed,
            features_path=features,
            utr_lengths_path=utr_tsv,
        )

        hit = lookup_global_cache(gtf_file)
        assert hit is not None
        assert hit["utr_lengths"] == utr_lengths
        assert Path(hit["endbed_path"]).exists()
        assert Path(hit["features_path"]).exists()

    def test_populate_returns_same_shape_as_lookup(
        self, gtf_file: Path, cache_root: Path, tmp_path: Path
    ) -> None:
        endbed, features, utr_tsv = self._make_artifacts(tmp_path)
        utr_lengths = {"ENSG001": 200}

        populated = populate_global_cache(
            gtf_path=gtf_file,
            utr_lengths=utr_lengths,
            endbed_path=endbed,
            features_path=features,
            utr_lengths_path=utr_tsv,
        )
        looked_up = lookup_global_cache(gtf_file)

        assert set(populated.keys()) == set(looked_up.keys())
        assert populated["utr_lengths"] == looked_up["utr_lengths"]

    def test_populate_idempotent_on_existing_entry(
        self, gtf_file: Path, cache_root: Path, tmp_path: Path
    ) -> None:
        endbed, features, utr_tsv = self._make_artifacts(tmp_path)
        utr_lengths = {"ENSG001": 200}

        populate_global_cache(gtf_file, utr_lengths, endbed, features, utr_tsv)
        # Second call must not raise even though entry dir already exists
        populate_global_cache(gtf_file, utr_lengths, endbed, features, utr_tsv)
        assert lookup_global_cache(gtf_file) is not None

    def test_no_stale_tmp_dirs_after_populate(
        self, gtf_file: Path, cache_root: Path, tmp_path: Path
    ) -> None:
        endbed, features, utr_tsv = self._make_artifacts(tmp_path)
        populate_global_cache(gtf_file, {"ENSG001": 200}, endbed, features, utr_tsv)

        tmp_dirs = list(cache_root.glob(".tmp_gtf_*"))
        assert tmp_dirs == [], f"Stale tmp dirs found: {tmp_dirs}"


# ---------------------------------------------------------------------------
# process_gtf_cached() — cache hierarchy
# ---------------------------------------------------------------------------

class TestProcessGtfCached:
    """Integration tests using a mock gtf_bed to avoid heavy parsing."""

    def _make_artifact_files(self, out_dir: Path) -> tuple[Path, Path, Path]:
        endbed = out_dir / "gene_ends.bed"
        features = out_dir / "raw_features.tsv"
        utr_tsv = out_dir / "utr_lengths.tsv"
        return endbed, features, utr_tsv

    def _mock_gtf_bed(self, endbed: Path, features: Path, utr_tsv: Path) -> MagicMock:
        """Return a mock that writes artifact stubs and returns utr_lengths."""
        def _side_effect(endbeddir, gtfdir, featuresdir, utr_lengths_dir=None, **kw):
            Path(endbeddir).write_text("bed_stub\n")
            Path(featuresdir).write_text("feat_stub\n")
            if utr_lengths_dir:
                Path(utr_lengths_dir).write_text("utr_stub\n")
            return {"ENSG001": 300}

        m = MagicMock(side_effect=_side_effect)
        return m

    def test_cold_cache_calls_gtf_bed(
        self, gtf_file: Path, cache_root: Path, output_dir: Path
    ) -> None:
        endbed, features, utr_tsv = self._make_artifact_files(output_dir)
        mock = self._mock_gtf_bed(endbed, features, utr_tsv)

        with patch("ema.annotate.gtf_cache._gtf_bed_fn", mock, create=True):
            result = process_gtf_cached(
                gtf_path=gtf_file,
                output_dir=output_dir,
                endbed_path=endbed,
                features_path=features,
                utr_lengths_path=utr_tsv,
            )

        assert mock.call_count == 1
        assert result == {"ENSG001": 300}

    def test_global_cache_hit_skips_gtf_bed(
        self, gtf_file: Path, cache_root: Path, output_dir: Path, tmp_path: Path
    ) -> None:
        # First: populate global cache via a separate artifacts dir
        art_dir = tmp_path / "arts"
        art_dir.mkdir()
        endbed_src = art_dir / "gene_ends.bed"
        features_src = art_dir / "raw_features.tsv"
        utr_src = art_dir / "utr_lengths.tsv"
        endbed_src.write_text("bed\n")
        features_src.write_text("feat\n")
        utr_src.write_text("utr\n")

        populate_global_cache(
            gtf_path=gtf_file,
            utr_lengths={"ENSG001": 400},
            endbed_path=endbed_src,
            features_path=features_src,
            utr_lengths_path=utr_src,
        )

        # Now call process_gtf_cached — must not call gtf_bed
        mock = MagicMock(return_value={"ENSG001": 999})
        endbed, features, utr_tsv = self._make_artifact_files(output_dir)

        with patch("ema.annotate.gtf_cache._gtf_bed_fn", mock, create=True):
            result = process_gtf_cached(
                gtf_path=gtf_file,
                output_dir=output_dir,
                endbed_path=endbed,
                features_path=features,
                utr_lengths_path=utr_tsv,
            )

        assert mock.call_count == 0
        assert result == {"ENSG001": 400}

    def test_after_cold_parse_global_cache_populated(
        self, gtf_file: Path, cache_root: Path, output_dir: Path
    ) -> None:
        endbed, features, utr_tsv = self._make_artifact_files(output_dir)
        mock = self._mock_gtf_bed(endbed, features, utr_tsv)

        with patch("ema.annotate.gtf_cache._gtf_bed_fn", mock, create=True):
            process_gtf_cached(
                gtf_path=gtf_file,
                output_dir=output_dir,
                endbed_path=endbed,
                features_path=features,
                utr_lengths_path=utr_tsv,
            )

        # Global cache must now be populated
        hit = lookup_global_cache(gtf_file)
        assert hit is not None
        assert hit["utr_lengths"] == {"ENSG001": 300}

    def test_second_call_same_output_dir_hits_local_cache(
        self, gtf_file: Path, cache_root: Path, output_dir: Path
    ) -> None:
        endbed, features, utr_tsv = self._make_artifact_files(output_dir)
        mock = self._mock_gtf_bed(endbed, features, utr_tsv)

        # First call — populates global + local
        with patch("ema.annotate.gtf_cache._gtf_bed_fn", mock, create=True):
            process_gtf_cached(
                gtf_path=gtf_file,
                output_dir=output_dir,
                endbed_path=endbed,
                features_path=features,
                utr_lengths_path=utr_tsv,
            )

        # Bust the global cache so only local remains
        fp = gtf_global_fingerprint(gtf_file)
        shutil.rmtree(cache_root / fp, ignore_errors=True)

        mock2 = MagicMock(return_value={"ENSG001": 999})
        with patch("ema.annotate.gtf_cache._gtf_bed_fn", mock2, create=True):
            result2 = process_gtf_cached(
                gtf_path=gtf_file,
                output_dir=output_dir,
                endbed_path=endbed,
                features_path=features,
                utr_lengths_path=utr_tsv,
            )

        assert mock2.call_count == 0
        assert result2 == {"ENSG001": 300}


# ---------------------------------------------------------------------------
# Speed: global cache hit is fast
# ---------------------------------------------------------------------------

class TestGlobalCacheSpeed:
    def test_global_hit_under_100ms(
        self, gtf_file: Path, cache_root: Path, tmp_path: Path
    ) -> None:
        art_dir = tmp_path / "arts"
        art_dir.mkdir()
        endbed = art_dir / "gene_ends.bed"
        features = art_dir / "raw_features.tsv"
        utr_tsv = art_dir / "utr_lengths.tsv"
        endbed.write_text("bed\n")
        features.write_text("feat\n")
        utr_tsv.write_text("utr\n")

        populate_global_cache(
            gtf_path=gtf_file,
            utr_lengths={f"ENSG{i:04d}": i * 100 for i in range(1000)},
            endbed_path=endbed,
            features_path=features,
            utr_lengths_path=utr_tsv,
        )

        t0 = time.monotonic()
        hit = lookup_global_cache(gtf_file)
        elapsed = time.monotonic() - t0

        assert hit is not None
        assert elapsed < 0.1, f"Global cache hit took {elapsed:.3f}s (> 100 ms)"
