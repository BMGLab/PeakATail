"""YAML loader contract tests."""
import logging
from pathlib import Path

import pytest
import yaml

from ema.cli.yaml_loader import load_run_yaml, RunYamlError


@pytest.fixture(autouse=True)
def _restore_ema_propagate():
    """Ensure ema logger propagates to root so caplog captures records.

    setup_logging() (called in test_logging_config.py) sets propagate=False
    on the ema logger. This fixture restores propagation before each test so
    caplog works correctly regardless of test execution order.
    """
    logger = logging.getLogger("ema")
    original = logger.propagate
    logger.propagate = True
    yield
    logger.propagate = original


def _write(tmp_path: Path, body: dict) -> Path:
    p = tmp_path / "config.yaml"
    p.write_text(yaml.safe_dump(body))
    return p


def test_minimal_yaml_returns_dict(tmp_path):
    p = _write(tmp_path, {
        "datasets": [{"id": "a", "merge_strategy": "none", "bams": ["a.bam"]}],
        "gtf": "g.gtf",
        "output_dir": "out",
    })
    cfg = load_run_yaml(p)
    assert cfg["datasets"][0]["id"] == "a"
    assert cfg["gtf"] == "g.gtf"


def test_dead_key_warns_but_loads(tmp_path, caplog):
    p = _write(tmp_path, {
        "datasets": [{"id": "a", "merge_strategy": "none", "bams": ["a.bam"]}],
        "gtf": "g.gtf",
        "pdui_method": "classic",   # DEAD in ema run
        "diff_method": "fisher",    # DEAD in ema run
    })
    with caplog.at_level(logging.WARNING):
        cfg = load_run_yaml(p)
    text = caplog.text
    assert "pdui_method" in text and "ema switch length" in text
    assert "diff_method" in text and "ema switch diff" in text
    # the keys ARE preserved in the cfg, callers can ignore them
    assert "pdui_method" in cfg


def test_unknown_key_warns_but_loads(tmp_path, caplog):
    p = _write(tmp_path, {
        "datasets": [{"id": "a", "merge_strategy": "none", "bams": ["a.bam"]}],
        "totally_made_up_key": 123,
    })
    with caplog.at_level(logging.WARNING):
        load_run_yaml(tmp_path / "config.yaml")
    assert "totally_made_up_key" in caplog.text


def test_missing_datasets_raises(tmp_path):
    p = _write(tmp_path, {"gtf": "g.gtf"})
    with pytest.raises(RunYamlError, match="datasets"):
        load_run_yaml(p)


def test_invalid_merge_strategy_raises(tmp_path):
    p = _write(tmp_path, {
        "datasets": [{"id": "a", "merge_strategy": "BOGUS", "bams": ["a.bam"]}],
    })
    with pytest.raises(RunYamlError, match="merge_strategy"):
        load_run_yaml(p)


def test_dataset_id_with_underscore_rejected(tmp_path):
    """Underscore in a dataset id breaks barcode namespacing → hard error."""
    import yaml
    from ema.cli.yaml_loader import load_run_yaml, RunYamlError
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump({"gtf": "x.gtf", "datasets": [
        {"id": "GSM1_StageI", "merge_strategy": "none", "bams": ["a.bam"]}]}))
    with pytest.raises(RunYamlError, match="underscore"):
        load_run_yaml(str(p))


def test_duplicate_dataset_id_rejected(tmp_path):
    """Reusing an id would collide two libraries' barcodes → hard error."""
    import yaml
    from ema.cli.yaml_loader import load_run_yaml, RunYamlError
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump({"gtf": "x.gtf", "datasets": [
        {"id": "GSM1-A", "bams": ["a.bam"]},
        {"id": "GSM1-A", "bams": ["b.bam"]}]}))
    with pytest.raises(RunYamlError, match="duplicate"):
        load_run_yaml(str(p))
