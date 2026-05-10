"""Wizard contract tests. Uses questionary's mock-friendly API."""
from unittest.mock import patch

from ema.cli.wizard import run as wizard_run, _ask_mode


def test_mode_picker_returns_string(monkeypatch):
    with patch("questionary.select") as mock_select:
        mock_select.return_value.ask.return_value = "run"
        result = _ask_mode()
    assert result == "run"


def test_wizard_aborts_on_none(monkeypatch):
    """If user hits Esc on the first prompt, wizard exits 1."""
    with patch("questionary.select") as mock_select:
        mock_select.return_value.ask.return_value = None
        rc = wizard_run()
    assert rc == 1


def test_wizard_dispatches_run(monkeypatch):
    """When mode=run is chosen and config is built, dispatch into run command."""
    called = {}
    def fake_run_dispatch(cfg):
        called["cfg"] = cfg
        return 0
    monkeypatch.setattr("ema.cli.wizard._dispatch_run", fake_run_dispatch)

    with patch("questionary.select") as mock_select, \
         patch("ema.cli.wizard._ask_run_config") as mock_run_cfg, \
         patch("ema.cli.wizard._confirm_proceed", return_value=True), \
         patch("ema.cli.wizard._confirm_save_yaml", return_value=False):
        mock_select.return_value.ask.return_value = "run"
        mock_run_cfg.return_value = {"datasets": [{"id": "x", "merge_strategy": "none", "bams": ["x.bam"]}]}
        rc = wizard_run()
    assert rc == 0
    assert "cfg" in called
