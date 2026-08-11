"""D9: atlas_mode / ip_filter_mode are selectable RunConfig fields, not
hardcoded, and default to "annotate" (never drop by default)."""
from __future__ import annotations

from dataclasses import fields

from ema.cli.config_schema import RunConfig, field_specs


def test_atlas_mode_field_defaults_to_annotate_and_is_bridged():
    rc = RunConfig()
    assert rc.atlas_mode == "annotate"

    specs = field_specs(RunConfig)
    spec = specs["atlas_mode"]
    assert spec.choice == ("annotate", "filter")
    assert spec.skip_legacy_bridge is False
    assert spec.legacy_args_attr == "atlas_mode"
    assert spec.cli_flag == "--atlas-mode"
    assert spec.yaml_key == "atlas_mode"


def test_ip_filter_mode_field_defaults_to_annotate_and_is_bridged():
    rc = RunConfig()
    assert rc.ip_filter_mode == "annotate"

    specs = field_specs(RunConfig)
    spec = specs["ip_filter_mode"]
    assert spec.choice == ("annotate", "filter")
    assert spec.skip_legacy_bridge is False
    assert spec.legacy_args_attr == "ip_filter_mode"
    assert spec.cli_flag == "--ip-filter-mode"
    assert spec.yaml_key == "ip_filter_mode"


def test_atlas_and_ip_filter_mode_are_declared_fields():
    field_names = {f.name for f in fields(RunConfig)}
    assert "atlas_mode" in field_names
    assert "ip_filter_mode" in field_names


def test_apply_to_legacy_globals_bridges_atlas_mode():
    from ema import config as _cfg

    rc = RunConfig()
    rc.atlas_mode = "filter"
    rc.ip_filter_mode = "filter"
    try:
        rc.apply_to_legacy_globals()
        assert _cfg.args.atlas_mode == "filter"
        assert _cfg.args.ip_filter_mode == "filter"
    finally:
        for attr in ("atlas_mode", "ip_filter_mode"):
            try:
                delattr(_cfg.args._get(), attr)
            except AttributeError:
                pass
