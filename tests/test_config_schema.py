"""Contract tests for the centralised RunConfig schema.

These tests pin the invariants Part B was designed to enforce:
* DEFAULTS, _LIVE_KEYS, and the legacy bridge are all derived from the
  schema -- modify the schema and the consumers update automatically.
* Every CLI flag has a YAML key (or is explicitly skip_legacy_bridge).
* No two fields share the same yaml_key, cli_flag, or legacy_alias.
* The schema itself can round-trip through its own constructors.
"""
from __future__ import annotations

import pytest

from ema.cli.config_schema import (
    RunConfig,
    FieldSpec,
    cli_defaults,
    click_options_from_schema,
    field_specs,
    legacy_alias_to_field_name,
    yaml_key_to_field_name,
    yaml_keys_from_schema,
)


def test_every_field_has_a_spec():
    """Every dataclass field must carry FieldSpec metadata."""
    specs = field_specs(RunConfig)
    from dataclasses import fields
    assert set(specs.keys()) == {f.name for f in fields(RunConfig)}, (
        "Some fields are missing FieldSpec metadata -- the auto-generators "
        "would silently skip them."
    )


def test_no_duplicate_cli_flags():
    flags = [s.cli_flag for s in field_specs(RunConfig).values() if s.cli_flag]
    assert len(flags) == len(set(flags)), (
        f"Duplicate CLI flags in schema: {sorted(set([f for f in flags if flags.count(f) > 1]))}"
    )


def test_no_duplicate_yaml_keys():
    keys = [s.yaml_key for s in field_specs(RunConfig).values() if s.yaml_key]
    assert len(keys) == len(set(keys)), (
        f"Duplicate YAML keys: {sorted(set([k for k in keys if keys.count(k) > 1]))}"
    )


def test_cli_defaults_table_covers_every_cli_field():
    """cli_defaults() must contain every field that has a cli_flag."""
    defaults = cli_defaults(RunConfig)
    for name, spec in field_specs(RunConfig).items():
        if spec.cli_flag is None:
            continue
        kebab = spec.cli_flag.lstrip("-")
        assert kebab in defaults, f"missing default for {kebab}"


def test_yaml_keys_from_schema_is_superset_of_live_keys():
    """The YAML loader's _LIVE_KEYS set must be a subset of (or equal to)
    the schema-derived set.  Drift here means the loader will warn on a
    valid key or accept an unknown one.
    """
    from ema.cli import yaml_loader
    schema_keys = yaml_keys_from_schema(RunConfig)
    # Every legacy _LIVE_KEYS entry must be in schema_keys -- if the
    # schema is missing one, that key would silently warn.
    missing_in_schema = yaml_loader._LIVE_KEYS - schema_keys
    # Allow a small set of legacy keys the schema doesn't model directly
    # (datasets is dict-shaped and validated by the loader; min_genes is
    # the legacy alias for min_pas_per_cell and gets handled in
    # apply_to_legacy_globals).
    allowed_legacy = {"datasets"}
    leftovers = missing_in_schema - allowed_legacy
    assert not leftovers, (
        f"_LIVE_KEYS contains keys missing from RunConfig schema: {sorted(leftovers)}.\n"
        "Add them as fields with a yaml_key, or move them to the allow-list."
    )


def test_defaults_module_derives_from_schema():
    """ema/cli/defaults.py::DEFAULTS must equal cli_defaults(RunConfig)
    on every overlapping key (no drift)."""
    from ema.cli.defaults import DEFAULTS
    schema_defaults = cli_defaults(RunConfig)
    for kebab, schema_default in schema_defaults.items():
        if kebab not in DEFAULTS:
            continue  # DEFAULTS may also include non-run subcommand keys
        assert DEFAULTS[kebab] == schema_default, (
            f"DEFAULTS[{kebab!r}] = {DEFAULTS[kebab]!r} but schema "
            f"default = {schema_default!r}; the two surfaces have drifted."
        )


def test_yaml_to_field_name_round_trip():
    """For every (field_name, spec.yaml_key) pair, the inverse map agrees."""
    inv = yaml_key_to_field_name(RunConfig)
    for name, spec in field_specs(RunConfig).items():
        if spec.yaml_key is not None:
            assert inv[spec.yaml_key] == name


def test_legacy_alias_resolves_to_field():
    """min_genes legacy alias should resolve to min_pas_per_cell."""
    aliases = legacy_alias_to_field_name(RunConfig)
    assert aliases.get("min_genes") == "min_pas_per_cell"


# ---------------------------------------------------------------------------
# Constructor / bridge tests
# ---------------------------------------------------------------------------

def test_from_yaml_dict_handles_legacy_alias():
    cfg = RunConfig.from_yaml_dict({
        "min_genes": 25,         # legacy alias
        "min_read": 800,
    })
    assert cfg.min_pas_per_cell == 25, (
        "legacy alias min_genes must be mapped to min_pas_per_cell"
    )
    assert cfg.min_read == 800


def test_from_yaml_dict_silently_ignores_unknown_keys():
    """The schema does not raise on unknown keys (loader handles warnings)."""
    cfg = RunConfig.from_yaml_dict({"totally_made_up": 99})
    # Should not raise; the field was simply skipped.
    assert isinstance(cfg, RunConfig)


def test_apply_to_legacy_globals_bridges_filter_config(monkeypatch):
    """apply_to_legacy_globals() must mutate filter_config attributes."""
    from ema import config as _cfg

    # Snapshot originals so other tests aren't polluted
    orig = (
        _cfg.filter_config.min_read,
        _cfg.filter_config.min_cells,
        _cfg.filter_config.min_pas_per_cell,
        _cfg.filter_config.min_genes,
    )
    try:
        cfg = RunConfig(min_read=999, min_cells=7, min_pas_per_cell=42)
        cfg.apply_to_legacy_globals()
        assert _cfg.filter_config.min_read == 999
        assert _cfg.filter_config.min_cells == 7
        assert _cfg.filter_config.min_pas_per_cell == 42
        # min_pas_per_cell also bridges to min_genes (preprocessing reads it)
        assert _cfg.filter_config.min_genes == 42
    finally:
        (_cfg.filter_config.min_read,
         _cfg.filter_config.min_cells,
         _cfg.filter_config.min_pas_per_cell,
         _cfg.filter_config.min_genes) = orig


def test_apply_to_legacy_globals_user_set_filter():
    """When user_set is given, only those fields are bridged."""
    from ema import config as _cfg

    orig_min_read = _cfg.filter_config.min_read
    try:
        cfg = RunConfig(min_read=4242, min_cells=99)
        # Only min_cells listed -> min_read should NOT be bridged
        cfg.apply_to_legacy_globals(user_set={"min_cells"})
        assert _cfg.filter_config.min_read == orig_min_read, (
            "min_read should not have been bridged because it wasn't in user_set"
        )
        assert _cfg.filter_config.min_cells == 99
    finally:
        _cfg.filter_config.min_read = orig_min_read


# ---------------------------------------------------------------------------
# click_options_from_schema generator
# ---------------------------------------------------------------------------

def test_click_options_from_schema_creates_options():
    """The decorator should attach a Click option for every cli_flag field."""
    import click
    from click.testing import CliRunner

    @click.command()
    @click_options_from_schema(RunConfig)
    def fake(**kwargs):
        click.echo(f"min_read={kwargs.get('min_read')}")

    runner = CliRunner()
    result = runner.invoke(fake, ["--help"])
    assert result.exit_code == 0
    # Sample a handful of flags from across the groups
    for flag in ("--min-read", "--peak-strategy", "--cluster-method",
                 "--match-method", "--threads", "--tile-size"):
        assert flag in result.output, f"{flag} missing from generated --help"


def test_click_options_from_schema_skip_omits_flags():
    """When skip=['threads'] the --threads option should be absent."""
    import click
    from click.testing import CliRunner

    @click.command()
    @click_options_from_schema(RunConfig, skip=["threads"])
    def fake(**kwargs):
        pass

    result = CliRunner().invoke(fake, ["--help"])
    assert "--threads" not in result.output


def test_click_option_default_matches_schema_default():
    """The Click default must be the schema field default verbatim."""
    import click
    from click.testing import CliRunner

    captured: dict = {}

    @click.command()
    @click_options_from_schema(RunConfig)
    def fake(**kwargs):
        captured.update(kwargs)

    CliRunner().invoke(fake, [])
    # Spot-check defaults
    assert captured["min_read"] == 1500
    assert captured["atlas_distance"] == 50
    assert captured["bam_threads"] == 4
    assert captured["lambda_fold_change"] == 2.0
    assert captured["pipeline"] is False
    assert captured["tiles"] is False
