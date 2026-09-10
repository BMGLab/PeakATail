"""Syntax + shape checks for the release-engineering assets.

These guard the tag-triggered PyPI workflow, the bioconda recipe skeleton, and
their agreement with pyproject.toml. They are pure static parsing -- no network,
no build -- so they run anywhere the repo is checked out.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

REPO_ROOT = Path(__file__).resolve().parents[1]
RELEASE_WF = REPO_ROOT / ".github" / "workflows" / "release.yml"
RECIPE = REPO_ROOT / "recipes" / "peakatail" / "meta.yaml"
PYPROJECT = REPO_ROOT / "pyproject.toml"
CONTRIBUTING = REPO_ROOT / "CONTRIBUTING.md"


def _load_pyproject() -> dict:
    try:
        import tomllib  # py311+
    except ModuleNotFoundError:  # pragma: no cover - local dev on <3.11
        tomllib = pytest.importorskip("tomli")
    return tomllib.loads(PYPROJECT.read_text())


def _render_recipe_jinja(text: str) -> str:
    """Render the conda meta.yaml Jinja so it parses as plain YAML.

    conda recipes use `{% set %}` and `{{ }}` plus build-time placeholders like
    {{ PYTHON }}. We substitute the recipe's own set-vars and blank the rest so
    yaml.safe_load sees valid YAML.
    """
    setvars: dict[str, str] = {}
    for name, value in re.findall(r"{%\s*set\s+(\w+)\s*=\s*\"([^\"]*)\"\s*%}", text):
        setvars[name] = value
    # drop the {% set %} statement lines
    text = re.sub(r"{%.*?%}", "", text)

    def repl(match: re.Match) -> str:
        expr = match.group(1).strip()
        # `name|lower`, `name[0]`, `name`, `version`
        base = re.split(r"[|\[]", expr)[0].strip()
        val = setvars.get(base)
        if val is None:
            return "PLACEHOLDER"
        if expr.endswith("|lower"):
            return val.lower()
        idx = re.search(r"\[(\d+)\]", expr)
        if idx:
            return val[int(idx.group(1))]
        return val

    return re.sub(r"{{\s*(.*?)\s*}}", repl, text)


def test_release_workflow_parses_and_is_tag_triggered() -> None:
    doc = yaml.safe_load(RELEASE_WF.read_text())
    assert doc["name"]
    # `on:` round-trips to the bool key True under YAML 1.1 -- accept either.
    on = doc.get("on", doc.get(True))
    assert "tags" in on["push"]
    assert any(t.startswith("v") for t in on["push"]["tags"])

    jobs = doc["jobs"]
    assert "build" in jobs and "publish" in jobs
    publish = jobs["publish"]
    assert publish.get("needs") == "build"
    # OIDC / trusted publishing must be wired: id-token write + the publish action
    assert publish["permissions"]["id-token"] == "write"
    steps_using = [s.get("uses", "") for s in publish["steps"]]
    assert any("gh-action-pypi-publish" in u for u in steps_using)
    # build job must actually build the dists
    build_run = " ".join(s.get("run", "") for s in jobs["build"]["steps"])
    assert "python -m build" in build_run


def test_recipe_parses_and_matches_pyproject() -> None:
    rendered = _render_recipe_jinja(RECIPE.read_text())
    doc = yaml.safe_load(rendered)

    proj = _load_pyproject()["project"]

    assert doc["package"]["name"] == proj["name"]
    assert str(doc["package"]["version"]) == str(proj["version"])

    # the recipe must declare EVERY console script pyproject declares, with the
    # same target -- including the deprecated `ema` alias, whose omission would
    # silently stop shipping it to conda users.
    eps = {e.split("=")[0].strip(): e.split("=", 1)[1].strip()
           for e in doc["build"]["entry_points"]}
    for script_name, target in proj["scripts"].items():
        assert script_name in eps, (
            f"pyproject declares the console script {script_name!r} but the conda "
            f"recipe does not; conda users would not get it"
        )
        assert eps[script_name] == target, (
            f"{script_name}: recipe points at {eps[script_name]!r}, "
            f"pyproject at {target!r}"
        )
    assert "peakatail" in eps, "the primary command must be `peakatail`"

    # every runtime dep from pyproject appears in the recipe run section, under
    # its conda channel name where that differs from the PyPI name.
    conda_name = {
        "igraph": "python-igraph",
        "pyyaml": "pyyaml",
        "kaleido": "python-kaleido",
        "matplotlib": "matplotlib-base",
    }
    run_names = {line.split()[0].lower() for line in doc["requirements"]["run"]}
    for dep in proj["dependencies"]:
        pypi_name = re.split(r"[<>=!~ ]", dep, maxsplit=1)[0].strip().lower()
        expected = conda_name.get(pypi_name, pypi_name)
        assert expected in run_names, f"{pypi_name} -> {expected} missing from recipe run deps"


def test_recipe_sha256_is_a_placeholder_todo() -> None:
    # The skeleton ships with a zero digest on purpose; make the TODO explicit so
    # nobody mistakes it for a real, verified hash.
    text = RECIPE.read_text()
    assert re.search(r"sha256:\s*0{64}", text)
    assert "TODO" in text


def test_contributing_documents_zenodo_release() -> None:
    text = CONTRIBUTING.read_text().lower()
    assert "zenodo" in text
    assert "doi" in text
    assert "github release" in text
