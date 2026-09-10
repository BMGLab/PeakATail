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


# ---------------------------------------------------------------------------
# The guards that protect a PERMANENT PyPI version and a permanent DOI.
#
# These were added after a review found that a `workflow_dispatch` run on a
# branch would upload to PyPI while the tag-guarded Release job was skipped --
# version burned, no Release, no DOI, every step green. A later review then
# found that deleting each of those guards left the entire suite passing: the
# only test parsing this workflow checked jobs/needs/id-token/the publish
# action, so the fixes themselves shipped unprotected. These pin them.
# ---------------------------------------------------------------------------

def _release_jobs() -> dict:
    return yaml.safe_load(RELEASE_WF.read_text())["jobs"]


def _step_names(job: dict) -> list[str]:
    return [s.get("name", "") for s in job["steps"]]


def test_publish_cannot_run_off_a_tag() -> None:
    """Without this, workflow_dispatch on a branch publishes to PyPI while
    `github_release` is skipped -- burning the version with no DOI."""
    guard = _release_jobs()["publish"].get("if", "")
    assert "refs/tags/" in guard, (
        "the `publish` job has no tag guard, so a non-tag run can upload to "
        f"PyPI: if={guard!r}"
    )


def test_github_release_job_exists_and_is_tag_guarded() -> None:
    """Zenodo mints a DOI from a published GitHub Release, not from a tag."""
    jobs = _release_jobs()
    assert "github_release" in jobs, (
        "no job creates a GitHub Release, so Zenodo would never fire and no "
        "DOI would be minted"
    )
    rel = jobs["github_release"]
    assert "refs/tags/" in rel.get("if", "")
    assert rel["permissions"]["contents"] == "write"
    assert "publish" in rel["needs"], (
        "the Release must not be created for a version that failed to reach PyPI"
    )


def test_everything_that_can_fail_runs_before_the_upload() -> None:
    """PyPI never lets a version be reused, so all checks live in `build`."""
    build = _release_jobs()["build"]
    names = " | ".join(_step_names(build)).lower()
    for needle, why in [
        ("tag matches", "tag-vs-pyproject check"),
        ("tag is on main", "tag ancestry check"),
        ("test suite", "the test suite on the tagged tree"),
        ("release notes", "release-note generation"),
    ]:
        assert needle in names, (
            f"the build job is missing {why}; it would run after the "
            f"irreversible PyPI upload, or not at all. Steps: {names}"
        )

    # ...and they must precede the build itself.
    order = _step_names(build)
    build_idx = next(i for i, n in enumerate(order) if "Build distributions" in n)
    for needle in ("tag matches", "tag is on main", "test suite", "release notes"):
        idx = next(i for i, n in enumerate(order) if needle in n.lower())
        assert idx < build_idx, f"{order[idx]!r} runs after the build step"


def test_release_is_serialised_and_survives_a_repushed_tag() -> None:
    doc = yaml.safe_load(RELEASE_WF.read_text())
    assert "concurrency" in doc, (
        "no concurrency group: a double tag push races two PyPI uploads"
    )
    assert doc["concurrency"].get("cancel-in-progress") is False, (
        "a release must never be cancelled mid-publish"
    )
    publish_steps = _release_jobs()["publish"]["steps"]
    pypi = next(s for s in publish_steps if "gh-action-pypi-publish" in s.get("uses", ""))
    assert pypi.get("with", {}).get("skip-existing") is True, (
        "without skip-existing a re-pushed tag fails on the duplicate upload "
        "and `github_release` (needs: publish) never runs -- package on PyPI, "
        "no Release, no DOI"
    )


def test_release_workflow_shell_steps_use_valid_git_flags() -> None:
    """A CI-only shell bug costs a whole release cycle.

    `git fetch --depth=0` is not "no limit" -- it is an error ("depth 0 is not
    a positive number"), and it failed the very ancestry guard it was meant to
    implement. Nothing was published, because that guard runs before the PyPI
    upload, but the tag had to be deleted and re-cut. Cheap to pin.
    """
    text = RELEASE_WF.read_text()
    assert "--depth=0" not in text, (
        "`--depth=0` is invalid for git fetch; use `fetch-depth: 0` on "
        "actions/checkout and a plain `git fetch` here"
    )
    assert "--depth 0" not in text, "same, spelled with a space"


def test_release_workflow_installs_the_same_system_deps_as_ci() -> None:
    """The release suite must run in the same environment CI validates.

    release.yml runs the full suite on the tagged tree, but it did not install
    bedtools -- which ci.yml does. Several tests shell out to bedtools and
    ERROR (not skip) without it, so the release run failed on code that was
    green in CI. Any apt package CI installs must be installed here too.
    """
    ci = (RELEASE_WF.parent / "ci.yml").read_text()
    rel = RELEASE_WF.read_text()

    def apt_packages(text: str) -> set[str]:
        """Packages from real `apt-get install` RUN lines, not prose.

        Matching the bare package name anywhere in the file is useless: the
        comment explaining why bedtools is needed also contains "bedtools", so
        deleting the install step left this test green.
        """
        found = set()
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith("#") or "apt-get install" not in stripped:
                continue
            tail = stripped.split("apt-get install", 1)[1]
            for token in tail.split():
                if not token.startswith("-") and token not in {"&&", "|"}:
                    found.add(token)
        return found

    packages = apt_packages(ci)
    assert packages, "could not parse CI's apt packages; update this test"
    missing = packages - apt_packages(rel)
    assert not missing, (
        f"ci.yml installs {sorted(missing)} but release.yml does not, so the "
        "test suite runs in a different environment during a release"
    )


# ---------------------------------------------------------------------------
# `requires-python` is a PROMISE. CI is what verifies it.
#
# The floor said >=3.11 while CI tested 3.11 and 3.12, so 3.10 was neither
# supported nor proven broken -- just untested. It in fact worked: the package
# has one `match` statement (3.10+) and no 3.11-only syntax, so the floor had
# overshot by one release. An unverified floor can be wrong in both directions:
# too high locks users out for no reason, too low promises a broken version.
# ---------------------------------------------------------------------------

def _declared_python_floor() -> tuple[int, int]:
    text = (RELEASE_WF.parents[2] / "pyproject.toml").read_text()
    match = re.search(r'requires-python\s*=\s*">=(\d+)\.(\d+)"', text)
    assert match, "pyproject.toml has no `requires-python = \">=X.Y\"`"
    return int(match.group(1)), int(match.group(2))


def _ci_matrix_versions() -> list[tuple[int, int]]:
    ci = yaml.safe_load((RELEASE_WF.parent / "ci.yml").read_text())
    versions = ci["jobs"]["test"]["strategy"]["matrix"]["python-version"]
    out = []
    for v in versions:
        parts = str(v).split(".")
        out.append((int(parts[0]), int(parts[1])))
    return out


def test_ci_tests_the_python_version_the_package_claims_to_support():
    """The declared floor must actually be exercised by CI."""
    floor = _declared_python_floor()
    matrix = _ci_matrix_versions()
    assert floor in matrix, (
        f"pyproject declares requires-python >={floor[0]}.{floor[1]} but the CI "
        f"matrix is {['.'.join(map(str, m)) for m in matrix]}. The lowest "
        "supported version must be tested, or the floor is an unverified claim."
    )


def test_ci_does_not_test_below_the_declared_floor():
    """The mirror image: CI must not imply support the package disclaims."""
    floor = _declared_python_floor()
    below = [m for m in _ci_matrix_versions() if m < floor]
    assert not below, (
        f"CI tests {['.'.join(map(str, m)) for m in below]}, below the declared "
        f"floor >={floor[0]}.{floor[1]}. Either lower requires-python or drop "
        "those jobs -- passing tests on an unsupported version is a promise "
        "pip will refuse to honour."
    )


def test_dependency_floors_are_verified_by_a_ci_job():
    """A lower bound nothing installs is a guess.

    CI installing only latest versions proves "works with today's releases".
    It said `pybedtools>=0.10` (no wheels, unbuildable) and `pyfaidx>=0.7`
    (fails to import -- `pkg_resources`), and both were invisible until a job
    resolved the declared minimums.
    """
    ci = yaml.safe_load((RELEASE_WF.parent / "ci.yml").read_text())
    jobs = ci["jobs"]
    lowest = [j for j in jobs.values()
              if any("lowest-direct" in str(s.get("run", "")) for s in j.get("steps", []))]
    assert lowest, (
        "no CI job installs with `--resolution lowest-direct`, so the declared "
        "dependency lower bounds are never exercised"
    )
    steps = " ".join(str(s.get("run", "")) for s in lowest[0]["steps"])
    assert "pytest" in steps, (
        "the lowest-direct job resolves the floors but never runs the tests "
        "against them, so it only proves they install"
    )


# ---------------------------------------------------------------------------
# Internal planning material must not reach the published branch.
# ---------------------------------------------------------------------------

def _internal_paths_from_docs() -> set[str]:
    """The canonical list, read from the table in docs/INTERNAL_ONLY_PATHS.md."""
    doc = RELEASE_WF.parents[2] / "docs" / "INTERNAL_ONLY_PATHS.md"
    text = doc.read_text()
    paths = set()
    for row in re.finditer(r"^\|\s*`([^`]+)`\s*\|", text, re.MULTILINE):
        paths.add(row.group(1).rstrip("/"))
    assert paths, "could not parse the internal-paths table"
    return paths


def test_release_refuses_to_publish_internal_only_paths():
    """Every path the doc calls internal must be checked before the upload.

    main is what the Zenodo archive contains, so this is checked on the tagged
    tree -- and, like every other release guard, before the irreversible step.
    """
    build = yaml.safe_load(RELEASE_WF.read_text())["jobs"]["build"]
    steps = build["steps"]
    names = [s.get("name", "") for s in steps]
    guard = [i for i, n in enumerate(names) if "internal-only" in n.lower()]
    assert guard, f"no internal-only path guard in the build job: {names}"

    body = steps[guard[0]].get("run", "")
    missing = {p for p in _internal_paths_from_docs() if p not in body}
    assert not missing, (
        f"docs/INTERNAL_ONLY_PATHS.md lists {sorted(missing)} but the release "
        "guard does not check them, so they could be published"
    )

    build_idx = next(i for i, n in enumerate(names) if "Build distributions" in n)
    assert guard[0] < build_idx, "the guard runs after the build"
