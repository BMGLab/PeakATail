# Contributing to PeakATail

Thanks for helping improve PeakATail. This file focuses on the **release
process**; for day-to-day development, install the test extras and run the
suite:

```bash
pip install -e '.[test]'
pytest -q
```

PRs target `develop`. CI (`.github/workflows/ci.yml`) runs the pytest suite on
every PR to `develop`/`main` across Python 3.11 and 3.12.

## Cutting a release

Releases are **tag-triggered**. The single source of truth for the version is
`version` in `pyproject.toml`.

1. Bump `version` in `pyproject.toml` and, in lockstep, `version` in
   `CITATION.cff` and `{% set version %}` in `recipes/peakatail/meta.yaml`.
2. Update `CHANGELOG.md`.
3. Merge to `main`, then tag and push:

   ```bash
   git tag -a v0.3.0 -m "PeakATail v0.3.0"
   git push origin v0.3.0
   ```

4. The tag push runs `.github/workflows/release.yml`, which does three things
   in order and stops at the first failure:

   1. **Verifies the tag against `pyproject.toml`.** A `vX.Y.Z` tag on a tree
      whose version says something else fails here, before anything is
      published — PyPI never lets a version number be reused, so a wrong upload
      is permanent.
   2. **Builds and publishes to [PyPI](https://pypi.org/project/peakatail/)**
      (`python -m build`, then `twine check`) via **trusted publishing
      (OIDC)** — no API token is stored in the repo. This must be configured
      once on PyPI (project → Settings → Publishing → add a GitHub publisher
      for `owner=BMGLab repo=PeakATail workflow=release.yml environment=pypi`).
   3. **Creates the GitHub Release** from the tag, with the matching
      `CHANGELOG.md` section as the body and the sdist + wheel attached. This
      is the step Zenodo watches (see below). Its notes come from
      `scripts/changelog_section.py`, which **fails the job** if `CHANGELOG.md`
      has no `## <version>` heading — so collapsing `## Unreleased` into
      `## X.Y.Z` in step 2 is not optional.

### bioconda (after the PyPI release exists)

`recipes/peakatail/meta.yaml` is a submission-ready skeleton. Once the sdist is
on PyPI, compute its `sha256`, drop it into the recipe, and open a PR against
[`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes).

## Zenodo DOI on release (for the paper's Availability statement)

A Q1 methods submission must cite an archived, versioned release. Zenodo's
GitHub integration mints a DOI automatically for each GitHub Release.

**One-time setup (repo admin):**

1. Sign in to [Zenodo](https://zenodo.org/) with GitHub (or link the accounts
   under Zenodo → *Account* → *GitHub*).
2. On the Zenodo *GitHub* page, flip the toggle **ON** for `BMGLab/PeakATail`.
   This installs a webhook; only releases created **after** the toggle is on
   are archived.
3. Zenodo immediately reserves a **concept DOI** (one stable DOI that always
   resolves to the latest version) plus a per-release **version DOI**.

**Each release:**

1. The `github_release` job in `release.yml` publishes the GitHub Release for
   you, and the Zenodo webhook fires on that publish. Nothing manual here.
2. Zenodo archives the repository's auto-generated **source zipball** for
   that tag (not the sdist/wheel attached to the Release) and issues the
   version DOI. Grab the badge from
   the Zenodo record and (optionally) add it to `README.md`.
3. Add the DOI to `CITATION.cff` (`doi:` and/or an `identifiers:` entry) and to
   the manuscript's *Data/Code Availability* statement. Cite the **concept DOI**
   for "PeakATail" in general and the **version DOI** for the exact release used
   in the paper.

> The GitHub **Release**, not the git tag, is what triggers Zenodo. Pushing a
> tag while `release.yml` is broken or disabled therefore mints no DOI, however
> healthy the tag looks. If a DOI is missing after a release, check that the
> `github_release` job ran and that the Release exists on the Releases page.

> Zenodo reads `CITATION.cff` for the archived record's authors and title.
> Whatever is in that file at tag time becomes the public citation metadata for
> that DOI, so settle the author list **before** tagging.
