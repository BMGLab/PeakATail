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

4. The tag push runs `.github/workflows/release.yml`, which builds the sdist +
   wheel (`python -m build`) and publishes them to
   [PyPI](https://pypi.org/project/peakatail/) via **trusted publishing
   (OIDC)** — no API token is stored in the repo. This must be configured once
   on PyPI (project → Settings → Publishing → add a GitHub publisher for
   `owner=BMGLab repo=PeakATail workflow=release.yml environment=pypi`).

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

1. After the tag/PyPI step above, create a **GitHub Release** from that tag
   (Releases → *Draft a new release* → pick the tag → publish). The Zenodo
   webhook fires on *publish*.
2. Zenodo archives the tarball and issues the version DOI. Grab the badge from
   the Zenodo record and (optionally) add it to `README.md`.
3. Add the DOI to `CITATION.cff` (`doi:` and/or an `identifiers:` entry) and to
   the manuscript's *Data/Code Availability* statement. Cite the **concept DOI**
   for "PeakATail" in general and the **version DOI** for the exact release used
   in the paper.

> The GitHub Release, not just the git tag, is what triggers Zenodo — remember
> to publish it.
