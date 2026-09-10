# Releasing PeakATail

## Versioning

[Semantic Versioning](https://semver.org). While the project is `0.x`:

| Change | Bump |
|---|---|
| Breaking (see below) | **minor** — `0.3.0` → `0.4.0` |
| Everything else | **patch** — `0.3.0` → `0.3.1` |

After `1.0.0`, normal SemVer (breaking → major).

### What counts as breaking here

PeakATail is an analysis tool, so "breaking" is **wider** than an API change. It is
**anything that changes the numbers a user gets from the same command on the same
data.** A user does not care that the old number was wrong — they care that it moved.

Breaking therefore includes:

- **A changed default.** `--ip-filter-default auto` (0.3.0) made the internal-priming
  veto run with no flags: ~17% fewer PAS called. `--marker-top-n` 200 → 0 changed which
  tests run at all.
- **A corrected result.** The spliced-3'UTR de-duplication fixed a silently inflated
  denominator (32 reads reported where 16 was correct). Still breaking: re-running
  gives different output.
- **A renamed command or flag**, even with an alias (`ema` → `peakatail`).

Not breaking: new opt-in flags that default off, docs, tests, performance work whose
output is byte-identical.

## The version lives in one place

`pyproject.toml` is the source of truth. Everything else mirrors it:
`CITATION.cff`, `recipes/peakatail/meta.yaml`, the README citation block,
`docs/tutorials/01-installation.md`, `docs/concepts/output-files.md`, and the example
commands in `CONTRIBUTING.md`.

Never edit those by hand:

```bash
python scripts/bump_version.py --check   # verify every site agrees
python scripts/bump_version.py 0.4.0     # rewrite every site
```

`tests/test_version_sync.py` runs `--check`, so drift fails CI rather than reaching a
release.

## Release procedure

1. **Land everything** on `develop`; CI green.
2. **Bump** — `python scripts/bump_version.py X.Y.Z`.
3. **CHANGELOG** — collapse the themed `## Unreleased — <topic>` sections into a single
   `## X.Y.Z — YYYY-MM-DD`. Lead with any BREAKING paragraph.
4. **Release PR** `develop` → `main`; merge once green.
5. **Tag on `main`**:
   ```bash
   git tag -a vX.Y.Z -m "PeakATail vX.Y.Z"
   git push origin vX.Y.Z
   ```
6. `.github/workflows/release.yml` fires. In order, stopping at the first
   failure: it **verifies the tag matches `pyproject.toml`**, **verifies the
   release notes can be generated** from `CHANGELOG.md`, builds sdist + wheel,
   publishes to PyPI via trusted publishing (OIDC — no token stored), and then
   **creates the GitHub Release** with both distributions attached. Everything
   that can fail runs *before* the PyPI upload, because a version number can
   never be reused.
7. **bioconda** — after the sdist is on PyPI, replace the placeholder `sha256` in
   `recipes/peakatail/meta.yaml` with the real digest and submit.

### One-time setup (already done)

- **PyPI**: a *pending publisher* at pypi.org → Account → Publishing, with
  project `peakatail`, owner `BMGLab`, repo `PeakATail`, workflow `release.yml`,
  environment `pypi`.
- **GitHub**: a repository environment named `pypi`.

## A published version can never be reused

PyPI permanently burns a version number on first upload. A wrong upload can be yanked,
but the number is gone. That is why `release.yml` refuses to build when the tag and
`pyproject.toml` disagree — treat that guard as load-bearing.

## Citing a release

Results depend on the version, so the manuscript must cite an **exact tag and DOI**.
The GitHub–Zenodo integration (see `CONTRIBUTING.md`) mints the DOI, and it
fires on a published **GitHub Release** — *not* on the tag. A tag alone mints
nothing; `release.yml`'s `github_release` job is what creates that Release, so
if a DOI is missing after a release, check that job ran. Quote the exact tag
and DOI in the Methods.
