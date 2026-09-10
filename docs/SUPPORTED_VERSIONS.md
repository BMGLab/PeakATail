# Supported Python and dependency versions

## The rule

**Every version this project claims to support is tested in CI.** A
`requires-python` floor or a dependency lower bound that no job exercises is a
guess, and it can be wrong in either direction — too high locks users out for
no reason, too low promises a version that is quietly broken.

Three tests in `tests/test_release_packaging.py` enforce this:

| Test | Prevents |
|---|---|
| `test_ci_tests_the_python_version_the_package_claims_to_support` | a floor no job runs |
| `test_ci_does_not_test_below_the_declared_floor` | CI implying support `pip` will refuse |
| `test_dependency_floors_are_verified_by_a_ci_job` | lower bounds nothing installs |

## Python

Currently **`>=3.10`**, tested on 3.10, 3.11 and 3.12.

`ema/annotate/gtftobed.py` uses a `match` statement, which is 3.10+. That is
the real constraint; the floor previously read `>=3.11`, one release higher
than anything in the code required.

### On SPEC 0

The scientific-Python ecosystem follows
[SPEC 0](https://scientific-python.org/specs/spec-0000/): drop Python versions
three years after release, core dependencies two years after release. By that
schedule 3.10 (released October 2021) is already out of support, and it reaches
upstream end-of-life in **October 2026**.

Supporting 3.10 past that line is a deliberate exception, not an oversight —
the project's own benchmark environment
(`reports/apabenchmark_runbook/environment.yml`) pins `python=3.10`. Revisit it
when a core dependency drops 3.10, which will force the issue regardless.

## Dependency lower bounds

CI installs *latest* versions, which only ever proves "works with today's
releases". The `lowest-direct` job installs what `pyproject.toml` actually
**promises**, via `uv pip install --resolution lowest-direct`, and runs the
full suite against it.

It found two bounds that were wrong at 0.3.0:

- **`pybedtools>=0.10`** — 0.10 and 0.11 ship no wheels for the supported
  interpreters, and their sdist build fails without setuptools present. The
  install simply fails. Raised to `>=0.12.1`.
- **`pyfaidx>=0.7`** — 0.7.0 imports `pkg_resources`, absent from modern
  environments, so it fails to import. This one was worse than an install
  error: the internal-priming filter caught the ImportError, logged, and
  **skipped**, so a user got a call set roughly 17 % larger than intended with
  no failure. Raised to `>=0.7.2.2`.

### Raising a bound

Set it to the lowest version you have actually verified, not the newest
release. `--resolution lowest-direct` is what verifies it.
