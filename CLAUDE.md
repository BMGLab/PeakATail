
## Branch hygiene — internal paths must not reach `main`

`main` is the published branch (GitHub landing page, release tarball, Zenodo
archive). Planning docs, internal reports and experiment scripts live on
`develop` only.

The list and the rules are in `docs/INTERNAL_ONLY_PATHS.md`. In short:

- **Never add** `PRIME_PLAN.md`, `HANDOFF.md`, `ROADMAP.md`, `CLAUDE.md`,
  `reports/`, `other_repos/`, `experiments/` or `test/` to `main`.
- **Never delete** them from `develop`.
- If a release merge reintroduces one on `main`, delete it again before
  tagging — do not restore the file to "resolve" the divergence.
