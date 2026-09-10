# Paths that belong on `develop`, never on `main`

`main` is the published branch: it is what people browse on GitHub, what the
release tarball and the Zenodo archive contain, and what a reviewer of the
manuscript sees first. Internal planning material does not belong there.

These paths are kept on `develop` and **deleted from `main`**:

| Path | Why it is internal |
|---|---|
| `PRIME_PLAN.md` | branch plan for `peakAtail-prime` |
| `HANDOFF.md` | engine-branch handoff notes |
| `ROADMAP.md` | development & publication roadmap |
| `CLAUDE.md` | AI agent configuration |
| `reports/` | internal analysis reports and figures |
| `other_repos/` | comparison notes on other tools |
| `experiments/` | Nextflow experiment scripts and configs |
| `test/` | the pre-restructure test directory; `testpaths = ["tests"]`, so nothing here ever runs |

## Rule for anyone — human or agent — working on this repo

1. **Do not add these paths to `main`.** If a release merge reintroduces one,
   delete it again on `main` before tagging; do not "fix" the divergence by
   restoring the file.
2. **Do not delete them from `develop`.** They are working material, not dead
   weight. `develop` is where they live.
3. **New internal material goes in one of the directories above**, not at the
   repository root, so it is stripped by the same rule rather than needing a
   new decision each time.
4. When adding a new internal path, add it to this table **and** to
   `INTERNAL_ONLY_PATHS` in `tests/test_internal_paths_absent_from_main.py`.

## Why the two branches differ here

`main` and `develop` are otherwise kept identical. This is the one deliberate
exception. A release merge (`develop` -> `main`) leaves these deletions in
place on its own, because `develop` does not modify the files relative to the
merge base. If `develop` *does* edit one, the merge raises a modify/delete
conflict: resolve it by **keeping the deletion**.
