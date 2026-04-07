(contributing)=
# Contributing to TSMP-PDAF

## Workflow

1. **Create a feature branch** from `tsmp-pdaf-patched`:
   ```
   git checkout -b tsmp-pdaf-patched-<short-description> tsmp-pdaf-patched
   ```
2. **Develop** on the feature branch. Keep commits focused; one logical
   change per commit.
3. **Open a pull request** against `tsmp-pdaf-patched` on GitHub when
   the work is ready for review.
4. **Address review feedback** with additional commits on the same branch.
5. **Squash-merge** into `tsmp-pdaf-patched` once approved (see below).

## Branch naming

Use the prefix `tsmp-pdaf-patched-` followed by a short, hyphenated
description, e.g. `tsmp-pdaf-patched-snow-da` or
`tsmp-pdaf-patched-fix-grace-localization`.

## Commits and merge strategy

- Write commit messages in the imperative mood and keep the subject line
  under 72 characters.
- **Always squash-merge** feature branches into `tsmp-pdaf-patched`.
  The merge commit message should summarise the entire change, not list
  the individual development commits. This keeps the history of the
  protected branch linear and readable.

## What a PR must include

| Concern | Requirement |
|---|---|
| **Tests** | Run the CLM3.5-PDAF, CLM3.5-ParFlow-PDAF, and eCLM-PDAF test suite and confirm it passes. |
| **Documentation** | Update `docs/users_guide/` for any user-facing change (new `enkfpf.par` parameter, changed behaviour). |
| **Interface changes** | Coordinate with the eCLM and ParFlow repository maintainers before merging anything that changes the API visible to component models. |
| **PDAF library changes** | Avoid changes in `src/`. If unavoidable, collaborate with repository administrators and document the change in {ref}`dev`. |

## Collaborating with repository administrators

Two activities require administrator involvement and must not be done
unilaterally:

- **PDAF version updates** — pulling a new upstream PDAF `master` into
  `tsmp-pdaf-patched`.
- **PDAF library changes** — any modification inside `src/`.

Open an issue or contact the administrators directly before starting
work that falls into these categories.
