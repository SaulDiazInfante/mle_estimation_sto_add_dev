# Contributing Guide

## Development workflow

1. Create a branch from `main`.
2. Write or update a test before changing production code whenever practical.
3. Run `make test` locally before opening a pull request.
4. Open a pull request that explains the motivation, the test coverage, and any numerical or performance impact.

## TDD expectations

- Prefer small unit tests under `tests/unit` for deterministic logic.
- Keep `tests/smoke_test.sh` fast and representative of the end-to-end workflow.
- When fixing a bug, add a test that would have failed before the fix.
- Avoid merging changes that only add features without test coverage unless there is a clear reason documented in the PR.

## CI/CD expectations

- Pull requests should pass the GitHub Actions `CI` workflow before merge.
- Versioned releases should be created by pushing a tag like `v0.1.0`.
- Tagged releases trigger the release workflow, which rebuilds the project, runs tests, and publishes a Linux artifact to GitHub Releases.

## Coding notes

- Keep one maintained Fortran module per file under `src/modules`.
- Put executable entry points under `src/app`.
- Generated data and plots belong in ignored output directories, not in version control.
- Preserve the `legacy` folder for reference implementations, but do not wire those files into the default build unless the change is intentional.
