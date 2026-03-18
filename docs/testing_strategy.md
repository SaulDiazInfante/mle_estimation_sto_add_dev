# Testing Strategy

## Goals

- Keep deterministic logic covered by small unit tests.
- Keep one end-to-end smoke test that exercises build, run, and plotting.
- Make test execution fast enough for every pull request.
- Prevent large-file drift from undermining the git and CI workflow.

## Current layers

- `tests/unit/test_core.f90`: deterministic unit coverage for core helper logic.
- `tests/run_unit_tests.sh`: compiles and runs the unit test binary against the built modules.
- `tests/smoke_test.sh`: small end-to-end run used in CI.
- `scripts/check_large_files.sh`: repository policy check for tracked file sizes.

## TDD workflow

1. Reproduce the desired behavior with a failing test.
2. Implement the smallest change needed to pass.
3. Run `make test`.
4. Refactor with tests still passing.

## Future improvements

- Add regression tests for estimator outputs on fixed seeds.
- Add golden-file checks for generated CSV structure when new outputs are introduced.
- Separate numerical correctness tests from performance-oriented tests so CI stays fast.
