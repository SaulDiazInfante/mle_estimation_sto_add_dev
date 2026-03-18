# Repository Design Notes

## Goals

- Keep maintained source code separate from generated artifacts.
- Support a clean git workflow by ignoring build and output files.
- Make room for incremental TDD by giving tests a stable home and a reproducible smoke target.
- Keep CI simple: build, run a lightweight end-to-end case, and validate the plotting path.

## Layout

- `src/app/main.f90`: production entry point.
- `src/modules/*.f90`: reusable simulation, estimation, validation, and output modules.
- `legacy/*.f90`: historical reference implementations not included in the default build.
- `data/input/`: future experiment inputs or configuration snapshots.
- `data/output/`: generated simulation and estimation artifacts.
- `visualization/scripts/`: plotting code.
- `visualization/plots/`: generated figures.
- `tests/`: unit and smoke tests, plus ignored test artifacts.
- `build/`: compiler outputs, module files, objects, and binaries.
- `docs/Doxyfile`: Doxygen configuration for the generated API/reference site.
- `.github/workflows/`: CI for pull requests and release automation for version tags.

## TDD and CI direction

- New features should start with a failing test under `tests/` whenever practical.
- Fast smoke coverage is preferred in CI; heavy numerical experiments should stay opt-in.
- Runtime parameters are controlled through environment variables so tests can be small and deterministic.
- Generated CSV and PNG artifacts use an ISO 8601 timestamp prefix so simulation runs stay traceable.
- The GitHub Actions workflows mirror local usage through `make build`, `make test`, `make docs`, and tagged release packaging.
