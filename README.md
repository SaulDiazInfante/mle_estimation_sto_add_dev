# Sargazo ML Estimation

This repository contains the maintained Fortran implementation for simulating the SDE system, estimating model parameters, and visualizing estimator trajectories.

## Repository layout

- `src/app`: executable entry points.
- `src/modules`: maintained Fortran modules, one module per file.
- `legacy`: preserved legacy Fortran sources kept out of the main build.
- `data/input`: reserved location for future input datasets and configuration files.
- `data/output`: generated CSV and text outputs from simulation and estimation runs.
- `visualization/scripts`: plotting utilities.
- `visualization/plots`: generated figures.
- `docs`: project documentation and engineering notes.
- `tests`: smoke tests and future regression tests.
- `.github/workflows`: CI definitions.

## Common commands

- `make build`: compile the maintained application into `build/bin/`.
- `make run`: run the default workflow and write timestamped outputs under `data/output/`.
- `make plot`: plot the latest generated estimator CSV already present in `data/output/`.
- `make docs`: build the Doxygen HTML site locally under `build/docs/doxygen/html`.
- `make test`: run the repository smoke test used by CI.
- `make test-unit`: run deterministic unit tests.
- `make test-smoke`: run the end-to-end smoke test only.
- `make check-large-files`: fail if tracked files exceed the repository size policy.
- `make setup-git-hooks`: enable the tracked pre-commit hook path for this clone.
- `make clean`: remove build products and timestamped CSV/PNG outputs.

Generated filenames use an ISO 8601-style timestamp prefix, for example
`2026-03-17T12:47:05_estimator_trajectory.csv`.

## Runtime configuration

The driver reads environment variables so tests and CI can run a smaller case without editing source code:

- `SARGAZO_N_OBSERVATIONS`
- `SARGAZO_MINIMUM_TRAJECTORY_OBSERVATIONS`
- `SARGAZO_REQUESTED_TRAJECTORY_POINTS`
- `SARGAZO_TIME_STEP`
- `SARGAZO_BETA`
- `SARGAZO_THETA`
- `SARGAZO_SIGMA`
- `SARGAZO_SEED`
- `SARGAZO_OUTPUT_TIMESTAMP`
- `SARGAZO_WRITE_STATE_HISTORY`
- `SARGAZO_STATE_HISTORY_FILE`
- `SARGAZO_ESTIMATOR_HISTORY_FILE`

You can also provide `TIMESTAMP=2026-03-17T12:47:05` to `make run` or `make plot`
to force a specific output prefix for `make run`.

`make plot` does not rerun the Fortran application. If there is no estimator CSV in
`data/output`, it stops with an error asking you to run `make run` first.

Long simulations print progress updates from the time-stepping loop at regular
checkpoints so you can monitor the run while it is executing.

`make test` runs both the unit-test layer and the fast smoke test used by CI.

## Git tracking policy

- Keep source code, docs, CI config, and small deterministic tests in git.
- Keep generated outputs in `data/output/`, `visualization/plots/`, and `tests/artifacts/`, which are ignored.
- Keep local or raw input data in the ignored paths under `data/input/`.
- The repository includes a size guard that rejects staged files larger than `50 MB` by default.

Run `make setup-git-hooks` once per clone to activate the local pre-commit size check.

## GitHub workflow

- CI runs on pull requests and pushes to `main`.
- Doxygen documentation is built on pull requests and published to GitHub Pages from `main`.
- Tagged releases like `v0.1.0` publish a Linux build artifact through GitHub Releases.
- Contribution templates and issue forms are included under `.github/`.

See [docs/github_publish.md](docs/github_publish.md) for the exact publish steps and [docs/testing_strategy.md](docs/testing_strategy.md) for the TDD-oriented test layout.
