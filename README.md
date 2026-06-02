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
- `make run-monte-carlo`: run the Monte Carlo study driver and write summary plus replicate CSV files under `data/output/`.
- `make run-monte-carlo-paper`: run the paper-style study preset with the selected sweep over `dt`, basis size, and observation count.
- `make run-snapshot-comparison`: reconstruct one stochastic path and its deterministic counterpart at selected times.
- `make plot`: plot the latest generated estimator CSV already present in `data/output/`.
- `make plot-snapshot-comparison`: plot the latest reconstructed snapshot comparison CSV.
- `make video-snapshot-comparison`: animate the latest reconstructed snapshot comparison CSV.
- `make docs`: build the Doxygen HTML site locally under `build/docs/doxygen/html`.
- `make test`: run the repository smoke test used by CI.
- `make test-unit`: run deterministic unit tests.
- `make test-smoke`: run the end-to-end smoke test only.
- `make test-monte-carlo`: run the Monte Carlo smoke test only.
- `make test-snapshot-comparison`: run the reconstructed snapshot smoke test only.
- `make check-large-files`: fail if tracked files exceed the repository size policy.
- `make setup-git-hooks`: enable the tracked pre-commit hook path for this clone.
- `make clean`: remove build products and timestamped CSV/PNG outputs.

Generated filenames use a filesystem-safe ISO 8601 basic timestamp prefix, for example
`20260317T124705_estimator_trajectory.csv`.

## Runtime configuration

The driver reads environment variables so tests and CI can run a smaller case without editing source code:

- `SARGAZO_N_OBSERVATIONS`
- `SARGAZO_MINIMUM_TRAJECTORY_OBSERVATIONS`
- `SARGAZO_REQUESTED_TRAJECTORY_POINTS`
- `SARGAZO_TIME_STEP`
- `SARGAZO_GRID_NX`
- `SARGAZO_GRID_NY`
- `SARGAZO_VELOCITY_MODE_X`
- `SARGAZO_VELOCITY_MODE_Y`
- `SARGAZO_LENGTH_X`
- `SARGAZO_LENGTH_Y`
- `SARGAZO_GAMMA`
- `SARGAZO_BETA`
- `SARGAZO_THETA`
- `SARGAZO_SIGMA`
- `SARGAZO_SEED`
- `SARGAZO_OUTPUT_TIMESTAMP`
- `SARGAZO_WRITE_STATE_HISTORY`
- `SARGAZO_STATE_HISTORY_FILE`
- `SARGAZO_ESTIMATOR_HISTORY_FILE`
- `SARGAZO_SNAPSHOT_TIMES`
- `SARGAZO_SNAPSHOT_FRAME_COUNT`
- `SARGAZO_SNAPSHOT_INITIAL_TIME`
- `SARGAZO_SNAPSHOT_FINAL_TIME`
- `SARGAZO_SNAPSHOT_GRID_NX`
- `SARGAZO_SNAPSHOT_GRID_NY`
- `SARGAZO_SNAPSHOT_COMPARISON_FILE`

You can also provide `TIMESTAMP=20260317T124705` to `make run`
to force a specific output prefix. If you pass an older extended form such as
`2026-03-17T12:47:05`, the application normalizes it before creating filenames.

`make plot` does not rerun the Fortran application. If there is no estimator CSV in
`data/output`, it stops with an error asking you to run `make run` first.

Long simulations print progress updates from the time-stepping loop at regular
checkpoints so you can monitor the run while it is executing.

`make test` runs both the unit-test layer and the fast smoke test used by CI.

## Monte Carlo study

The Monte Carlo driver uses the same simulation and MLE pipeline as `make run`,
but sweeps over combinations of trajectory length, integration step size, and
square spectral basis resolution. The study writes two CSV files:

- a summary table with per-combination means, standard deviations, and average runtimes;
- a replicate table with one row per simulated dataset, suitable for histograms.

The study driver reads these additional environment variables:

- `SARGAZO_MC_REPLICATES`
- `SARGAZO_MC_BASIS_LEVELS`
- `SARGAZO_MC_N_OBSERVATIONS`
- `SARGAZO_MC_TIME_STEPS`
- `SARGAZO_MC_SUMMARY_FILE`
- `SARGAZO_MC_REPLICATE_FILE`

List-valued study variables use comma-separated values. The selected paper-style
preset can be launched directly with:

```bash
make run-monte-carlo-paper
```

The preset uses:

- `SARGAZO_MC_REPLICATES=1000`
- `SARGAZO_MC_BASIS_LEVELS=20`
- `SARGAZO_MC_N_OBSERVATIONS=1000,5000,20000,50000`
- `SARGAZO_MC_TIME_STEPS=1.0e-7,1.0e-6,1.0e-5,1.0e-4`

You can still override any preset component from the command line. For example:

```bash
make run-monte-carlo-paper PAPER_MC_REPLICATES=250
```

## Snapshot comparison

The snapshot-comparison driver reconstructs the physical-space solution from
the spectral coefficients for one stochastic realization and the matching
deterministic trajectory with `sigma = 0`. It writes a long-form CSV with both
fields at the requested output times, then the plotting utility turns that data
into either a comparison figure or a time-evolution video, both with the
stationary velocity field overlaid as vectors.

By default the snapshots are taken at:

- `t = 0`
- `t = 0.5`
- `t = 2`

Run it with:

```bash
make run-snapshot-comparison
make plot-snapshot-comparison
make video-snapshot-comparison
```

The targets use different snapshot-selection rules:

- `make run-snapshot-comparison` and `make plot-snapshot-comparison` use `SARGAZO_SNAPSHOT_TIMES`.
- `make video-snapshot-comparison` regenerates the snapshot CSV using `SARGAZO_SNAPSHOT_FRAME_COUNT`, `SARGAZO_SNAPSHOT_INITIAL_TIME`, and `SARGAZO_SNAPSHOT_FINAL_TIME`.

You can override the snapshot-specific settings with:

- `SARGAZO_SNAPSHOT_TIMES`
- `SARGAZO_SNAPSHOT_FRAME_COUNT`
- `SARGAZO_SNAPSHOT_INITIAL_TIME`
- `SARGAZO_SNAPSHOT_FINAL_TIME`
- `SARGAZO_SNAPSHOT_GRID_NX`
- `SARGAZO_SNAPSHOT_GRID_NY`
- `SARGAZO_SNAPSHOT_COMPARISON_FILE`

`SARGAZO_SNAPSHOT_GRID_NX` and `SARGAZO_SNAPSHOT_GRID_NY` control only the
physical plotting grid used for reconstructed fields. They default to
`SARGAZO_GRID_NX/Y`, but can be larger or smaller than the spectral truncation.

For videos, use a frame count instead of listing every time manually. For example:

```bash
env \
  SARGAZO_SNAPSHOT_TIMES=0,0.5,2.0 \
  SARGAZO_SNAPSHOT_FRAME_COUNT=121 \
  SARGAZO_SNAPSHOT_FINAL_TIME=2.0 \
  make video-snapshot-comparison
```

`SARGAZO_SNAPSHOT_FRAME_COUNT` generates evenly spaced, time-grid-aligned
snapshots between `SARGAZO_SNAPSHOT_INITIAL_TIME` and
`SARGAZO_SNAPSHOT_FINAL_TIME`. This is only activated by
`make video-snapshot-comparison`; it does not replace
`SARGAZO_SNAPSHOT_TIMES` for the static figure target.

A tracked shell environment file is available at
[`data/input/snapshot_comparison.env`](/home/saul/Desktop/2026_SargazoMLDE/mle_estimation_sto_add_dev/data/input/snapshot_comparison.env:1).
Edit the values there and run:

```bash
source data/input/snapshot_comparison.env
make run-snapshot-comparison
```

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
