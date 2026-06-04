# Stochastic Convection-Diffusion MLE

This repository contains the maintained Fortran implementation for simulating
the spectral Galerkin system associated with a stochastic convection-diffusion
equation, estimating model parameters, and visualizing estimator trajectories.

## Repository layout

- `src/app`: executable entry points.
- `src/modules`: maintained Fortran modules, one module per file.
- `legacy`: preserved legacy Fortran sources kept out of the main build.
- `data/input`: reserved location for future input datasets and configuration files.
- `data/output`: generated CSV and text outputs from simulation and estimation runs.
- `visualization/scripts`: plotting utilities.
- `visualization/plots`: generated figures.
- `visualization/paraview`: generated ParaView VTK exports.
- `docs`: project documentation and engineering notes.
- `tests`: smoke tests and future regression tests.
- `.github/workflows`: CI definitions.

## Common commands

- `make build`: compile the maintained application into `build/bin/`.
- `make run`: run the default workflow and write timestamped outputs under `data/output/`.
- `make run-monte-carlo`: run the Monte Carlo study driver and write summary plus replicate CSV files under `data/output/`.
- `make run-monte-carlo-paper`: run the paper-style study preset with the selected sweep over `dt`, basis size, and observation count.
- `make run-snapshot-comparison`: reconstruct one stochastic path and its deterministic counterpart at selected times.
- `make run-with-log`: same as `make run` but captures all output to a timestamped log file.
- `make run-snapshot-comparison-with-log`: same as `make run-snapshot-comparison` with logging.
- `make run-monte-carlo-with-log`: same as `make run-monte-carlo` with logging.
- `make run-monte-carlo-paper-with-log`: same as `make run-monte-carlo-paper` with logging (recommended for paper studies).
- `make plot`: plot the latest generated estimator CSV already present in `data/output/`.
- `make plot-snapshot-comparison`: plot the latest reconstructed snapshot comparison CSV.
- `make video-snapshot-comparison`: animate the latest reconstructed snapshot comparison CSV.
- `make export-snapshot-paraview`: export the latest reconstructed snapshot CSV as ParaView `.pvd` and `.vts` files.
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

- `SPDE_N_OBSERVATIONS`
- `SPDE_MINIMUM_TRAJECTORY_OBSERVATIONS`
- `SPDE_REQUESTED_TRAJECTORY_POINTS`
- `SPDE_TIME_STEP`
- `SPDE_GRID_NX`
- `SPDE_GRID_NY`
- `SPDE_VELOCITY_MODE_X`
- `SPDE_VELOCITY_MODE_Y`
- `SPDE_LENGTH_X`
- `SPDE_LENGTH_Y`
- `SPDE_GAMMA`
- `SPDE_BETA`
- `SPDE_THETA`
- `SPDE_SIGMA`
- `SPDE_SEED`
- `SPDE_OUTPUT_TIMESTAMP`
- `SPDE_WRITE_STATE_HISTORY`
- `SPDE_STATE_HISTORY_FILE`
- `SPDE_ESTIMATOR_HISTORY_FILE`
- `SPDE_SNAPSHOT_TIMES`
- `SPDE_SNAPSHOT_FRAME_COUNT`
- `SPDE_SNAPSHOT_INITIAL_TIME`
- `SPDE_SNAPSHOT_FINAL_TIME`
- `SPDE_SNAPSHOT_GRID_NX`
- `SPDE_SNAPSHOT_GRID_NY`
- `SPDE_SNAPSHOT_COMPARISON_FILE`

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

- `SPDE_MC_REPLICATES`
- `SPDE_MC_BASIS_LEVELS`
- `SPDE_MC_N_OBSERVATIONS`
- `SPDE_MC_TIME_STEPS`
- `SPDE_MC_SUMMARY_FILE`
- `SPDE_MC_REPLICATE_FILE`

List-valued study variables use comma-separated values. The selected paper-style
preset can be launched directly with:

```bash
make run-monte-carlo-paper
```

The preset uses:

- `SPDE_MC_REPLICATES=1000`
- `SPDE_MC_BASIS_LEVELS=20`
- `SPDE_MC_N_OBSERVATIONS=1000,5000,20000,50000`
- `SPDE_MC_TIME_STEPS=1.0e-7,1.0e-6,1.0e-5,1.0e-4`

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
make export-snapshot-paraview
```

The targets use different snapshot-selection rules:

- `make run-snapshot-comparison` and `make plot-snapshot-comparison` use `SPDE_SNAPSHOT_TIMES`.
- `make video-snapshot-comparison` regenerates the snapshot CSV using `SPDE_SNAPSHOT_FRAME_COUNT`, `SPDE_SNAPSHOT_INITIAL_TIME`, and `SPDE_SNAPSHOT_FINAL_TIME`.

You can override the snapshot-specific settings with:

- `SPDE_SNAPSHOT_TIMES`
- `SPDE_SNAPSHOT_FRAME_COUNT`
- `SPDE_SNAPSHOT_INITIAL_TIME`
- `SPDE_SNAPSHOT_FINAL_TIME`
- `SPDE_SNAPSHOT_GRID_NX`
- `SPDE_SNAPSHOT_GRID_NY`
- `SPDE_SNAPSHOT_COMPARISON_FILE`

`SPDE_SNAPSHOT_GRID_NX` and `SPDE_SNAPSHOT_GRID_NY` control only the
physical plotting grid used for reconstructed fields. They default to
`SPDE_GRID_NX/Y`, but can be larger or smaller than the spectral truncation.

For videos, use a frame count instead of listing every time manually. For example:

```bash
env \
  SPDE_SNAPSHOT_TIMES=0,0.5,2.0 \
  SPDE_SNAPSHOT_FRAME_COUNT=121 \
  SPDE_SNAPSHOT_FINAL_TIME=2.0 \
  make video-snapshot-comparison
```

`SPDE_SNAPSHOT_FRAME_COUNT` generates evenly spaced, time-grid-aligned
snapshots between `SPDE_SNAPSHOT_INITIAL_TIME` and
`SPDE_SNAPSHOT_FINAL_TIME`. This is only activated by
`make video-snapshot-comparison`; it does not replace
`SPDE_SNAPSHOT_TIMES` for the static figure target.

A tracked shell environment file is available at
[`data/input/snapshot_comparison.env`](data/input/snapshot_comparison.env).
`make run-snapshot-comparison` now automatically sources this file when it is present.
You can still override any runtime setting manually before running the target.

To export the same animation data for ParaView, run:

```bash
make export-snapshot-paraview \
  TIMESTAMP=20260317T124705 \
  SNAPSHOT_PARAVIEW_INPUT=data/output/spde_animation_solution_snapshot_comparison.csv
```

This creates a timestamped directory under `visualization/paraview/`, for example:

```text
visualization/paraview/20260317T124705_snapshot_paraview/
```

Open `solution_snapshots.pvd` in ParaView. The `frames/` subdirectory contains
one `.vts` structured-grid file per snapshot time, with point-data arrays named
`deterministic` and `stochastic`.

### ParaView figure and animation recommendations

For manuscript-quality figures, export a high-resolution PNG and use the
`stochastic` surface preset when the goal is to highlight uncertainty patterns.
Use a size such as `1920 1080` or larger for Elsevier/CNSNS-quality figures.

For side-by-side comparison of deterministic and stochastic fields, use the
`side-by-side-comparison` preset and export distinct outputs for each field.
Add `--frame-label` when creating video output to overlay frame numbering like
`25/50,000` in the top-right corner.

To export static 3-D initial/final comparisons, run:

```bash
make paraview-initial-final-images PVPYTHON=/opt/paraview/bin/pvpython
```

This writes `solution_initial_final_deterministic.png` and
`solution_initial_final_stochastic.png` under `visualization/paraview/`.
Each image contains the initial and final simulation states with a warped
surface and edge overlay for a clearer 3-D reading.

If the stochastic panel looks over-warped, reduce the vertical scale with
`--stochastic-warp-scale 0.8` or use `--warp-scale 1.0` to keep the surface
visible without producing a wall-like distortion.

Example commands:

```bash
pvpython visualization/scripts/paraview_solution_snapshot_visualization.py \
  --preset stochastic-surface \
  --output-image visualization/paraview/solution_snapshots_stochastic.png \
  --output-video visualization/paraview/solution_snapshots_stochastic.mp4 \
  --fps 8 \
  --size 1920 1080 \
  --frame-label \
  --show-legend
```

```bash
pvpython visualization/scripts/paraview_solution_snapshot_visualization.py \
  --preset side-by-side-comparison \
  --output-image visualization/paraview/solution_snapshots.png \
  --output-video visualization/paraview/solution_snapshots.mp4 \
  --frame-label
```

You can also use the new Plotly exporter for an interactive HTML version of the
same comparison. It uses shared color and value-axis ranges across all panels,
a golden 3-D camera/aspect button, and a tracker marker at the point of largest
final stochastic-minus-deterministic difference:

```bash
python visualization/scripts/plotly_solution_snapshot_comparison.py \
  --output visualization/plots/solution_snapshot_comparison_plotly.html
```

## Progress monitoring and logging

### Estimated remaining time (ETA)

All long-running simulations display a progress bar with an estimated remaining time. The ETA format automatically adapts to the remaining duration:

```
Simulation progress [########--] 80.00% (8000/10000) eta 2m 15s
Simulation progress [##########] 100.00% (10000/10000) eta 0s
```

The ETA is calculated based on the current completion rate and updates continuously. It provides a useful estimate even though actual runtime may vary.

### Comprehensive simulation logging

For production runs (especially `make run-monte-carlo-paper`), use the logging variants to capture all output with timestamps and environment configuration:

```bash
# Run with logging enabled (output captured to timestamped log file)
make run-snapshot-comparison-with-log
make run-monte-carlo-with-log
make run-monte-carlo-paper-with-log
```

Log files are written to `data/output/{TIMESTAMP}_{TARGET_NAME}.log` and include:

- **Simulation metadata**: start/end times, elapsed wall time, target name
- **Environment configuration**: all `SPDE_*` variables used by the run
- **Fortran compiler version**
- **Complete simulation output**: progress bars, console messages, any warnings or errors
- **Summary**: exit code and final status (SUCCESS or FAILED)

Example log file location:

```
data/output/20260603T143015_run_monte_carlo_paper.log
```

For details, see [docs/logging_and_eta.md](docs/logging_and_eta.md).

## Git tracking policy

- Keep source code, docs, CI config, and small deterministic tests in git.
- Keep generated outputs in `data/output/`, `visualization/plots/`,
  `visualization/paraview/`, and `tests/artifacts/`, which are ignored.
- Keep local or raw input data in the ignored paths under `data/input/`.
- The repository includes a size guard that rejects staged files larger than `50 MB` by default.

Run `make setup-git-hooks` once per clone to activate the local pre-commit size check.

## GitHub workflow

- CI runs on pull requests and pushes to `main`.
- Doxygen documentation is built on pull requests and published to GitHub Pages from `main`.
- Tagged releases like `v0.1.0` publish a Linux build artifact through GitHub Releases.
- Contribution templates and issue forms are included under `.github/`.

See [docs/github_publish.md](docs/github_publish.md) for the exact publish steps and [docs/testing_strategy.md](docs/testing_strategy.md) for the TDD-oriented test layout.
