# Sargazo ML Estimation Documentation

## Overview

This documentation site is generated with Doxygen from the maintained Fortran
source tree and selected project documentation.

The codebase centers on four responsibilities:

- defining shared grid, operator, and estimator data types
- assembling the spectral operators for the SDE model
- simulating state trajectories and estimating model parameters
- exporting CSV outputs and plotting estimator trajectories

## Where to start

- Browse the file list for the high-level source layout.
- Start with `src/app/main.f90` for the end-to-end execution flow.
- Read the modules under `src/modules/` for the maintained API.
- Use the `legacy/` files only as historical references; they are not part of the default build.

## Local build

Generate the HTML documentation locally with:

```bash
make docs
```

The generated site is written to `build/docs/doxygen/html`.

## GitHub integration

- Pull requests build the documentation as part of the `Docs` workflow.
- Pushes to `main` publish the generated Doxygen site to GitHub Pages.
- See `docs/github_publish.md` for the repository and Pages setup steps.
