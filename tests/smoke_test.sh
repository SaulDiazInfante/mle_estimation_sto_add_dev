#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ARTIFACT_DIR="$ROOT_DIR/tests/artifacts"
TIMESTAMP="20260317T124705"
LEGACY_TIMESTAMP="2026-03-17T12:47:05"
TRAJECTORY_FILE="$ARTIFACT_DIR/${TIMESTAMP}_estimator_trajectory_smoke.csv"
STATE_HISTORY_FILE="$ARTIFACT_DIR/${TIMESTAMP}_state_history_smoke.csv"
PLOT_FILE="$ARTIFACT_DIR/${TIMESTAMP}_estimator_trajectory_smoke.png"
NORMALIZED_DEFAULT_FILE="$ROOT_DIR/data/output/${TIMESTAMP}_estimator_trajectory.csv"

mkdir -p "$ARTIFACT_DIR"
rm -f "$TRAJECTORY_FILE" "$STATE_HISTORY_FILE" "$PLOT_FILE" "$NORMALIZED_DEFAULT_FILE"

make -C "$ROOT_DIR" build

env \
    SARGAZO_N_OBSERVATIONS=25 \
    SARGAZO_MINIMUM_TRAJECTORY_OBSERVATIONS=5 \
    SARGAZO_OUTPUT_TIMESTAMP="$TIMESTAMP" \
    SARGAZO_REQUESTED_TRAJECTORY_POINTS=4 \
    SARGAZO_WRITE_STATE_HISTORY=0 \
    SARGAZO_ESTIMATOR_HISTORY_FILE="$TRAJECTORY_FILE" \
    SARGAZO_STATE_HISTORY_FILE="$STATE_HISTORY_FILE" \
    "$ROOT_DIR/build/bin/multivariate_modular"

test -f "$TRAJECTORY_FILE"
test ! -f "$STATE_HISTORY_FILE"
case "$TRAJECTORY_FILE" in
    *"${TIMESTAMP}"*) ;;
    *) exit 1 ;;
esac
case "$TRAJECTORY_FILE" in
    *:*) exit 1 ;;
esac
grep -q '^n_obs,time,sigma_hat,sigma_true,beta_hat,beta_true,theta_hat,theta_true$' "$TRAJECTORY_FILE"

env \
    SARGAZO_N_OBSERVATIONS=12 \
    SARGAZO_MINIMUM_TRAJECTORY_OBSERVATIONS=4 \
    SARGAZO_OUTPUT_TIMESTAMP="$LEGACY_TIMESTAMP" \
    SARGAZO_REQUESTED_TRAJECTORY_POINTS=3 \
    SARGAZO_WRITE_STATE_HISTORY=0 \
    "$ROOT_DIR/build/bin/multivariate_modular"

test -f "$NORMALIZED_DEFAULT_FILE"
rm -f "$NORMALIZED_DEFAULT_FILE"

"${PYTHON:-python3}" "$ROOT_DIR/visualization/scripts/plot_estimator_trajectory.py" \
    --input "$TRAJECTORY_FILE" \
    --output "$PLOT_FILE"

test -f "$PLOT_FILE"
