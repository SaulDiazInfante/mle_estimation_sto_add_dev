#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ARTIFACT_DIR="$ROOT_DIR/tests/artifacts"
MPL_CONFIG_DIR="$ARTIFACT_DIR/matplotlib-cache"
PLOT_TIMESTAMP="20260515T130000"
VIDEO_TIMESTAMP="20260515T130001"
PLOT_SNAPSHOT_FILE="$ROOT_DIR/data/output/${PLOT_TIMESTAMP}_solution_snapshot_comparison.csv"
VIDEO_SNAPSHOT_FILE="$ROOT_DIR/data/output/${VIDEO_TIMESTAMP}_solution_snapshot_comparison.csv"
PLOT_FILE="$ARTIFACT_DIR/${PLOT_TIMESTAMP}_solution_snapshot_comparison_smoke.png"
VIDEO_FILE="$ARTIFACT_DIR/${VIDEO_TIMESTAMP}_solution_snapshot_comparison_smoke.gif"

mkdir -p "$ARTIFACT_DIR" "$MPL_CONFIG_DIR"
rm -f "$PLOT_SNAPSHOT_FILE" "$VIDEO_SNAPSHOT_FILE" "$PLOT_FILE" "$VIDEO_FILE"

make -C "$ROOT_DIR" build

env \
    SPDE_GRID_NX=6 \
    SPDE_GRID_NY=6 \
    SPDE_N_OBSERVATIONS=201 \
    SPDE_SEED=12345 \
    SPDE_SNAPSHOT_TIMES=0,0.5,2.0 \
    SPDE_TIME_STEP=1.0e-2 \
    make -C "$ROOT_DIR" run-snapshot-comparison \
    TIMESTAMP="$PLOT_TIMESTAMP"

env \
    MPLCONFIGDIR="$MPL_CONFIG_DIR" \
    make -C "$ROOT_DIR" plot-snapshot-comparison \
    SNAPSHOT_PLOT_ARGS="--input $PLOT_SNAPSHOT_FILE --output $PLOT_FILE"

test -f "$PLOT_SNAPSHOT_FILE"
grep -q '^# length_x=' "$PLOT_SNAPSHOT_FILE"
grep -q '^# velocity_mode_x=' "$PLOT_SNAPSHOT_FILE"
grep -q '^solution_kind,time,x_index,y_index,x,y,value$' "$PLOT_SNAPSHOT_FILE"
test "$(wc -l < "$PLOT_SNAPSHOT_FILE")" -eq 221
test -f "$PLOT_FILE"

env \
    SPDE_GRID_NX=6 \
    SPDE_GRID_NY=6 \
    SPDE_N_OBSERVATIONS=201 \
    SPDE_SEED=12345 \
    SPDE_SNAPSHOT_TIMES=0,0.5,2.0 \
    SPDE_SNAPSHOT_FINAL_TIME=2.0 \
    SPDE_TIME_STEP=1.0e-2 \
    make -C "$ROOT_DIR" run-snapshot-comparison \
    TIMESTAMP="$VIDEO_TIMESTAMP" \
    SNAPSHOT_FRAME_COUNT=5

env \
    MPLCONFIGDIR="$MPL_CONFIG_DIR" \
    make -C "$ROOT_DIR" video-snapshot-comparison \
    SNAPSHOT_VIDEO_ARGS="--input $VIDEO_SNAPSHOT_FILE --writer pillow --output $VIDEO_FILE --fps 2"

test -f "$VIDEO_SNAPSHOT_FILE"
test "$(wc -l < "$VIDEO_SNAPSHOT_FILE")" -eq 365
test -f "$VIDEO_FILE"
