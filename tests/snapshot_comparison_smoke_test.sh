#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ARTIFACT_DIR="$ROOT_DIR/tests/artifacts"
MPL_CONFIG_DIR="$ARTIFACT_DIR/matplotlib-cache"
PLOT_TIMESTAMP="20260515T130000"
VIDEO_TIMESTAMP="20260515T130001"
PLOT_DET_FILE="$ROOT_DIR/data/output/${PLOT_TIMESTAMP}_deterministic_path.csv"
PLOT_STO_FILE="$ROOT_DIR/data/output/${PLOT_TIMESTAMP}_stochastic_path.csv"
VIDEO_DET_FILE="$ROOT_DIR/data/output/${VIDEO_TIMESTAMP}_deterministic_path.csv"
VIDEO_STO_FILE="$ROOT_DIR/data/output/${VIDEO_TIMESTAMP}_stochastic_path.csv"
PLOT_FILE="$ARTIFACT_DIR/${PLOT_TIMESTAMP}_solution_snapshot_comparison_smoke.png"
VIDEO_FILE="$ARTIFACT_DIR/${VIDEO_TIMESTAMP}_solution_snapshot_comparison_smoke.gif"
PARAVIEW_DIR="$ARTIFACT_DIR/${VIDEO_TIMESTAMP}_paraview_smoke"

mkdir -p "$ARTIFACT_DIR" "$MPL_CONFIG_DIR"
rm -f "$PLOT_DET_FILE" "$PLOT_STO_FILE" "$VIDEO_DET_FILE" "$VIDEO_STO_FILE" "$PLOT_FILE" "$VIDEO_FILE"
rm -rf "$PARAVIEW_DIR"

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
    SNAPSHOT_PLOT_ARGS="--input $PLOT_DET_FILE --output $PLOT_FILE"

test -f "$PLOT_DET_FILE"
test -f "$PLOT_STO_FILE"
grep -q '^# length_x=' "$PLOT_DET_FILE"
grep -q '^# velocity_mode_x=' "$PLOT_DET_FILE"
grep -q '^time,x_index,y_index,x,y,value$' "$PLOT_DET_FILE"
grep -q '^time,x_index,y_index,x,y,value$' "$PLOT_STO_FILE"
test "$(wc -l < "$PLOT_DET_FILE")" -eq 113
test "$(wc -l < "$PLOT_STO_FILE")" -eq 113
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
    SNAPSHOT_VIDEO_ARGS="--input $VIDEO_DET_FILE --writer pillow --output $VIDEO_FILE --fps 2"

test -f "$VIDEO_DET_FILE"
test -f "$VIDEO_STO_FILE"
test "$(wc -l < "$VIDEO_DET_FILE")" -eq 185
test "$(wc -l < "$VIDEO_STO_FILE")" -eq 185
test -f "$VIDEO_FILE"

make -C "$ROOT_DIR" export-snapshot-paraview \
    TIMESTAMP="$VIDEO_TIMESTAMP" \
    SNAPSHOT_PARAVIEW_INPUT="$VIDEO_DET_FILE" \
    SNAPSHOT_PARAVIEW_DIR="$PARAVIEW_DIR"

test -f "$PARAVIEW_DIR/solution_snapshots.pvd"
test "$(find "$PARAVIEW_DIR/frames" -maxdepth 1 -type f -name '*.vts' | wc -l)" -eq 5
grep -q 'Name="deterministic"' "$PARAVIEW_DIR/frames/solution_00000.vts"
grep -q 'Name="stochastic"' "$PARAVIEW_DIR/frames/solution_00000.vts"
