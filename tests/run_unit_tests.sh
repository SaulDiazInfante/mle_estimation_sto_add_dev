#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$ROOT_DIR/build"
TEST_BIN_DIR="$BUILD_DIR/tests"
TEST_BIN="$TEST_BIN_DIR/test_core"

mkdir -p "$TEST_BIN_DIR"

make -C "$ROOT_DIR" build

gfortran \
    -std=f2008 \
    -Wall -Wextra -Wimplicit-interface \
    -fopenmp \
    -I "$BUILD_DIR/mod" \
    "$ROOT_DIR/tests/unit/test_core.f90" \
    "$BUILD_DIR/obj/model_types_mod.o" \
    "$BUILD_DIR/obj/parameter_ml_estimation_mod.o" \
    -o "$TEST_BIN"

"$TEST_BIN"
