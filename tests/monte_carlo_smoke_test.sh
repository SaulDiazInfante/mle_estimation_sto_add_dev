#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ARTIFACT_DIR="$ROOT_DIR/tests/artifacts"
TIMESTAMP="20260515T120000"
SUMMARY_FILE="$ARTIFACT_DIR/${TIMESTAMP}_monte_carlo_summary_smoke.csv"
REPLICATE_FILE="$ARTIFACT_DIR/${TIMESTAMP}_monte_carlo_replicates_smoke.csv"
PAPER_TIMESTAMP="20260515T120001"
PAPER_SUMMARY_FILE="$ARTIFACT_DIR/${PAPER_TIMESTAMP}_monte_carlo_summary_paper_smoke.csv"
PAPER_REPLICATE_FILE="$ARTIFACT_DIR/${PAPER_TIMESTAMP}_monte_carlo_replicates_paper_smoke.csv"

mkdir -p "$ARTIFACT_DIR"
rm -f "$SUMMARY_FILE" "$REPLICATE_FILE" "$PAPER_SUMMARY_FILE" "$PAPER_REPLICATE_FILE"

make -C "$ROOT_DIR" build

env \
    SARGAZO_OUTPUT_TIMESTAMP="$TIMESTAMP" \
    SARGAZO_SEED=31415 \
    SARGAZO_BETA=0.1 \
    SARGAZO_THETA=0.1 \
    SARGAZO_SIGMA=1.0 \
    SARGAZO_MC_REPLICATES=2 \
    SARGAZO_MC_BASIS_LEVELS=4,6 \
    SARGAZO_MC_N_OBSERVATIONS=20,25 \
    SARGAZO_MC_TIME_STEPS=1.0e-4,5.0e-5 \
    SARGAZO_MC_SUMMARY_FILE="$SUMMARY_FILE" \
    SARGAZO_MC_REPLICATE_FILE="$REPLICATE_FILE" \
    "$ROOT_DIR/build/bin/monte_carlo_study"

test -f "$SUMMARY_FILE"
test -f "$REPLICATE_FILE"

grep -q '^basis_level,nx,ny,state_dimension,n_obs,time_step,n_replicates,sigma_mean,sigma_sd,sigma_true,beta_mean,beta_sd,beta_true,theta_mean,theta_sd,theta_true,average_setup_time,average_estimation_time,average_total_time$' "$SUMMARY_FILE"
grep -q '^basis_level,nx,ny,state_dimension,n_obs,time_step,replicate,sigma_hat,sigma_true,beta_hat,beta_true,theta_hat,theta_true,setup_time,estimation_time,total_time$' "$REPLICATE_FILE"

test "$(wc -l < "$SUMMARY_FILE")" -eq 9
test "$(wc -l < "$REPLICATE_FILE")" -eq 17

env \
    SARGAZO_SEED=27182 \
    SARGAZO_BETA=0.1 \
    SARGAZO_THETA=0.1 \
    SARGAZO_SIGMA=1.0 \
    SARGAZO_MC_SUMMARY_FILE="$PAPER_SUMMARY_FILE" \
    SARGAZO_MC_REPLICATE_FILE="$PAPER_REPLICATE_FILE" \
    make -C "$ROOT_DIR" run-monte-carlo-paper \
    TIMESTAMP="$PAPER_TIMESTAMP" \
    PAPER_MC_REPLICATES=1 \
    PAPER_MC_BASIS_LEVELS=4 \
    PAPER_MC_N_OBSERVATIONS=12 \
    PAPER_MC_TIME_STEPS=1.0e-4

test -f "$PAPER_SUMMARY_FILE"
test -f "$PAPER_REPLICATE_FILE"
