#!/bin/bash
# run_with_logging.sh
# Wrapper script to execute a simulation target with comprehensive logging.
# Usage: run_with_logging.sh <target_name> <command_to_run> [args...]

set -u

TARGET_NAME="${1:-}"
shift || true

if [ -z "$TARGET_NAME" ]; then
    echo "Error: target name required" >&2
    exit 1
fi

TIMESTAMP=$(date '+%Y%m%dT%H%M%S')
COMMAND_NAME=$(echo "$TARGET_NAME" | sed 's/-/_/g')
LOG_FILE="data/output/${TIMESTAMP}_${COMMAND_NAME}.log"

# Ensure output directory exists
mkdir -p data/output

# Function to log with timestamp
log_msg() {
    local msg="$1"
    local ts
    ts=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${ts}] ${msg}" | tee -a "$LOG_FILE"
}

# Start logging
{
    echo "============================================================"
    echo "Simulation Log: ${TARGET_NAME}"
    echo "============================================================"
    echo ""
    log_msg "Target: ${TARGET_NAME}"
    log_msg "Log file: ${LOG_FILE}"
    log_msg "Timestamp: ${TIMESTAMP}"
    echo ""
    
    echo "--- Environment Configuration ---"
    log_msg "Working directory: $(pwd)"
    log_msg "Fortran compiler: $(which gfortran 2>/dev/null || echo 'not found')"
    log_msg "Python: $(which python3 2>/dev/null || echo 'not found')"
    log_msg "ParaView: $(which pvpython 2>/dev/null || echo 'not found')"
    echo ""
    
    echo "--- Simulation Parameters ---"
    if [ -n "${SPDE_N_OBSERVATIONS:-}" ]; then
        log_msg "SPDE_N_OBSERVATIONS: ${SPDE_N_OBSERVATIONS}"
    fi
    if [ -n "${SPDE_TIME_STEP:-}" ]; then
        log_msg "SPDE_TIME_STEP: ${SPDE_TIME_STEP}"
    fi
    if [ -n "${SPDE_GRID_NX:-}" ]; then
        log_msg "SPDE_GRID_NX: ${SPDE_GRID_NX}"
    fi
    if [ -n "${SPDE_GRID_NY:-}" ]; then
        log_msg "SPDE_GRID_NY: ${SPDE_GRID_NY}"
    fi
    if [ -n "${SPDE_BETA:-}" ]; then
        log_msg "SPDE_BETA: ${SPDE_BETA}"
    fi
    if [ -n "${SPDE_THETA:-}" ]; then
        log_msg "SPDE_THETA: ${SPDE_THETA}"
    fi
    if [ -n "${SPDE_SIGMA:-}" ]; then
        log_msg "SPDE_SIGMA: ${SPDE_SIGMA}"
    fi
    if [ -n "${SPDE_SEED:-}" ]; then
        log_msg "SPDE_SEED: ${SPDE_SEED}"
    fi
    if [ -n "${SPDE_MC_REPLICATES:-}" ]; then
        log_msg "SPDE_MC_REPLICATES: ${SPDE_MC_REPLICATES}"
    fi
    if [ -n "${SPDE_MC_BASIS_LEVELS:-}" ]; then
        log_msg "SPDE_MC_BASIS_LEVELS: ${SPDE_MC_BASIS_LEVELS}"
    fi
    if [ -n "${SPDE_MC_N_OBSERVATIONS:-}" ]; then
        log_msg "SPDE_MC_N_OBSERVATIONS: ${SPDE_MC_N_OBSERVATIONS}"
    fi
    if [ -n "${SPDE_MC_TIME_STEPS:-}" ]; then
        log_msg "SPDE_MC_TIME_STEPS: ${SPDE_MC_TIME_STEPS}"
    fi
    echo ""
    
    echo "--- Full Environment (all SPDE_* variables) ---"
    env | grep '^SPDE_' | sort | tee -a "$LOG_FILE"
    echo ""
    
    echo "--- Execution ---"
    START_TIME=$(date '+%Y-%m-%d %H:%M:%S')
    START_EPOCH=$(date '+%s')
    log_msg "Start time: ${START_TIME}"
    log_msg "Command: $*"
    echo ""
    
    # Execute the command with all output logged
    if eval "$@" >> "$LOG_FILE" 2>&1; then
        EXIT_CODE=0
        STATUS="SUCCESS"
    else
        EXIT_CODE=$?
        STATUS="FAILED"
    fi
    
    END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
    END_EPOCH=$(date '+%s')
    ELAPSED=$((END_EPOCH - START_EPOCH))
    
    echo ""
    echo "--- Execution Summary ---"
    log_msg "End time: ${END_TIME}"
    HOURS=$((ELAPSED / 3600))
    MINUTES=$(((ELAPSED % 3600) / 60))
    SECONDS=$((ELAPSED % 60))
    if [ $HOURS -gt 0 ]; then
        log_msg "Elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED} seconds)"
    elif [ $MINUTES -gt 0 ]; then
        log_msg "Elapsed time: ${MINUTES}m ${SECONDS}s (${ELAPSED} seconds)"
    else
        log_msg "Elapsed time: ${SECONDS}s"
    fi
    log_msg "Exit code: ${EXIT_CODE}"
    log_msg "Status: ${STATUS}"
    
    if [ -d data/output ]; then
        echo ""
        echo "--- Generated output files ---"
        log_msg "Output files in data/output/:"
        find data/output -maxdepth 1 -type f -newermt "$(date -d @$START_EPOCH '+%Y-%m-%d %H:%M:%S')" 2>/dev/null | \
            while read -r file; do
                size=$(du -h "$file" 2>/dev/null | cut -f1)
                log_msg "  $(basename "$file") (${size})"
            done
    fi
    
    echo ""
    echo "============================================================"
    log_msg "Log complete"
    echo "============================================================"
    
    exit $EXIT_CODE
} | tee "$LOG_FILE"

