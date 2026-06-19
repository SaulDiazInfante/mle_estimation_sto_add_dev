#!/bin/bash
# scripts/sweep_gamma.sh
# Automates sweeping the SPDE_GAMMA parameter, generating manuscript panels,
# and bundling everything into tar archives.

set -e

# Parameter values to sweep
GAMMAS=(0.5 0.6 0.7 0.8 0.9 1.0)
ENV_FILE="data/input/snapshot_comparison.env"
RESULTS_DIR="sweep_results"
SWEEP_LOG="$RESULTS_DIR/sweep.log"

mkdir -p "$RESULTS_DIR"
> "$SWEEP_LOG" # clear/create the log file

# Ensure the .env file exists and has SPDE_GAMMA exported
if ! grep -q "^export SPDE_GAMMA=" "$ENV_FILE"; then
    echo "Error: Could not find 'export SPDE_GAMMA=' in $ENV_FILE"
    exit 1
fi

TOTAL=${#GAMMAS[@]}
CURRENT=0

# Function to draw a tqdm-like progress bar
draw_progress_bar() {
    local current=$1
    local total=$2
    local prefix=$3
    local bar_width=40
    
    local percent=$((current * 100 / total))
    local completed=$((bar_width * current / total))
    local remaining=$((bar_width - completed))
    
    # Build the bar string
    local bar=""
    for ((i=0; i<completed; i++)); do bar="${bar}█"; done
    for ((i=0; i<remaining; i++)); do bar="${bar}-"; done
    
    # \e[K clears the line from the cursor to the end
    printf "\r%s |%s| %d%% [%d/%d]\e[K" "$prefix" "$bar" "$percent" "$current" "$total"
}

echo "Starting parameter sweep. Detailed logs are being saved to $SWEEP_LOG"
draw_progress_bar $CURRENT $TOTAL "Sweeping gamma"

for gamma in "${GAMMAS[@]}"; do
    echo "============================================================" >> "$SWEEP_LOG"
    echo "Running parameter sweep for SPDE_GAMMA = $gamma" >> "$SWEEP_LOG"
    echo "============================================================" >> "$SWEEP_LOG"
    
    # 1. Update the parameter in the environment file
    sed -i "s/^export SPDE_GAMMA=.*/export SPDE_GAMMA=$gamma/" "$ENV_FILE"
    
    # 2. Run simulation to generate the data and the parameter txt backup
    make run-snapshot-comparison-with-log >> "$SWEEP_LOG" 2>&1
    
    # Find the timestamp prefix of the run we just completed
    LATEST_LOG=$(ls data/output/*_run_snapshot_comparison.log 2>/dev/null | sort | tail -n 1)
    TIMESTAMP=$(basename "$LATEST_LOG" | sed "s/_run_snapshot_comparison.log//")
    
    # 3. Generate the manuscript panels (2x2 plots matching the theta_09_gamma_05 format)
    make export-snapshot-manuscript-panels >> "$SWEEP_LOG" 2>&1
    
    # 4. Copy the generated plots into the data/output directory.
    PLOT_SOURCE_DIR="visualization/plots/${TIMESTAMP}_manuscript_panels"
    PLOT_DEST_DIR="data/output/${TIMESTAMP}_manuscript_panels"
    
    if [ -d "$PLOT_SOURCE_DIR" ]; then
        cp -r "$PLOT_SOURCE_DIR" "$PLOT_DEST_DIR" >> "$SWEEP_LOG" 2>&1
    fi
    
    # 5. Create the tarball backup (includes CSVs, param txt, and the plots folder)
    make package-outputs >> "$SWEEP_LOG" 2>&1
    
    # 6. Move the tar file to the results directory and rename it to indicate the gamma value
    TAG=$(grep "Start Time:" "$LATEST_LOG" | sed "s/Start Time: //" | sed "s/-/_/g; s/ /-/" | tr -d "\r")
    TAR_FILE="${TAG}_sto_adv_diff_outputs.tar"
    
    if [ -f "$TAR_FILE" ]; then
        FINAL_TAR="$RESULTS_DIR/gamma_${gamma}_${TAR_FILE}"
        mv "$TAR_FILE" "$FINAL_TAR" >> "$SWEEP_LOG" 2>&1
    fi

    # Update progress bar
    CURRENT=$((CURRENT + 1))
    draw_progress_bar $CURRENT $TOTAL "Sweeping gamma"
done

echo "" # New line after the progress bar finishes
echo "Parameter sweep complete! All tar files are located in the '$RESULTS_DIR/' folder."
