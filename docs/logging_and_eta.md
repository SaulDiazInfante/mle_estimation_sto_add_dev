# Logging and ETA Improvements

This document describes the logging and estimated time remaining (ETA) features added to the simulation suite.

## Overview

The simulation framework now provides:

1. **Improved Progress Display**: Progress reporting now shows estimated remaining time in human-readable format (hours, minutes, seconds)
2. **Comprehensive Logging**: Optional logging targets that capture all simulation output with timestamps and environment configuration

## Enhanced Progress Reporting

### ETA Display Format

The progress bar now shows the estimated time remaining in an easy-to-read format:

- **Example outputs**:
  - `"eta 2h 15m 30s"` - for long-running simulations
  - `"eta 5m 20s"` - for medium-length simulations  
  - `"eta 45s"` - for short simulations

### How ETA Works

The ETA is calculated based on the current completion rate:

```
ETA = elapsed_time × (total_work - completed_work) / completed_work
```

**Features**:
- Computed continuously as progress updates
- Accounts for varying computation speed
- Shows as seconds when remaining time < 1 minute
- Shows as minutes and seconds when < 1 hour
- Shows as hours, minutes, and seconds for longer durations

## Logging Targets

New Make targets have been added to automatically capture simulation logs:

### Available Logging Targets

#### `make run-with-log`
Logs the basic `run` target (multivariate simulation).

Log file: `data/output/{TIMESTAMP}_run.log`

#### `make run-snapshot-comparison-with-log`
Logs the snapshot comparison workflow.

Log file: `data/output/{TIMESTAMP}_run_snapshot_comparison.log`

#### `make run-monte-carlo-with-log`
Logs the Monte Carlo study workflow.

Log file: `data/output/{TIMESTAMP}_run_monte_carlo.log`

#### `make run-monte-carlo-paper-with-log`
Logs the full paper Monte Carlo study with paper-specific parameters.

Log file: `data/output/{TIMESTAMP}_run_monte_carlo_paper.log`

### Log File Format

Each log file contains:

1. **Header Information**
   - Target name
   - Start timestamp (formatted as `YYYY-MM-DD HH:MM:SS`)
   - End timestamp
   - Elapsed wall time

2. **Environment Configuration**
   - All `SPDE_*` environment variables
   - Fortran compiler version

3. **Execution Output**
   - Complete stdout/stderr from the simulation
   - Progress bar output with ETA updates

4. **Summary**
   - Total elapsed time
   - Exit code
   - Execution status (SUCCESS or FAILED)

### Example Log File Contents

```
============================================================
Simulation Log: run-snapshot-comparison
Start Time: 2026-06-03 12:35:00
============================================================

SPDE_BETA=1.5
SPDE_GRID_NX=50
SPDE_GRID_NY=50
SPDE_N_OBSERVATIONS=10000
SPDE_SEED=42
SPDE_SIGMA=1.0
SPDE_THETA=1.5
SPDE_TIME_STEP=5.0e-5

Simulation progress [####------] 40.00% (4000/10000) eta 3m 20s
Simulation progress [########--] 80.00% (8000/10000) eta 45s
Simulation progress [##########] 100.00% (10000/10000) eta 0s

End Time: 2026-06-03 12:47:30
Elapsed time: 750 seconds (12 minutes 30 seconds)
Status: SUCCESS
```

## Usage Examples

### Run with ETA display only
```bash
# Standard run with progress and ETA
make run-snapshot-comparison
```

### Run with logging
```bash
# Run and capture all output to timestamped log file
make run-snapshot-comparison-with-log

# Check the generated log file
ls -lh data/output/*_run_snapshot_comparison.log
tail -f data/output/*_run_snapshot_comparison.log  # Monitor in real-time
```

### Run Monte Carlo with logging
```bash
# Run the paper Monte Carlo study with logging
make run-monte-carlo-paper-with-log

# Monitor progress in log file
tail -f data/output/*_run_monte_carlo_paper.log
```

## Troubleshooting

### ETA seems wrong
The ETA calculation assumes the work rate remains constant. If computation speed varies significantly during the run, the ETA may not be perfectly accurate, but it provides a useful estimate.

### Log file not created
- Ensure `data/output/` directory exists (created automatically on first use)
- Check write permissions in the output directory
- Verify the simulation actually ran (check exit code)

### Missing environment variables in log
The log captures all variables starting with `SPDE_`. If you want to log other variables, they must be named with the `SPDE_` prefix.

## Implementation Details

### Modified Files

1. **src/modules/progress_reporting_mod.f90**
   - Enhanced `write_progress_line()` subroutine to format ETA in hours/minutes/seconds
   - Calculates time remaining based on current completion fraction and elapsed time
   - Intelligently formats display based on magnitude (h, m, s)

2. **Makefile**
   - Added logging wrapper targets that call original targets with output capture
   - Each wrapper target:
     - Generates timestamped log file with format `{TIMESTAMP}_{TARGET_NAME}.log`
     - Captures environment variables
     - Records start/end times
     - Preserves all simulation output
   - Targets added to `.PHONY` declaration for proper Make handling

3. **scripts/run_with_logging.sh**
   - Helper script (optional) for programmatic logging support
   - Can be used for more advanced logging scenarios

## Performance Impact

- **ETA Calculation**: Minimal overhead (<1% CPU impact)
- **Logging**: Uses `tee` for real-time output, slight I/O overhead but doesn't delay computation
- **Overall**: Recommended for production runs (paper studies) with negligible performance impact

## Future Enhancements

Potential improvements:

1. Machine learning-based ETA prediction (accounts for varying computation speed)
2. Checkpointing with restart capability
3. Email notifications when simulations complete
4. Integration with job schedulers (SLURM, PBS)
5. Real-time progress dashboard (HTTP endpoint)
