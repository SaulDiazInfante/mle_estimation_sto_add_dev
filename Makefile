FC = gfortran
FFLAGS ?= -std=f2008 -O3 -Wall -Wextra -Wimplicit-interface
OMPFLAGS ?= -fopenmp
PYTHON ?= python3
DOXYGEN ?= doxygen
TIMESTAMP ?= $(shell date '+%Y%m%dT%H%M%S')
LATEST_ESTIMATOR_PATTERN := *_estimator_trajectory.csv
LATEST_SNAPSHOT_PATTERN := *_deterministic_path.csv
# Find latest snapshot data (matches both old and new naming schemes)
# Old: *_solution_snapshot_comparison.csv
# New: *_deterministic_path.csv or *_stochastic_path.csv
FIND_LATEST_SNAPSHOT_CMD = find $(DATA_OUTPUT_DIR) -maxdepth 1 -type f \( -name '*_deterministic_path.csv' -o -name '*_stochastic_path.csv' -o -name '*_solution_snapshot_comparison.csv' \) | sort | tail -n 1
FIND_LATEST_LOG_CMD = ls $(DATA_OUTPUT_DIR)/*_run*.log 2>/dev/null | sort | tail -n 1
PAPER_MC_REPLICATES ?= 1000
PAPER_MC_BASIS_LEVELS ?= 20
PAPER_MC_N_OBSERVATIONS ?= 100,500
PAPER_MC_TIME_STEPS ?= 1.0e-5,1.0e-4
#PAPER_MC_N_OBSERVATIONS ?= 1000,5000,20000,50000
#PAPER_MC_TIME_STEPS ?= 1.0e-7,1.0e-6,1.0e-5,1.0e-4

SRC_DIR := src
APP_DIR := $(SRC_DIR)/app
MODULE_DIR := $(SRC_DIR)/modules
BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
MOD_DIR := $(BUILD_DIR)/mod
BIN_DIR := $(BUILD_DIR)/bin
RUN_TARGET := $(BIN_DIR)/multivariate_modular
MONTE_CARLO_TARGET := $(BIN_DIR)/monte_carlo_study
SNAPSHOT_COMPARISON_TARGET := $(BIN_DIR)/snapshot_comparison

DATA_OUTPUT_DIR := data/output
PLOT_DIR := visualization/plots
PARAVIEW_DIR := visualization/paraview
PVPYTHON ?= pvpython
PLOT_SCRIPT := visualization/scripts/plot_estimator_trajectory.py
SNAPSHOT_PLOT_SCRIPT := visualization/scripts/plot_solution_snapshot_comparison.py
SNAPSHOT_VIDEO_SCRIPT := visualization/scripts/animate_solution_snapshot_comparison.py
SNAPSHOT_3D_SCRIPT := visualization/scripts/plot_3d_rugosity_comparison.py
SNAPSHOT_MANUSCRIPT_SCRIPT := visualization/scripts/export_snapshot_manuscript_panels.py
SNAPSHOT_PARAVIEW_SCRIPT := visualization/scripts/export_snapshot_paraview.py
PARAVIEW_VIS_SCRIPT := visualization/scripts/paraview_solution_snapshot_visualization.py
PARAVIEW_IMAGE ?= $(PARAVIEW_DIR)/solution_snapshots.png
PARAVIEW_VIDEO ?= $(PARAVIEW_DIR)/solution_snapshots.mp4
PARAVIEW_PRESET ?= side-by-side-comparison
SNAPSHOT_PLOT_ARGS ?=
SNAPSHOT_VIDEO_ARGS ?=
SNAPSHOT_3D_ARGS ?=
SNAPSHOT_MANUSCRIPT_ARGS ?=
SNAPSHOT_PARAVIEW_ARGS ?=
SNAPSHOT_MANUSCRIPT_INPUT ?=
SNAPSHOT_PARAVIEW_INPUT ?=
SNAPSHOT_PARAVIEW_DIR ?= $(PARAVIEW_DIR)/$(TIMESTAMP)_snapshot_paraview
SNAPSHOT_FRAME_COUNT ?=
DOXYFILE := docs/Doxyfile
DOCS_DIR := $(BUILD_DIR)/docs/doxygen
DOCS_HTML_DIR := $(DOCS_DIR)/html

MODULE_SRC := \
	$(MODULE_DIR)/model_types_mod.f90 \
	$(MODULE_DIR)/driver_support_mod.f90 \
	$(MODULE_DIR)/validation_mod.f90 \
	$(MODULE_DIR)/csv_output_mod.f90 \
	$(MODULE_DIR)/spectral_operators_mod.f90 \
	$(MODULE_DIR)/solution_reconstruction_mod.f90 \
	$(MODULE_DIR)/progress_reporting_mod.f90 \
	$(MODULE_DIR)/sde_simulation_mod.f90 \
	$(MODULE_DIR)/parameter_ml_estimation_mod.f90 \
	$(MODULE_DIR)/workflow_mod.f90 \
	$(MODULE_DIR)/monte_carlo_study_mod.f90
COMMON_MODULE_OBJ := \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/driver_support_mod.o \
	$(OBJ_DIR)/validation_mod.o \
	$(OBJ_DIR)/csv_output_mod.o \
	$(OBJ_DIR)/spectral_operators_mod.o \
	$(OBJ_DIR)/solution_reconstruction_mod.o \
	$(OBJ_DIR)/progress_reporting_mod.o \
	$(OBJ_DIR)/sde_simulation_mod.o \
	$(OBJ_DIR)/parameter_ml_estimation_mod.o \
	$(OBJ_DIR)/workflow_mod.o \
	$(OBJ_DIR)/monte_carlo_study_mod.o
APP_OBJ := \
	$(OBJ_DIR)/main.o
APP_OBJ += $(OBJ_DIR)/monte_carlo_study.o
APP_OBJ += $(OBJ_DIR)/snapshot_comparison.o
OBJ := $(COMMON_MODULE_OBJ) $(APP_OBJ)

.PHONY: all build run run-monte-carlo run-monte-carlo-paper run-snapshot-comparison run-with-log run-snapshot-comparison-with-log run-monte-carlo-with-log run-monte-carlo-paper-with-log plot plot-snapshot-comparison video-snapshot-comparison plot-3d-rugosity export-snapshot-manuscript-panels export-snapshot-paraview paraview-figure-video docs test test-smoke test-monte-carlo test-snapshot-comparison test-unit check-large-files setup-git-hooks clean distclean package-outputs

all: build

build: $(RUN_TARGET) $(MONTE_CARLO_TARGET) $(SNAPSHOT_COMPARISON_TARGET)

run: build | $(DATA_OUTPUT_DIR)
	SPDE_OUTPUT_TIMESTAMP='$(TIMESTAMP)' $(RUN_TARGET)

run-monte-carlo: build | $(DATA_OUTPUT_DIR)
	SPDE_OUTPUT_TIMESTAMP='$(TIMESTAMP)' $(MONTE_CARLO_TARGET)

run-monte-carlo-paper: build | $(DATA_OUTPUT_DIR)
	SPDE_OUTPUT_TIMESTAMP='$(TIMESTAMP)' \
	SPDE_MC_REPLICATES='$(PAPER_MC_REPLICATES)' \
	SPDE_MC_BASIS_LEVELS='$(PAPER_MC_BASIS_LEVELS)' \
	SPDE_MC_N_OBSERVATIONS='$(PAPER_MC_N_OBSERVATIONS)' \
	SPDE_MC_TIME_STEPS='$(PAPER_MC_TIME_STEPS)' \
	$(MONTE_CARLO_TARGET)

run-snapshot-comparison: build | $(DATA_OUTPUT_DIR)
	@if [ -f data/input/snapshot_comparison.env ]; then . data/input/snapshot_comparison.env; fi; \
	SPDE_OUTPUT_TIMESTAMP='$(TIMESTAMP)' \
	$(if $(SNAPSHOT_FRAME_COUNT),SPDE_SNAPSHOT_USE_FRAME_COUNT=1 SPDE_SNAPSHOT_FRAME_COUNT='$(SNAPSHOT_FRAME_COUNT)') \
	$(SNAPSHOT_COMPARISON_TARGET)

# Logging wrapper targets
# These targets execute the simulation targets with comprehensive logging
run-with-log: build | $(DATA_OUTPUT_DIR)
	@target_log="data/output/$(TIMESTAMP)_run.log"; \
	echo "Logging simulation to $$target_log"; \
	echo "============================================================" > "$$target_log"; \
	echo "Simulation Log: run" >> "$$target_log"; \
	echo "Target: run" >> "$$target_log"; \
	echo "Start Time: $$(date '+%Y-%m-%d %H:%M:%S')" >> "$$target_log"; \
	echo "============================================================" >> "$$target_log"; \
	echo "" >> "$$target_log"; \
	env | grep '^SPDE_' | sort >> "$$target_log"; \
	echo "" >> "$$target_log"; \
	($(MAKE) run 2>&1 | tee -a "$$target_log") && \
	echo "End Time: $$(date '+%Y-%m-%d %H:%M:%S')" >> "$$target_log" && \
	echo "Status: SUCCESS" >> "$$target_log" || \
	echo "Status: FAILED" >> "$$target_log"

run-snapshot-comparison-with-log: build | $(DATA_OUTPUT_DIR)
	@target_log="data/output/$(TIMESTAMP)_run_snapshot_comparison.log"; \
	target_params="data/output/$(TIMESTAMP)_snapshot_comparison_params.txt"; \
	echo "Logging simulation to $$target_log"; \
	if [ -f data/input/snapshot_comparison.env ]; then \
		cp data/input/snapshot_comparison.env "$$target_params"; \
		echo "Saved environment parameters to $$target_params"; \
	fi; \
	echo "============================================================" > "$$target_log"; \
	echo "Simulation Log: run-snapshot-comparison" >> "$$target_log"; \
	echo "Start Time: $$(date '+%Y-%m-%d %H:%M:%S')" >> "$$target_log"; \
	echo "============================================================" >> "$$target_log"; \
	echo "" >> "$$target_log"; \
	env | grep '^SPDE_' | sort >> "$$target_log"; \
	echo "" >> "$$target_log"; \
	($(MAKE) run-snapshot-comparison 2>&1 | tee -a "$$target_log") && \
	echo "End Time: $$(date '+%Y-%m-%d %H:%M:%S')" >> "$$target_log" && \
	echo "Status: SUCCESS" >> "$$target_log" || \
	echo "Status: FAILED" >> "$$target_log"

run-monte-carlo-with-log: build | $(DATA_OUTPUT_DIR)
	@target_log="data/output/$(TIMESTAMP)_run_monte_carlo.log"; \
	echo "Logging simulation to $$target_log"; \
	echo "============================================================" > "$$target_log"; \
	echo "Simulation Log: run-monte-carlo" >> "$$target_log"; \
	echo "Start Time: $$(date '+%Y-%m-%d %H:%M:%S')" >> "$$target_log"; \
	echo "============================================================" >> "$$target_log"; \
	echo "" >> "$$target_log"; \
	env | grep '^SPDE_' | sort >> "$$target_log"; \
	echo "" >> "$$target_log"; \
	($(MAKE) run-monte-carlo 2>&1 | tee -a "$$target_log") && \
	echo "End Time: $$(date '+%Y-%m-%d %H:%M:%S')" >> "$$target_log" && \
	echo "Status: SUCCESS" >> "$$target_log" || \
	echo "Status: FAILED" >> "$$target_log"

run-monte-carlo-paper-with-log: build | $(DATA_OUTPUT_DIR)
	@target_log="data/output/$(TIMESTAMP)_run_monte_carlo_paper.log"; \
	echo "Logging simulation to $$target_log"; \
	echo "============================================================" > "$$target_log"; \
	echo "Simulation Log: run-monte-carlo-paper" >> "$$target_log"; \
	echo "Start Time: $$(date '+%Y-%m-%d %H:%M:%S')" >> "$$target_log"; \
	echo "Compiler: $$(gfortran --version | head -n 1)" >> "$$target_log"; \
	echo "============================================================" >> "$$target_log"; \
	echo "" >> "$$target_log"; \
	env | grep '^SPDE_' | sort >> "$$target_log"; \
	echo "" >> "$$target_log"; \
	($(MAKE) run-monte-carlo-paper 2>&1 | tee -a "$$target_log") && \
	echo "End Time: $$(date '+%Y-%m-%d %H:%M:%S')" >> "$$target_log" && \
	echo "Status: SUCCESS" >> "$$target_log" || \
	echo "Status: FAILED" >> "$$target_log"

plot: | $(PLOT_DIR)
	@latest_file="$$(find $(DATA_OUTPUT_DIR) -maxdepth 1 -type f -name '$(LATEST_ESTIMATOR_PATTERN)' | sort | tail -n 1)"; \
	if [ -z "$$latest_file" ]; then \
		echo "Error: no generated estimator data was found in $(DATA_OUTPUT_DIR). First run the application with 'make run'." >&2; \
		exit 1; \
	fi; \
	$(PYTHON) $(PLOT_SCRIPT) --input "$$latest_file"

plot-snapshot-comparison: | $(DATA_OUTPUT_DIR) $(PLOT_DIR)
	@latest_file="$$( $(FIND_LATEST_SNAPSHOT_CMD) )"; \
	if [ -z "$$latest_file" ]; then \
		echo "Error: no snapshot comparison data found in $(DATA_OUTPUT_DIR). First run 'make run-snapshot-comparison'." >&2; \
		exit 1; \
	fi; \
	$(PYTHON) $(SNAPSHOT_PLOT_SCRIPT) --input "$$latest_file" $(SNAPSHOT_PLOT_ARGS)

video-snapshot-comparison: | $(DATA_OUTPUT_DIR) $(PLOT_DIR)
	@latest_file="$$( $(FIND_LATEST_SNAPSHOT_CMD) )"; \
	if [ -z "$$latest_file" ]; then \
		echo "Error: no snapshot comparison data found in $(DATA_OUTPUT_DIR). First run 'make run-snapshot-comparison'." >&2; \
		exit 1; \
	fi; \
	$(PYTHON) $(SNAPSHOT_VIDEO_SCRIPT) --input "$$latest_file" $(SNAPSHOT_VIDEO_ARGS)

plot-3d-rugosity: | $(DATA_OUTPUT_DIR) $(PLOT_DIR)
	@latest_file="$$( $(FIND_LATEST_SNAPSHOT_CMD) )"; \
	if [ -z "$$latest_file" ]; then \
		echo "Error: no snapshot comparison data found in $(DATA_OUTPUT_DIR). First run 'make run-snapshot-comparison'." >&2; \
		exit 1; \
	fi; \
	$(PYTHON) $(SNAPSHOT_3D_SCRIPT) --input "$$latest_file" $(SNAPSHOT_3D_ARGS)

export-snapshot-manuscript-panels: | $(DATA_OUTPUT_DIR) $(PLOT_DIR)
	@if [ -n "$(SNAPSHOT_MANUSCRIPT_INPUT)" ]; then \
		snapshot_file="$(SNAPSHOT_MANUSCRIPT_INPUT)"; \
	else \
		snapshot_file="$$( $(FIND_LATEST_SNAPSHOT_CMD) )"; \
	fi; \
	if [ -z "$$snapshot_file" ]; then \
		echo "Error: no snapshot comparison data found in $(DATA_OUTPUT_DIR). First run 'make run-snapshot-comparison'." >&2; \
		exit 1; \
	fi; \
	$(PYTHON) $(SNAPSHOT_MANUSCRIPT_SCRIPT) --input "$$snapshot_file" $(SNAPSHOT_MANUSCRIPT_ARGS)

export-snapshot-paraview: | $(DATA_OUTPUT_DIR) $(PARAVIEW_DIR)
	@if [ -n "$(SNAPSHOT_PARAVIEW_INPUT)" ]; then \
		snapshot_file="$(SNAPSHOT_PARAVIEW_INPUT)"; \
	else \
		snapshot_file="$$( $(FIND_LATEST_SNAPSHOT_CMD) )"; \
	fi; \
	if [ -z "$$snapshot_file" ]; then \
		echo "Error: no snapshot comparison data found in $(DATA_OUTPUT_DIR). First run 'make run-snapshot-comparison'." >&2; \
		exit 1; \
	fi; \
	$(PYTHON) $(SNAPSHOT_PARAVIEW_SCRIPT) --input "$$snapshot_file" --output-dir "$(SNAPSHOT_PARAVIEW_DIR)" $(SNAPSHOT_PARAVIEW_ARGS)
package-outputs: | $(DATA_OUTPUT_DIR)
	@if [ -n "$(LOG_FILE)" ]; then log_file="$(LOG_FILE)"; else log_file=$$($(FIND_LATEST_LOG_CMD)); fi; \
	if [ -z "$$log_file" ] || [ ! -f "$$log_file" ]; then \
		echo "Error: log file not found." >&2; \
		exit 1; \
	fi; \
	raw_time=$$(grep "Start Time:" "$$log_file" | sed "s/Start Time: //"); \
	if [ -z "$$raw_time" ]; then \
		echo "Error: could not find Start Time in $$log_file." >&2; \
		exit 1; \
	fi; \
	tag=$$(echo "$$raw_time" | sed "s/-/_/g; s/ /-/" | tr -d "\r"); \
	prefix=$$(basename "$$log_file" | sed "s/_run.*//"); \
	tar_name="$${tag}_sto_adv_diff_outputs.tar"; \
	echo "Packaging output files with prefix $$prefix into $$tar_name..."; \
	tar --force-local -cvf "$$tar_name" -C $(DATA_OUTPUT_DIR) $$(ls $(DATA_OUTPUT_DIR)/$${prefix}* | sed "s|.*/||"); \
	echo "Archive created: $$tar_name"
paraview-figure-video: | $(PARAVIEW_DIR)
	@latest_pvd="$$(find $(PARAVIEW_DIR) -maxdepth 2 -type f -name 'solution_snapshots.pvd' | sort | tail -n 1)"; \
	if [ -z "$$latest_pvd" ]; then \
		echo "Error: no solution_snapshots.pvd found in $(PARAVIEW_DIR). Run 'make export-snapshot-paraview' first." >&2; \
		exit 1; \
	fi; \
	start_time=$$(date +%s); \
	echo "Starting paraview-figure-video at $$(date '+%Y-%m-%d %H:%M:%S')"; \
	$(PVPYTHON) $(PARAVIEW_VIS_SCRIPT) \
		--input "$$latest_pvd" \
		--preset $(PARAVIEW_PRESET) \
		--output-image "$(PARAVIEW_IMAGE)" \
		--output-video "$(PARAVIEW_VIDEO)" \
		--frame-label \
		--show-legend \
		--size 1920 1080; \
	ret=$$?; \
	end_time=$$(date +%s); \
	elapsed=$$((end_time - start_time)); \
	echo "Finished paraview-figure-video at $$(date '+%Y-%m-%d %H:%M:%S') in $$elapsed seconds"; \
	exit $$ret

docs:
	mkdir -p $(DOCS_DIR)
	$(DOXYGEN) $(DOXYFILE)
	@test -d $(DOCS_HTML_DIR) || { echo "Error: Doxygen did not produce $(DOCS_HTML_DIR)"; exit 1; }
	touch $(DOCS_HTML_DIR)/.nojekyll

test: build
	tests/run_unit_tests.sh
	tests/smoke_test.sh
	tests/monte_carlo_smoke_test.sh
	tests/snapshot_comparison_smoke_test.sh

test-unit: build
	tests/run_unit_tests.sh

test-smoke: build
	tests/smoke_test.sh

test-monte-carlo: build
	tests/monte_carlo_smoke_test.sh

test-snapshot-comparison: build
	tests/snapshot_comparison_smoke_test.sh

check-large-files:
	scripts/check_large_files.sh tracked

setup-git-hooks:
	git config core.hooksPath .githooks

$(RUN_TARGET): $(COMMON_MODULE_OBJ) $(OBJ_DIR)/main.o | $(BIN_DIR)
	$(FC) $(FFLAGS) $(OMPFLAGS) -J $(MOD_DIR) -I $(MOD_DIR) -o $@ $^

$(MONTE_CARLO_TARGET): $(COMMON_MODULE_OBJ) $(OBJ_DIR)/monte_carlo_study.o | $(BIN_DIR)
	$(FC) $(FFLAGS) $(OMPFLAGS) -J $(MOD_DIR) -I $(MOD_DIR) -o $@ $^

$(SNAPSHOT_COMPARISON_TARGET): $(COMMON_MODULE_OBJ) $(OBJ_DIR)/snapshot_comparison.o | $(BIN_DIR)
	$(FC) $(FFLAGS) $(OMPFLAGS) -J $(MOD_DIR) -I $(MOD_DIR) -o $@ $^

$(OBJ_DIR)/%.o: $(APP_DIR)/%.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OMPFLAGS) -J $(MOD_DIR) -I $(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/%.o: $(MODULE_DIR)/%.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) $(OMPFLAGS) -J $(MOD_DIR) -I $(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/main.o: \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/driver_support_mod.o \
	$(OBJ_DIR)/validation_mod.o \
	$(OBJ_DIR)/csv_output_mod.o \
	$(OBJ_DIR)/spectral_operators_mod.o \
	$(OBJ_DIR)/progress_reporting_mod.o \
	$(OBJ_DIR)/sde_simulation_mod.o \
	$(OBJ_DIR)/parameter_ml_estimation_mod.o \
	$(OBJ_DIR)/workflow_mod.o
$(OBJ_DIR)/monte_carlo_study.o: \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/driver_support_mod.o \
	$(OBJ_DIR)/csv_output_mod.o \
	$(OBJ_DIR)/progress_reporting_mod.o \
	$(OBJ_DIR)/workflow_mod.o \
	$(OBJ_DIR)/monte_carlo_study_mod.o
$(OBJ_DIR)/snapshot_comparison.o: \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/driver_support_mod.o \
	$(OBJ_DIR)/csv_output_mod.o \
	$(OBJ_DIR)/workflow_mod.o
$(OBJ_DIR)/driver_support_mod.o: $(OBJ_DIR)/model_types_mod.o
$(OBJ_DIR)/validation_mod.o: $(OBJ_DIR)/model_types_mod.o
$(OBJ_DIR)/csv_output_mod.o: $(OBJ_DIR)/model_types_mod.o
$(OBJ_DIR)/spectral_operators_mod.o: $(OBJ_DIR)/model_types_mod.o
$(OBJ_DIR)/solution_reconstruction_mod.o: \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/spectral_operators_mod.o
$(OBJ_DIR)/sde_simulation_mod.o: \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/progress_reporting_mod.o
$(OBJ_DIR)/parameter_ml_estimation_mod.o: \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/progress_reporting_mod.o
$(OBJ_DIR)/workflow_mod.o: \
	$(OBJ_DIR)/driver_support_mod.o \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/parameter_ml_estimation_mod.o \
	$(OBJ_DIR)/sde_simulation_mod.o \
	$(OBJ_DIR)/solution_reconstruction_mod.o \
	$(OBJ_DIR)/spectral_operators_mod.o \
	$(OBJ_DIR)/validation_mod.o
$(OBJ_DIR)/monte_carlo_study_mod.o: \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/progress_reporting_mod.o \
	$(OBJ_DIR)/workflow_mod.o

$(OBJ_DIR) $(MOD_DIR) $(BIN_DIR) $(DATA_OUTPUT_DIR) $(PLOT_DIR) $(PARAVIEW_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR)
	find $(DATA_OUTPUT_DIR) -maxdepth 1 -type f \
		\( -name 'estimator_trajectory.csv' -o -name '*_estimator_trajectory.csv' -o \
		-name 'state_history.csv' -o -name '*_state_history.csv' -o \
		-name 'monte_carlo_summary.csv' -o -name '*_monte_carlo_summary.csv' -o \
		-name 'monte_carlo_replicates.csv' -o -name '*_monte_carlo_replicates.csv' -o \
		-name 'deterministic_path.csv' -o -name '*_deterministic_path.csv' -o \
		-name 'stochastic_path.csv' -o -name '*_stochastic_path.csv' \) -delete
	find $(PLOT_DIR) -mindepth 1 -maxdepth 1 -type d -name '*_manuscript_panels' -exec rm -rf {} +
	find $(PLOT_DIR) -maxdepth 1 -type f \
		\( -name 'estimator_trajectory.png' -o -name '*_estimator_trajectory.png' -o \
		-name 'solution_snapshot_comparison.png' -o -name '*_solution_snapshot_comparison.png' -o \
		-name 'solution_snapshot_comparison.gif' -o -name '*_solution_snapshot_comparison.gif' -o \
		-name 'solution_snapshot_comparison.mp4' -o -name '*_solution_snapshot_comparison.mp4' \) -delete
	find $(PARAVIEW_DIR) -mindepth 1 ! -name '.gitkeep' -exec rm -rf {} +
	rm -f *.o *.mod multivariate_modular monte_carlo_study snapshot_comparison

distclean: clean
	find tests/artifacts -mindepth 1 ! -name '.gitkeep' -delete
