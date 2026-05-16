FC = gfortran
FFLAGS ?= -std=f2008 -O3 -Wall -Wextra -Wimplicit-interface
OMPFLAGS ?= -fopenmp
PYTHON ?= python3
DOXYGEN ?= doxygen
TIMESTAMP ?= $(shell date '+%Y%m%dT%H%M%S')
LATEST_ESTIMATOR_PATTERN := *_estimator_trajectory.csv
LATEST_SNAPSHOT_PATTERN := *_solution_snapshot_comparison.csv
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
PLOT_SCRIPT := visualization/scripts/plot_estimator_trajectory.py
SNAPSHOT_PLOT_SCRIPT := visualization/scripts/plot_solution_snapshot_comparison.py
SNAPSHOT_VIDEO_SCRIPT := visualization/scripts/animate_solution_snapshot_comparison.py
SNAPSHOT_3D_SCRIPT := visualization/scripts/plot_3d_rugosity_comparison.py
SNAPSHOT_PLOT_ARGS ?=
SNAPSHOT_VIDEO_ARGS ?=
SNAPSHOT_3D_ARGS ?=
DOXYFILE := docs/Doxyfile
DOCS_DIR := $(BUILD_DIR)/docs/doxygen
DOCS_HTML_DIR := $(DOCS_DIR)/html

MODULE_SRC := \
	$(MODULE_DIR)/model_types_mod.f90 \
	$(MODULE_DIR)/driver_support_mod.f90 \
	$(MODULE_DIR)/validation_mod.f90 \
	$(MODULE_DIR)/csv_output_mod.f90 \
	$(MODULE_DIR)/spectral_operators_mod.f90 \
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

.PHONY: all build run run-monte-carlo run-monte-carlo-paper run-snapshot-comparison plot plot-snapshot-comparison video-snapshot-comparison plot-3d-rugosity docs test test-smoke test-monte-carlo test-snapshot-comparison test-unit check-large-files setup-git-hooks clean distclean

all: build

build: $(RUN_TARGET) $(MONTE_CARLO_TARGET) $(SNAPSHOT_COMPARISON_TARGET)

run: build | $(DATA_OUTPUT_DIR)
	SARGAZO_OUTPUT_TIMESTAMP='$(TIMESTAMP)' $(RUN_TARGET)

run-monte-carlo: build | $(DATA_OUTPUT_DIR)
	SARGAZO_OUTPUT_TIMESTAMP='$(TIMESTAMP)' $(MONTE_CARLO_TARGET)

run-monte-carlo-paper: build | $(DATA_OUTPUT_DIR)
	SARGAZO_OUTPUT_TIMESTAMP='$(TIMESTAMP)' \
	SARGAZO_MC_REPLICATES='$(PAPER_MC_REPLICATES)' \
	SARGAZO_MC_BASIS_LEVELS='$(PAPER_MC_BASIS_LEVELS)' \
	SARGAZO_MC_N_OBSERVATIONS='$(PAPER_MC_N_OBSERVATIONS)' \
	SARGAZO_MC_TIME_STEPS='$(PAPER_MC_TIME_STEPS)' \
	$(MONTE_CARLO_TARGET)

run-snapshot-comparison: build | $(DATA_OUTPUT_DIR)
	SARGAZO_OUTPUT_TIMESTAMP='$(TIMESTAMP)' $(SNAPSHOT_COMPARISON_TARGET)

plot: | $(PLOT_DIR)
	@latest_file="$$(find $(DATA_OUTPUT_DIR) -maxdepth 1 -type f -name '$(LATEST_ESTIMATOR_PATTERN)' | sort | tail -n 1)"; \
	if [ -z "$$latest_file" ]; then \
		echo "Error: no generated estimator data was found in $(DATA_OUTPUT_DIR). First run the application with 'make run'." >&2; \
		exit 1; \
	fi; \
	$(PYTHON) $(PLOT_SCRIPT) --input "$$latest_file"

plot-snapshot-comparison: | $(DATA_OUTPUT_DIR) $(PLOT_DIR)
	@latest_file="$$(find $(DATA_OUTPUT_DIR) -maxdepth 1 -type f -name '$(LATEST_SNAPSHOT_PATTERN)' | sort | tail -n 1)"; \
	if [ -z "$$latest_file" ]; then \
		echo "Error: no snapshot comparison data found in $(DATA_OUTPUT_DIR). First run 'make run-snapshot-comparison'." >&2; \
		exit 1; \
	fi; \
	$(PYTHON) $(SNAPSHOT_PLOT_SCRIPT) --input "$$latest_file" $(SNAPSHOT_PLOT_ARGS)

video-snapshot-comparison: | $(DATA_OUTPUT_DIR) $(PLOT_DIR)
	@latest_file="$$(find $(DATA_OUTPUT_DIR) -maxdepth 1 -type f -name '$(LATEST_SNAPSHOT_PATTERN)' | sort | tail -n 1)"; \
	if [ -z "$$latest_file" ]; then \
		echo "Error: no snapshot comparison data found in $(DATA_OUTPUT_DIR). First run 'make run-snapshot-comparison'." >&2; \
		exit 1; \
	fi; \
	$(PYTHON) $(SNAPSHOT_VIDEO_SCRIPT) --input "$$latest_file" $(SNAPSHOT_VIDEO_ARGS)

plot-3d-rugosity: | $(DATA_OUTPUT_DIR) $(PLOT_DIR)
	@latest_file="$$(find $(DATA_OUTPUT_DIR) -maxdepth 1 -type f -name '$(LATEST_SNAPSHOT_PATTERN)' | sort | tail -n 1)"; \
	if [ -z "$$latest_file" ]; then \
		echo "Error: no snapshot comparison data found in $(DATA_OUTPUT_DIR). First run 'make run-snapshot-comparison'." >&2; \
		exit 1; \
	fi; \
	$(PYTHON) $(SNAPSHOT_3D_SCRIPT) --input "$$latest_file" $(SNAPSHOT_3D_ARGS)

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
	$(OBJ_DIR)/spectral_operators_mod.o \
	$(OBJ_DIR)/validation_mod.o
$(OBJ_DIR)/monte_carlo_study_mod.o: \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/progress_reporting_mod.o \
	$(OBJ_DIR)/workflow_mod.o

$(OBJ_DIR) $(MOD_DIR) $(BIN_DIR) $(DATA_OUTPUT_DIR) $(PLOT_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR)
	find $(DATA_OUTPUT_DIR) -maxdepth 1 -type f \
		\( -name 'estimator_trajectory.csv' -o -name '*_estimator_trajectory.csv' -o \
		-name 'state_history.csv' -o -name '*_state_history.csv' -o \
		-name 'monte_carlo_summary.csv' -o -name '*_monte_carlo_summary.csv' -o \
		-name 'monte_carlo_replicates.csv' -o -name '*_monte_carlo_replicates.csv' -o \
		-name 'solution_snapshot_comparison.csv' -o -name '*_solution_snapshot_comparison.csv' \) -delete
	find $(PLOT_DIR) -maxdepth 1 -type f \
		\( -name 'estimator_trajectory.png' -o -name '*_estimator_trajectory.png' -o \
		-name 'solution_snapshot_comparison.png' -o -name '*_solution_snapshot_comparison.png' -o \
		-name 'solution_snapshot_comparison.gif' -o -name '*_solution_snapshot_comparison.gif' -o \
		-name 'solution_snapshot_comparison.mp4' -o -name '*_solution_snapshot_comparison.mp4' \) -delete
	rm -f *.o *.mod multivariate_modular monte_carlo_study snapshot_comparison

distclean: clean
	find tests/artifacts -mindepth 1 ! -name '.gitkeep' -delete
