FC = gfortran
FFLAGS ?= -std=f2008 -O3 -Wall -Wextra -Wimplicit-interface
PYTHON ?= python3
DOXYGEN ?= doxygen
TIMESTAMP ?= $(shell date '+%Y%m%dT%H%M%S')
LATEST_ESTIMATOR_PATTERN := *_estimator_trajectory.csv

SRC_DIR := src
APP_DIR := $(SRC_DIR)/app
MODULE_DIR := $(SRC_DIR)/modules
BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
MOD_DIR := $(BUILD_DIR)/mod
BIN_DIR := $(BUILD_DIR)/bin
TARGET := $(BIN_DIR)/multivariate_modular

DATA_OUTPUT_DIR := data/output
PLOT_DIR := visualization/plots
PLOT_SCRIPT := visualization/scripts/plot_estimator_trajectory.py
DOXYFILE := docs/Doxyfile
DOCS_DIR := $(BUILD_DIR)/docs/doxygen
DOCS_HTML_DIR := $(DOCS_DIR)/html

MODULE_SRC := \
	$(MODULE_DIR)/model_types_mod.f90 \
	$(MODULE_DIR)/validation_mod.f90 \
	$(MODULE_DIR)/csv_output_mod.f90 \
	$(MODULE_DIR)/spectral_operators_mod.f90 \
	$(MODULE_DIR)/progress_reporting_mod.f90 \
	$(MODULE_DIR)/sde_simulation_mod.f90 \
	$(MODULE_DIR)/parameter_ml_estimation_mod.f90
MAIN_SRC := $(APP_DIR)/main.f90
SRC := $(MODULE_SRC) $(MAIN_SRC)
OBJ := \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/validation_mod.o \
	$(OBJ_DIR)/csv_output_mod.o \
	$(OBJ_DIR)/spectral_operators_mod.o \
	$(OBJ_DIR)/progress_reporting_mod.o \
	$(OBJ_DIR)/sde_simulation_mod.o \
	$(OBJ_DIR)/parameter_ml_estimation_mod.o \
	$(OBJ_DIR)/main.o

.PHONY: all build run plot docs test test-smoke test-unit check-large-files setup-git-hooks clean distclean

all: build

build: $(TARGET)

run: build | $(DATA_OUTPUT_DIR)
	SARGAZO_OUTPUT_TIMESTAMP='$(TIMESTAMP)' $(TARGET)

plot: | $(PLOT_DIR)
	@latest_file="$$(find $(DATA_OUTPUT_DIR) -maxdepth 1 -type f -name '$(LATEST_ESTIMATOR_PATTERN)' | sort | tail -n 1)"; \
	if [ -z "$$latest_file" ]; then \
		echo "Error: no generated estimator data was found in $(DATA_OUTPUT_DIR). First run the application with 'make run'." >&2; \
		exit 1; \
	fi; \
	$(PYTHON) $(PLOT_SCRIPT) --input "$$latest_file"

docs:
	mkdir -p $(DOCS_DIR)
	$(DOXYGEN) $(DOXYFILE)
	@test -d $(DOCS_HTML_DIR) || { echo "Error: Doxygen did not produce $(DOCS_HTML_DIR)"; exit 1; }
	touch $(DOCS_HTML_DIR)/.nojekyll

test: build
	tests/run_unit_tests.sh
	tests/smoke_test.sh

test-unit: build
	tests/run_unit_tests.sh

test-smoke: build
	tests/smoke_test.sh

check-large-files:
	scripts/check_large_files.sh tracked

setup-git-hooks:
	git config core.hooksPath .githooks

$(TARGET): $(OBJ) | $(BIN_DIR)
	$(FC) $(FFLAGS) -J $(MOD_DIR) -I $(MOD_DIR) -o $@ $(OBJ)

$(OBJ_DIR)/main.o: $(MAIN_SRC) | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J $(MOD_DIR) -I $(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/%.o: $(MODULE_DIR)/%.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J $(MOD_DIR) -I $(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/main.o: \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/validation_mod.o \
	$(OBJ_DIR)/csv_output_mod.o \
	$(OBJ_DIR)/spectral_operators_mod.o \
	$(OBJ_DIR)/progress_reporting_mod.o \
	$(OBJ_DIR)/sde_simulation_mod.o \
	$(OBJ_DIR)/parameter_ml_estimation_mod.o
$(OBJ_DIR)/validation_mod.o: $(OBJ_DIR)/model_types_mod.o
$(OBJ_DIR)/csv_output_mod.o: $(OBJ_DIR)/model_types_mod.o
$(OBJ_DIR)/spectral_operators_mod.o: $(OBJ_DIR)/model_types_mod.o
$(OBJ_DIR)/sde_simulation_mod.o: \
	$(OBJ_DIR)/model_types_mod.o \
	$(OBJ_DIR)/progress_reporting_mod.o
$(OBJ_DIR)/parameter_ml_estimation_mod.o: $(OBJ_DIR)/model_types_mod.o

$(OBJ_DIR) $(MOD_DIR) $(BIN_DIR) $(DATA_OUTPUT_DIR) $(PLOT_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR)
	find $(DATA_OUTPUT_DIR) -maxdepth 1 -type f \
		\( -name 'estimator_trajectory.csv' -o -name '*_estimator_trajectory.csv' -o \
		-name 'state_history.csv' -o -name '*_state_history.csv' \) -delete
	find $(PLOT_DIR) -maxdepth 1 -type f \
		\( -name 'estimator_trajectory.png' -o -name '*_estimator_trajectory.png' \) -delete
	rm -f *.o *.mod multivariate_modular

distclean: clean
	find tests/artifacts -mindepth 1 ! -name '.gitkeep' -delete
