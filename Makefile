.PHONY: all clean task1 task2

R := Rscript --vanilla

# Auto-detect script locations (root or R/)
UTILS   := $(firstword $(wildcard R/00_utils.R 00_utils.R))
TASK1   := $(firstword $(wildcard R/10_task1_gap_hypercube.R 10_task1_gap_hypercube.R))
TASK2   := $(firstword $(wildcard R/20_task2_spectral_shells.R 20_task2_spectral_shells.R))

# Sanity checks
ifeq ($(strip $(UTILS)),)
  $(error Could not find 00_utils.R in R/ or project root)
endif
ifeq ($(strip $(TASK1)),)
  $(error Could not find 10_task1_gap_hypercube.R in R/ or project root)
endif
ifeq ($(strip $(TASK2)),)
  $(error Could not find 20_task2_spectral_shells.R in R/ or project root)
endif

all: task1 task2

# Task 1: Gap statistic on hypercube clusters
task1: output/figures/task1_gap_by_dim.png output/tables/task1_gap_thresholds.csv

output/figures/task1_gap_by_dim.png output/tables/task1_gap_thresholds.csv: $(TASK1) $(UTILS)
	@mkdir -p output/figures output/tables
	$(R) $(TASK1)

# Task 2: Spectral clustering on concentric shells
task2: output/figures/task2_gap_vs_radius.png output/tables/task2_gap_thresholds.csv output/figures/task2_shells_preview.html

output/figures/task2_gap_vs_radius.png output/tables/task2_gap_thresholds.csv output/figures/task2_shells_preview.html: $(TASK2) $(UTILS)
	@mkdir -p output/figures output/tables
	$(R) $(TASK2)

clean:
	rm -rf output/figures/*.png output/figures/*.html output/tables/*.csv
