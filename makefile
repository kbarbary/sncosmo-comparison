DATA_FILES := $(shell ls jla_light_curves/lc-*)
SNE := $(patsubst jla_light_curves/lc-%.list,%,$(DATA_FILES))

.PHONY: all

all: $(SNE)

$(SNE) : % : results_snfit/result-%.dat

results_snfit/result-%.dat: jla_light_curves/lc-%.list
	time bin/snfit $< -o $@ 2>&1 | tee results_snfit/result-$*.log
