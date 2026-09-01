SHELL := bash

.PHONY: all clean

all: \
	../output/lambda_r1_pctile.csv \
	../output/year_fe.csv \
	../output/lambda_pre1992.csv \
	../output/alpha_additive_fe.csv \
	../output/model_comparison_cv.csv \
	../output/cohort_heterogeneity.csv \
	../output/structural_breaks.csv \
	../output/peer_effects_colocated.csv

# --- upstream inputs, symlinked in from earlier pipeline stages ------------
# (mirrors doctor-tasks' `input/%: ../tasks/X/output/% | input/` convention:
# build the upstream stage first if its output is missing, then link it in)

../input/publication_history.csv: | ../input/
	$(MAKE) -C ../../02_publication_history/code
	ln -sf ../../02_publication_history/output/publication_history.csv $@

../input/citation_10y.csv: | ../input/
	$(MAKE) -C ../../03_citation_10y/code
	ln -sf ../../03_citation_10y/output/citation_10y.csv $@

../input/t6_authors.csv: | ../input/
	$(MAKE) -C ../../02_publication_history/code ../input/t6_authors.csv
	ln -sf ../../02_publication_history/input/t6_authors.csv $@

# ../input/nobel_priority_authors.csv is hand-curated, not derived -- no rule.

COMMON_INPUTS = ../input/publication_history.csv ../input/citation_10y.csv

# --- 02: interactive fixed-effects model (lambda_i * F_t) ------------------
../output/lambda_r1_pctile.csv ../output/Ft_r1_pctile.csv ../output/nobel_validation_pctile.csv: \
		02_ife_cohort_percentile.R $(COMMON_INPUTS) ../input/nobel_priority_authors.csv
	Rscript 02_ife_cohort_percentile.R

# --- 03: explicit calendar-year fixed effects (gamma_y) --------------------
../output/lambda_r1_yearfe.csv ../output/Ft_r1_yearfe.csv ../output/year_fe.csv: \
		03_ife_year_fe.R $(COMMON_INPUTS)
	Rscript 03_ife_year_fe.R

# --- 04: pre/post-1992 T6-entry cohort comparison (interactive lambda_i) ---
../output/lambda_pre1992.csv ../output/lambda_post1992.csv ../output/Ft_cohort_compare.csv: \
		04_cohort_compare_1992.R $(COMMON_INPUTS) ../input/t6_authors.csv
	Rscript 04_cohort_compare_1992.R

# --- 05: additive two-way FE (alpha_i + gamma_t) -- PRIMARY estimator ------
../output/alpha_additive_fe.csv ../output/gamma_additive_fe.csv ../output/nobel_validation_additive.csv: \
		05_additive_fe_pctile.R $(COMMON_INPUTS) ../input/nobel_priority_authors.csv
	Rscript 05_additive_fe_pctile.R

# --- 06: out-of-sample model comparison + placebo validation ---------------
../output/model_comparison_cv.csv: 06_model_comparison_cv.R $(COMMON_INPUTS)
	Rscript 06_model_comparison_cv.R

# --- 07: cohort x collaboration-structure heterogeneity --------------------
# depends on 05's alpha_additive_fe.csv
../output/cohort_heterogeneity.csv: \
		07_cohort_heterogeneity.R $(COMMON_INPUTS) ../input/t6_authors.csv ../output/alpha_additive_fe.csv
	Rscript 07_cohort_heterogeneity.R

# --- 08: calendar-year structural-break scan --------------------------------
# depends on 03's year_fe.csv
../output/structural_breaks.csv: 08_structural_break_scan.R ../output/year_fe.csv ../input/publication_history.csv
	Rscript 08_structural_break_scan.R

# --- 09: peer-environment design (coauthors vs. co-located non-coauthors) --
# depends on 05's alpha_additive_fe.csv
../output/peer_effects_coauthors.csv ../output/peer_effects_colocated.csv: \
		09_peer_effects.R $(COMMON_INPUTS) ../output/alpha_additive_fe.csv
	Rscript 09_peer_effects.R

clean:
	rm -f ../output/*.csv
