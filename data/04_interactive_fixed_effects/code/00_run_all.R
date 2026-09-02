#!/usr/bin/env Rscript
# Stage 04 orchestrator -- every analysis behind draft.qmd ("Who Benefits
# from Technological Progress"), run in the order the paper needs them.
#
#   02  Interactive fixed-effects model (Bai 2009), lambda_i * F_t
#         -> Section VI.C (lambda_i robustness), motivates model selection
#   03  Explicit calendar-year fixed effects; writes year_fe.csv
#         -> input to 08's structural-break scan (Section VI.E)
#   04  Pre/post-1992 T6-entry cohort comparison, interactive lambda_i
#         -> Section VI.C
#   05  Additive two-way FE (alpha_i + gamma_t) -- PRIMARY human-capital
#       estimator selected in Section VI.B; writes alpha_additive_fe.csv,
#       used as an input by 06, 07, and 09 below
#         -> Section VI.A (alpha_i distribution), Section VII.B (Nobel check)
#   06  Out-of-sample cross-validation + placebo validation of 02 vs 05
#         -> Section VI.B, Section VII.A
#   07  Cohort heterogeneity by collaboration structure (team size,
#       international collaboration) + team-size return interaction test
#         -> Section VI.D
#   08  Quandt-Andrews structural-break scan of year_fe.csv and auxiliary
#       corpus-coverage series
#         -> Section VI.E
#   09  Peer-environment design: actual coauthors vs co-located non-
#       coauthors, with and without institution fixed effects
#         -> Section VI.C
#   10  How trustworthy is alpha_i: split-half reliability, and whether it
#       is actually time-invariant (early- vs late-career alpha)
#   11  Does 09's congestion result survive measuring peers at the time of
#       co-location rather than by their career-average alpha?
#   12  Monte Carlo confirming the asymptotic rates behind 10's reliability:
#       lambda_i converges at sqrt(T_i) and not at all in N (simulation only,
#       reads no data)
#
# 06, 07, and 09 depend on 05's output (alpha_additive_fe.csv); 08 depends
# on 03's output (year_fe.csv). The order below respects both.
#
# code/_archive/ holds two earlier exploratory scripts not used by
# draft.qmd: a raw (non-percentile) baseline model and a rank-2 factor
# experiment.
#
# All steps are pure local computation (no external API calls). A full run
# takes several minutes: 02/05 are dominated by SVD/demeaning iterations
# over the ~10k-author panel, 06 repeats that five times per model for
# cross-validation, and 09's co-located-peer step loops over ~13k focal
# authors.

code_dir <- dirname(this.path::this.path())
steps <- c(
  "02_ife_cohort_percentile.R",
  "03_ife_year_fe.R",
  "04_cohort_compare_1992.R",
  "05_additive_fe_pctile.R",
  "06_model_comparison_cv.R",
  "07_cohort_heterogeneity.R",
  "08_structural_break_scan.R",
  "09_peer_effects.R",
  "10_alpha_reliability.R",
  "11_peer_effects_contemporaneous.R",
  "12_ife_consistency_mc.R"
)

for (step in steps) {
  cat(sprintf("\n=== %s ===\n", step))
  status <- system2("Rscript", file.path(code_dir, step))
  if (status != 0) {
    cat(sprintf("[stop] %s exited with %d\n", step, status))
    quit(status = status)
  }
}
