#!/usr/bin/env Rscript
# Stage 04g -- Heterogeneity of the pre/post-1992 T6-entry cohort comparison
# by collaboration structure (draft.qmd Section VI.D "Heterogeneity by
# Collaboration Structure").
#
# Extends 04_cohort_compare_1992.R's pre/post-1992 T6-entry split (same
# t6_first_year cutoff, same academic-age window 0-25) from the interactive
# lambda_i to the additive alpha_i (fit SEPARATELY within each cohort, as
# 04_cohort_compare_1992.R does for lambda_i), then asks whether the
# widening of the alpha_i distribution documented there is concentrated in
# a particular collaboration subgroup:
#   - international collaboration: share of an author's papers with
#     coauthors affiliated with more than one country (career-long)
#   - team size: mean number of coauthors per paper (career-long)
#
# Also runs the direct interaction regression alpha_i ~ mean_team * post1993
# on the STANDARD (not per-cohort-split) alpha_i from 05_additive_fe_pctile.R,
# testing whether the *return* to team size changed across cohorts (as
# opposed to team size itself, or the unconditional dispersion of alpha_i,
# which the subgroup table already speaks to).
#
# Input : ../input/publication_history.csv, ../input/citation_10y.csv,
#         ../input/t6_authors.csv, ../output/alpha_additive_fe.csv
# Output: ../output/cohort_heterogeneity.csv (SD by subgroup x cohort)

suppressPackageStartupMessages(library(data.table))

ROOT <- dirname(dirname(this.path::this.path()))
AGE_MIN <- 0L; AGE_MAX <- 25L  # matches 04_cohort_compare_1992.R's window
MIN_OBS_PER_AUTHOR <- 2L
CUTOFF_YEAR <- 1992L          # matches 04_cohort_compare_1992.R (Hamermesh 2013)
RETURN_TEST_CUTOFF <- 1993L   # cutoff for the team-size RETURN interaction test below

fit_twoway_fe <- function(Y, mask_obs, n_iter = 200L, tol = 1e-10) {
  global_mean <- mean(Y[mask_obs]); alpha <- rep(0, nrow(Y)); gamma <- rep(0, ncol(Y))
  for (it in seq_len(n_iter)) {
    resid <- Y - outer(alpha, gamma, "+") - global_mean; resid[!mask_obs] <- NA
    row_adj <- rowMeans(resid, na.rm = TRUE); row_adj[is.na(row_adj)] <- 0
    alpha <- alpha + row_adj
    resid <- Y - outer(alpha, gamma, "+") - global_mean; resid[!mask_obs] <- NA
    col_adj <- colMeans(resid, na.rm = TRUE); col_adj[is.na(col_adj)] <- 0
    gamma <- gamma + col_adj
    if (max(abs(row_adj)) < tol && max(abs(col_adj)) < tol) break
  }
  alpha
}

build_alpha <- function(sub_pub) {
  valid <- sub_pub[age >= AGE_MIN & age <= AGE_MAX & !is.na(pctile)]
  cell <- valid[, .(pctile = mean(pctile)), by = .(author_id, age)]
  n_obs <- cell[, .N, by = author_id]
  keep <- n_obs[N >= MIN_OBS_PER_AUTHOR, author_id]
  cell <- cell[author_id %in% keep]
  authors <- sort(unique(cell$author_id)); ages <- AGE_MIN:AGE_MAX
  a_idx <- setNames(seq_along(authors), authors); t_idx <- setNames(seq_along(ages), ages)
  Y <- matrix(NA_real_, length(authors), length(ages))
  Y[cbind(a_idx[cell$author_id], t_idx[as.character(cell$age)])] <- cell$pctile
  mask_obs <- !is.na(Y)
  data.table(author_id = authors, alpha = fit_twoway_fe(Y, mask_obs))
}

main <- function() {
  t6 <- fread(file.path(ROOT, "input", "t6_authors.csv"), colClasses = "character")
  t6[, t6_first_year := as.numeric(first_year)]
  cohort_map <- setNames(t6$t6_first_year, t6$author_id)

  pub <- fread(file.path(ROOT, "input", "publication_history.csv"),
               select = c("author_id", "eid", "publication_year", "affiliation_country", "author_count"),
               colClasses = "character")
  pub[, publication_year := as.numeric(publication_year)]
  pub[, author_count := as.numeric(author_count)]
  c10 <- fread(file.path(ROOT, "input", "citation_10y.csv"), colClasses = list(character = "eid"))
  pub <- merge(pub, c10, by = "eid", all.x = TRUE)

  cohort <- unique(pub, by = "eid")
  cohort <- cohort[cited_by_10y_complete == 1 & !is.na(cited_by_10y)]
  cohort[, pctile := frank(cited_by_10y, ties.method = "average") / .N, by = publication_year]
  pub <- merge(pub, cohort[, .(eid, pctile)], by = "eid", all.x = TRUE)
  pub[, first_year := min(publication_year, na.rm = TRUE), by = author_id]
  pub[, age := publication_year - first_year]
  pub[, t6_first_year := cohort_map[author_id]]

  # career-long collaboration profile (independent of the age window above)
  pub[, n_countries := sapply(strsplit(affiliation_country, ";"), function(x) length(unique(trimws(x))))]
  pub[, is_intl := as.integer(n_countries >= 2)]
  author_profile <- pub[, .(prop_intl = mean(is_intl, na.rm = TRUE),
                             mean_team = mean(author_count, na.rm = TRUE)), by = author_id]

  # per-cohort additive alpha (mirrors 04_cohort_compare_1992.R's per-cohort lambda_i)
  pre  <- pub[t6_first_year < CUTOFF_YEAR]
  post <- pub[t6_first_year >= CUTOFF_YEAR]
  alpha_pre  <- merge(build_alpha(pre),  author_profile, by = "author_id")
  alpha_post <- merge(build_alpha(post), author_profile, by = "author_id")
  alpha_pre[, cohort := "pre"]; alpha_post[, cohort := "post"]
  d <- rbind(alpha_pre, alpha_post)
  cat(sprintf("[cohorts] pre N=%d, post N=%d\n", nrow(alpha_pre), nrow(alpha_post)))

  d[, g_intl := ifelse(prop_intl >= median(prop_intl), "high_intl", "low_intl")]
  d[, g_team := ifelse(mean_team >= median(mean_team), "large_team", "small_team")]

  het <- d[, .(N = .N, sd_alpha = sd(alpha)), by = .(subgroup = g_intl, cohort)]
  het <- rbind(het, d[, .(N = .N, sd_alpha = sd(alpha)), by = .(subgroup = g_team, cohort)])
  het_wide <- dcast(het, subgroup ~ cohort, value.var = c("N", "sd_alpha"))
  het_wide[, pct_increase := 100 * (sd_alpha_post / sd_alpha_pre - 1)]
  cat("\n=== SD(alpha_i) by collaboration subgroup x cohort ===\n")
  print(het_wide)
  fwrite(het_wide, file.path(ROOT, "output", "cohort_heterogeneity.csv"))

  # direct test: did the RETURN to team size change across cohorts?
  # uses the standard (non-cohort-split) alpha_i from 05_additive_fe_pctile.R
  alpha_std <- fread(file.path(ROOT, "output", "alpha_additive_fe.csv"), colClasses = list(character = "author_id"))
  d2 <- merge(alpha_std, author_profile, by = "author_id")
  d2 <- merge(d2, t6[, .(author_id, t6_first_year)], by = "author_id")
  d2[, post := as.integer(t6_first_year >= RETURN_TEST_CUTOFF)]
  m <- lm(alpha ~ mean_team * post, data = d2)
  cat(sprintf("\n=== Return to team size: alpha ~ mean_team * post%d ===\n", RETURN_TEST_CUTOFF))
  print(summary(m)$coefficients)
  cat("[done] cohort_heterogeneity.csv\n")
}

main()
