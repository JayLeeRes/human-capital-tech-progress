#!/usr/bin/env Rscript
# Stage 04e -- Additive two-way fixed-effects model on cohort-percentile
# citations (PRIMARY human-capital estimator used in draft.qmd).
#
#     Y_it = alpha_i + gamma_t + eps_it
#     i = author, t = academic age (years since author's first-ever publication)
#
# Same panel construction as 02_ife_cohort_percentile.R (cohort-percentile of
# cited_by_10y, academic ages 0-30, authors with >=2 observed ages), but fit
# with the additive two-way fixed-effects model instead of the interactive
# (Bai 2009) rank-1 factor model. This is the model draft.qmd's Section III.C
# selects as the primary specification, on the strength of the out-of-sample
# model-comparison exercise in draft.qmd's Table "Model comparison via
# five-fold out-of-sample cross-validation" (interactive and additive tie at
# R2=0.132 out-of-sample; additive is preferred for its simplicity and
# absence of rotation/scale indeterminacy -- see 02_ife_cohort_percentile.R's
# header for that indeterminacy).
#
# Input : ../input/publication_history.csv, ../input/citation_10y.csv,
#         ../input/nobel_priority_authors.csv
# Output: ../output/alpha_additive_fe.csv, ../output/gamma_additive_fe.csv,
#         ../output/nobel_validation_additive.csv

suppressPackageStartupMessages(library(data.table))

ROOT <- dirname(dirname(this.path::this.path()))
AGE_MIN <- 0L; AGE_MAX <- 30L
MIN_OBS_PER_AUTHOR <- 2L

# Fits Y_it = alpha_i + gamma_t + eps_it by iterative demeaning: repeatedly
# subtract the current row means (author effect) then the current column
# means (age effect) from the observed-cell residuals, until both stop
# moving. Equivalent to the within estimator for a two-way fixed-effects
# panel with missing cells.
fit_twoway_fe <- function(Y, mask_obs, n_iter = 200L, tol = 1e-10) {
  global_mean <- mean(Y[mask_obs])
  alpha <- rep(0, nrow(Y)); gamma <- rep(0, ncol(Y))
  for (it in seq_len(n_iter)) {
    resid <- Y - outer(alpha, gamma, "+") - global_mean; resid[!mask_obs] <- NA
    row_adj <- rowMeans(resid, na.rm = TRUE); row_adj[is.na(row_adj)] <- 0
    alpha <- alpha + row_adj
    resid <- Y - outer(alpha, gamma, "+") - global_mean; resid[!mask_obs] <- NA
    col_adj <- colMeans(resid, na.rm = TRUE); col_adj[is.na(col_adj)] <- 0
    gamma <- gamma + col_adj
    if (max(abs(row_adj)) < tol && max(abs(col_adj)) < tol) break
  }
  Yhat <- outer(alpha, gamma, "+") + global_mean
  resid_obs <- Y[mask_obs] - Yhat[mask_obs]
  r2 <- 1 - sum(resid_obs^2) / sum((Y[mask_obs] - global_mean)^2)
  list(alpha = alpha, gamma = gamma, r2 = r2, iters = it)
}

main <- function() {
  pub <- fread(file.path(ROOT, "input", "publication_history.csv"),
               select = c("author_id", "eid", "publication_year"), colClasses = "character")
  pub[, publication_year := as.numeric(publication_year)]
  c10 <- fread(file.path(ROOT, "input", "citation_10y.csv"), colClasses = list(character = "eid"))
  pub <- merge(pub, c10, by = "eid", all.x = TRUE)

  cohort <- unique(pub, by = "eid")
  cohort <- cohort[cited_by_10y_complete == 1 & !is.na(cited_by_10y)]
  cohort[, pctile := frank(cited_by_10y, ties.method = "average") / .N, by = publication_year]
  pub <- merge(pub, cohort[, .(eid, pctile)], by = "eid", all.x = TRUE)

  pub[, first_year := min(publication_year, na.rm = TRUE), by = author_id]
  pub[, age := publication_year - first_year]

  valid <- pub[age >= AGE_MIN & age <= AGE_MAX & !is.na(pctile)]
  cell <- valid[, .(pctile = mean(pctile)), by = .(author_id, age)]
  n_obs <- cell[, .N, by = author_id]
  keep <- n_obs[N >= MIN_OBS_PER_AUTHOR, author_id]
  cell <- cell[author_id %in% keep]

  authors <- sort(unique(cell$author_id))
  ages <- AGE_MIN:AGE_MAX
  a_idx <- setNames(seq_along(authors), authors)
  t_idx <- setNames(seq_along(ages), ages)
  N <- length(authors); T <- length(ages)

  Y <- matrix(NA_real_, N, T)
  Y[cbind(a_idx[cell$author_id], t_idx[as.character(cell$age)])] <- cell$pctile
  mask_obs <- !is.na(Y)
  cat(sprintf("[panel] N=%s authors x T=%d ages, observed=%s/%s (%.1f%%)\n",
              format(N, big.mark = ","), T, format(sum(mask_obs), big.mark = ","),
              format(N * T, big.mark = ","), 100 * mean(mask_obs)))

  fit <- fit_twoway_fe(Y, mask_obs)
  cat(sprintf("in-sample R2 = %.4f (converged in %d iters)\n", fit$r2, fit$iters))

  alpha <- fit$alpha
  cat(sprintf("\n[alpha distribution] mean=%.4f median=%.4f sd=%.4f min=%.4f max=%.4f\n",
              mean(alpha), median(alpha), sd(alpha), min(alpha), max(alpha)))
  if (requireNamespace("moments", quietly = TRUE)) {
    set.seed(1)
    samp <- if (length(alpha) > 5000) sample(alpha, 5000) else alpha
    sw <- shapiro.test(samp)
    cat(sprintf("[alpha distribution] skewness=%.4f excess kurtosis=%.4f Shapiro-Wilk W=%.4f (n=%d)\n",
                moments::skewness(alpha), moments::kurtosis(alpha) - 3, sw$statistic, length(samp)))
  } else {
    cat("[alpha distribution] 'moments' package not installed; skipping skewness/kurtosis/Shapiro-Wilk\n")
  }

  alpha_dt <- data.table(author_id = authors, alpha = alpha)
  alpha_dt[, pctile_of_alpha := frank(alpha, ties.method = "average") / .N * 100]
  fwrite(alpha_dt, file.path(ROOT, "output", "alpha_additive_fe.csv"))
  fwrite(data.table(age = ages, gamma_t = fit$gamma), file.path(ROOT, "output", "gamma_additive_fe.csv"))

  # sanity check: known Nobel laureates should skew toward high alpha percentiles
  pri <- fread(file.path(ROOT, "input", "nobel_priority_authors.csv"), colClasses = "character")
  nobel <- pri[priority_group == "nobel"]
  merged <- merge(nobel, alpha_dt, by = "author_id")
  setorder(merged, -pctile_of_alpha)
  fwrite(merged[, .(author_id, label, alpha, pctile_of_alpha)],
         file.path(ROOT, "output", "nobel_validation_additive.csv"))
  cat(sprintf("\n[validation] %d Nobel-laureate profiles found in panel; mean alpha percentile=%.1f (50=random)\n",
              nrow(merged), mean(merged$pctile_of_alpha)))
  cat("[done] alpha_additive_fe.csv, gamma_additive_fe.csv, nobel_validation_additive.csv\n")
}

main()
