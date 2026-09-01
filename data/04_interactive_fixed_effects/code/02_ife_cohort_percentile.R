#!/usr/bin/env Rscript
# Stage 04b -- Interactive Fixed Effects (Bai 2009), cohort-percentile-normalized.
#
#     Y_it = lambda_i' F_t + eps_it
#     i = author, t = academic age (years since author's first-ever publication)
#
# Y_it is the paper's PERCENTILE RANK of cited_by_10y within its own
# publication-year cohort (0-1), instead of raw log(1+cited_by_10y). This
# removes the calendar-era citation-scale confound documented in
# code/_archive/01_ife_raw_baseline.R (not part of the draft.qmd analysis):
# every era is put on the same 0-1 "how well did this paper do relative to
# its peers published the same year" scale, so lambda_i reflects relative
# standing rather than which decade someone published in.
#
# NOTE: draft.qmd's out-of-sample model comparison (Section VI.B) finds this
# interactive model statistically indistinguishable from the much simpler
# additive alternative in 05_additive_fe_pctile.R, which is preferred there
# as the primary human-capital estimator. This script is kept because its
# lambda_i is still used for the interactive-model robustness check in
# Section VI.C and the r=1..5 R2 comparison motivating that model-selection
# exercise.
#
# Input : ../input/publication_history.csv, ../input/citation_10y.csv
# Output: ../output/lambda_r1_pctile.csv, ../output/Ft_r1_pctile.csv
#         ../output/nobel_validation_pctile.csv (sanity check vs known laureates)

suppressPackageStartupMessages(library(data.table))

ROOT <- dirname(dirname(this.path::this.path()))
AGE_MIN <- 0L; AGE_MAX <- 30L
MIN_OBS_PER_AUTHOR <- 2L
R_MAX <- 5L

# Fits Y_it = lambda_i' F_t + eps_it via "hard-impute" matrix completion:
# repeatedly (1) do a rank-r SVD on the current best guess of the full
# N x T matrix, (2) use that low-rank reconstruction to refill only the
# cells we never observed, (3) repeat until the reconstruction stops
# moving. Observed cells (mask_obs == TRUE) are always reset back to
# their true value at the start of each loop, so the SVD is never
# allowed to "smooth over" real data -- only to guess the missing cells.
fit_ife <- function(Y, mask_obs, r, n_iter = 200L, tol = 1e-6) {
  global_mean <- mean(Y[mask_obs])
  Yc <- Y
  Yc[!mask_obs] <- global_mean  # initial guess for missing cells: the overall mean
  prev <- NULL
  Yhat <- NULL
  for (it in seq_len(n_iter)) {
    sv <- svd(Yc)
    idx <- seq_len(r)
    # keep only the top r singular vectors/values -> best rank-r approximation of Yc
    Yhat <- sv$u[, idx, drop = FALSE] %*% diag(sv$d[idx], r, r) %*% t(sv$v[, idx, drop = FALSE])
    Yc <- ifelse(mask_obs, Y, Yhat)  # real data stays real; only missing cells get the new guess
    if (!is.null(prev) && max(abs(Yhat - prev)) < tol) break  # stop once the guess barely changes
    prev <- Yhat
  }
  lam <- sv$u[, idx, drop = FALSE] %*% diag(sv$d[idx], r, r)  # author loadings (lambda_i), N x r
  Ft <- sv$v[, idx, drop = FALSE]                             # age factors (F_t), T x r
  resid <- Y[mask_obs] - Yhat[mask_obs]
  r2 <- 1 - sum(resid^2) / sum((Y[mask_obs] - global_mean)^2)  # fit quality on observed cells only
  list(lam = lam, Ft = Ft, r2 = r2, iters = it)
}

main <- function() {
  pub <- fread(file.path(ROOT, "input", "publication_history.csv"),
               select = c("author_id", "eid", "publication_year"), colClasses = "character")
  pub[, publication_year := as.numeric(publication_year)]
  c10 <- fread(file.path(ROOT, "input", "citation_10y.csv"), colClasses = list(character = "eid"))
  pub <- merge(pub, c10, by = "eid", all.x = TRUE)

  # Step 1: cohort percentile computed on the full corpus (one row per unique
  # eid), so it's the true "same publication-year peer set" regardless of how
  # many T6 co-authors a paper has. e.g. pctile = 0.9 means "this paper got
  # more 10y citations than 90% of papers published the same year."
  cohort <- unique(pub, by = "eid")
  cohort <- cohort[cited_by_10y_complete == 1 & !is.na(cited_by_10y)]
  cohort[, pctile := frank(cited_by_10y, ties.method = "average") / .N, by = publication_year]
  pub <- merge(pub, cohort[, .(eid, pctile)], by = "eid", all.x = TRUE)

  # Step 2: academic age = years since each author's first publication (0, 1, 2, ...)
  pub[, first_year := min(publication_year, na.rm = TRUE), by = author_id]
  pub[, age := publication_year - first_year]

  # Step 3: collapse to one row per (author, age) -- if an author published
  # several papers at the same age, average their percentiles into one cell
  valid <- pub[age >= AGE_MIN & age <= AGE_MAX & !is.na(pctile)]
  cell <- valid[, .(pctile = mean(pctile)), by = .(author_id, age)]
  n_obs <- cell[, .N, by = author_id]
  keep <- n_obs[N >= MIN_OBS_PER_AUTHOR, author_id]  # drop authors with too few (author, age) cells to fit
  cell <- cell[author_id %in% keep]

  # Step 4: reshape the long (author, age, pctile) table into an
  # authors x ages matrix Y. Most cells are NA (that author had no
  # publication at that age) -- mask_obs marks which cells are real data.
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

  # Fit rank r = 1..R_MAX just to see how much extra variance higher ranks
  # explain (printed as R2), but only r=1 is actually used/saved below --
  # r=1 gives a single scalar lambda_i per author, i.e. "one number that
  # summarizes how good this author is," which is what we want.
  fit1 <- NULL
  for (r in 1:R_MAX) {
    fit <- fit_ife(Y, mask_obs, r)
    cat(sprintf("r=%d: R2=%.4f (converged in %d iters)\n", r, fit$r2, fit$iters))
    if (r == 1) fit1 <- fit
  }

  # SVD sign is arbitrary (u,v could both flip sign and represent the same
  # fit) -- pin it down so higher lambda always means higher F_t/citations,
  # not the reverse, by forcing the age-profile F_t to have positive mean.
  Ft1 <- fit1$Ft[, 1]; lam1 <- fit1$lam[, 1]
  if (mean(Ft1) < 0) { Ft1 <- -Ft1; lam1 <- -lam1 }

  cat("\n=== r=1: F_t over academic age ===\n")
  for (i in seq_along(ages)) cat(sprintf("  age %2d  F_t=%+.3f\n", ages[i], Ft1[i]))

  # Convert lambda (an arbitrary-scale number) into a 0-100 percentile
  # ranking of authors, which is easier to interpret and compare.
  lam_dt <- data.table(author_id = authors, lambda_r1 = lam1)
  lam_dt[, pctile_of_lambda := frank(lambda_r1, ties.method = "average") / .N * 100]
  fwrite(lam_dt, file.path(ROOT, "output", "lambda_r1_pctile.csv"))
  fwrite(data.table(age = ages, F_t_r1 = Ft1), file.path(ROOT, "output", "Ft_r1_pctile.csv"))

  # sanity check: known Nobel laureates should skew toward high lambda percentiles
  pri <- fread(file.path(ROOT, "input", "nobel_priority_authors.csv"), colClasses = "character")
  nobel <- pri[priority_group == "nobel"]
  merged <- merge(nobel, lam_dt, by = "author_id")
  setorder(merged, -pctile_of_lambda)
  fwrite(merged[, .(author_id, label, lambda_r1, pctile_of_lambda)],
         file.path(ROOT, "output", "nobel_validation_pctile.csv"))
  cat(sprintf("\n[validation] %d Nobel-laureate profiles found in panel; mean lambda percentile=%.1f (50=random)\n",
              nrow(merged), mean(merged$pctile_of_lambda)))
  cat("[done] lambda_r1_pctile.csv, Ft_r1_pctile.csv, nobel_validation_pctile.csv\n")
}

main()
