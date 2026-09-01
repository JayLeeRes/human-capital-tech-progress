#!/usr/bin/env Rscript
# Stage 04f -- Out-of-sample model comparison and placebo validation
# (draft.qmd Section VI.B "Model Selection", Section VII.A "Out-of-Sample
# and Placebo Validation").
#
# Compares, via 5-fold cross-validation on observed author-age cells of the
# standard cohort-percentile panel (same construction as
# 02_ife_cohort_percentile.R / 05_additive_fe_pctile.R):
#   - two uninformative baselines (author-mean only, academic-age-mean only)
#   - the additive two-way FE model (alpha_i + gamma_t)
#   - the interactive rank-1 factor model (lambda_i * F_t, Bai 2009)
#   - a placebo: the additive model refit on a version of the panel where
#     observed cohort-percentile values are randomly permuted across cells,
#     holding the missingness pattern fixed
#
# The placebo's in-sample R2 can be substantial (many free per-author
# parameters overfit noise) but its out-of-sample R2 should collapse to
# roughly zero or below; a real, generalizable model should not.
#
# Input : ../input/publication_history.csv, ../input/citation_10y.csv
# Output: ../output/model_comparison_cv.csv (per-fold R2 by model)

suppressPackageStartupMessages(library(data.table))

ROOT <- dirname(dirname(this.path::this.path()))
AGE_MIN <- 0L; AGE_MAX <- 30L
MIN_OBS_PER_AUTHOR <- 2L
K_FOLDS <- 5L
SEED <- 42L

fit_ife_r1 <- function(Y, mask_obs, n_iter = 200L, tol = 1e-6) {
  global_mean <- mean(Y[mask_obs]); Yc <- Y; Yc[!mask_obs] <- global_mean
  prev <- NULL; Yhat <- NULL
  for (it in seq_len(n_iter)) {
    sv <- svd(Yc)
    Yhat <- sv$u[, 1, drop = FALSE] %*% sv$d[1] %*% t(sv$v[, 1, drop = FALSE])
    Yc <- ifelse(mask_obs, Y, Yhat)
    if (!is.null(prev) && max(abs(Yhat - prev)) < tol) break
    prev <- Yhat
  }
  Yhat
}

fit_twoway_fe <- function(Y, mask_obs, n_iter = 100L, tol = 1e-8) {
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
  outer(alpha, gamma, "+") + global_mean
}

build_panel <- function(pub) {
  valid <- pub[age >= AGE_MIN & age <= AGE_MAX & !is.na(pctile)]
  cell <- valid[, .(pctile = mean(pctile)), by = .(author_id, age)]
  n_obs <- cell[, .N, by = author_id]
  keep <- n_obs[N >= MIN_OBS_PER_AUTHOR, author_id]
  cell <- cell[author_id %in% keep]
  authors <- sort(unique(cell$author_id)); ages <- AGE_MIN:AGE_MAX
  a_idx <- setNames(seq_along(authors), authors); t_idx <- setNames(seq_along(ages), ages)
  Y <- matrix(NA_real_, length(authors), length(ages))
  Y[cbind(a_idx[cell$author_id], t_idx[as.character(cell$age)])] <- cell$pctile
  list(Y = Y, mask = !is.na(Y))
}

main <- function() {
  set.seed(SEED)

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

  panel <- build_panel(pub)
  Y_full <- panel$Y; mask_full <- panel$mask
  cat(sprintf("[panel] N=%d authors x T=%d ages, observed=%d/%d (%.1f%%)\n",
              nrow(Y_full), ncol(Y_full), sum(mask_full), length(Y_full), 100 * mean(mask_full)))

  obs_idx <- which(mask_full)
  folds <- sample(rep(1:K_FOLDS, length.out = length(obs_idx)))

  Y_placebo <- Y_full
  Y_placebo[mask_full] <- sample(Y_full[mask_full])

  results <- data.table(model = character(), fold = integer(), r2 = numeric())
  for (k in 1:K_FOLDS) {
    test_idx <- obs_idx[folds == k]
    mask_train <- mask_full; mask_train[test_idx] <- FALSE
    train_mean <- mean(Y_full[mask_train])

    row_mean <- rowMeans(ifelse(mask_train, Y_full, NA), na.rm = TRUE)
    pred_author <- row_mean[row(Y_full)[test_idx]]; pred_author[is.na(pred_author)] <- train_mean
    r2_author <- 1 - sum((Y_full[test_idx] - pred_author)^2) / sum((Y_full[test_idx] - train_mean)^2)

    col_mean <- colMeans(ifelse(mask_train, Y_full, NA), na.rm = TRUE)
    pred_age <- col_mean[col(Y_full)[test_idx]]; pred_age[is.na(pred_age)] <- train_mean
    r2_age <- 1 - sum((Y_full[test_idx] - pred_age)^2) / sum((Y_full[test_idx] - train_mean)^2)

    Yhat_add <- fit_twoway_fe(Y_full, mask_train)
    r2_add <- 1 - sum((Y_full[test_idx] - Yhat_add[test_idx])^2) / sum((Y_full[test_idx] - train_mean)^2)

    Yhat_ife <- fit_ife_r1(Y_full, mask_train)
    r2_ife <- 1 - sum((Y_full[test_idx] - Yhat_ife[test_idx])^2) / sum((Y_full[test_idx] - train_mean)^2)

    placebo_train_mean <- mean(Y_placebo[mask_train])
    Yhat_placebo <- fit_twoway_fe(Y_placebo, mask_train)
    r2_placebo <- 1 - sum((Y_placebo[test_idx] - Yhat_placebo[test_idx])^2) /
      sum((Y_placebo[test_idx] - placebo_train_mean)^2)

    results <- rbind(results, data.table(
      model = c("author_mean_only", "age_mean_only", "additive_fe", "interactive_r1", "placebo_additive"),
      fold = k, r2 = c(r2_author, r2_age, r2_add, r2_ife, r2_placebo)))
  }

  summary_tbl <- results[, .(mean_r2 = mean(r2), sd_r2 = sd(r2)), by = model][order(-mean_r2)]
  cat("\n=== out-of-sample R2 (5-fold mean) ===\n")
  print(summary_tbl)

  # in-sample R2 for the two real models and the placebo, for comparison
  Yhat_add_full <- fit_twoway_fe(Y_full, mask_full)
  r2_add_in <- 1 - sum((Y_full[mask_full] - Yhat_add_full[mask_full])^2) / sum((Y_full[mask_full] - mean(Y_full[mask_full]))^2)
  Yhat_ife_full <- fit_ife_r1(Y_full, mask_full)
  r2_ife_in <- 1 - sum((Y_full[mask_full] - Yhat_ife_full[mask_full])^2) / sum((Y_full[mask_full] - mean(Y_full[mask_full]))^2)
  Yhat_placebo_full <- fit_twoway_fe(Y_placebo, mask_full)
  r2_placebo_in <- 1 - sum((Y_placebo[mask_full] - Yhat_placebo_full[mask_full])^2) / sum((Y_placebo[mask_full] - mean(Y_placebo[mask_full]))^2)
  cat(sprintf("\n[in-sample] additive=%.4f  interactive_r1=%.4f  placebo=%.4f\n", r2_add_in, r2_ife_in, r2_placebo_in))

  fwrite(results, file.path(ROOT, "output", "model_comparison_cv.csv"))
  cat("[done] model_comparison_cv.csv\n")
}

main()
