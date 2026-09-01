#!/usr/bin/env Rscript
# Stage 04c -- Interactive Fixed Effects with an EXPLICIT calendar-year
# fixed effect (alternative fix to the vintage confound; compare against
# 02_ife_cohort_percentile.R's implicit percentile-normalization fix).
#
#     Y_it = gamma_{year(i,t)} + lambda_i' F_t + eps_it
#
#   i = author, t = academic age, Y_it = raw log(1+cited_by_10y)
#   year(i,t) = author i's first-publication calendar year + academic age t
#
# Algorithm: alternate between (a) OLS of the year-effect-residualized
# outcome on year dummies given current (lambda, F), and (b) low-rank SVD
# completion of the residual to update (lambda, F), until beta converges
# (Bai 2009 with covariates).
#
# Also reports the year fixed effects gamma_y themselves, useful for
# spotting any sudden calendar-year jump/break in citation levels (there
# isn't one in this data -- see ../output/year_fe.csv: it's a smooth
# secular increase from 1927 to 2016, no discontinuity around any
# particular year).
#
# Input : ../input/publication_history.csv, ../input/citation_10y.csv
# Output: ../output/lambda_r1_yearfe.csv, ../output/Ft_r1_yearfe.csv,
#         ../output/year_fe.csv

suppressPackageStartupMessages(library(data.table))

ROOT <- dirname(dirname(this.path::this.path()))
AGE_MIN <- 0L; AGE_MAX <- 30L
MIN_OBS_PER_AUTHOR <- 2L
R <- 1L

main <- function() {
  pub <- fread(file.path(ROOT, "input", "publication_history.csv"),
               select = c("author_id", "eid", "publication_year"), colClasses = "character")
  pub[, publication_year := as.numeric(publication_year)]
  c10 <- fread(file.path(ROOT, "input", "citation_10y.csv"), colClasses = list(character = "eid"))
  pub <- merge(pub, c10, by = "eid", all.x = TRUE)

  pub[, first_year := min(publication_year, na.rm = TRUE), by = author_id]
  pub[, age := publication_year - first_year]

  valid <- pub[cited_by_10y_complete == 1 & age >= AGE_MIN & age <= AGE_MAX & !is.na(cited_by_10y)]
  valid[, logC10 := log1p(cited_by_10y)]

  cell <- valid[, .(logC10 = mean(logC10), year = round(mean(publication_year))), by = .(author_id, age)]
  n_obs <- cell[, .N, by = author_id]
  keep <- n_obs[N >= MIN_OBS_PER_AUTHOR, author_id]
  cell <- cell[author_id %in% keep]

  authors <- sort(unique(cell$author_id))
  ages <- AGE_MIN:AGE_MAX
  a_idx <- setNames(seq_along(authors), authors)
  t_idx <- setNames(seq_along(ages), ages)
  N <- length(authors); T <- length(ages)

  Y <- matrix(NA_real_, N, T)
  Y[cbind(a_idx[cell$author_id], t_idx[as.character(cell$age)])] <- cell$logC10
  mask_obs <- !is.na(Y)

  # calendar year is defined for EVERY (i,t) cell (= first_year_i + t),
  # observed or not
  fy <- setNames(pub[!duplicated(author_id), first_year], pub[!duplicated(author_id), author_id])
  YR <- outer(fy[authors], ages, "+")

  cat(sprintf("[panel] N=%s authors x T=%d ages, observed=%s/%s (%.1f%%)\n",
              format(N, big.mark = ","), T, format(sum(mask_obs), big.mark = ","),
              format(N * T, big.mark = ","), 100 * mean(mask_obs)))

  years_all <- sort(unique(as.integer(YR[mask_obs])))
  cat(sprintf("[years] %d-%d, %d distinct calendar years\n", min(years_all), max(years_all), length(years_all)))

  global_mean <- mean(Y[mask_obs])

  obs_idx <- which(mask_obs, arr.ind = TRUE)
  obs_i <- obs_idx[, 1]; obs_t <- obs_idx[, 2]
  y_obs <- Y[mask_obs]

  year_dummies <- function(rows_i, rows_t) {
    yrs <- as.integer(YR[cbind(rows_i, rows_t)])
    X <- matrix(0, length(yrs), length(years_all))
    X[cbind(seq_along(yrs), match(yrs, years_all))] <- 1
    X
  }
  X_obs <- year_dummies(obs_i, obs_t)

  beta <- as.numeric(.lm.fit(X_obs, y_obs - global_mean)$coefficients)

  lam <- NULL; Ft <- NULL
  for (outer_it in 1:30) {
    yhat_year_obs <- as.numeric(X_obs %*% beta) + global_mean
    resid_obs <- y_obs - yhat_year_obs

    R_mat <- matrix(0, N, T)
    R_mat[mask_obs] <- resid_obs
    Rc <- R_mat
    prev <- NULL
    Rhat <- NULL
    for (it in 1:150) {
      sv <- svd(Rc)
      Rhat <- sv$u[, 1:R, drop = FALSE] %*% diag(sv$d[1:R], R, R) %*% t(sv$v[, 1:R, drop = FALSE])
      Rc <- ifelse(mask_obs, R_mat, Rhat)
      if (!is.null(prev) && max(abs(Rhat - prev)) < 1e-7) break
      prev <- Rhat
    }
    lam <- sv$u[, 1:R, drop = FALSE] %*% diag(sv$d[1:R], R, R)
    Ft <- sv$v[, 1:R, drop = FALSE]

    factor_obs <- rowSums(lam[obs_i, , drop = FALSE] * Ft[obs_t, , drop = FALSE])
    target <- y_obs - global_mean - factor_obs
    beta_new <- as.numeric(.lm.fit(X_obs, target)$coefficients)

    delta <- max(abs(beta_new - beta))
    beta <- beta_new
    if (delta < 1e-6) {
      cat(sprintf("[converged] outer iter %d, beta delta=%.2e\n", outer_it, delta))
      break
    }
  }

  lam1 <- lam[, 1]; Ft1 <- Ft[, 1]
  if (mean(Ft1) < 0) { Ft1 <- -Ft1; lam1 <- -lam1 }

  yhat_full <- global_mean + as.numeric(X_obs %*% beta) + (lam[obs_i, 1] * Ft[obs_t, 1])
  r2 <- 1 - sum((y_obs - yhat_full)^2) / sum((y_obs - global_mean)^2)
  cat(sprintf("\nR2 (year FE + r=1 interactive age factor) = %.4f\n", r2))

  lam_dt <- data.table(author_id = authors, lambda_r1 = lam1)
  lam_dt[, pctile_of_lambda := frank(lambda_r1, ties.method = "average") / .N * 100]
  fwrite(lam_dt, file.path(ROOT, "output", "lambda_r1_yearfe.csv"))
  fwrite(data.table(age = ages, F_t_r1 = Ft1), file.path(ROOT, "output", "Ft_r1_yearfe.csv"))

  n_obs_per_year <- table(as.integer(YR[mask_obs]))
  year_fe <- data.table(year = years_all, gamma_y = beta)
  setorder(year_fe, year)
  year_fe[, n_obs := as.integer(n_obs_per_year[as.character(year)])]
  year_fe[is.na(n_obs), n_obs := 0L]
  year_fe[, delta_yoy := gamma_y - shift(gamma_y)]
  fwrite(year_fe, file.path(ROOT, "output", "year_fe.csv"))

  cat("[done] lambda_r1_yearfe.csv, Ft_r1_yearfe.csv, year_fe.csv\n")
}

main()
