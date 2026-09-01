#!/usr/bin/env Rscript
# Stage 04d -- Split authors into two T6-ENTRY cohorts (first T6 publication
# year <1992 vs >=1992; 1992 follows Hamermesh 2013's early-1990s empirical
# turn used elsewhere in this project as a cohort demarcation point), refit
# the cohort-percentile IFE(r=1) model (see 02_ife_cohort_percentile.R)
# separately in each cohort, and compare:
#   (a) the common academic-age factor F_t shape
#   (b) the lambda_i distribution (mean/median/spread, KS test, t-test)
#
# Input : ../input/publication_history.csv, ../input/citation_10y.csv,
#         ../input/t6_authors.csv (first_year = first T6 publication year)
# Output: ../output/lambda_pre1992.csv, ../output/lambda_post1992.csv,
#         ../output/Ft_cohort_compare.csv

suppressPackageStartupMessages(library(data.table))

ROOT <- dirname(dirname(this.path::this.path()))
AGE_MIN <- 0L; AGE_MAX <- 25L  # shorter than the other stages so the (younger) post-1992 cohort has enough data
MIN_OBS_PER_AUTHOR <- 2L
CUTOFF_YEAR <- 1992

fit_ife <- function(Y, mask_obs, r = 1L, n_iter = 200L, tol = 1e-6) {
  global_mean <- mean(Y[mask_obs])
  Yc <- Y
  Yc[!mask_obs] <- global_mean
  prev <- NULL
  Yhat <- NULL
  for (it in seq_len(n_iter)) {
    sv <- svd(Yc)
    idx <- seq_len(r)
    Yhat <- sv$u[, idx, drop = FALSE] %*% diag(sv$d[idx], r, r) %*% t(sv$v[, idx, drop = FALSE])
    Yc <- ifelse(mask_obs, Y, Yhat)
    if (!is.null(prev) && max(abs(Yhat - prev)) < tol) break
    prev <- Yhat
  }
  lam <- sv$u[, idx, drop = FALSE] %*% diag(sv$d[idx], r, r)
  Ft <- sv$v[, idx, drop = FALSE]
  resid <- Y[mask_obs] - Yhat[mask_obs]
  r2 <- 1 - sum(resid^2) / sum((Y[mask_obs] - global_mean)^2)
  list(lam = lam[, 1], Ft = Ft[, 1], r2 = r2)
}

build_and_fit <- function(sub_pub, label) {
  valid <- sub_pub[age >= AGE_MIN & age <= AGE_MAX & !is.na(pctile)]
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
  cat(sprintf("[%s] N=%s authors, observed=%s/%s (%.1f%%)\n",
              label, format(N, big.mark = ","), format(sum(mask_obs), big.mark = ","),
              format(N * T, big.mark = ","), 100 * mean(mask_obs)))

  fit <- fit_ife(Y, mask_obs, r = 1L)
  lam <- fit$lam; Ft <- fit$Ft
  if (mean(Ft) < 0) { Ft <- -Ft; lam <- -lam }
  cat(sprintf("[%s] R2=%.4f\n", label, fit$r2))
  list(authors = authors, ages = ages, lam = lam, Ft = Ft, r2 = fit$r2, N = N)
}

main <- function() {
  t6 <- fread(file.path(ROOT, "input", "t6_authors.csv"), colClasses = "character")
  t6[, t6_first_year := as.numeric(first_year)]
  cohort_map <- setNames(t6$t6_first_year, t6$author_id)

  pub <- fread(file.path(ROOT, "input", "publication_history.csv"),
               select = c("author_id", "eid", "publication_year"), colClasses = "character")
  pub[, publication_year := as.numeric(publication_year)]
  c10 <- fread(file.path(ROOT, "input", "citation_10y.csv"), colClasses = list(character = "eid"))
  pub <- merge(pub, c10, by = "eid", all.x = TRUE)

  cohort_full <- unique(pub, by = "eid")
  cohort_full <- cohort_full[cited_by_10y_complete == 1 & !is.na(cited_by_10y)]
  cohort_full[, pctile := frank(cited_by_10y, ties.method = "average") / .N, by = publication_year]
  pub <- merge(pub, cohort_full[, .(eid, pctile)], by = "eid", all.x = TRUE)

  pub[, first_year := min(publication_year, na.rm = TRUE), by = author_id]
  pub[, age := publication_year - first_year]
  pub[, t6_first_year := cohort_map[author_id]]

  pre <- pub[t6_first_year < CUTOFF_YEAR]
  post <- pub[t6_first_year >= CUTOFF_YEAR]

  fit_pre <- build_and_fit(pre, sprintf("pre-%d T6 entry", CUTOFF_YEAR))
  fit_post <- build_and_fit(post, sprintf("post-%d T6 entry", CUTOFF_YEAR))

  cat(sprintf("\n=== F_t (academic-age career curve): pre-%d vs post-%d T6 entry ===\n", CUTOFF_YEAR, CUTOFF_YEAR))
  for (i in seq_along(fit_pre$ages)) {
    cat(sprintf("  age %2d  pre=%+.3f  post=%+.3f\n", fit_pre$ages[i], fit_pre$Ft[i], fit_post$Ft[i]))
  }

  cat("\n=== lambda_i distribution comparison ===\n")
  for (nm in c("pre", "post")) {
    fit <- if (nm == "pre") fit_pre else fit_post
    l <- fit$lam
    cat(sprintf("%s: N=%s  mean=%.3f  median=%.3f  std=%.3f  p10=%.3f  p90=%.3f\n",
                nm, format(fit$N, big.mark = ","), mean(l), median(l), sd(l),
                quantile(l, 0.10), quantile(l, 0.90)))
  }

  ks <- ks.test(fit_pre$lam, fit_post$lam)
  tt <- t.test(fit_pre$lam, fit_post$lam)  # Welch's t-test is R's default (var.equal=FALSE)
  cat(sprintf("\nKS test: D=%.4f  p=%.4g\n", ks$statistic, ks$p.value))
  cat(sprintf("Welch t-test: t=%.3f  p=%.4g\n", tt$statistic, tt$p.value))

  fwrite(data.table(author_id = fit_pre$authors, lambda = fit_pre$lam),
         file.path(ROOT, "output", "lambda_pre1992.csv"))
  fwrite(data.table(author_id = fit_post$authors, lambda = fit_post$lam),
         file.path(ROOT, "output", "lambda_post1992.csv"))
  fwrite(data.table(age = fit_pre$ages, Ft_pre = fit_pre$Ft, Ft_post = fit_post$Ft),
         file.path(ROOT, "output", "Ft_cohort_compare.csv"))
  cat("\n[done] lambda_pre1992.csv, lambda_post1992.csv, Ft_cohort_compare.csv\n")
}

main()
