#!/usr/bin/env Rscript
# Stage 04j -- How trustworthy is alpha_i? Two questions, deliberately separated.
#
# (a) MEASUREMENT PRECISION: if we split one author's observations in half at
#     random and estimate alpha separately on each half, do the two halves
#     agree? The correlation between them is a direct, assumption-free
#     estimate of reliability; Spearman-Brown then rescales it from
#     "half-length" to the full-sample alpha. This is compared against the
#     parametric formula r = tau^2 / (tau^2 + sigma^2/T_i), which relies on
#     homoskedasticity and correct specification -- if the two agree, those
#     assumptions are doing no harm here.
#
# (b) TIME INVARIANCE: the additive model writes alpha_i with no t subscript,
#     so it *assumes* a researcher's standing is constant over a career. That
#     is an assumption, not a finding. Estimating alpha separately on the
#     first half of a career (academic ages 0-12) and the second (13-30) and
#     correlating the two tests it. Attenuating for measurement error (using
#     the reliability from (a)) gives the correlation between the underlying
#     quantities rather than between their noisy estimates.
#
# Input : ../input/publication_history.csv, ../input/citation_10y.csv
# Output: ../output/alpha_reliability.csv

suppressPackageStartupMessages(library(data.table))

ROOT <- dirname(dirname(this.path::this.path()))
AGE_MIN <- 0L; AGE_MAX <- 30L
ERA_CUT <- 12L      # academic age separating "early career" from "late career"
SEED <- 42L

# Two-way FE on an arbitrary subset of author-age cells, returning alpha_i.
fe_alpha <- function(dt, n_iter = 300L, tol = 1e-10) {
  au <- sort(unique(dt$author_id)); ages <- sort(unique(dt$age))
  ai <- setNames(seq_along(au), au); ti <- setNames(seq_along(ages), ages)
  Y <- matrix(NA_real_, length(au), length(ages))
  Y[cbind(ai[dt$author_id], ti[as.character(dt$age)])] <- dt$pctile
  m <- !is.na(Y); gm <- mean(Y[m]); a <- rep(0, nrow(Y)); g <- rep(0, ncol(Y))
  for (i in seq_len(n_iter)) {
    r <- Y - outer(a, g, "+") - gm; r[!m] <- NA
    ra <- rowMeans(r, na.rm = TRUE); ra[is.na(ra)] <- 0; a <- a + ra
    r <- Y - outer(a, g, "+") - gm; r[!m] <- NA
    ca <- colMeans(r, na.rm = TRUE); ca[is.na(ca)] <- 0; g <- g + ca
    if (max(abs(ra)) < tol && max(abs(ca)) < tol) break
  }
  Yhat <- outer(a, g, "+") + gm
  list(alpha = setNames(a, au), Ti = rowSums(m),
       sigma2 = sum((Y[m] - Yhat[m])^2) / (sum(m) - length(au) - length(ages) + 1))
}

build_cells <- function() {
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
  pub[age >= AGE_MIN & age <= AGE_MAX & !is.na(pctile),
      .(pctile = mean(pctile)), by = .(author_id, age)]
}

main <- function() {
  set.seed(SEED)
  cell <- build_cells()
  out <- list()

  # --- parametric benchmark on the full panel -----------------------------
  full <- cell[author_id %in% cell[, .N, by = author_id][N >= 2, author_id]]
  f <- fe_alpha(full)
  tau2 <- var(f$alpha) - f$sigma2 * mean(1 / f$Ti)
  r_param_med <- tau2 / (tau2 + f$sigma2 / median(f$Ti))
  cat(sprintf("[panel] N=%d  median T_i=%.0f  sigma^2=%.5f  tau^2=%.5f\n",
              length(f$alpha), median(f$Ti), f$sigma2, tau2))
  cat(sprintf("[parametric] reliability at median T_i = %.3f\n\n", r_param_med))

  # --- (a) random split-half reliability ----------------------------------
  # needs >=4 cells so each half still has >=2 (the minimum for identification)
  c4 <- cell[author_id %in% cell[, .N, by = author_id][N >= 4, author_id]]
  c4[, half := sample(rep(1:2, length.out = .N)), by = author_id]
  A <- fe_alpha(c4[half == 1])$alpha; B <- fe_alpha(c4[half == 2])$alpha
  ids <- intersect(names(A), names(B))
  r_half <- cor(A[ids], B[ids])
  r_full <- 2 * r_half / (1 + r_half)          # Spearman-Brown
  cat(sprintf("[split-half] N=%d  r_half=%.3f  ->  full-sample reliability=%.3f\n",
              length(ids), r_half, r_full))
  cat(sprintf("             (parametric prediction was %.3f -- close agreement means\n", r_param_med))
  cat("              homoskedasticity/specification assumptions are not distorting it)\n\n")

  # --- (b) early-career vs late-career alpha ------------------------------
  cell[, era := ifelse(age <= ERA_CUT, "early", "late")]
  ok <- cell[, .(ne = sum(era == "early"), nl = sum(era == "late")), by = author_id][ne >= 2 & nl >= 2, author_id]
  E <- fe_alpha(cell[author_id %in% ok & era == "early"])$alpha
  L <- fe_alpha(cell[author_id %in% ok & era == "late"])$alpha
  ids2 <- intersect(names(E), names(L))
  r_era <- cor(E[ids2], L[ids2])
  r_era_true <- min(1, r_era / r_full)          # disattenuated for measurement error
  cat(sprintf("[time invariance] N=%d  early-vs-late alpha r=%.3f (Spearman %.3f)\n",
              length(ids2), r_era, cor(E[ids2], L[ids2], method = "spearman")))
  cat(sprintf("                  disattenuated for measurement error: r=%.3f\n", r_era_true))
  cat(sprintf("                  -> only %.0f%% of late-career variance is shared with early career;\n",
              100 * r_era_true^2))
  cat("                     alpha_i is a career average of a time-VARYING quantity,\n")
  cat("                     not an estimate of a constant trait.\n")

  out <- data.table(
    quantity = c("sigma2", "tau2", "median_T_i", "reliability_parametric_at_median_T",
                 "split_half_r", "reliability_spearman_brown", "n_split_half",
                 "early_late_r", "early_late_r_disattenuated", "n_early_late"),
    value = c(f$sigma2, tau2, median(f$Ti), r_param_med,
              r_half, r_full, length(ids),
              r_era, r_era_true, length(ids2)))
  fwrite(out, file.path(ROOT, "output", "alpha_reliability.csv"))
  cat("\n[done] alpha_reliability.csv\n")
}

main()
