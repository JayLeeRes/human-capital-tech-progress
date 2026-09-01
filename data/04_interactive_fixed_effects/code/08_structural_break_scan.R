#!/usr/bin/env Rscript
# Stage 04h -- Quandt-Andrews-style structural break scan for a discrete
# calendar-year ICT effect (draft.qmd Section VI.E "No Evidence of a
# Discrete Calendar-Year ICT Effect").
#
# THE IDEA: if ICT caused a one-time jump in how fast citation levels grew,
# there should be some calendar year where the trend visibly "kinks" --
# growing at one rate before that year, a different rate after. We don't
# know which year in advance, so we try EVERY plausible year as a
# candidate, measure how strong the kink would be at each one, and see
# which candidate (if any) stands out.
#
# THE REGRESSION AT ONE CANDIDATE YEAR tau:
#     yv = b0 + b1*(year - tau) + b2*(year - tau)*1[year >= tau] + error
# Before tau, 1[year>=tau]=0, so the fitted line has slope b1 ("slope_before").
# From tau onward, the fitted line has slope b1+b2 ("slope_after"). So b2 is
# exactly "how much the slope changes at tau", and testing b2=0 (an F-test
# comparing this model to one without the interaction term) tells us how
# surprising that particular kink is.
#
# THE SCAN: repeat that regression for every candidate tau, and note which
# one gives the biggest F-statistic -- that's the year the data itself
# points to as "most likely break point", the same logic as a
# Quandt-Andrews test for an unknown structural break.
#
# APPLIED TO:
#   (a) gamma_y, the calendar-year fixed effect from 03_ife_year_fe.R
#       (the primary test: a break here in the 1991-2000 ICT-diffusion
#       window would be direct evidence for the returns channel)
#   (b)-(d) three auxiliary series constructed from publication_history.csv
#       and year_fe.csv, used to check whether any break in (a) is better
#       explained by a change in Scopus's own indexing coverage than by
#       researcher behavior: number of author-age cells observed, number
#       of distinct journals, and mean coauthors per paper
#
# Input : ../output/year_fe.csv (from 03_ife_year_fe.R),
#         ../input/publication_history.csv
# Output: ../output/structural_breaks.csv (every candidate year, every series)

suppressPackageStartupMessages(library(data.table))

MIN_N <- 50L        # a candidate year needs at least this many underlying observations to be scanned
EDGE_TRIM <- 10L    # exclude candidates within this many years of either end of the sample
                    # (too close to the edge -> "before" or "after" has almost no data to fit a slope to)

# Fits the "slope kink at tau" regression described above for ONE candidate
# year, and returns the F-statistic for whether the slope actually changed.
fit_break_at_tau <- function(d, tau, weight_col = NULL) {
  d <- copy(d)
  d[, post := as.integer(year >= tau)]
  d[, year_c := year - tau]        # recentered so the intercept is "at tau", not at year 0
  d[, year_post := year_c * post]  # this is the b2 term: the extra slope kicking in from tau onward
  if (!is.null(weight_col)) {
    full <- lm(yv ~ year_c + year_post, data = d, weights = d[[weight_col]])
    null <- lm(yv ~ year_c, data = d, weights = d[[weight_col]])          # same model, minus the kink term
  } else {
    full <- lm(yv ~ year_c + year_post, data = d)
    null <- lm(yv ~ year_c, data = d)
  }
  f_stat <- anova(null, full)$F[2]   # F-test of H0: no slope change (b2 = 0)
  b <- coef(full)
  list(f_stat = f_stat, slope_before = b["year_c"], slope_after = b["year_c"] + b["year_post"])
}

# Runs fit_break_at_tau() for every plausible candidate year and ranks them
# by F-statistic, largest (= most convincing break) first.
scan_breaks <- function(dt, year_col, y_col, weight_col = NULL, label) {
  d <- copy(dt)
  setnames(d, c(year_col, y_col), c("year", "yv"))
  yr_range <- range(d$year)
  candidates <- d$year[d$year > yr_range[1] + EDGE_TRIM & d$year < yr_range[2] - EDGE_TRIM]

  res <- data.table(series = character(), tau = integer(), f_stat = numeric(),
                     slope_before = numeric(), slope_after = numeric())
  for (tau in candidates) {
    fit <- fit_break_at_tau(d, tau, weight_col)
    res <- rbind(res, data.table(series = label, tau = tau, f_stat = fit$f_stat,
                                  slope_before = fit$slope_before, slope_after = fit$slope_after))
  }
  setorder(res, -f_stat)
  res
}

main <- function() {
  ROOT <- dirname(dirname(this.path::this.path()))

  # ---- (a) gamma_y: the calendar-year citation-level fixed effect --------
  # year_fe.csv comes from 03_ife_year_fe.R; n_obs there is how many
  # author-age cells fell in that calendar year, used both to drop
  # sparse/noisy years and as a regression weight (a year built from more
  # data should count for more).
  yf <- fread(file.path(ROOT, "output", "year_fe.csv"))
  yf <- yf[n_obs >= MIN_N]
  res_gamma <- scan_breaks(yf, "year", "gamma_y", "n_obs", "gamma_y (citation level)")
  cat(sprintf("[gamma_y] scanned %d candidate years (%d-%d, n_obs>=%d)\n",
              nrow(res_gamma), min(yf$year), max(yf$year), MIN_N))
  cat("Top 5 candidates (ranked by F-statistic):\n"); print(res_gamma[1:5])
  cat("\nCandidates specifically within the 1991-2000 ICT-diffusion window:\n")
  print(res_gamma[tau %in% 1991:2000][order(tau)])

  # ---- (b) does the break in (a) just track how much data each year has? -
  # If gamma_y's break coincides with a break in raw data volume, that's a
  # sign the "break" might be a data-coverage artifact rather than a real
  # change in citation behavior.
  yf2 <- copy(yf); yf2[, log_n_obs := log(n_obs)]
  res_nobs <- scan_breaks(yf2, "year", "log_n_obs", label = "log(author-age cells observed)")
  cat("\n[log author-age cells observed] top candidate:\n"); print(res_nobs[1])

  # ---- (c)-(d) two more coverage/behavior checks from the raw corpus -----
  pub <- fread(file.path(ROOT, "input", "publication_history.csv"),
               select = c("publication_year", "source_id", "author_count"))
  pub[, publication_year := as.numeric(publication_year)]
  yearly <- pub[publication_year >= 1960 & publication_year <= 2016, .(
    n_journals = uniqueN(source_id),               # how many distinct journals contributed a paper this year
    mean_authors = mean(author_count, na.rm = TRUE) # mean coauthors per paper this year
  ), by = publication_year][order(publication_year)]
  yearly[, log_journals := log(n_journals)]

  res_journals <- scan_breaks(yearly, "publication_year", "log_journals", label = "log(distinct journals)")
  res_authors  <- scan_breaks(yearly, "publication_year", "mean_authors", label = "mean coauthors/paper")

  cat("\n[log distinct journals] top candidate:\n"); print(res_journals[1])
  cat("\n[mean coauthors/paper] top candidate:\n"); print(res_authors[1])
  cat("\n[mean coauthors/paper] candidates in the 1993-1997 window:\n")
  print(res_authors[tau %in% 1993:1997][order(tau)])

  all_res <- rbind(res_gamma, res_nobs, res_journals, res_authors)
  fwrite(all_res, file.path(ROOT, "output", "structural_breaks.csv"))
  cat("\n[done] structural_breaks.csv\n")
}

main()
