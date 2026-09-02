#!/usr/bin/env Rscript
# Stage 04k -- Does the congestion result survive measuring peers at the time
# of exposure rather than over their whole career?
#
# THE PROBLEM: 09_peer_effects.R measures a peer's quality by their
# career-average alpha. That is the right measure only if alpha is a fixed
# trait. 10_alpha_reliability.R shows it is not: early- and late-career
# alpha correlate about 0.49 once measurement error is removed. So the
# career-average measure answers "how good did this person turn out to be
# over 40 years?" when the design needs "how strong was this person at the
# moment we were in the same building."
#
# THE FIX: alpha_i is, by construction, the average over an author's cells of
#     s_it = Y_it - gamma_t - global_mean
# (relative standing at academic age t, net of the common age profile). The
# contemporaneous analogue simply averages s_it over only those cells whose
# CALENDAR year falls inside the co-location window, rather than over the
# peer's whole career.
#
# Everything else -- the co-located, never-coauthored peer pool, the outcome,
# the controls, the institution fixed effect -- is held identical to
# 09_peer_effects.R, so the two coefficients are directly comparable.
#
# Input : ../input/publication_history.csv, ../input/citation_10y.csv
# Output: ../output/peer_effects_contemporaneous.csv

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

ROOT <- dirname(dirname(this.path::this.path()))
AGE_MIN <- 0L; AGE_MAX <- 30L
EARLY_MAX_AGE <- 3L
LATE_MIN_AGE <- 15L; LATE_MAX_AGE <- 25L
WINDOW <- 1L        # co-location window: same calendar year +/- 1

main <- function() {
  pub <- fread(file.path(ROOT, "input", "publication_history.csv"),
               select = c("author_id", "eid", "publication_year", "author_ids", "own_affiliation"),
               colClasses = "character")
  pub[, publication_year := as.numeric(publication_year)]
  c10 <- fread(file.path(ROOT, "input", "citation_10y.csv"), colClasses = list(character = "eid"))
  pub <- merge(pub, c10, by = "eid", all.x = TRUE)
  cohort <- unique(pub, by = "eid")
  cohort <- cohort[cited_by_10y_complete == 1 & !is.na(cited_by_10y)]
  cohort[, pctile := frank(cited_by_10y, ties.method = "average") / .N, by = publication_year]
  pub <- merge(pub, cohort[, .(eid, pctile)], by = "eid", all.x = TRUE)
  pub[, first_year := min(publication_year, na.rm = TRUE), by = author_id]
  pub[, age := publication_year - first_year]

  # --- fit the age profile once, then form per-cell relative standing s_it --
  cell <- pub[age >= AGE_MIN & age <= AGE_MAX & !is.na(pctile),
              .(pctile = mean(pctile)), by = .(author_id, age)]
  cell <- cell[author_id %in% cell[, .N, by = author_id][N >= 2, author_id]]
  au <- sort(unique(cell$author_id)); ages <- AGE_MIN:AGE_MAX
  ai <- setNames(seq_along(au), au); ti <- setNames(seq_along(ages), ages)
  Y <- matrix(NA_real_, length(au), length(ages))
  Y[cbind(ai[cell$author_id], ti[as.character(cell$age)])] <- cell$pctile
  m <- !is.na(Y); gm <- mean(Y[m]); a <- rep(0, nrow(Y)); g <- rep(0, ncol(Y))
  for (i in 1:200) {
    r <- Y - outer(a, g, "+") - gm; r[!m] <- NA
    ra <- rowMeans(r, na.rm = TRUE); ra[is.na(ra)] <- 0; a <- a + ra
    r <- Y - outer(a, g, "+") - gm; r[!m] <- NA
    ca <- colMeans(r, na.rm = TRUE); ca[is.na(ca)] <- 0; g <- g + ca
    if (max(abs(ra)) < 1e-10 && max(abs(ca)) < 1e-10) break
  }
  alpha_map <- setNames(a, au)
  gamma <- setNames(g, as.character(ages))

  # s_it, tagged with the calendar year the cell falls in
  fy <- unique(pub[, .(author_id, first_year)])
  s <- merge(cell, fy, by = "author_id")
  s[, s_it := pctile - gamma[as.character(age)] - gm]
  s[, cyear := first_year + age]
  s <- s[, .(author_id, cyear, s_it)]
  setkey(s, author_id, cyear)

  # --- co-located, never-coauthored peer pool (identical to 09) ------------
  # every (author, ever-coauthor) pair, kept as a table so peers can be
  # excluded by anti-join rather than a per-row lookup (there are ~1.4M
  # focal-peer-year rows; a row-wise check here dominates the runtime)
  ever_co <- unique(pub[, .(cid = unlist(strsplit(author_ids, ";"))),
                        by = .(author_id, eid)][cid != author_id, .(author_id, peer = cid)])

  inst_year <- pub[!is.na(own_affiliation) & own_affiliation != "",
                   .(peer = unique(author_id)), by = .(own_affiliation, publication_year)]
  pe <- pub[age >= 0 & age <= EARLY_MAX_AGE & !is.na(own_affiliation) & own_affiliation != ""]
  focal_iy <- unique(pe[, .(author_id, own_affiliation, publication_year)])

  # expand each focal institution-year over the +/- WINDOW, then join to find
  # everyone else publishing from that institution in that window
  enc <- focal_iy[, .(author_id, own_affiliation,
                      yr = publication_year,
                      publication_year = publication_year + rep(-WINDOW:WINDOW, each = .N))]
  enc <- merge(enc, inst_year, by = c("own_affiliation", "publication_year"), allow.cartesian = TRUE)
  enc <- enc[peer != author_id]
  enc <- enc[!ever_co, on = .(author_id, peer)]             # drop actual coauthors
  cat(sprintf("[encounters] %s focal-peer-year rows over %d focal authors\n",
              format(nrow(enc), big.mark = ","), uniqueN(enc$author_id)))

  # --- contemporaneous peer standing: average s over the window ------------
  win <- unique(enc[, .(author_id, peer, yr)])
  win <- win[, .(cyear = yr + (-WINDOW:WINDOW)), by = .(author_id, peer, yr)]
  win <- merge(win, s, by.x = c("peer", "cyear"), by.y = c("author_id", "cyear"))
  peer_now <- win[, .(s_peer = mean(s_it)), by = .(author_id, peer, yr)][
    , .(mean_peer_now = mean(s_peer), n_peers = uniqueN(peer)), by = author_id]

  # --- career-average measure on the same rows, for a like-for-like check --
  peer_career <- unique(enc[, .(author_id, peer)])
  peer_career[, pa := alpha_map[peer]]
  peer_career <- peer_career[!is.na(pa)][, .(mean_peer_career = mean(pa)), by = author_id]

  # --- outcome, control, institution --------------------------------------
  early <- pub[age >= 0 & age <= EARLY_MAX_AGE & !is.na(pctile), .(ep = mean(pctile)), by = author_id]
  late  <- pub[age >= LATE_MIN_AGE & age <= LATE_MAX_AGE & !is.na(pctile), .(lp = mean(pctile)), by = author_id]
  inst  <- pe[, .N, by = .(author_id, own_affiliation)][order(author_id, -N)][
    , .SD[1], by = author_id][, .(author_id, inst = own_affiliation)]

  d <- Reduce(function(x, y) merge(x, y, by = "author_id"),
              list(early, late, peer_now, peer_career, inst))
  d <- merge(d, d[, .N, by = inst], by = "inst")[N >= 2]
  cat(sprintf("[sample] N=%d authors, %d institutions\n\n", nrow(d), uniqueN(d$inst)))

  m_career <- feols(lp ~ ep + mean_peer_career | inst, data = d)
  m_now    <- feols(lp ~ ep + mean_peer_now    | inst, data = d)
  cat("--- peer quality = career-average alpha (as in 09_peer_effects.R) ---\n")
  print(coeftable(m_career))
  cat("\n--- peer quality = contemporaneous standing during co-location ---\n")
  print(coeftable(m_now))

  cat(sprintf("\ncor(career-average, contemporaneous) = %.3f\n",
              cor(d$mean_peer_career, d$mean_peer_now)))
  cat(sprintf("SD: career-average %.4f, contemporaneous %.4f\n",
              sd(d$mean_peer_career), sd(d$mean_peer_now)))
  d[, `:=`(c_dm = mean_peer_career - mean(mean_peer_career),
           n_dm = mean_peer_now - mean(mean_peer_now)), by = inst]
  cat(sprintf("1 within-institution SD effect on late percentile: career %+.2f pp, contemporaneous %+.2f pp\n",
              100 * coef(m_career)["mean_peer_career"] * sd(d$c_dm),
              100 * coef(m_now)["mean_peer_now"] * sd(d$n_dm)))

  fwrite(d, file.path(ROOT, "output", "peer_effects_contemporaneous.csv"))
  cat("\n[done] peer_effects_contemporaneous.csv\n")
}

main()
