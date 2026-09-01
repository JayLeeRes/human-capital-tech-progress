#!/usr/bin/env Rscript
# Stage 04i -- Peer-environment design separating ability from resource
# congestion (draft.qmd Section VI.F "Peer Environment and the Origins of
# alpha_i").
#
# THE QUESTION: does an author's estimated ability (alpha_i) partly reflect
# who happened to be around them early in their career, rather than a pure,
# portable trait of the individual?
#
# THE PROBLEM WITH A NAIVE TEST: researchers choose their coauthors, and
# they tend to choose coauthors of similar ability to their own (sorting).
# So "authors with strong coauthors do well later" could just mean "able
# authors attract able coauthors AND are able themselves" -- no causal
# peer effect required.
#
# THE FIX: compare three increasingly strict versions of the same
# regression
#     pctile_late_i = b0 + b1 * pctile_early_i + b2 * peer_quality_i + [institution FE] + u_i
# where pctile_early/late are an author's mean cohort-percentile citation
# rank over academic ages 0-3 ("early career") and 15-25 ("late career"):
#
#   (1) peer_quality = mean alpha_i of the author's ACTUAL coauthors during
#       academic age 0-3. Endogenous: authors choose coauthors, so this
#       partly just re-measures the author's own ability (sorting).
#   (2) peer_quality = mean alpha_i of researchers who published from the
#       SAME institution in the same (or an adjacent) calendar year as the
#       focal author's early papers, but who the focal author NEVER
#       coauthored with, at any point in either career. Much less sorted
#       on ability (you don't choose who else happens to work down the
#       hall), so closer to a genuinely exogenous "exposure".
#   (3) same regression as (2), but with a fixed effect for the author's
#       own early-career institution. This throws away the comparison
#       "strong-peer-pool institutions vs weak-peer-pool institutions"
#       (which could just reflect that better institutions are better in
#       every way) and keeps only "within the same institution, did this
#       author happen to arrive during an unusually strong cohort".
#
# Reuses alpha_i from 05_additive_fe_pctile.R, both as (a) the peer-quality
# measure (mean alpha_i of coauthors / co-located researchers) and (b) via
# early/late cohort-percentile computed directly in this script.
#
# Input : ../input/publication_history.csv, ../input/citation_10y.csv,
#         ../output/alpha_additive_fe.csv
# Output: ../output/peer_effects_coauthors.csv, ../output/peer_effects_colocated.csv

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

ROOT <- dirname(dirname(this.path::this.path()))
EARLY_MAX_AGE <- 3L   # "early career" = academic age 0-3
LATE_MIN_AGE <- 15L; LATE_MAX_AGE <- 25L  # "late career" = academic age 15-25

# ---------------------------------------------------------------------------
# Step 1: rebuild the same author x paper panel used everywhere else in this
# pipeline (cohort-percentile of 10-year citations, academic age), plus the
# two extra raw columns (author_ids, own_affiliation) this script needs.
# ---------------------------------------------------------------------------
load_panel <- function() {
  pub <- fread(file.path(ROOT, "input", "publication_history.csv"),
               select = c("author_id", "eid", "publication_year", "author_ids", "own_affiliation"),
               colClasses = "character")
  pub[, publication_year := as.numeric(publication_year)]
  c10 <- fread(file.path(ROOT, "input", "citation_10y.csv"), colClasses = list(character = "eid"))
  pub <- merge(pub, c10, by = "eid", all.x = TRUE)

  # same cohort-percentile construction as 02_ife_cohort_percentile.R
  cohort <- unique(pub, by = "eid")
  cohort <- cohort[cited_by_10y_complete == 1 & !is.na(cited_by_10y)]
  cohort[, pctile := frank(cited_by_10y, ties.method = "average") / .N, by = publication_year]
  pub <- merge(pub, cohort[, .(eid, pctile)], by = "eid", all.x = TRUE)

  pub[, first_year := min(publication_year, na.rm = TRUE), by = author_id]
  pub[, age := publication_year - first_year]
  pub
}

# each author's own early- and late-career mean percentile (the "own_early"
# control and the outcome variable in every spec below)
own_pctiles <- function(pub) {
  early <- pub[age >= 0 & age <= EARLY_MAX_AGE & !is.na(pctile), .(early_pctile = mean(pctile)), by = author_id]
  late  <- pub[age >= LATE_MIN_AGE & age <= LATE_MAX_AGE & !is.na(pctile), .(late_pctile = mean(pctile)), by = author_id]
  merge(early, late, by = "author_id")
}

# ---------------------------------------------------------------------------
# Step 2: for every author, the set of everyone they EVER coauthored with
# (used in two places: to build spec (1)'s coauthor-quality measure, and to
# EXCLUDE actual coauthors from spec (2)-(3)'s "never coauthored" pool).
# ---------------------------------------------------------------------------
build_ever_coauthor_map <- function(pub) {
  # one row per (author_id, eid, coauthor_id) -- unnest the semicolon list
  pairs <- pub[, .(coauthor_id = unlist(strsplit(author_ids, ";"))), by = .(author_id, eid)]
  pairs <- pairs[coauthor_id != author_id]  # drop the author appearing as their own "coauthor"
  ever <- pairs[, .(coauthors = list(unique(coauthor_id))), by = author_id]
  list(pairs = pairs, map = setNames(ever$coauthors, ever$author_id))
}

# ---------------------------------------------------------------------------
# Spec (1): actual early-career coauthors
# ---------------------------------------------------------------------------
run_spec1_actual_coauthors <- function(pub, own, coauthor_pairs, alpha_map) {
  early_papers <- pub[age >= 0 & age <= EARLY_MAX_AGE, .(author_id, eid)]
  co <- merge(early_papers, coauthor_pairs, by = c("author_id", "eid"))  # keep only early-career papers' coauthors
  co[, coauthor_alpha := alpha_map[coauthor_id]]
  co <- co[!is.na(coauthor_alpha)]  # coauthor must also be in the alpha-estimation panel
  coauthor_quality <- co[, .(mean_coauthor_alpha = mean(coauthor_alpha),
                              n_coauthors = uniqueN(coauthor_id)), by = author_id]

  d <- merge(own, coauthor_quality, by = "author_id")
  cat(sprintf("[spec 1: actual coauthors] N=%d\n", nrow(d)))
  cat(sprintf("  cor(early, late)=%.3f  cor(late, coauthor_alpha)=%.3f  cor(early, coauthor_alpha)=%.3f (sorting)\n",
              cor(d$early_pctile, d$late_pctile), cor(d$late_pctile, d$mean_coauthor_alpha),
              cor(d$early_pctile, d$mean_coauthor_alpha)))

  m0 <- feols(late_pctile ~ early_pctile, data = d)                          # own early performance alone
  m1 <- feols(late_pctile ~ early_pctile + mean_coauthor_alpha, data = d)    # + coauthor quality
  cat("  (unconditional own-effect, no coauthor control):\n"); print(coeftable(m0))
  cat("  (with coauthor quality):\n"); print(coeftable(m1))
  d
}

# ---------------------------------------------------------------------------
# Specs (2)-(3): co-located, never-coauthored peers
# ---------------------------------------------------------------------------

# for a given (institution, year), who published from there that year, the
# year before, or the year after? (a +-1 year window, since "who was around
# at the same time" is inherently a bit fuzzy)
build_colocated_pool_lookup <- function(pub) {
  inst_year_authors <- pub[!is.na(own_affiliation) & own_affiliation != "",
                            .(author_id = unique(author_id)), by = .(own_affiliation, publication_year)]
  setkey(inst_year_authors, own_affiliation, publication_year)
  function(inst, yr) {
    hits <- inst_year_authors[.(inst, (yr - 1):(yr + 1)), on = .(own_affiliation, publication_year), nomatch = 0]
    unique(hits$author_id)
  }
}

# for each focal author: union the co-located pools across all their early
# (institution, year) combinations, then remove (a) themselves and (b)
# anyone they ever actually coauthored with -- what remains is "people who
# were physically nearby but never became collaborators"
build_peer_quality <- function(pub, ever_co_map, alpha_map, get_colocated_pool) {
  pub_early <- pub[age >= 0 & age <= EARLY_MAX_AGE & !is.na(own_affiliation) & own_affiliation != ""]
  early_inst_years <- unique(pub_early[, .(author_id, own_affiliation, publication_year)])
  focal_authors <- unique(early_inst_years$author_id)

  cat(sprintf("[spec 2-3: co-located peers] computing pools for %d focal authors (this takes a few minutes)...\n",
              length(focal_authors)))
  out <- vector("list", length(focal_authors))
  for (i in seq_along(focal_authors)) {
    a <- focal_authors[i]
    a_rows <- early_inst_years[author_id == a]
    pool <- unique(unlist(mapply(get_colocated_pool, a_rows$own_affiliation, a_rows$publication_year, SIMPLIFY = FALSE)))
    pool <- setdiff(pool, a)                       # (a) drop self
    actual_co <- ever_co_map[[a]]
    if (!is.null(actual_co)) pool <- setdiff(pool, actual_co)  # (b) drop actual coauthors
    pool_alpha <- alpha_map[pool]
    pool_alpha <- pool_alpha[!is.na(pool_alpha)]
    out[[i]] <- data.table(author_id = a,
                            mean_peer_alpha = if (length(pool_alpha) > 0) mean(pool_alpha) else NA_real_,
                            n_peers = length(pool_alpha))
  }
  peer_dt <- rbindlist(out)
  peer_dt[!is.na(mean_peer_alpha) & n_peers >= 1]
}

# each author's most common early-career institution (used as the fixed
# effect in spec 3, and to check how many focal authors share it)
build_primary_institution <- function(pub) {
  pub_early <- pub[age >= 0 & age <= EARLY_MAX_AGE & !is.na(own_affiliation) & own_affiliation != ""]
  by_inst <- pub_early[, .N, by = .(author_id, own_affiliation)][order(author_id, -N)]
  by_inst[, .SD[1], by = author_id][, .(author_id, primary_affiliation = own_affiliation)]
}

run_specs2_3_colocated <- function(pub, own, ever_co_map, alpha_map, d1) {
  get_colocated_pool <- build_colocated_pool_lookup(pub)
  peer_dt <- build_peer_quality(pub, ever_co_map, alpha_map, get_colocated_pool)
  primary_inst <- build_primary_institution(pub)

  d <- merge(own, peer_dt, by = "author_id")
  d <- merge(d, primary_inst, by = "author_id")
  inst_n <- d[, .N, by = primary_affiliation]
  d <- merge(d, inst_n, by = "primary_affiliation")
  d_fe <- d[N >= 2]  # keep only institutions with >=2 focal authors, so a within-institution
                      # fixed effect actually has variation to work with in spec (3)

  cat(sprintf("[spec 2-3] N=%d authors, %d institutions with >=2 authors\n",
              nrow(d_fe), uniqueN(d_fe$primary_affiliation)))
  cat(sprintf("  cor(early, peer_alpha)=%.3f (sorting; compare to %.3f for actual coauthors)\n",
              cor(d_fe$early_pctile, d_fe$mean_peer_alpha), cor(d1$early_pctile, d1$mean_coauthor_alpha)))

  m_pooled <- feols(late_pctile ~ early_pctile + mean_peer_alpha, data = d_fe)                         # spec (2)
  m_fe     <- feols(late_pctile ~ early_pctile + mean_peer_alpha | primary_affiliation, data = d_fe)   # spec (3)
  cat("\n  spec (2) pooled, no institution FE:\n"); print(m_pooled)
  cat("\n  spec (3) institution FE:\n"); print(m_fe)
  d_fe
}

main <- function() {
  pub <- load_panel()
  own <- own_pctiles(pub)

  alpha_dt <- fread(file.path(ROOT, "output", "alpha_additive_fe.csv"), colClasses = list(character = "author_id"))
  alpha_map <- setNames(alpha_dt$alpha, alpha_dt$author_id)

  ever <- build_ever_coauthor_map(pub)

  d1 <- run_spec1_actual_coauthors(pub, own, ever$pairs, alpha_map)
  fwrite(d1, file.path(ROOT, "output", "peer_effects_coauthors.csv"))

  cat("\n")
  d2_fe <- run_specs2_3_colocated(pub, own, ever$map, alpha_map, d1)
  fwrite(d2_fe, file.path(ROOT, "output", "peer_effects_colocated.csv"))

  cat("\n[done] peer_effects_coauthors.csv, peer_effects_colocated.csv\n")
}

main()
