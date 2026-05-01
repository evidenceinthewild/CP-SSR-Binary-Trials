# R/analysis.R
# ----------------------------------------------------------------------------
# Helper functions for the manuscript:
#   "Conditional Power Promising Zone Sample Size Re-estimation Inflates Type I
#    Error in Single-Arm Binary Trials"
#
# This file is sourced from CP_PromisingZone_SSR_TypeIError.qmd. It defines:
#   - critical-count helpers: x_crit_fn (Bayesian), x_crit_z_fn (frequentist)
#   - posterior + predictive probability: post_prob, pred_prob_fn
#   - conditional power: cp_fn
#   - SSR rules: ssr_cp_fn (Mehta-Pocock), ssr_bayes_fn (PP-zone)
#   - simulators: sim_trial_cp, sim_trial_bayes, run_cp_scenario, run_bayes_scenario
#   - exact OC computation: cp_oc_exact, bayes_oc_exact, cp_oc_post_final
#   - utilities: fixed_n_t1e, discrete_quantile
# ----------------------------------------------------------------------------

# Required packages (the manuscript loads these in its setup chunk):
#   tidyverse, glue
# Loaded but used only by the manuscript:
#   gt, scales, patchwork, latex2exp

# run_cp_scenario() and run_bayes_scenario() use the global `N_SIM` as
# their default replication count. The manuscript sets `N_SIM` explicitly
# in its setup chunk; for standalone use of analysis.R we provide a
# fallback so that calls without an explicit `n_sim` argument do not
# raise "object 'N_SIM' not found". The manuscript's later assignment
# overrides this default.
if (!exists("N_SIM", envir = .GlobalEnv, inherits = FALSE)) {
  assign("N_SIM", 10000L, envir = .GlobalEnv)
}

# =============================================================================
# Critical-count helpers
# =============================================================================

#' Min integer x such that the Beta(ap+x, bp+n-x) posterior puts >= thresh
#' probability above p0 (used by the Bayesian final analysis and pred_prob_fn).
x_crit_fn <- function(n, p0, ap, bp, thresh) {
  for (x in 0L:n) {
    if (1 - pbeta(p0, ap + x, bp + n - x) >= thresh) return(x)
  }
  return(n + 1L)
}

#' Min integer x such that the one-sided one-proportion z-test rejects
#' H0: p <= p0 at level `alpha` (normal approximation).
#' MUST match the final-analysis decision rule used in sim_trial_cp().
x_crit_z_fn <- function(n, p0, alpha) {
  thresh <- n * p0 + qnorm(1 - alpha) * sqrt(n * p0 * (1 - p0))
  as.integer(floor(thresh) + 1L)
}

#' Posterior P(p > p0 | x successes in n, Beta(ap, bp) prior).
post_prob <- function(x, n, p0, ap, bp) {
  1 - pbeta(p0, ap + x, bp + n - x)
}

#' Beta-Binomial predictive probability of eventual success at the final
#' analysis (declared when posterior >= gamma_fin), given the interim count.
pred_prob_fn <- function(x_int, n_int, n_final, p0, ap, bp, gamma_fin) {
  a_post <- ap + x_int
  b_post <- bp + n_int - x_int
  n_rem  <- n_final - n_int
  xc     <- x_crit_fn(n_final, p0, ap, bp, gamma_fin)
  if (xc > n_final) return(0)
  x_need <- max(0L, xc - x_int)
  if (x_need > n_rem) return(0)
  sum(vapply(x_need:n_rem, function(xr) {
    exp(
      lchoose(n_rem, xr) +
      lbeta(a_post + xr, b_post + n_rem - xr) -
      lbeta(a_post, b_post)
    )
  }, numeric(1)))
}

# =============================================================================
# Conditional power and SSR rules
# =============================================================================

#' Frequentist conditional power: projects post-interim data under p1 against
#' the same z-test critical count that the final analysis uses.
cp_fn <- function(x_int, n_int, n_final, p0, p1, alpha) {
  xc     <- x_crit_z_fn(n_final, p0, alpha)
  x_need <- max(0L, xc - x_int)
  n_rem  <- n_final - n_int
  if (x_need > n_rem) return(0)
  1 - pbinom(x_need - 1L, n_rem, p1)
}

#' Mehta-Pocock CP-promising-zone SSR: search upward from n_init for the
#' smallest n_final at which CP rises back to cp_upper, capped at n_max.
#' Returns an integer scalar so downstream `vapply(..., integer(1))` calls
#' work even when callers pass plain numerics for n_init / n_max.
ssr_cp_fn <- function(x_int, n_int, p0, p1, alpha,
                      cp_upper, n_init, n_max) {
  n_init <- as.integer(n_init)
  n_max  <- as.integer(n_max)
  for (n_cand in seq.int(n_init, n_max)) {
    cp <- cp_fn(x_int, n_int, n_cand, p0, p1, alpha)
    if (cp >= cp_upper) return(as.integer(n_cand))
  }
  return(n_max)
}

#' Bayesian PP-promising-zone SSR: search upward from n_init for the smallest
#' n_final at which predictive probability rises to pp_upper, capped at n_max.
#' Returns an integer scalar (see ssr_cp_fn note).
ssr_bayes_fn <- function(x_int, n_int, p0, ap, bp, gamma_fin,
                         pp_upper, n_init, n_max) {
  n_init <- as.integer(n_init)
  n_max  <- as.integer(n_max)
  for (n_cand in seq.int(n_init, n_max)) {
    pp <- pred_prob_fn(x_int, n_int, n_cand, p0, ap, bp, gamma_fin)
    if (pp >= pp_upper) return(as.integer(n_cand))
  }
  return(n_max)
}

# =============================================================================
# Trial simulators (single-replicate)
# =============================================================================

#' Simulate one single-arm CP-promising-zone trial. Returns success/n/zone/n_cap.
#' Final analysis: one-sided one-proportion z-test at level alpha.
sim_trial_cp <- function(p_true, cp_lower, cp_upper = 0.80,
                         cp_futility = 0.10, n_init = 84L,
                         n_int = 65L, n_max = 200L,
                         p0 = 0.23, p1 = 0.35, alpha = 0.05,
                         ap = 0.5, bp = 0.5) {
  x_int  <- rbinom(1L, n_int, p_true)
  cp_obs <- cp_fn(x_int, n_int, n_init, p0, p1, alpha)

  # Futility
  if (cp_obs < cp_futility)
    return(c(success = 0L, n = n_int, zone = 0L, n_cap = 0L))

  # Favorable: no extension
  if (cp_obs >= cp_upper) {
    n_final <- n_init; zone <- 3L
  # Promising: Mehta-Pocock SSR
  } else if (cp_obs >= cp_lower) {
    n_final <- ssr_cp_fn(x_int, n_int, p0, p1, alpha,
                         cp_upper, n_init, n_max)
    zone    <- 2L
  # Unfavorable: continue at planned N
  } else {
    n_final <- n_init; zone <- 1L
  }

  x_rem   <- rbinom(1L, n_final - n_int, p_true)
  x_total <- x_int + x_rem
  phat    <- x_total / n_final
  z       <- (phat - p0) / sqrt(p0 * (1 - p0) / n_final)
  pval    <- 1 - pnorm(z)

  c(success = as.integer(pval < alpha),
    n       = n_final,
    zone    = zone,
    n_cap   = as.integer(n_final == n_max))
}

#' Run n_sim independent CP-zone trials and return a 1-row OC tibble.
run_cp_scenario <- function(p_true, cp_lower, n_sim = N_SIM, ...) {
  one_trial <- function() sim_trial_cp(p_true, cp_lower, ...)
  mat <- t(replicate(n_sim, one_trial()))
  tibble::tibble(
    p_true   = p_true,
    cp_lower = cp_lower,
    power    = mean(mat[, "success"]),
    exp_n    = mean(mat[, "n"]),
    pr_fut   = mean(mat[, "zone"] == 0),
    pr_unf   = mean(mat[, "zone"] == 1),
    pr_prom  = mean(mat[, "zone"] == 2),
    pr_fav   = mean(mat[, "zone"] == 3),
    pr_cap   = mean(mat[, "n_cap"]),
    n_p10    = quantile(mat[, "n"], 0.10),
    n_p50    = quantile(mat[, "n"], 0.50),
    n_p90    = quantile(mat[, "n"], 0.90)
  )
}

#' Simulate one single-arm Bayesian PP-promising-zone trial.
#' Final analysis: posterior P(p > p0 | x_total, n_final) >= gamma_fin.
sim_trial_bayes <- function(p_true, gamma_fin, gamma_int = 0.99,
                             pp_fut = 0.05, pp_upper = 0.50,
                             n_init = 84L, n_int = 65L, n_max = 200L,
                             p0 = 0.23, ap = 0.5, bp = 0.5) {
  x_int <- rbinom(1L, n_int, p_true)

  # Interim efficacy stop
  post_int <- post_prob(x_int, n_int, p0, ap, bp)
  if (post_int >= gamma_int) {
    return(c(success = 1L, n = n_int, zone = 3L,
             pr_eff = 1L, stop_fut = 0L, n_cap = 0L))
  }

  pp <- pred_prob_fn(x_int, n_int, n_init, p0, ap, bp, gamma_fin)

  if (pp <= pp_fut)
    return(c(success = 0L, n = n_int, zone = 0L,
             pr_eff = 0L, stop_fut = 1L, n_cap = 0L))

  if (pp >= pp_upper) {
    n_final <- n_init; zone <- 2L
  } else {
    n_final <- ssr_bayes_fn(x_int, n_int, p0, ap, bp, gamma_fin,
                            pp_upper, n_init, n_max)
    zone    <- 1L
  }

  x_rem    <- rbinom(1L, n_final - n_int, p_true)
  x_total  <- x_int + x_rem
  post_fin <- post_prob(x_total, n_final, p0, ap, bp)

  c(success  = as.integer(post_fin >= gamma_fin),
    n        = n_final,
    zone     = zone,
    pr_eff   = 0L,
    stop_fut = 0L,
    n_cap    = as.integer(n_final == n_max))
}

#' Run n_sim independent Bayesian-PP trials and return a 1-row OC tibble.
run_bayes_scenario <- function(p_true, gamma_fin, n_sim = N_SIM, ...) {
  one_trial <- function() sim_trial_bayes(p_true, gamma_fin, ...)
  mat <- t(replicate(n_sim, one_trial()))
  tibble::tibble(
    p_true    = p_true,
    gamma_fin = gamma_fin,
    power     = mean(mat[, "success"]),
    exp_n     = mean(mat[, "n"]),
    pr_eff    = mean(mat[, "pr_eff"]),
    pr_fut    = mean(mat[, "stop_fut"]),
    pr_prom   = mean(mat[, "zone"] == 1),
    pr_fav    = mean(mat[, "zone"] == 2),
    pr_cap    = mean(mat[, "n_cap"]),
    n_p10     = quantile(mat[, "n"], 0.10),
    n_p50     = quantile(mat[, "n"], 0.50),
    n_p90     = quantile(mat[, "n"], 0.90)
  )
}

# =============================================================================
# Exact operating characteristics (no Monte Carlo)
# =============================================================================

#' Fixed-N z-test Type I error: rejection probability of the one-sided
#' one-proportion z-test at sample size n under H0: p = p0.
fixed_n_t1e <- function(n, p0 = 0.23, alpha = 0.05) {
  xc <- x_crit_z_fn(n, p0, alpha)
  if (xc > n) return(0)
  1 - pbinom(xc - 1L, n, p0)
}

#' Quantile of a discrete distribution given (values, probs).
discrete_quantile <- function(values, probs, q) {
  ord <- order(values)
  v   <- values[ord]
  cp  <- cumsum(probs[ord])
  v[which(cp >= q)[1]]
}

#' Exact OC for the CP-promising-zone design via enumeration of x_int.
#' Returns the same columns as run_cp_scenario but deterministically.
cp_oc_exact <- function(p_true, cp_lower,
                        cp_upper = 0.80, cp_futility = 0.10,
                        n_init = 84L, n_int = 65L, n_max = 200L,
                        p0 = 0.23, p1 = 0.35, alpha = 0.05,
                        ap = 0.5, bp = 0.5) {
  n_init <- as.integer(n_init)
  n_int  <- as.integer(n_int)
  n_max  <- as.integer(n_max)
  x_grid    <- 0L:n_int
  pmf_x_int <- dbinom(x_grid, n_int, p_true)
  cp_obs    <- vapply(x_grid, function(x) cp_fn(x, n_int, n_init, p0, p1, alpha),
                      numeric(1))

  zone <- integer(length(x_grid))
  zone[cp_obs <  cp_futility]                       <- 0L
  zone[cp_obs >= cp_upper]                          <- 3L
  zone[cp_obs >= cp_lower & cp_obs < cp_upper]      <- 2L
  zone[cp_obs >= cp_futility & cp_obs < cp_lower]   <- 1L

  n_final <- ifelse(zone == 0L, n_int,
              ifelse(zone == 2L,
                vapply(x_grid, function(x) ssr_cp_fn(x, n_int, p0, p1, alpha,
                                                      cp_upper, n_init, n_max),
                       integer(1)),
                n_init))

  p_succ <- vapply(seq_along(x_grid), function(i) {
    if (zone[i] == 0L) return(0)
    xc     <- x_crit_z_fn(n_final[i], p0, alpha)
    x_need <- max(0L, xc - x_grid[i])
    n_rem  <- n_final[i] - n_int
    if (x_need > n_rem) return(0)
    1 - pbinom(x_need - 1L, n_rem, p_true)
  }, numeric(1))

  tibble::tibble(
    p_true   = p_true,
    cp_lower = cp_lower,
    power    = sum(pmf_x_int * p_succ),
    exp_n    = sum(pmf_x_int * n_final),
    pr_fut   = sum(pmf_x_int * (zone == 0L)),
    pr_unf   = sum(pmf_x_int * (zone == 1L)),
    pr_prom  = sum(pmf_x_int * (zone == 2L)),
    pr_fav   = sum(pmf_x_int * (zone == 3L)),
    pr_cap   = sum(pmf_x_int * (n_final == n_max)),
    n_p10    = discrete_quantile(n_final, pmf_x_int, 0.10),
    n_p50    = discrete_quantile(n_final, pmf_x_int, 0.50),
    n_p90    = discrete_quantile(n_final, pmf_x_int, 0.90)
  )
}

#' Exact OC for the Bayesian PP design via enumeration of x_int.
bayes_oc_exact <- function(p_true, gamma_fin, gamma_int = 0.99,
                            pp_fut = 0.05, pp_upper = 0.50,
                            n_init = 84L, n_int = 65L, n_max = 200L,
                            p0 = 0.23, ap = 0.5, bp = 0.5) {
  n_init <- as.integer(n_init)
  n_int  <- as.integer(n_int)
  n_max  <- as.integer(n_max)
  x_grid    <- 0L:n_int
  pmf_x_int <- dbinom(x_grid, n_int, p_true)

  post_int <- vapply(x_grid, function(x) post_prob(x, n_int, p0, ap, bp),
                     numeric(1))
  pp <- vapply(x_grid, function(x)
                pred_prob_fn(x, n_int, n_init, p0, ap, bp, gamma_fin),
                numeric(1))

  # Zones: 3=interim efficacy, 0=futility, 2=favorable, 1=promising
  zone <- ifelse(post_int >= gamma_int, 3L,
           ifelse(pp <= pp_fut,         0L,
            ifelse(pp >= pp_upper,      2L, 1L)))

  n_final <- ifelse(zone %in% c(0L, 3L), n_int,
              ifelse(zone == 2L, n_init,
                vapply(x_grid, function(x) ssr_bayes_fn(x, n_int, p0, ap, bp,
                                                         gamma_fin, pp_upper,
                                                         n_init, n_max),
                       integer(1))))

  p_succ <- vapply(seq_along(x_grid), function(i) {
    if (zone[i] == 0L) return(0)
    if (zone[i] == 3L) return(1)
    xc     <- x_crit_fn(n_final[i], p0, ap, bp, gamma_fin)
    x_need <- max(0L, xc - x_grid[i])
    n_rem  <- n_final[i] - n_int
    if (x_need > n_rem) return(0)
    1 - pbinom(x_need - 1L, n_rem, p_true)
  }, numeric(1))

  tibble::tibble(
    p_true    = p_true,
    gamma_fin = gamma_fin,
    power     = sum(pmf_x_int * p_succ),
    exp_n     = sum(pmf_x_int * n_final),
    pr_eff    = sum(pmf_x_int * (zone == 3L)),
    pr_fut    = sum(pmf_x_int * (zone == 0L)),
    pr_prom   = sum(pmf_x_int * (zone == 1L)),
    pr_fav    = sum(pmf_x_int * (zone == 2L)),
    pr_cap    = sum(pmf_x_int * (n_final == n_max)),
    n_p10     = discrete_quantile(n_final, pmf_x_int, 0.10),
    n_p50     = discrete_quantile(n_final, pmf_x_int, 0.50),
    n_p90     = discrete_quantile(n_final, pmf_x_int, 0.90)
  )
}

#' Symmetric variant: CP-zone SSR with a Bayesian posterior final analysis.
#' Same interim/zone/SSR logic as cp_oc_exact, but the final rule declares
#' success when post_prob(x_total, n_final) >= gamma_fin (not z-test).
cp_oc_post_final <- function(p_true, cp_lower, gamma_fin = 0.955,
                              cp_upper = 0.80, cp_futility = 0.10,
                              n_init = 84L, n_int = 65L, n_max = 200L,
                              p0 = 0.23, p1 = 0.35, alpha = 0.05,
                              ap = 0.5, bp = 0.5) {
  n_init <- as.integer(n_init)
  n_int  <- as.integer(n_int)
  n_max  <- as.integer(n_max)
  x_grid    <- 0L:n_int
  pmf_x_int <- dbinom(x_grid, n_int, p_true)
  cp_obs    <- vapply(x_grid, function(x) cp_fn(x, n_int, n_init, p0, p1, alpha),
                      numeric(1))

  zone <- integer(length(x_grid))
  zone[cp_obs <  cp_futility]                     <- 0L
  zone[cp_obs >= cp_upper]                        <- 3L
  zone[cp_obs >= cp_lower & cp_obs < cp_upper]    <- 2L
  zone[cp_obs >= cp_futility & cp_obs < cp_lower] <- 1L

  n_final <- ifelse(zone == 0L, n_int,
              ifelse(zone == 2L,
                vapply(x_grid, function(x) ssr_cp_fn(x, n_int, p0, p1, alpha,
                                                      cp_upper, n_init, n_max),
                       integer(1)),
                n_init))

  p_succ <- vapply(seq_along(x_grid), function(i) {
    if (zone[i] == 0L) return(0)
    xc <- x_crit_fn(n_final[i], p0, ap, bp, gamma_fin)
    if (xc > n_final[i]) return(0)
    x_need <- max(0L, xc - x_grid[i])
    n_rem  <- n_final[i] - n_int
    if (x_need > n_rem) return(0)
    1 - pbinom(x_need - 1L, n_rem, p_true)
  }, numeric(1))

  tibble::tibble(
    p_true    = p_true,
    cp_lower  = cp_lower,
    gamma_fin = gamma_fin,
    power     = sum(pmf_x_int * p_succ),
    exp_n     = sum(pmf_x_int * n_final),
    pr_fut    = sum(pmf_x_int * (zone == 0L)),
    pr_unf    = sum(pmf_x_int * (zone == 1L)),
    pr_prom   = sum(pmf_x_int * (zone == 2L)),
    pr_fav    = sum(pmf_x_int * (zone == 3L)),
    pr_cap    = sum(pmf_x_int * (n_final == n_max))
  )
}
