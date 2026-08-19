# R/analysis.R
# ----------------------------------------------------------------------------
# Helper functions for the manuscript:
#   "Final-Analysis Discreteness Drives Type I Error Miscalibration in
#    Single-Arm Binary Promising-Zone Designs"
#
# This file is sourced from CP_PromisingZone_SSR_TypeIError.qmd. It defines:
#   - critical-count helpers: x_crit_fn, x_crit_z_fn, x_crit_exact_fn
#   - posterior + predictive probability: post_prob, pred_prob_fn
#   - conditional power: cp_fn
#   - SSR rules: ssr_cp_fn (Mehta-Pocock), ssr_bayes_fn (PP-zone)
#   - simulators: sim_trial_cp, sim_trial_bayes, run_cp_scenario, run_bayes_scenario
#   - exact OC computation: cp_oc_exact, bayes_oc_exact, cp_oc_post_final,
#     cp_oc_exact_final, bayes_oc_exact_final, bayes_oc_ztest_final,
#     comb_oc_exact
#   - utilities: fixed_n_t1e, fixed_n_t1e_exact, least_favorable_null
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
#'
#' Formulation: SCORE (null-variance) z-test. The standard error is computed
#' under H0 as sqrt(n * p0 * (1 - p0)), NOT from the sample proportion
#' (Wald variance). The test rejects when
#'     (x/n - p0) / sqrt(p0 (1 - p0) / n) > z_{1-alpha},
#' i.e. when x > n*p0 + z_{1-alpha} * sqrt(n p0 (1-p0)).
#'
#' `correct = TRUE` applies the usual +1/2 continuity correction, which
#' shifts the critical count upward and makes the test more conservative.
#' The manuscript's primary analysis uses correct = FALSE.
x_crit_z_fn <- function(n, p0, alpha, correct = FALSE) {
  thresh <- n * p0 + qnorm(1 - alpha) * sqrt(n * p0 * (1 - p0))
  if (correct) thresh <- thresh + 0.5
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

#' Projection proportion used by conditional power.
#'
#' `cp_mode = "design"` projects under the designed alternative p1 (the
#' convention used in the originally submitted manuscript).
#'
#' `cp_mode = "estimated"` projects under the interim-estimated response
#' proportion p-hat = x_int / n_int. This is the convention of Mehta and
#' Pocock (2011, Section 3), whose promising-zone construction evaluates CP
#' at the *observed* interim effect rather than the design alternative.
#'
#' Optional shrinkage pulls p-hat toward `shrink_to` with weight `shrink_w`
#' (shrink_w = 0 is raw p-hat, shrink_w = 1 recovers the design mode when
#' shrink_to = p1). Optional bounding clips the projection proportion into
#' [p_min, p_max], which stabilises CP at extreme interim counts.
cp_proj_p <- function(x_int, n_int, p1,
                      cp_mode = c("design", "estimated"),
                      shrink_to = p1, shrink_w = 0,
                      p_min = NULL, p_max = NULL) {
  cp_mode <- match.arg(cp_mode)
  p_proj <- if (cp_mode == "design") p1 else x_int / n_int
  if (cp_mode == "estimated" && shrink_w > 0) {
    p_proj <- (1 - shrink_w) * p_proj + shrink_w * shrink_to
  }
  if (!is.null(p_min)) p_proj <- max(p_proj, p_min)
  if (!is.null(p_max)) p_proj <- min(p_proj, p_max)
  p_proj
}

#' Frequentist conditional power against the same z-test critical count that
#' the final analysis uses. `alpha` here is the FINAL-ANALYSIS level: CP must
#' be computed against the boundary the trial will actually be judged by.
#'
#' See cp_proj_p() for the `cp_mode` / shrinkage / bounding arguments. With
#' cp_mode = "design" (the default) this is identical to the function used in
#' the originally submitted manuscript.
cp_fn <- function(x_int, n_int, n_final, p0, p1, alpha,
                  cp_mode = c("design", "estimated"),
                  shrink_to = p1, shrink_w = 0,
                  p_min = NULL, p_max = NULL,
                  correct = FALSE) {
  cp_mode <- match.arg(cp_mode)
  p_proj <- cp_proj_p(x_int, n_int, p1, cp_mode, shrink_to, shrink_w,
                      p_min, p_max)
  xc     <- x_crit_z_fn(n_final, p0, alpha, correct = correct)
  x_need <- max(0L, xc - x_int)
  n_rem  <- n_final - n_int
  if (x_need > n_rem) return(0)
  if (n_rem <= 0L) return(as.numeric(x_need <= 0L))
  1 - pbinom(x_need - 1L, n_rem, p_proj)
}

#' Mehta-Pocock CP-promising-zone SSR: search upward from n_init for the
#' smallest n_final at which CP rises back to cp_upper, capped at n_max.
#' Returns an integer scalar so downstream `vapply(..., integer(1))` calls
#' work even when callers pass plain numerics for n_init / n_max.
ssr_cp_fn <- function(x_int, n_int, p0, p1, alpha,
                      cp_upper, n_init, n_max,
                      cp_mode = c("design", "estimated"),
                      shrink_to = p1, shrink_w = 0,
                      p_min = NULL, p_max = NULL,
                      correct = FALSE) {
  cp_mode <- match.arg(cp_mode)
  n_init  <- as.integer(n_init)
  n_max   <- as.integer(n_max)
  for (n_cand in seq.int(n_init, n_max)) {
    cp <- cp_fn(x_int, n_int, n_cand, p0, p1, alpha,
                cp_mode, shrink_to, shrink_w, p_min, p_max, correct)
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
                         ap = 0.5, bp = 0.5,
                         cp_mode = c("design", "estimated"),
                         shrink_to = p1, shrink_w = 0,
                         p_min = NULL, p_max = NULL,
                         correct = FALSE) {
  cp_mode <- match.arg(cp_mode)
  x_int <- rbinom(1L, n_int, p_true)
  cp_obs <- cp_fn(
    x_int, n_int, n_init, p0, p1, alpha,
    cp_mode = cp_mode, shrink_to = shrink_to, shrink_w = shrink_w,
    p_min = p_min, p_max = p_max, correct = correct
  )

  # Futility
  if (cp_obs < cp_futility)
    return(c(success = 0L, n = n_int, zone = 0L, n_cap = 0L))

  # Favorable: no extension
  if (cp_obs >= cp_upper) {
    n_final <- n_init
    zone <- 3L
  # Promising: Mehta-Pocock SSR
  } else if (cp_obs >= cp_lower) {
    n_final <- ssr_cp_fn(
      x_int, n_int, p0, p1, alpha, cp_upper, n_init, n_max,
      cp_mode = cp_mode, shrink_to = shrink_to, shrink_w = shrink_w,
      p_min = p_min, p_max = p_max, correct = correct
    )
    zone <- 2L
  # Unfavorable: continue at planned N
  } else {
    n_final <- n_init
    zone <- 1L
  }

  x_rem <- rbinom(1L, n_final - n_int, p_true)
  x_total <- x_int + x_rem
  x_crit <- x_crit_z_fn(n_final, p0, alpha, correct)

  c(success = as.integer(x_total >= x_crit),
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
#' `cp_mode` selects the CP projection effect (see cp_proj_p): "design"
#' reproduces the originally submitted results, "estimated" implements the
#' Mehta-Pocock (2011) convention of projecting under the observed interim
#' proportion.
#'
#' `alpha_fin` is the final-analysis level used to set the z-test critical
#' count. It defaults to `alpha` (the two coincide in the manuscript's primary
#' analysis). Because conditional power and the SSR target are evaluated
#' against that same boundary, sweeping `alpha_fin` jointly recalibrates the
#' final rejection rule and the interim adaptation mapping.
#'
#' Futility can be disabled with `cp_futility = -Inf`.
cp_oc_exact <- function(p_true, cp_lower,
                        cp_upper = 0.80, cp_futility = 0.10,
                        n_init = 84L, n_int = 65L, n_max = 200L,
                        p0 = 0.23, p1 = 0.35, alpha = 0.05,
                        ap = 0.5, bp = 0.5,
                        cp_mode = c("design", "estimated"),
                        shrink_to = p1, shrink_w = 0,
                        p_min = NULL, p_max = NULL,
                        alpha_fin = alpha, correct = FALSE) {
  cp_mode <- match.arg(cp_mode)
  n_init <- as.integer(n_init)
  n_int  <- as.integer(n_int)
  n_max  <- as.integer(n_max)
  x_grid    <- 0L:n_int
  pmf_x_int <- dbinom(x_grid, n_int, p_true)
  cp_obs    <- vapply(x_grid, function(x)
                 cp_fn(x, n_int, n_init, p0, p1, alpha_fin,
                       cp_mode, shrink_to, shrink_w, p_min, p_max, correct),
                 numeric(1))

  # Zone boundaries must be ordered for the assignments below to partition the
  # support. Without this guard, cp_lower < cp_futility would let the later
  # "promising" assignment silently overwrite an earlier "futility" one.
  stopifnot(cp_futility <= cp_lower, cp_lower <= cp_upper)

  zone <- integer(length(x_grid))
  zone[cp_obs <  cp_futility]                       <- 0L
  zone[cp_obs >= cp_upper]                          <- 3L
  zone[cp_obs >= cp_lower & cp_obs < cp_upper]      <- 2L
  zone[cp_obs >= cp_futility & cp_obs < cp_lower]   <- 1L

  n_final <- ifelse(zone == 0L, n_int,
              ifelse(zone == 2L,
                vapply(x_grid, function(x)
                  ssr_cp_fn(x, n_int, p0, p1, alpha_fin,
                            cp_upper, n_init, n_max,
                            cp_mode, shrink_to, shrink_w, p_min, p_max,
                            correct),
                  integer(1)),
                n_init))

  p_succ <- vapply(seq_along(x_grid), function(i) {
    if (zone[i] == 0L) return(0)
    xc     <- x_crit_z_fn(n_final[i], p0, alpha_fin, correct = correct)
    x_need <- max(0L, xc - x_grid[i])
    n_rem  <- n_final[i] - n_int
    if (x_need > n_rem) return(0)
    if (n_rem <= 0L) return(as.numeric(x_need <= 0L))
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

#' Interim-count-to-zone/sample-size mapping for Bayesian PP SSR.
#'
#' The final posterior threshold is part of the predictive-probability
#' calculation, so changing `gamma_fin` changes this mapping as well as the
#' final critical count. Returning the mapping explicitly lets final-rule
#' contrasts hold it fixed.
bayes_pp_mapping <- function(gamma_fin, gamma_int = 0.99,
                             pp_fut = 0.05, pp_upper = 0.50,
                             n_init = 84L, n_int = 65L, n_max = 200L,
                             p0 = 0.23, ap = 0.5, bp = 0.5) {
  n_init <- as.integer(n_init)
  n_int  <- as.integer(n_int)
  n_max  <- as.integer(n_max)
  x_grid <- 0L:n_int

  post_int <- vapply(x_grid, function(x) post_prob(x, n_int, p0, ap, bp),
                     numeric(1))
  pp <- vapply(x_grid, function(x)
                pred_prob_fn(x, n_int, n_init, p0, ap, bp, gamma_fin),
                numeric(1))

  # Zones: 3=interim efficacy, 0=futility, 2=favorable, 1=promising
  zone <- ifelse(post_int >= gamma_int, 3L,
           ifelse(pp <= pp_fut,         0L,
            ifelse(pp >= pp_upper,      2L, 1L)))

  ssr_n <- vapply(x_grid, function(x)
                    ssr_bayes_fn(x, n_int, p0, ap, bp, gamma_fin, pp_upper,
                                 n_init, n_max),
                  integer(1))
  n_final <- as.integer(ifelse(zone %in% c(0L, 3L), n_int,
                        ifelse(zone == 2L, n_init, ssr_n)))

  tibble::tibble(
    x_int = x_grid,
    post_int = post_int,
    pp = pp,
    zone = as.integer(zone),
    n_final = n_final
  )
}

# Evaluate a posterior-PP mapping under either the posterior or exact-binomial
# final boundary. Interim efficacy and futility decisions remain frozen.
.bayes_oc_from_mapping <- function(p_true, mapping,
                                   final_rule = c("posterior", "exact"),
                                   gamma_fin, p0, ap, bp, alpha,
                                   n_int, n_max) {
  final_rule <- match.arg(final_rule)
  x_grid <- mapping$x_int
  pmf_x_int <- dbinom(x_grid, n_int, p_true)

  p_succ <- vapply(seq_along(x_grid), function(i) {
    if (mapping$zone[i] == 0L) return(0)
    if (mapping$zone[i] == 3L) return(1)
    xc <- if (final_rule == "posterior") {
      x_crit_fn(mapping$n_final[i], p0, ap, bp, gamma_fin)
    } else {
      x_crit_exact_fn(mapping$n_final[i], p0, alpha)
    }
    if (xc > mapping$n_final[i]) return(0)
    x_need <- max(0L, xc - x_grid[i])
    n_rem  <- mapping$n_final[i] - n_int
    if (n_rem <= 0L) return(as.numeric(x_need <= 0L))
    if (x_need > n_rem) return(0)
    1 - pbinom(x_need - 1L, n_rem, p_true)
  }, numeric(1))

  tibble::tibble(
    p_true    = p_true,
    gamma_fin = gamma_fin,
    power     = sum(pmf_x_int * p_succ),
    exp_n     = sum(pmf_x_int * mapping$n_final),
    pr_eff    = sum(pmf_x_int * (mapping$zone == 3L)),
    pr_fut    = sum(pmf_x_int * (mapping$zone == 0L)),
    pr_prom   = sum(pmf_x_int * (mapping$zone == 1L)),
    pr_fav    = sum(pmf_x_int * (mapping$zone == 2L)),
    pr_cap    = sum(pmf_x_int * (mapping$n_final == n_max)),
    n_p10     = discrete_quantile(mapping$n_final, pmf_x_int, 0.10),
    n_p50     = discrete_quantile(mapping$n_final, pmf_x_int, 0.50),
    n_p90     = discrete_quantile(mapping$n_final, pmf_x_int, 0.90)
  )
}

#' Exact OC for the Bayesian PP design via enumeration of x_int.
bayes_oc_exact <- function(p_true, gamma_fin, gamma_int = 0.99,
                           pp_fut = 0.05, pp_upper = 0.50,
                           n_init = 84L, n_int = 65L, n_max = 200L,
                           p0 = 0.23, ap = 0.5, bp = 0.5) {
  mapping <- bayes_pp_mapping(
    gamma_fin, gamma_int, pp_fut, pp_upper,
    n_init, n_int, n_max, p0, ap, bp
  )
  .bayes_oc_from_mapping(
    p_true, mapping, "posterior", gamma_fin, p0, ap, bp,
    alpha = NA_real_, n_int, n_max
  )
}

#' Bayesian PP adaptation mapping with an exact-binomial final analysis.
#'
#' The mapping is constructed using the posterior threshold `gamma_fin`, then
#' frozen while the final boundary is replaced by the level-`alpha` exact
#' binomial critical count. This directly tests whether posterior and exact
#' control are numerically equivalent in the realized design.
bayes_oc_exact_final <- function(p_true, gamma_fin, gamma_int = 0.99,
                                 pp_fut = 0.05, pp_upper = 0.50,
                                 n_init = 84L, n_int = 65L, n_max = 200L,
                                 p0 = 0.23, alpha = 0.05,
                                 ap = 0.5, bp = 0.5) {
  mapping <- bayes_pp_mapping(
    gamma_fin, gamma_int, pp_fut, pp_upper,
    n_init, n_int, n_max, p0, ap, bp
  )
  .bayes_oc_from_mapping(
    p_true, mapping, "exact", gamma_fin, p0, ap, bp,
    alpha, n_int, n_max
  )
}

#' Symmetric variant: CP-zone SSR with a Bayesian posterior final analysis.
#' Same interim/zone/SSR logic as cp_oc_exact, but the final rule declares
#' success when post_prob(x_total, n_final) >= gamma_fin (not z-test).
cp_oc_post_final <- function(p_true, cp_lower, gamma_fin = 0.955,
                              cp_upper = 0.80, cp_futility = 0.10,
                              n_init = 84L, n_int = 65L, n_max = 200L,
                              p0 = 0.23, p1 = 0.35, alpha = 0.05,
                              ap = 0.5, bp = 0.5,
                              cp_mode = c("design", "estimated"),
                              shrink_to = p1, shrink_w = 0,
                              p_min = NULL, p_max = NULL,
                              correct = FALSE) {
  cp_mode <- match.arg(cp_mode)
  n_init <- as.integer(n_init)
  n_int  <- as.integer(n_int)
  n_max  <- as.integer(n_max)
  x_grid    <- 0L:n_int
  pmf_x_int <- dbinom(x_grid, n_int, p_true)
  cp_obs    <- vapply(x_grid, function(x)
                 cp_fn(x, n_int, n_init, p0, p1, alpha,
                       cp_mode, shrink_to, shrink_w, p_min, p_max, correct),
                 numeric(1))

  # See cp_oc_exact(): boundaries must be ordered for these assignments to
  # partition the support.
  stopifnot(cp_futility <= cp_lower, cp_lower <= cp_upper)

  zone <- integer(length(x_grid))
  zone[cp_obs <  cp_futility]                     <- 0L
  zone[cp_obs >= cp_upper]                        <- 3L
  zone[cp_obs >= cp_lower & cp_obs < cp_upper]    <- 2L
  zone[cp_obs >= cp_futility & cp_obs < cp_lower] <- 1L

  n_final <- ifelse(zone == 0L, n_int,
              ifelse(zone == 2L,
                vapply(x_grid, function(x)
                  ssr_cp_fn(x, n_int, p0, p1, alpha,
                            cp_upper, n_init, n_max,
                            cp_mode, shrink_to, shrink_w, p_min, p_max,
                            correct),
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

# =============================================================================
# Revision additions (SBR-26-068, R1 + R2)
# =============================================================================

#' Critical count under any of the four final decision rules, on a common
#' scale so they can be compared directly.
#'
#'   "z"         uncorrected score-form normal approximation
#'   "z_cc"      same, with the +1/2 continuity correction
#'   "exact"     exact binomial (Clopper-Pearson-equivalent)
#'   "posterior" Beta posterior threshold at gamma_fin
#'
#' Every one of these designs rejects iff x_total >= crit_count(...). Two rules
#' with the same critical count at a given n are THE SAME TEST at that n,
#' whatever their motivation. This function exists so that claim can be
#' checked rather than assumed.
crit_count <- function(n, p0 = 0.23, alpha = 0.05,
                       rule = c("z", "z_cc", "exact", "posterior"),
                       ap = 0.5, bp = 0.5, gamma_fin = 0.955) {
  rule <- match.arg(rule)
  switch(rule,
    z         = x_crit_z_fn(n, p0, alpha, correct = FALSE),
    z_cc      = x_crit_z_fn(n, p0, alpha, correct = TRUE),
    exact     = x_crit_exact_fn(n, p0, alpha),
    posterior = x_crit_fn(n, p0, ap, bp, gamma_fin)
  )
}

#' Critical count for every rule across a range of n, with the attained
#' fixed-n Type I error of each. Use this to establish where rules coincide
#' and where they genuinely differ.
crit_count_table <- function(n_seq, p0 = 0.23, alpha = 0.05,
                             ap = 0.5, bp = 0.5, gamma_fin = 0.955) {
  rules <- c("z", "z_cc", "exact", "posterior")
  out <- lapply(rules, function(r) {
    xc <- vapply(n_seq, function(n)
      crit_count(n, p0, alpha, r, ap, bp, gamma_fin), integer(1))
    t1e <- vapply(seq_along(n_seq), function(i) {
      if (xc[i] > n_seq[i]) return(0)
      1 - pbinom(xc[i] - 1L, n_seq[i], p0)
    }, numeric(1))
    tibble::tibble(n = n_seq, rule = r, crit = xc, t1e = t1e)
  })
  do.call(rbind, out)
}

#' Fixed-n operating characteristics with NO interim analysis and NO sample
#' size re-estimation, under any final rule.
#'
#' This supplies the no-adaptation reference needed to distinguish
#' adaptation-versus-fixed contrasts from CP-versus-PP mapping contrasts.
fixed_n_oc <- function(p_true, n, p0 = 0.23, alpha = 0.05,
                       rule = c("z", "z_cc", "exact", "posterior"),
                       ap = 0.5, bp = 0.5, gamma_fin = 0.955) {
  rule <- match.arg(rule)
  n  <- as.integer(n)
  xc <- crit_count(n, p0, alpha, rule, ap, bp, gamma_fin)
  pw <- if (xc > n) 0 else 1 - pbinom(xc - 1L, n, p_true)
  tibble::tibble(
    p_true = p_true, n = n, rule = rule, crit = xc,
    power = pw, exp_n = as.numeric(n)
  )
}

#' Critical count for the one-sided EXACT binomial test of H0: p <= p0.
#'
#' Returns the smallest integer x with P(X >= x | X ~ Bin(n, p0)) <= alpha.
#' This is the Clopper-Pearson-equivalent decision rule: the exact test
#' rejects precisely when the Clopper-Pearson lower confidence bound at level
#' 1 - alpha exceeds p0. Unlike the normal-approximation z-test it is
#' guaranteed non-anticonservative at every n, at the cost of being strictly
#' conservative wherever the binomial atom straddles alpha.
x_crit_exact_fn <- function(n, p0, alpha) {
  n <- as.integer(n)
  for (x in 0L:n) {
    if (1 - pbinom(x - 1L, n, p0) <= alpha) return(as.integer(x))
  }
  as.integer(n + 1L)
}

#' Fixed-N Type I error of the exact binomial test (companion to fixed_n_t1e,
#' which reports the same quantity for the z-test).
fixed_n_t1e_exact <- function(n, p0 = 0.23, alpha = 0.05) {
  xc <- x_crit_exact_fn(n, p0, alpha)
  if (xc > n) return(0)
  1 - pbinom(xc - 1L, n, p0)
}

#' EXACT-FINAL COMPARATOR (Referee 2, minor point 2).
#'
#' CP-promising-zone SSR paired with an exact binomial (Clopper-Pearson-type)
#' final analysis. The interim zone logic and the SSR mapping are identical to
#' cp_oc_exact() and cp_oc_post_final(): conditional power is evaluated against
#' the z-test boundary, so the adaptation is held fixed and only the final rule
#' changes. That makes the three CP-row designs directly comparable.
#'
#' The Discussion of the submitted manuscript asserted that an exact final test
#' "would push Type I error uniformly below alpha". This function evaluates
#' that assertion rather than assuming it.
cp_oc_exact_final <- function(p_true, cp_lower,
                              cp_upper = 0.80, cp_futility = 0.10,
                              n_init = 84L, n_int = 65L, n_max = 200L,
                              p0 = 0.23, p1 = 0.35, alpha = 0.05,
                              ap = 0.5, bp = 0.5,
                              cp_mode = c("design", "estimated"),
                              shrink_to = p1, shrink_w = 0,
                              p_min = NULL, p_max = NULL,
                              alpha_fin = alpha, correct = FALSE) {
  cp_mode <- match.arg(cp_mode)
  n_init  <- as.integer(n_init)
  n_int   <- as.integer(n_int)
  n_max   <- as.integer(n_max)

  x_grid    <- 0L:n_int
  pmf_x_int <- dbinom(x_grid, n_int, p_true)
  cp_obs    <- vapply(x_grid, function(x)
                 cp_fn(x, n_int, n_init, p0, p1, alpha_fin,
                       cp_mode, shrink_to, shrink_w, p_min, p_max, correct),
                 numeric(1))

  # See cp_oc_exact(): boundaries must be ordered for these assignments to
  # partition the support.
  stopifnot(cp_futility <= cp_lower, cp_lower <= cp_upper)

  zone <- integer(length(x_grid))
  zone[cp_obs <  cp_futility]                     <- 0L
  zone[cp_obs >= cp_upper]                        <- 3L
  zone[cp_obs >= cp_lower & cp_obs < cp_upper]    <- 2L
  zone[cp_obs >= cp_futility & cp_obs < cp_lower] <- 1L

  n_final <- ifelse(zone == 0L, n_int,
              ifelse(zone == 2L,
                vapply(x_grid, function(x)
                  ssr_cp_fn(x, n_int, p0, p1, alpha_fin,
                            cp_upper, n_init, n_max,
                            cp_mode, shrink_to, shrink_w, p_min, p_max,
                            correct),
                  integer(1)),
                n_init))

  # FINAL RULE: exact binomial test at level alpha.
  p_succ <- vapply(seq_along(x_grid), function(i) {
    if (zone[i] == 0L) return(0)
    xc <- x_crit_exact_fn(n_final[i], p0, alpha_fin)
    if (xc > n_final[i]) return(0)
    x_need <- max(0L, xc - x_grid[i])
    n_rem  <- n_final[i] - n_int
    if (n_rem <= 0L) return(as.numeric(x_need <= 0L))
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
    pr_cap   = sum(pmf_x_int * (n_final == n_max))
  )
}

#' Beta-Binomial predictive probability that the final count reaches an
#' arbitrary critical count `xc`. Generalises pred_prob_fn(), which is the
#' special case xc = x_crit_fn(...) for a posterior-threshold final rule.
pred_prob_generic <- function(x_int, n_int, n_final, ap, bp, xc) {
  a_post <- ap + x_int
  b_post <- bp + n_int - x_int
  n_rem  <- n_final - n_int
  if (xc > n_final) return(0)
  x_need <- max(0L, xc - x_int)
  if (n_rem <= 0L) return(as.numeric(x_need <= 0L))
  if (x_need > n_rem) return(0)
  sum(vapply(x_need:n_rem, function(xr) {
    exp(
      lchoose(n_rem, xr) +
      lbeta(a_post + xr, b_post + n_rem - xr) -
      lbeta(a_post, b_post)
    )
  }, numeric(1)))
}

#' Predictive probability of eventual success when "success" is defined by the
#' frequentist z-test critical count rather than a posterior threshold. This is
#' the coherent PP rule to pair with a z-test final analysis: the SSR targets
#' the boundary the trial will actually be judged by.
pred_prob_z_fn <- function(x_int, n_int, n_final, p0, ap, bp, alpha,
                           correct = FALSE) {
  xc <- x_crit_z_fn(n_final, p0, alpha, correct = correct)
  pred_prob_generic(x_int, n_int, n_final, ap, bp, xc)
}

#' PP-promising-zone SSR targeting the z-test boundary.
ssr_bayes_z_fn <- function(x_int, n_int, p0, ap, bp, alpha,
                           pp_upper, n_init, n_max, correct = FALSE) {
  n_init <- as.integer(n_init)
  n_max  <- as.integer(n_max)
  for (n_cand in seq.int(n_init, n_max)) {
    pp <- pred_prob_z_fn(x_int, n_int, n_cand, p0, ap, bp, alpha, correct)
    if (pp >= pp_upper) return(as.integer(n_cand))
  }
  return(n_max)
}

#' PP-driven SSR paired with a z-test final analysis.
#'
#' Completes the four-way comparison begun by cp_oc_exact (CP SSR + z final),
#' cp_oc_post_final (CP SSR + posterior final) and bayes_oc_exact
#' (PP SSR + posterior final). Together these four functions compare frozen
#' adaptation mappings and final-analysis rules. A fixed-N row is needed to
#' estimate adaptation-versus-no-adaptation contrasts.
#'
#' `pp_target = "ztest"` (default) computes the predictive probability against
#' the z-test critical count, so the SSR rule and the final rule are coherent.
#' `pp_target = "posterior"` leaves the PP machinery exactly as in
#' bayes_oc_exact and swaps only the final rule, which isolates the final-rule
#' contrast in a descriptive frozen-mapping comparison.
#'
#' Interim efficacy stopping is disabled by setting `gamma_int` above 1
#' (the posterior probability is strictly < 1, so the branch never fires).
#' Futility is disabled with `pp_fut = -Inf`.
bayes_oc_ztest_final <- function(p_true, gamma_int = 0.99,
                                 pp_fut = 0.05, pp_upper = 0.50,
                                 n_init = 84L, n_int = 65L, n_max = 200L,
                                 p0 = 0.23, alpha = 0.05,
                                 ap = 0.5, bp = 0.5,
                                 gamma_fin = 0.955,
                                 pp_target = c("ztest", "posterior"),
                                 correct = FALSE) {
  pp_target <- match.arg(pp_target)
  n_init <- as.integer(n_init)
  n_int  <- as.integer(n_int)
  n_max  <- as.integer(n_max)

  x_grid    <- 0L:n_int
  pmf_x_int <- dbinom(x_grid, n_int, p_true)

  post_int <- vapply(x_grid, function(x) post_prob(x, n_int, p0, ap, bp),
                     numeric(1))

  pp <- if (pp_target == "ztest") {
    vapply(x_grid, function(x)
      pred_prob_z_fn(x, n_int, n_init, p0, ap, bp, alpha, correct), numeric(1))
  } else {
    vapply(x_grid, function(x)
      pred_prob_fn(x, n_int, n_init, p0, ap, bp, gamma_fin), numeric(1))
  }

  # Zones: 3 = interim efficacy, 0 = futility, 2 = favorable, 1 = promising
  zone <- ifelse(post_int >= gamma_int, 3L,
           ifelse(pp <= pp_fut,         0L,
            ifelse(pp >= pp_upper,      2L, 1L)))

  ssr_n <- if (pp_target == "ztest") {
    vapply(x_grid, function(x)
      ssr_bayes_z_fn(x, n_int, p0, ap, bp, alpha, pp_upper,
                     n_init, n_max, correct), integer(1))
  } else {
    vapply(x_grid, function(x)
      ssr_bayes_fn(x, n_int, p0, ap, bp, gamma_fin, pp_upper,
                   n_init, n_max), integer(1))
  }

  n_final <- ifelse(zone %in% c(0L, 3L), n_int,
              ifelse(zone == 2L, n_init, ssr_n))

  # FINAL RULE: one-sided one-proportion z-test (the swap under study).
  p_succ <- vapply(seq_along(x_grid), function(i) {
    if (zone[i] == 0L) return(0)
    if (zone[i] == 3L) return(1)   # interim efficacy stop declares success
    xc     <- x_crit_z_fn(n_final[i], p0, alpha, correct = correct)
    x_need <- max(0L, xc - x_grid[i])
    n_rem  <- n_final[i] - n_int
    if (n_rem <= 0L) return(as.numeric(x_need <= 0L))
    if (x_need > n_rem) return(0)
    1 - pbinom(x_need - 1L, n_rem, p_true)
  }, numeric(1))

  tibble::tibble(
    p_true    = p_true,
    pp_target = pp_target,
    power     = sum(pmf_x_int * p_succ),
    exp_n     = sum(pmf_x_int * n_final),
    pr_eff    = sum(pmf_x_int * (zone == 3L)),
    pr_fut    = sum(pmf_x_int * (zone == 0L)),
    pr_prom   = sum(pmf_x_int * (zone == 1L)),
    pr_fav    = sum(pmf_x_int * (zone == 2L)),
    pr_cap    = sum(pmf_x_int * (n_final == n_max))
  )
}

#' One-sided stagewise p-value for X ~ Bin(n, p0), valid under H0.
#'
#' "exact" returns P(X >= x), which is stochastically larger than U(0,1) on
#' the discrete support and hence conservative. "midp" returns
#' P(X > x) + 0.5 P(X = x), which is less conservative but not guaranteed
#' valid, so it can re-introduce inflation.
stagewise_p <- function(x, n, p0, type = c("exact", "midp")) {
  type <- match.arg(type)
  if (n <= 0L) return(1)
  if (type == "exact") {
    pbinom(x - 1L, n, p0, lower.tail = FALSE)
  } else {
    pbinom(x, n, p0, lower.tail = FALSE) + 0.5 * dbinom(x, n, p0)
  }
}

#' COMBINATION-TEST COMPARATOR (Referee 2, major point 2).
#'
#' Inverse-normal combination of exact binomial stagewise p-values with
#' weights fixed at the design stage:
#'
#'   Z_comb = w1 * qnorm(1 - p_1) + w2 * qnorm(1 - p_2),
#'   w1 = sqrt(n_int / N_init),  w2 = sqrt(1 - n_int / N_init),  w1^2 + w2^2 = 1
#'
#' reject iff Z_comb > z_{1-alpha}. Because the weights do not depend on the
#' re-estimated N, the Lehmacher-Wassmer / Cui-Hung-Wang invariance argument
#' applies to any pre-specified choice of n_final based only on stage-1 data.
#' Discreteness makes the exact binomial p-values stochastically larger than
#' uniform, so control is expected to be conservative; `pval = "midp"`
#' quantifies the middle ground.
#'
#' SSR logic is the CP promising-zone rule, so this is the direct comparator
#' for cp_oc_exact: identical adaptation, different final analysis.
comb_oc_exact <- function(p_true, cp_lower,
                          cp_upper = 0.80, cp_futility = 0.10,
                          n_init = 84L, n_int = 65L, n_max = 200L,
                          p0 = 0.23, p1 = 0.35, alpha = 0.05,
                          pval = c("exact", "midp"),
                          cp_mode = c("design", "estimated"),
                          shrink_to = p1, shrink_w = 0,
                          p_min = NULL, p_max = NULL,
                          alpha_fin = alpha, correct = FALSE) {
  pval    <- match.arg(pval)
  cp_mode <- match.arg(cp_mode)
  n_init  <- as.integer(n_init)
  n_int   <- as.integer(n_int)
  n_max   <- as.integer(n_max)

  # Design-stage weights: fixed BEFORE any interim data are seen.
  w1 <- sqrt(n_int / n_init)
  w2 <- sqrt(1 - n_int / n_init)
  z_crit <- qnorm(1 - alpha_fin)

  x_grid    <- 0L:n_int
  pmf_x_int <- dbinom(x_grid, n_int, p_true)

  cp_obs <- vapply(x_grid, function(x)
              cp_fn(x, n_int, n_init, p0, p1, alpha_fin,
                    cp_mode, shrink_to, shrink_w, p_min, p_max, correct),
              numeric(1))

  # See cp_oc_exact(): boundaries must be ordered for these assignments to
  # partition the support.
  stopifnot(cp_futility <= cp_lower, cp_lower <= cp_upper)

  zone <- integer(length(x_grid))
  zone[cp_obs <  cp_futility]                     <- 0L
  zone[cp_obs >= cp_upper]                        <- 3L
  zone[cp_obs >= cp_lower & cp_obs < cp_upper]    <- 2L
  zone[cp_obs >= cp_futility & cp_obs < cp_lower] <- 1L

  n_final <- ifelse(zone == 0L, n_int,
              ifelse(zone == 2L,
                vapply(x_grid, function(x)
                  ssr_cp_fn(x, n_int, p0, p1, alpha_fin,
                            cp_upper, n_init, n_max,
                            cp_mode, shrink_to, shrink_w, p_min, p_max,
                            correct),
                  integer(1)),
                n_init))

  # For each interim count, enumerate the stage-2 count exactly.
  p_succ <- vapply(seq_along(x_grid), function(i) {
    if (zone[i] == 0L) return(0)             # futility stop, no rejection
    n_rem <- n_final[i] - n_int
    if (n_rem <= 0L) return(0)
    p_1 <- stagewise_p(x_grid[i], n_int, p0, pval)
    z_1 <- qnorm(p_1, lower.tail = FALSE)
    xr  <- 0L:n_rem
    p_2 <- vapply(xr, function(x) stagewise_p(x, n_rem, p0, pval),
                  numeric(1))
    z_2 <- qnorm(p_2, lower.tail = FALSE)
    rej <- (w1 * z_1 + w2 * z_2) > z_crit
    sum(dbinom(xr, n_rem, p_true)[rej])
  }, numeric(1))

  tibble::tibble(
    p_true   = p_true,
    cp_lower = cp_lower,
    pval     = pval,
    power    = sum(pmf_x_int * p_succ),
    exp_n    = sum(pmf_x_int * n_final),
    pr_fut   = sum(pmf_x_int * (zone == 0L)),
    pr_unf   = sum(pmf_x_int * (zone == 1L)),
    pr_prom  = sum(pmf_x_int * (zone == 2L)),
    pr_fav   = sum(pmf_x_int * (zone == 3L)),
    pr_cap   = sum(pmf_x_int * (n_final == n_max))
  )
}

#' LEAST-FAVOURABLE NULL CHECK (Referee 2, minor point 3).
#'
#' Adaptation and discreteness can break monotonicity of the rejection
#' probability in p, so the largest Type I error over the composite null
#' p <= p0 need not occur at p0. Evaluates any OC function on a numerical grid
#' and reports the largest value on that grid; it does not claim an analytic
#' supremum over the continuum.
#'
#' `oc_fn` must accept `p_true` as its first argument and return a one-row
#' tibble containing a `power` column.
#'
#' The null boundary argument is named `p_null` (not `p0`) so that `p0` can be
#' passed through `...` to the OC function without a partial-matching clash.
#' The grid spans the FULL composite null [0, p_null] by default. An earlier
#' version truncated to [p_null - 0.15, p_null]; independent recomputation over
#' the full range confirmed the supremum is attained at p_null for every design
#' evaluated here, but the untruncated default removes the question entirely.
least_favorable_null <- function(oc_fn, p_null = 0.23, n_grid = 201L,
                                 lower = NULL, ...) {
  lower <- if (is.null(lower)) 0 else lower
  p_seq <- seq(lower, p_null, length.out = n_grid)
  t1e   <- vapply(p_seq, function(p) oc_fn(p_true = p, ...)$power, numeric(1))
  i_max <- which.max(t1e)
  list(
    grid        = tibble::tibble(p = p_seq, t1e = t1e),
    p_worst     = p_seq[i_max],
    t1e_worst   = t1e[i_max],
    t1e_at_p0   = t1e[length(t1e)],
    at_boundary = isTRUE(all.equal(p_seq[i_max], p_null))
  )
}
