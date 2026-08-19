# Focused regression checks for the exact-enumeration engine.
#
# Run from the repository root:
#   Rscript --vanilla R/regression_checks.R

.find_analysis <- function() {
  candidates <- c("R/analysis.R", "analysis.R", "../R/analysis.R")
  own <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA)
  if (!is.na(own)) {
    candidates <- c(candidates, file.path(dirname(own), "analysis.R"))
  }
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) stop("Could not locate R/analysis.R")
  hit[[1L]]
}

source(.find_analysis())

p0 <- 0.23
p1 <- 0.35
alpha <- 0.05
gamma_grid <- seq(0.95, 0.99, by = 0.005)
gamma_cal <- 0.955

# Published posterior-PP operating characteristics must remain reproducible.
post_t1e <- bayes_oc_exact(p0, gamma_cal)$power
post_power <- bayes_oc_exact(p1, gamma_cal)$power
stopifnot(
  isTRUE(all.equal(post_t1e, 0.046889133272, tolerance = 1e-11)),
  isTRUE(all.equal(post_power, 0.827493919619, tolerance = 1e-11))
)

# Posterior and exact-binomial boundaries are not equivalent. The frozen PP
# mapping reaches two sample sizes where their critical counts differ.
mapping <- bayes_pp_mapping(gamma_cal)
candidate_n <- 65L:200L
exact_crit <- vapply(candidate_n, x_crit_exact_fn, integer(1), p0, alpha)
post_crit <- vapply(
  candidate_n, x_crit_fn, integer(1),
  p0 = p0, ap = 0.5, bp = 0.5, thresh = gamma_cal
)
mismatch_n <- candidate_n[exact_crit != post_crit]
stopifnot(
  length(mismatch_n) == 33L,
  identical(
    intersect(
      sort(unique(mapping$n_final[mapping$zone %in% c(1L, 2L)])),
      mismatch_n
    ),
    c(87L, 125L)
  )
)

exact_t1e <- bayes_oc_exact_final(p0, gamma_cal)$power
exact_power <- bayes_oc_exact_final(p1, gamma_cal)$power
stopifnot(
  isTRUE(all.equal(exact_t1e, 0.0396897338242, tolerance = 1e-11)),
  isTRUE(all.equal(exact_power, 0.808112133288, tolerance = 1e-11)),
  exact_t1e < post_t1e,
  exact_power < post_power
)

# The selected posterior threshold controls; the immediately preceding
# prespecified grid value does not.
gamma_t1e <- vapply(
  gamma_grid,
  function(g) bayes_oc_exact(p0, g)$power,
  numeric(1)
)
selected <- min(gamma_grid[gamma_t1e <= alpha])
selected_i <- match(selected, gamma_grid)
stopifnot(
  identical(selected, gamma_cal),
  gamma_t1e[selected_i] <= alpha,
  gamma_t1e[selected_i - 1L] > alpha
)

# The standalone timing analysis must preserve the current 65/84 result,
# evaluate all four prespecified looks, and recalibrate PP prospectively at
# each timing.
source(file.path(dirname(.find_analysis()), "interim_timing_sensitivity.R"))
stopifnot(
  identical(interim_timing_cp_summary$n_int, c(34L, 42L, 55L, 65L)),
  isTRUE(all.equal(
    interim_timing_cp_summary$z_t1e,
    c(
      0.0507802106938570,
      0.0533984597420976,
      0.0544466463138578,
      0.0554621610159405
    ),
    tolerance = 1e-11
  )),
  identical(
    interim_timing_pp_calibrated$gamma_fin,
    c(0.960, 0.955, 0.955, 0.955)
  ),
  all(interim_timing_cp$exact_t1e <= alpha),
  max(abs(interim_timing_cp$mapping_only_delta)) <
    min(interim_timing_cp$final_rule_delta)
)

# Practice-oriented sensitivities must preserve coherent calibration and the
# published post-adaptation diagnostics. Sourcing does not run the command-line
# Monte Carlo check; that check is exercised by invoking the script directly.
source(file.path(dirname(.find_analysis()), "standard_practice_sensitivity.R"))
stopifnot(
  nrow(continuity_corrected_timing) == 44L,
  all(continuity_corrected_timing$t1e <= alpha),
  identical(pp_fine_selected$gamma_fin, gamma_cal),
  isTRUE(all.equal(
    historical_drift_results |>
      dplyr::filter(p_actual == p0, Design == "CP + z") |>
      dplyr::pull(positive_probability),
    0.0554621610159405,
    tolerance = 1e-11
  )),
  isTRUE(all.equal(
    estimation_sensitivity$bias,
    c(-0.004628246680, -0.005857515794,
      -0.000821441845, 0.005929880484),
    tolerance = 1e-9
  )),
  all(estimation_sensitivity$exact_coverage > 0.95),
  all(estimation_sensitivity$wilson_coverage > 0.93)
)

# The compact repeated-look sensitivity uses an exact pathwise recursion and
# must preserve the reported ordering across both schedules and the full grid.
source(file.path(dirname(.find_analysis()), "two_look_sensitivity.R"))
stopifnot(
  nrow(two_look_cp_results) == 22L,
  isTRUE(all.equal(
    two_look_cp_summary$z_t1e,
    c(0.0491690131935, 0.0511295275615),
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    two_look_cp_summary$z_power,
    c(0.792346351127, 0.823144337287),
    tolerance = 1e-10
  )),
  all(two_look_cp_results$exact_t1e <= alpha),
  max(abs(two_look_cp_results$mapping_only_delta)) <
    min(two_look_cp_results$final_rule_delta)
)

message("All statistical regression checks passed.")
