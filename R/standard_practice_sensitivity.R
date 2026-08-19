# R/standard_practice_sensitivity.R
# -----------------------------------------------------------------------------
# Practice-oriented sensitivity analyses for the one-look binary design:
#   1. continuity-corrected score-z operating characteristics;
#   2. finer-grid posterior-threshold calibration;
#   3. historical-benchmark drift (decision risk, not formal Type I error);
#   4. exact bias, RMSE and naive interval coverage after adaptation; and
#   5. a seeded Monte Carlo cross-check when this file is run directly.
#
# Run from the repository root:
#   Rscript --vanilla R/standard_practice_sensitivity.R
# -----------------------------------------------------------------------------

.find_repo_file_sp <- function(relative_path) {
  candidates <- c(
    relative_path,
    file.path("..", relative_path),
    file.path(dirname(tryCatch(
      normalizePath(sys.frame(1)$ofile),
      error = function(e) "R/standard_practice_sensitivity.R"
    )), basename(relative_path))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("Could not locate ", relative_path, ". Run from the repository root.")
  }
  hit[[1L]]
}

if (!exists("cp_oc_exact", mode = "function")) {
  source(.find_repo_file_sp("R/analysis.R"))
}
if (!exists("PRIMARY_CP_MODE", inherits = TRUE)) {
  source(.find_repo_file_sp("R/params.R"))
}
if (!exists("interim_timing_grid", inherits = TRUE)) {
  source(.find_repo_file_sp("R/interim_timing_sensitivity.R"))
}

.scalar_sp <- function(x, column) {
  value <- x[[column]]
  stopifnot(length(value) == 1L)
  value
}

# -----------------------------------------------------------------------------
# Continuity-corrected score-z design
# -----------------------------------------------------------------------------
# This is a coherent design alternative: conditional power and the SSR target
# are recalculated against the corrected final boundary. It is therefore not a
# frozen-mapping final-rule contrast.

.cc_timing_cell <- function(n_int_timing, information_fraction,
                            timing_label, cp_lower) {
  common <- list(
    cp_lower = cp_lower,
    cp_upper = cp_upper,
    cp_futility = cp_futility,
    n_init = n_init,
    n_int = n_int_timing,
    n_max = n_max,
    p0 = p0,
    p1 = p1,
    alpha = alpha,
    ap = ap,
    bp = bp,
    cp_mode = PRIMARY_CP_MODE,
    correct = TRUE
  )
  null <- do.call(cp_oc_exact, c(list(p_true = p0), common))
  alt <- do.call(cp_oc_exact, c(list(p_true = p1), common))
  tibble::tibble(
    n_int = n_int_timing,
    information_fraction = information_fraction,
    timing_label = timing_label,
    cp_lower = cp_lower,
    t1e = .scalar_sp(null, "power"),
    power = .scalar_sp(alt, "power"),
    exp_n_null = .scalar_sp(null, "exp_n"),
    pr_fut_null = .scalar_sp(null, "pr_fut"),
    pr_prom_null = .scalar_sp(null, "pr_prom")
  )
}

continuity_corrected_timing <- tidyr::crossing(
  interim_timing_grid,
  cp_lower = cp_lower_grid
) |>
  dplyr::select(-target_fraction) |>
  purrr::pmap_dfr(.cc_timing_cell)

continuity_corrected_base <- continuity_corrected_timing |>
  dplyr::filter(n_int == !!n_int) |>
  dplyr::arrange(cp_lower)

stopifnot(all(continuity_corrected_timing$t1e <= alpha + 1e-12))

# -----------------------------------------------------------------------------
# Finer posterior-threshold calibration grid
# -----------------------------------------------------------------------------

gamma_fin_fine_grid <- seq(0.950, 0.990, by = 0.001)
pp_fine_calibration <- purrr::map_dfr(gamma_fin_fine_grid, function(gamma) {
  null <- bayes_oc_exact(
    p_true = p0, gamma_fin = gamma, gamma_int = gamma_int,
    pp_fut = pp_fut, pp_upper = pp_upper,
    n_init = n_init, n_int = n_int, n_max = n_max,
    p0 = p0, ap = ap, bp = bp
  )
  alt <- bayes_oc_exact(
    p_true = p1, gamma_fin = gamma, gamma_int = gamma_int,
    pp_fut = pp_fut, pp_upper = pp_upper,
    n_init = n_init, n_int = n_int, n_max = n_max,
    p0 = p0, ap = ap, bp = bp
  )
  tibble::tibble(
    gamma_fin = gamma,
    t1e = .scalar_sp(null, "power"),
    power = .scalar_sp(alt, "power"),
    exp_n_null = .scalar_sp(null, "exp_n")
  )
})

pp_fine_selected <- pp_fine_calibration |>
  dplyr::filter(t1e <= alpha) |>
  dplyr::slice_min(gamma_fin, n = 1L, with_ties = FALSE)

stopifnot(
  nrow(pp_fine_selected) == 1L,
  pp_fine_selected$t1e <= alpha
)

# -----------------------------------------------------------------------------
# Historical-benchmark drift
# -----------------------------------------------------------------------------
# These are positive-decision probabilities when the actual untreated
# response rate differs from the design benchmark. Values above p0 are not
# Type I error under the formal hypothesis; they quantify application risk
# from using a misspecified historical benchmark.

historical_rate_grid <- p0 + c(-0.04, -0.02, 0, 0.02, 0.04)

historical_drift_results <- purrr::map_dfr(
  historical_rate_grid,
  function(p_actual) {
    dplyr::bind_rows(
      fixed_n_oc(
        p_actual, n_init, p0, alpha, rule = "z",
        ap = ap, bp = bp, gamma_fin = gamma_fin_cal
      ) |>
        dplyr::transmute(
          Design = "Fixed-N z",
          positive_probability = power
        ),
      cp_oc_exact(
        p_actual, 0.50, cp_upper = cp_upper,
        cp_futility = cp_futility, n_init = n_init, n_int = n_int,
        n_max = n_max, p0 = p0, p1 = p1, alpha = alpha,
        ap = ap, bp = bp, cp_mode = PRIMARY_CP_MODE
      ) |>
        dplyr::transmute(
          Design = "CP + z",
          positive_probability = power
        ),
      cp_oc_exact_final(
        p_actual, 0.50, cp_upper = cp_upper,
        cp_futility = cp_futility, n_init = n_init, n_int = n_int,
        n_max = n_max, p0 = p0, p1 = p1, alpha = alpha,
        ap = ap, bp = bp, cp_mode = PRIMARY_CP_MODE
      ) |>
        dplyr::transmute(
          Design = "CP + exact final",
          positive_probability = power
        ),
      bayes_oc_exact(
        p_actual, gamma_fin_cal, gamma_int = gamma_int,
        pp_fut = pp_fut, pp_upper = pp_upper,
        n_init = n_init, n_int = n_int, n_max = n_max,
        p0 = p0, ap = ap, bp = bp
      ) |>
        dplyr::transmute(
          Design = "PP + posterior",
          positive_probability = power
        )
    ) |>
      dplyr::mutate(
        p_actual = p_actual,
        drift = p_actual - p0,
        .before = 1L
      )
  }
)

# -----------------------------------------------------------------------------
# Exact estimation and naive interval coverage
# -----------------------------------------------------------------------------

.cp_mapping_sp <- function(cp_lower = 0.50) {
  x_int_grid <- 0L:n_int
  cp_obs <- vapply(
    x_int_grid,
    function(x) cp_fn(
      x, n_int, n_init, p0, p1, alpha,
      cp_mode = PRIMARY_CP_MODE
    ),
    numeric(1)
  )
  zone <- integer(length(x_int_grid))
  zone[cp_obs < cp_futility] <- 0L
  zone[cp_obs >= cp_upper] <- 3L
  zone[cp_obs >= cp_lower & cp_obs < cp_upper] <- 2L
  zone[cp_obs >= cp_futility & cp_obs < cp_lower] <- 1L
  n_final <- ifelse(
    zone == 0L,
    n_int,
    ifelse(
      zone == 2L,
      vapply(
        x_int_grid,
        function(x) ssr_cp_fn(
          x, n_int, p0, p1, alpha, cp_upper, n_init, n_max,
          cp_mode = PRIMARY_CP_MODE
        ),
        integer(1)
      ),
      n_init
    )
  )
  tibble::tibble(
    x_int = x_int_grid,
    zone = zone,
    n_final = as.integer(n_final)
  )
}

.binom_intervals_sp <- function(x, n, conf = 0.95) {
  tail <- (1 - conf) / 2
  exact_lower <- ifelse(x == 0L, 0, qbeta(tail, x, n - x + 1L))
  exact_upper <- ifelse(x == n, 1, qbeta(1 - tail, x + 1L, n - x))

  z <- qnorm(1 - tail)
  phat <- x / n
  denom <- 1 + z^2 / n
  wilson_center <- (phat + z^2 / (2 * n)) / denom
  wilson_half <- z * sqrt(
    phat * (1 - phat) / n + z^2 / (4 * n^2)
  ) / denom

  tibble::tibble(
    exact_lower = exact_lower,
    exact_upper = exact_upper,
    wilson_lower = pmax(0, wilson_center - wilson_half),
    wilson_upper = pmin(1, wilson_center + wilson_half)
  )
}

.estimation_oc_sp <- function(p_true, mapping, design) {
  pieces <- purrr::map_dfr(seq_len(nrow(mapping)), function(i) {
    x_int_value <- mapping$x_int[i]
    n_final_value <- mapping$n_final[i]
    n_rem <- n_final_value - n_int
    x_rem <- 0L:n_rem
    probability <- dbinom(x_int_value, n_int, p_true) *
      dbinom(x_rem, n_rem, p_true)
    x_total <- x_int_value + x_rem
    intervals <- .binom_intervals_sp(x_total, n_final_value)
    tibble::tibble(
      probability = probability,
      x_total = x_total,
      n_final = n_final_value,
      phat = x_total / n_final_value
    ) |>
      dplyr::bind_cols(intervals)
  })

  stopifnot(abs(sum(pieces$probability) - 1) < 1e-10)
  tibble::tibble(
    Design = design,
    p_true = p_true,
    bias = sum(pieces$probability * (pieces$phat - p_true)),
    rmse = sqrt(sum(pieces$probability * (pieces$phat - p_true)^2)),
    exact_coverage = sum(
      pieces$probability *
        (pieces$exact_lower <= p_true & pieces$exact_upper >= p_true)
    ),
    wilson_coverage = sum(
      pieces$probability *
        (pieces$wilson_lower <= p_true & pieces$wilson_upper >= p_true)
    ),
    exp_n = sum(pieces$probability * pieces$n_final)
  )
}

cp_estimation_mapping <- .cp_mapping_sp(0.50)
pp_estimation_mapping <- bayes_pp_mapping(
  gamma_fin_cal, gamma_int, pp_fut, pp_upper,
  n_init, n_int, n_max, p0, ap, bp
)

estimation_sensitivity <- purrr::map_dfr(c(p0, p1), function(p_true) {
  dplyr::bind_rows(
    .estimation_oc_sp(p_true, cp_estimation_mapping, "CP mapping"),
    .estimation_oc_sp(p_true, pp_estimation_mapping, "PP mapping")
  )
})

# -----------------------------------------------------------------------------
# Seeded stochastic validation, run only from the command line
# -----------------------------------------------------------------------------

.script_file_sp <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
)

if (length(.script_file_sp) == 1L &&
    basename(.script_file_sp) == "standard_practice_sensitivity.R") {
  n_sim <- as.integer(Sys.getenv("CP_SSR_MC_N", unset = "25000"))
  set.seed(20260812)

  exact_cp <- cp_oc_exact(
    p0, 0.50, cp_upper = cp_upper, cp_futility = cp_futility,
    n_init = n_init, n_int = n_int, n_max = n_max,
    p0 = p0, p1 = p1, alpha = alpha, ap = ap, bp = bp,
    cp_mode = PRIMARY_CP_MODE
  )$power
  sim_cp <- run_cp_scenario(
    p0, 0.50, n_sim = n_sim,
    cp_upper = cp_upper, cp_futility = cp_futility,
    n_init = n_init, n_int = n_int, n_max = n_max,
    p0 = p0, p1 = p1, alpha = alpha, ap = ap, bp = bp,
    cp_mode = PRIMARY_CP_MODE
  )$power

  exact_pp <- bayes_oc_exact(
    p0, gamma_fin_cal, gamma_int = gamma_int,
    pp_fut = pp_fut, pp_upper = pp_upper,
    n_init = n_init, n_int = n_int, n_max = n_max,
    p0 = p0, ap = ap, bp = bp
  )$power
  sim_pp <- run_bayes_scenario(
    p0, gamma_fin_cal, n_sim = n_sim, gamma_int = gamma_int,
    pp_fut = pp_fut, pp_upper = pp_upper,
    n_init = n_init, n_int = n_int, n_max = n_max,
    p0 = p0, ap = ap, bp = bp
  )$power

  validation <- tibble::tibble(
    Design = c("CP + z", "PP + posterior"),
    exact = c(exact_cp, exact_pp),
    simulated = c(sim_cp, sim_pp)
  ) |>
    dplyr::mutate(
      mcse = sqrt(exact * (1 - exact) / n_sim),
      standardized_difference = (simulated - exact) / mcse
    )

  stopifnot(all(abs(validation$standardized_difference) <= 4))
  cat("Seeded Monte Carlo cross-check (n =", n_sim, ")\n")
  print(validation, n = Inf, width = Inf)
  cat("\nContinuity-corrected score-z range across timings:\n")
  print(
    continuity_corrected_timing |>
      dplyr::group_by(timing_label) |>
      dplyr::summarise(
        t1e_min = min(t1e),
        t1e_max = max(t1e),
        power_min = min(power),
        power_max = max(power),
        .groups = "drop"
      ),
    n = Inf,
    width = Inf
  )
  cat("\nFine-grid PP calibration:\n")
  print(pp_fine_selected, n = Inf, width = Inf)
  cat("\nExact estimation diagnostics:\n")
  print(estimation_sensitivity, n = Inf, width = Inf)
}
