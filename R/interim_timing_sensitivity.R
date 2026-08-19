# R/interim_timing_sensitivity.R
# -----------------------------------------------------------------------------
# Exact operating-characteristic sensitivity analysis for the timing of a
# single interim look. The script can be sourced by the manuscript or run
# directly from the repository root:
#
#   Rscript --vanilla R/interim_timing_sensitivity.R
#
# It deliberately retains one interim analysis. Adding repeated looks would
# require a different design specification for repeated stopping, re-estimation
# and combination-test weights rather than a sensitivity analysis of this
# manuscript's one-look design.
# -----------------------------------------------------------------------------

.find_repo_file <- function(relative_path) {
  candidates <- c(
    relative_path,
    file.path("..", relative_path),
    file.path(dirname(tryCatch(
      normalizePath(sys.frame(1)$ofile),
      error = function(e) "R/interim_timing_sensitivity.R"
    )), basename(relative_path))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("Could not locate ", relative_path, ". Run from the repository root.")
  }
  hit[[1L]]
}

if (!exists("cp_oc_exact", mode = "function")) {
  source(.find_repo_file("R/analysis.R"))
}
if (!exists("PRIMARY_CP_MODE", inherits = TRUE)) {
  source(.find_repo_file("R/params.R"))
}

# Prespecified information fractions. Rounding converts them to feasible
# integer patient counts; the current 65/84 look is retained exactly.
interim_timing_targets <- c(0.40, 0.50, 0.65, n_int / n_init)
interim_timing_grid <- tibble::tibble(
  target_fraction = interim_timing_targets,
  n_int_timing = as.integer(round(n_init * interim_timing_targets))
) |>
  dplyr::distinct(n_int_timing, .keep_all = TRUE) |>
  dplyr::mutate(
    information_fraction = n_int_timing / n_init,
    timing_label = sprintf(
      "%d%% (n=%d)",
      round(100 * information_fraction),
      n_int_timing
    )
  )

stopifnot(
  all(interim_timing_grid$n_int_timing > 0L),
  all(interim_timing_grid$n_int_timing < n_init),
  n_int %in% interim_timing_grid$n_int_timing
)

.scalar <- function(x, column) {
  value <- x[[column]]
  stopifnot(length(value) == 1L)
  value
}

fixed_z_null <- .scalar(
  fixed_n_oc(p0, n_init, p0, alpha, rule = "z"),
  "power"
)
fixed_z_power <- .scalar(
  fixed_n_oc(p1, n_init, p0, alpha, rule = "z"),
  "power"
)

.cp_timing_cell <- function(n_int_timing, information_fraction,
                            timing_label, cp_lower) {
  common <- list(
    cp_lower = cp_lower,
    cp_upper = cp_upper,
    n_init = n_init,
    n_int = n_int_timing,
    n_max = n_max,
    p0 = p0,
    p1 = p1,
    alpha = alpha,
    ap = ap,
    bp = bp,
    cp_mode = PRIMARY_CP_MODE
  )

  z_null <- do.call(
    cp_oc_exact,
    c(list(p_true = p0, cp_futility = cp_futility), common)
  )
  z_alt <- do.call(
    cp_oc_exact,
    c(list(p_true = p1, cp_futility = cp_futility), common)
  )
  z_null_no_futility <- do.call(
    cp_oc_exact,
    c(list(p_true = p0, cp_futility = -Inf), common)
  )
  exact_null <- do.call(
    cp_oc_exact_final,
    c(list(p_true = p0, cp_futility = cp_futility), common)
  )
  exact_alt <- do.call(
    cp_oc_exact_final,
    c(list(p_true = p1, cp_futility = cp_futility), common)
  )

  tibble::tibble(
    n_int = n_int_timing,
    information_fraction = information_fraction,
    timing_label = timing_label,
    cp_lower = cp_lower,
    z_t1e = .scalar(z_null, "power"),
    z_power = .scalar(z_alt, "power"),
    exact_t1e = .scalar(exact_null, "power"),
    exact_power = .scalar(exact_alt, "power"),
    exp_n_null = .scalar(z_null, "exp_n"),
    pr_fut_null = .scalar(z_null, "pr_fut"),
    pr_prom_null = .scalar(z_null, "pr_prom"),
    pr_cap_null = .scalar(z_null, "pr_cap"),
    total_adaptation_delta = .scalar(z_null, "power") - fixed_z_null,
    mapping_only_delta =
      .scalar(z_null_no_futility, "power") - fixed_z_null,
    final_rule_delta =
      .scalar(z_null, "power") - .scalar(exact_null, "power")
  )
}

interim_timing_cp <- tidyr::crossing(
  interim_timing_grid,
  cp_lower = cp_lower_grid
) |>
  dplyr::select(-target_fraction) |>
  purrr::pmap_dfr(.cp_timing_cell)

# A compact operating point for prose and the supplementary table. The full
# CP_L sweep remains available in interim_timing_cp.
interim_timing_cp_summary <- interim_timing_cp |>
  dplyr::filter(abs(cp_lower - 0.50) < 1e-12) |>
  dplyr::arrange(information_fraction)

.pp_timing_cell <- function(n_int_timing, information_fraction,
                            timing_label, gamma_fin) {
  common <- list(
    gamma_fin = gamma_fin,
    gamma_int = gamma_int,
    pp_fut = pp_fut,
    pp_upper = pp_upper,
    n_init = n_init,
    n_int = n_int_timing,
    n_max = n_max,
    p0 = p0,
    ap = ap,
    bp = bp
  )
  null <- do.call(bayes_oc_exact, c(list(p_true = p0), common))
  alt <- do.call(bayes_oc_exact, c(list(p_true = p1), common))

  tibble::tibble(
    n_int = n_int_timing,
    information_fraction = information_fraction,
    timing_label = timing_label,
    gamma_fin = gamma_fin,
    t1e = .scalar(null, "power"),
    power = .scalar(alt, "power"),
    exp_n_null = .scalar(null, "exp_n"),
    pr_eff_null = .scalar(null, "pr_eff"),
    pr_fut_null = .scalar(null, "pr_fut"),
    pr_prom_null = .scalar(null, "pr_prom"),
    pr_cap_null = .scalar(null, "pr_cap")
  )
}

interim_timing_pp <- tidyr::crossing(
  interim_timing_grid,
  gamma_fin = gamma_fin_grid
) |>
  dplyr::select(-target_fraction) |>
  purrr::pmap_dfr(.pp_timing_cell)

# Recalibrate prospectively at each timing: choose the smallest prespecified
# posterior threshold whose exact null rejection probability is <= alpha.
interim_timing_pp_calibrated <- interim_timing_pp |>
  dplyr::filter(t1e <= alpha) |>
  dplyr::group_by(n_int, information_fraction, timing_label) |>
  dplyr::slice_min(gamma_fin, n = 1L, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::arrange(information_fraction)

stopifnot(
  nrow(interim_timing_pp_calibrated) == nrow(interim_timing_grid),
  all(interim_timing_cp$exact_t1e <= alpha + 1e-12),
  max(abs(interim_timing_cp$mapping_only_delta)) <
    min(interim_timing_cp$final_rule_delta),
  interim_timing_pp_calibrated$gamma_fin[
    interim_timing_pp_calibrated$n_int == n_int
  ] == gamma_fin_cal
)

.script_file <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
)
if (length(.script_file) == 1L &&
    basename(.script_file) == "interim_timing_sensitivity.R") {
  cat("Interim-timing sensitivity: CP design at CP_L = 0.50\n")
  print(interim_timing_cp_summary, n = Inf, width = Inf)
  cat("\nInterim-timing sensitivity: prospectively calibrated PP design\n")
  print(interim_timing_pp_calibrated, n = Inf, width = Inf)
}
