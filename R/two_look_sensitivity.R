# Exact two-look sensitivity for CP promising-zone re-estimation.
#
# The design has no interim efficacy rejection. At each prespecified look it
# may stop for futility, retain its current target, or increase that target.
# A later look may never reduce a target chosen earlier. The CP mapping always
# uses the score-z boundary; final-rule comparisons freeze that mapping and
# replace only the final rejection count.

.find_repo_file_tl <- function(relative_path) {
  candidates <- c(
    relative_path,
    file.path("..", relative_path),
    file.path(dirname(tryCatch(
      normalizePath(sys.frame(1)$ofile),
      error = function(e) "R/two_look_sensitivity.R"
    )), basename(relative_path))
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Could not locate ", relative_path, ". Run from the repository root.")
  }
  hit[[1L]]
}

if (!exists("cp_fn", mode = "function")) {
  source(.find_repo_file_tl("R/analysis.R"))
}
if (!exists("PRIMARY_CP_MODE", inherits = TRUE)) {
  source(.find_repo_file_tl("R/params.R"))
}

two_look_schedules <- tibble::tribble(
  ~schedule, ~n_look_1, ~n_look_2,
  "40% and 65% (n=34, 55)", 34L, 55L,
  "50% and 77% (n=42, 65)", 42L, 65L
)

.cp_target_tl <- function(x, n_seen, current_target,
                          cp_upper_value, corrected = FALSE) {
  candidates <- seq.int(current_target, n_max)
  values <- vapply(candidates, function(n_candidate) {
    cp_fn(
      x, n_seen, n_candidate, p0, p1, alpha,
      cp_mode = PRIMARY_CP_MODE, correct = corrected
    )
  }, numeric(1))
  hit <- candidates[values >= cp_upper_value]
  as.integer(if (length(hit)) hit[[1L]] else n_max)
}

.cp_action_tl <- function(x, n_seen, current_target, cp_lower_value,
                          cp_futility_value, cp_upper_value,
                          corrected = FALSE) {
  cp_now <- cp_fn(
    x, n_seen, current_target, p0, p1, alpha,
    cp_mode = PRIMARY_CP_MODE, correct = corrected
  )
  if (cp_now < cp_futility_value) {
    return(list(stop = TRUE, zone = "futility", target = n_seen, cp = cp_now))
  }
  if (cp_now >= cp_lower_value && cp_now < cp_upper_value) {
    increased_target <- .cp_target_tl(
      x, n_seen, current_target, cp_upper_value, corrected
    )
    return(list(
      stop = FALSE,
      zone = "promising",
      target = max(as.integer(current_target), increased_target),
      cp = cp_now
    ))
  }
  list(
    stop = FALSE,
    zone = if (cp_now >= cp_upper_value) "favorable" else "unfavorable",
    target = as.integer(current_target),
    cp = cp_now
  )
}

two_look_cp_oc_exact <- function(
    p_true, cp_lower_value, n_look_1, n_look_2,
    cp_futility_value = cp_futility,
    cp_upper_value = cp_upper,
    final_rule = c("z", "exact"),
    corrected = FALSE) {
  final_rule <- match.arg(final_rule)
  stopifnot(
    n_look_1 < n_look_2,
    n_look_2 < n_init,
    cp_futility_value <= cp_lower_value,
    cp_lower_value <= cp_upper_value
  )

  success_probability <- 0
  expected_n <- 0
  futility_probability <- 0
  increase_probability <- 0
  cap_probability <- 0
  terminal_probability <- 0

  for (x1 in 0L:n_look_1) {
    mass1 <- dbinom(x1, n_look_1, p_true)
    action1 <- .cp_action_tl(
      x1, n_look_1, n_init, cp_lower_value,
      cp_futility_value, cp_upper_value, corrected
    )
    if (action1$stop) {
      expected_n <- expected_n + mass1 * n_look_1
      futility_probability <- futility_probability + mass1
      terminal_probability <- terminal_probability + mass1
      next
    }

    stage2_n <- n_look_2 - n_look_1
    for (x_stage2 in 0L:stage2_n) {
      x2 <- x1 + x_stage2
      path_mass <- mass1 * dbinom(x_stage2, stage2_n, p_true)
      action2 <- .cp_action_tl(
        x2, n_look_2, action1$target, cp_lower_value,
        cp_futility_value, cp_upper_value, corrected
      )
      if (action2$stop) {
        expected_n <- expected_n + path_mass * n_look_2
        futility_probability <- futility_probability + path_mass
        terminal_probability <- terminal_probability + path_mass
        next
      }

      final_n <- action2$target
      critical_count <- if (final_rule == "z") {
        x_crit_z_fn(final_n, p0, alpha, correct = corrected)
      } else {
        x_crit_exact_fn(final_n, p0, alpha)
      }
      needed <- critical_count - x2
      remaining <- final_n - n_look_2
      conditional_success <- if (needed <= 0L) {
        1
      } else if (needed > remaining) {
        0
      } else {
        pbinom(needed - 1L, remaining, p_true, lower.tail = FALSE)
      }

      success_probability <- success_probability +
        path_mass * conditional_success
      expected_n <- expected_n + path_mass * final_n
      increase_probability <- increase_probability +
        path_mass * as.numeric(final_n > n_init)
      cap_probability <- cap_probability +
        path_mass * as.numeric(final_n == n_max)
      terminal_probability <- terminal_probability + path_mass
    }
  }

  stopifnot(
    success_probability >= 0,
    success_probability <= 1,
    abs(terminal_probability - 1) < 1e-12
  )

  tibble::tibble(
    p_true = p_true,
    cp_lower = cp_lower_value,
    n_look_1 = n_look_1,
    n_look_2 = n_look_2,
    final_rule = final_rule,
    power = success_probability,
    exp_n = expected_n,
    pr_fut = futility_probability,
    pr_increase = increase_probability,
    pr_cap = cap_probability
  )
}

.two_look_cell <- function(schedule, n_look_1, n_look_2, cp_lower) {
  null_z <- two_look_cp_oc_exact(
    p0, cp_lower, n_look_1, n_look_2, final_rule = "z"
  )
  alt_z <- two_look_cp_oc_exact(
    p1, cp_lower, n_look_1, n_look_2, final_rule = "z"
  )
  null_exact <- two_look_cp_oc_exact(
    p0, cp_lower, n_look_1, n_look_2, final_rule = "exact"
  )
  alt_exact <- two_look_cp_oc_exact(
    p1, cp_lower, n_look_1, n_look_2, final_rule = "exact"
  )
  null_mapping_only <- two_look_cp_oc_exact(
    p0, cp_lower, n_look_1, n_look_2,
    cp_futility_value = -Inf, final_rule = "z"
  )
  fixed_z <- fixed_n_oc(p0, n_init, p0, alpha, rule = "z")$power

  tibble::tibble(
    schedule = schedule,
    n_look_1 = n_look_1,
    n_look_2 = n_look_2,
    cp_lower = cp_lower,
    z_t1e = null_z$power,
    z_power = alt_z$power,
    exact_t1e = null_exact$power,
    exact_power = alt_exact$power,
    exp_n_null = null_z$exp_n,
    pr_fut_null = null_z$pr_fut,
    pr_increase_null = null_z$pr_increase,
    mapping_only_delta = null_mapping_only$power - fixed_z,
    final_rule_delta = null_z$power - null_exact$power
  )
}

two_look_cp_results <- tidyr::crossing(
  two_look_schedules,
  cp_lower = cp_lower_grid
) |>
  purrr::pmap_dfr(.two_look_cell)

two_look_cp_summary <- two_look_cp_results |>
  dplyr::filter(abs(cp_lower - 0.50) < 1e-12) |>
  dplyr::arrange(n_look_1)

stopifnot(
  all(two_look_cp_results$exact_t1e <= alpha + 1e-12),
  max(abs(two_look_cp_results$mapping_only_delta)) <
    min(two_look_cp_results$final_rule_delta)
)

.script_file_tl <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
)
if (length(.script_file_tl) == 1L &&
    basename(.script_file_tl) == "two_look_sensitivity.R") {
  cat("Exact two-look sensitivity at CP_L = 0.50\n")
  print(two_look_cp_summary, n = Inf, width = Inf)
  cat("\nAcross both schedules and the CP_L grid:\n")
  print(two_look_cp_results |>
    dplyr::summarise(
      z_t1e_min = min(z_t1e),
      z_t1e_max = max(z_t1e),
      exact_t1e_min = min(exact_t1e),
      exact_t1e_max = max(exact_t1e),
      mapping_abs_max = max(abs(mapping_only_delta)),
      final_rule_min = min(final_rule_delta),
      final_rule_max = max(final_rule_delta)
    ), n = Inf, width = Inf)
}
