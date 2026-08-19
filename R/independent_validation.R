# Independent base-R validation of selected published operating characteristics.
#
# This file deliberately does not source analysis.R or params.R and does not
# use any manuscript decision helper. It reimplements the base design from
# first principles so a shared helper error cannot make both calculations
# agree. Run from any directory with:
#
#   Rscript --vanilla R/independent_validation.R

p0_v <- 0.23
p1_v <- 0.35
alpha_v <- 0.05
n_plan_v <- 84L
n_look_v <- 65L
n_cap_v <- 200L
cp_lower_v <- 0.50
cp_upper_v <- 0.80
cp_futility_v <- 0.10
gamma_final_v <- 0.955
gamma_interim_v <- 0.99
pp_futility_v <- 0.05
pp_upper_v <- 0.50
prior_a_v <- 0.5
prior_b_v <- 0.5

z_count_v <- function(n, corrected = FALSE) {
  as.integer(floor(
    n * p0_v + qnorm(1 - alpha_v) * sqrt(n * p0_v * (1 - p0_v)) +
      if (corrected) 0.5 else 0
  ) + 1L)
}

exact_count_v <- function(n) {
  candidates <- 0L:n
  tails <- pbinom(candidates - 1L, n, p0_v, lower.tail = FALSE)
  hit <- candidates[tails <= alpha_v]
  if (length(hit)) hit[[1L]] else as.integer(n + 1L)
}

posterior_probability_v <- function(x, n) {
  pbeta(
    p0_v,
    prior_a_v + x,
    prior_b_v + n - x,
    lower.tail = FALSE
  )
}

posterior_count_v <- function(n, threshold = gamma_final_v) {
  candidates <- 0L:n
  hit <- candidates[vapply(
    candidates,
    posterior_probability_v,
    numeric(1),
    n = n
  ) >= threshold]
  if (length(hit)) hit[[1L]] else as.integer(n + 1L)
}

binomial_success_v <- function(x_seen, n_seen, n_final, p_true, count) {
  remaining <- n_final - n_seen
  needed <- count - x_seen
  if (needed <= 0L) return(1)
  if (needed > remaining) return(0)
  pbinom(needed - 1L, remaining, p_true, lower.tail = FALSE)
}

conditional_power_v <- function(x, n_seen, n_final, corrected = FALSE) {
  p_projected <- x / n_seen
  binomial_success_v(
    x,
    n_seen,
    n_final,
    p_projected,
    z_count_v(n_final, corrected)
  )
}

cp_target_v <- function(x, n_seen, corrected = FALSE) {
  candidates <- n_plan_v:n_cap_v
  cp_values <- vapply(
    candidates,
    function(n) conditional_power_v(x, n_seen, n, corrected),
    numeric(1)
  )
  hit <- candidates[cp_values >= cp_upper_v]
  if (length(hit)) hit[[1L]] else n_cap_v
}

cp_mapping_v <- function(corrected = FALSE) {
  x <- 0L:n_look_v
  cp_at_plan <- vapply(
    x,
    conditional_power_v,
    numeric(1),
    n_seen = n_look_v,
    n_final = n_plan_v,
    corrected = corrected
  )
  zone <- ifelse(
    cp_at_plan < cp_futility_v,
    "futility",
    ifelse(
      cp_at_plan < cp_lower_v,
      "unfavorable",
      ifelse(cp_at_plan < cp_upper_v, "promising", "favorable")
    )
  )
  target <- rep.int(n_plan_v, length(x))
  target[zone == "futility"] <- n_look_v
  target[zone == "promising"] <- vapply(
    x[zone == "promising"],
    cp_target_v,
    integer(1),
    n_seen = n_look_v,
    corrected = corrected
  )
  data.frame(x = x, zone = zone, target = as.integer(target))
}

cp_probability_v <- function(p_true, final_rule = c("z", "exact")) {
  final_rule <- match.arg(final_rule)
  mapping <- cp_mapping_v(FALSE)
  interim_mass <- dbinom(mapping$x, n_look_v, p_true)
  conditional_success <- vapply(seq_len(nrow(mapping)), function(i) {
    if (mapping$zone[[i]] == "futility") return(0)
    count <- if (final_rule == "z") {
      z_count_v(mapping$target[[i]])
    } else {
      exact_count_v(mapping$target[[i]])
    }
    binomial_success_v(
      mapping$x[[i]], n_look_v, mapping$target[[i]], p_true, count
    )
  }, numeric(1))
  sum(interim_mass * conditional_success)
}

beta_binomial_mass_v <- function(x, size, shape1, shape2) {
  exp(lchoose(size, x) + lbeta(x + shape1, size - x + shape2) -
        lbeta(shape1, shape2))
}

predictive_probability_v <- function(x, n_seen, n_final) {
  remaining <- n_final - n_seen
  future <- 0L:remaining
  final_count <- posterior_count_v(n_final)
  success <- x + future >= final_count
  shape1 <- prior_a_v + x
  shape2 <- prior_b_v + n_seen - x
  sum(beta_binomial_mass_v(future, remaining, shape1, shape2) * success)
}

pp_target_v <- function(x, n_seen) {
  candidates <- n_plan_v:n_cap_v
  values <- vapply(
    candidates,
    predictive_probability_v,
    numeric(1),
    x = x,
    n_seen = n_seen
  )
  hit <- candidates[values >= pp_upper_v]
  if (length(hit)) hit[[1L]] else n_cap_v
}

pp_mapping_v <- function() {
  x <- 0L:n_look_v
  posterior_at_look <- vapply(
    x, posterior_probability_v, numeric(1), n = n_look_v
  )
  pp_at_plan <- vapply(
    x,
    predictive_probability_v,
    numeric(1),
    n_seen = n_look_v,
    n_final = n_plan_v
  )
  zone <- ifelse(
    posterior_at_look >= gamma_interim_v,
    "efficacy",
    ifelse(
      pp_at_plan <= pp_futility_v,
      "futility",
      ifelse(pp_at_plan >= pp_upper_v, "favorable", "promising")
    )
  )
  target <- rep.int(n_plan_v, length(x))
  target[zone %in% c("efficacy", "futility")] <- n_look_v
  target[zone == "promising"] <- vapply(
    x[zone == "promising"],
    pp_target_v,
    integer(1),
    n_seen = n_look_v
  )
  data.frame(x = x, zone = zone, target = as.integer(target))
}

pp_probability_v <- function(p_true) {
  mapping <- pp_mapping_v()
  interim_mass <- dbinom(mapping$x, n_look_v, p_true)
  conditional_success <- vapply(seq_len(nrow(mapping)), function(i) {
    if (mapping$zone[[i]] == "futility") return(0)
    if (mapping$zone[[i]] == "efficacy") return(1)
    binomial_success_v(
      mapping$x[[i]],
      n_look_v,
      mapping$target[[i]],
      p_true,
      posterior_count_v(mapping$target[[i]])
    )
  }, numeric(1))
  sum(interim_mass * conditional_success)
}

validated <- c(
  fixed_z_t1e = pbinom(
    z_count_v(n_plan_v) - 1L,
    n_plan_v,
    p0_v,
    lower.tail = FALSE
  ),
  cp_z_t1e = cp_probability_v(p0_v, "z"),
  cp_z_power = cp_probability_v(p1_v, "z"),
  cp_exact_t1e = cp_probability_v(p0_v, "exact"),
  pp_posterior_t1e = pp_probability_v(p0_v),
  pp_posterior_power = pp_probability_v(p1_v)
)

published <- c(
  fixed_z_t1e = 0.0579856428370138,
  cp_z_t1e = 0.0554621610159405,
  cp_z_power = 0.827760454813420,
  cp_exact_t1e = 0.0364134942985834,
  pp_posterior_t1e = 0.0468891332720135,
  pp_posterior_power = 0.827493919619187
)

print(data.frame(
  estimand = names(validated),
  independent = unname(validated),
  published = unname(published),
  absolute_difference = abs(unname(validated - published)),
  row.names = NULL
), digits = 12)

stopifnot(
  identical(z_count_v(n_plan_v), 26L),
  identical(exact_count_v(n_plan_v), 27L),
  max(abs(validated - published)) < 1e-12
)

cat("Independent base-R validation passed.\n")
