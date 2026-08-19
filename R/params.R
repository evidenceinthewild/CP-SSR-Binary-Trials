# R/params.R
# ----------------------------------------------------------------------------
# Design parameters shared by the manuscript and the supplementary material.
#
# Both documents source this file so that a change to the design cannot leave
# the two out of step. Nothing here computes operating characteristics; it only
# fixes the setting. Sourced AFTER R/analysis.R.
# ----------------------------------------------------------------------------

# ── Fixed design parameters (representative oncology setting) ─────────────
p0      <- 0.23   # null response rate (historical control)
p1      <- 0.35   # target response rate
alpha   <- 0.05   # one-sided Type I error target
power   <- 0.80   # target power
n_init  <- 84L    # planned N (one-sided alpha = 0.05 hypothesis test)
n_int   <- 65L    # interim analysis (absolute N)
n_max   <- 200L   # hard budget cap

# ── Bayesian prior ────────────────────────────────────────────────────────
ap <- 0.5; bp <- 0.5   # Jeffreys Beta(0.5, 0.5)

# ── CP zone parameters ────────────────────────────────────────────────────
cp_upper    <- 0.80   # fixed throughout; favorable zone above this
cp_futility <- 0.10   # fixed; stop for futility below this

# ── Bayesian thresholds ───────────────────────────────────────────────────
gamma_int   <- 0.99   # interim efficacy stop bar, kept high so the interim
                      # branch contributes only a small share of T1E. At
                      # n_int = 65, P(post >= 0.95 | p = p0) is already ~5.5%,
                      # which by itself would consume more than alpha; 0.95
                      # is too lenient once the interim stop is active.
pp_fut      <- 0.05   # PP futility boundary
pp_upper    <- 0.50   # PP favorable boundary (no increase at or above)

# ── Scenario grid: null + insurance region + alternative ─────────────────
delta     <- p1 - p0
scenarios <- round(c(
  p0,
  p0 + 0.25 * delta,
  p0 + 0.50 * delta,
  p0 + 0.75 * delta,
  p1,
  p1 + 0.25 * delta
), 4)

scen_labels <- c(
  "p0",
  "p0 + 0.25 delta",
  "p0 + 0.50 delta",
  "p0 + 0.75 delta",
  "p1",
  "p1 + 0.25 delta"
)

# Lookup table for joining labels, avoids fragile rep() calls throughout
scen_label_tbl <- tibble::tibble(p_true = scenarios, Label = scen_labels)

# ── Swept grids ───────────────────────────────────────────────────────────
cp_lower_grid  <- seq(0.30, 0.80, by = 0.05)
gamma_fin_grid <- seq(0.95, 0.99, by = 0.005)

# Primary analysis follows Mehta-Pocock by projecting under the observed
# interim response rate. The design-alternative projection is retained as
# a sensitivity analysis in R/revision_analyses.R.
PRIMARY_CP_MODE <- "estimated"

# Calibrated posterior threshold: the smallest value on the prespecified grid
# whose exact Type I error does not exceed alpha. Defined here so that the
# manuscript, the supplement and the regression checks all use one value.
gamma_fin_cal <- {
  g_t1e <- vapply(gamma_fin_grid, function(g) bayes_oc_exact(
    p_true = p0, gamma_fin = g, gamma_int = gamma_int, pp_fut = pp_fut,
    pp_upper = pp_upper, n_init = n_init, n_int = n_int, n_max = n_max,
    p0 = p0, ap = ap, bp = bp)$power, numeric(1))
  min(gamma_fin_grid[g_t1e <= alpha])
}
