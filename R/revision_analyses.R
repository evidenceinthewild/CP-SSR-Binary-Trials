# R/revision_analyses.R
# ----------------------------------------------------------------------------
# Driver for the SBR-26-068 major revision (Referees 1 and 2).
#
# Run from the repository root:
#     source("R/revision_analyses.R")
#
# Produces, in order:
#   [0] Reproduction check against the submitted manuscript
#   [1] Mehta-Pocock fidelity: CP under interim-estimated effect vs design p1
#   [2] Joint calibration: CP over the final z-test boundary (vs gamma_final)
#   [3] No-futility sensitivity
#   [4] Frozen adaptation mapping x final rule, plus fixed-N reference
#   [5] Combination-test comparator (exact and mid-p)
#   [6] Least-favourable-null check over the composite null p <= p0
#
# Every number is exact enumeration. No Monte Carlo, no seed dependence.
# ----------------------------------------------------------------------------

# Locate analysis.R regardless of the session's working directory. Tries the
# repo root first, then this script's own directory (which works when the file
# is sourced by path, e.g. source("~/.../R/revision_analyses.R")).
.find_analysis <- function() {
  cands <- c("R/analysis.R", "analysis.R", "../R/analysis.R")
  own <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA)
  if (!is.na(own)) cands <- c(cands, file.path(dirname(own), "analysis.R"))
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0L) {
    stop("Could not find analysis.R. setwd() to the CP-SSR-Binary-Trials ",
         "repo root and re-run, or source this file by its full path.")
  }
  hit[[1L]]
}
source(.find_analysis())
suppressPackageStartupMessages(library(tidyverse))

# ---- Base design (as submitted) --------------------------------------------
P0        <- 0.23
P1        <- 0.35
ALPHA     <- 0.05
N_INIT    <- 84L
N_INT     <- 65L
N_MAX     <- 200L
AP        <- 0.5
BP        <- 0.5
CP_UPPER  <- 0.80
CP_FUT    <- 0.10
GAMMA_INT <- 0.99
PP_UPPER  <- 0.50
CP_L_GRID <- seq(0.30, 0.80, by = 0.05)
PRIMARY_CP_MODE <- "estimated"

# Derive the calibrated posterior threshold from the prespecified grid rather
# than duplicating the manuscript's selected literal.
GAMMA_GRID <- seq(0.95, 0.99, by = 0.005)
gamma_t1e <- vapply(GAMMA_GRID, function(g) {
  bayes_oc_exact(
    p_true = P0, gamma_fin = g, gamma_int = GAMMA_INT,
    pp_fut = 0.05, pp_upper = PP_UPPER,
    n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
    p0 = P0, ap = AP, bp = BP
  )$power
}, numeric(1))
GAMMA_CAL <- min(GAMMA_GRID[gamma_t1e <= ALPHA])
gamma_cal_idx <- match(GAMMA_CAL, GAMMA_GRID)
stopifnot(gamma_t1e[gamma_cal_idx] <= ALPHA)
if (gamma_cal_idx > 1L) stopifnot(gamma_t1e[gamma_cal_idx - 1L] > ALPHA)

pct <- function(x) sprintf("%.2f%%", 100 * x)
rule <- function(s) cat("\n", strrep("=", 74), "\n", s, "\n",
                        strrep("=", 74), "\n", sep = "")

# ============================================================================
# [0] Reproduction check
# ============================================================================
rule("[0] REPRODUCTION CHECK (must match the submitted manuscript)")

fixed_t1e <- fixed_n_t1e(N_INIT, P0, ALPHA)
cat("Fixed-N z-test T1E at N_init =", N_INIT, ":", pct(fixed_t1e),
    "  (bias vs alpha:", sprintf("%+.2f pp", 100 * (fixed_t1e - ALPHA)), ")\n")
cat("This is score-z critical-count miscalibration with NO interim and NO SSR.\n")
cat("Referee 2's central observation: the adaptive design must be judged\n",
    "against THIS baseline, not against the nominal 5%.\n", sep = "")

base_cp <- map_dfr(CP_L_GRID, ~cp_oc_exact(
  p_true = P0, cp_lower = .x, cp_upper = CP_UPPER, cp_futility = CP_FUT,
  n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
  p0 = P0, p1 = P1, alpha = ALPHA, ap = AP, bp = BP))

cat("\nCP SSR T1E across the CP_L grid (design-alternative CP, as submitted):\n")
print(base_cp |> transmute(CP_L = cp_lower, T1E = pct(power)), n = Inf)
cat("\nRange:", pct(min(base_cp$power)), "to", pct(max(base_cp$power)),
    "| fixed-N baseline:", pct(fixed_t1e), "\n")

# ============================================================================
# [1] Mehta-Pocock fidelity (Referee 1)
# ============================================================================
rule("[1] CP UNDER INTERIM-ESTIMATED EFFECT vs DESIGN ALTERNATIVE")

cat("Mehta-Pocock (2011, Sec. 3) evaluate CP at the OBSERVED interim effect.\n",
    "The submitted manuscript evaluated it at the design alternative p1.\n",
    "If these differ materially, the submitted design must be renamed\n",
    "'design-alternative CP promising-zone SSR' or the analysis rerun.\n",
    sep = "")

fidelity <- map_dfr(CP_L_GRID, function(pl) {
  d <- cp_oc_exact(p_true = P0, cp_lower = pl, cp_upper = CP_UPPER,
                   cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                   n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                   ap = AP, bp = BP, cp_mode = "design")
  e <- cp_oc_exact(p_true = P0, cp_lower = pl, cp_upper = CP_UPPER,
                   cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                   n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                   ap = AP, bp = BP, cp_mode = "estimated")
  tibble(CP_L = pl,
         T1E_design = d$power, T1E_estimated = e$power,
         EN_design  = d$exp_n, EN_estimated  = e$exp_n,
         PrProm_design = d$pr_prom, PrProm_estimated = e$pr_prom)
})

print(fidelity |> transmute(
  CP_L,
  `T1E design`    = pct(T1E_design),
  `T1E estimated` = pct(T1E_estimated),
  `diff (pp)`     = sprintf("%+.2f", 100 * (T1E_estimated - T1E_design)),
  `E[N] design`   = round(EN_design, 1),
  `E[N] est`      = round(EN_estimated, 1)), n = Inf)

cat("\nPower at p1 under each convention (CP_L = 0.5):\n")
pw <- bind_rows(
  cp_oc_exact(p_true = P1, cp_lower = 0.50, cp_upper = CP_UPPER,
              cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
              n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
              ap = AP, bp = BP, cp_mode = "design") |> mutate(mode = "design"),
  cp_oc_exact(p_true = P1, cp_lower = 0.50, cp_upper = CP_UPPER,
              cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
              n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
              ap = AP, bp = BP, cp_mode = "estimated") |> mutate(mode = "estimated"))
print(pw |> transmute(mode, Power = pct(power), `E[N]` = round(exp_n, 1)))

cat("\nAll remaining revision analyses use cp_mode = '", PRIMARY_CP_MODE,
    "' so the primary results follow the observed-effect Mehta-Pocock ",
    "convention.\n", sep = "")

# ============================================================================
# [2] Joint boundary-and-mapping calibration
# ============================================================================
rule("[2] JOINT CALIBRATION: CP OVER ITS FINAL BOUNDARY (vs gamma_final)")

cat("The submitted comparison swept CP over CP_L but Bayesian over\n",
    "gamma_final. The coherent frequentist tuning parameter is the final\n",
    "z-test level alpha_fin, not CP_L. As with gamma_final in the PP rule,\n",
    "changing alpha_fin also updates the boundary targeted by SSR and can\n",
    "change the interim mapping; this is joint calibration, not a frozen-\n",
    "mapping final-rule contrast.\n",
    "Holding the promising zone fixed at CP_L = 0.5:\n\n", sep = "")

cp_boundary_oc <- function(p_true, alpha_fin) {
  cp_oc_exact(p_true = p_true, cp_lower = 0.50, cp_upper = CP_UPPER,
              cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
              n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
              ap = AP, bp = BP, cp_mode = PRIMARY_CP_MODE,
              alpha_fin = alpha_fin)
}

alpha_fin_grid <- c(0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02)
joint_calibration <- map_dfr(alpha_fin_grid, function(af) {
  t1e <- cp_boundary_oc(P0, af)
  pwr <- cp_boundary_oc(P1, af)
  tibble(alpha_fin = af, T1E = t1e$power, Power = pwr$power,
         EN = pwr$exp_n)
})
# NOTE: compute `controls` from the NUMERIC column before formatting it,
# otherwise the comparison is made against the formatted character string.
print(joint_calibration |> transmute(alpha_fin,
                           controls = ifelse(T1E <= ALPHA, "YES", "no"),
                           T1E = pct(T1E), Power = pct(Power),
                           `E[N] at p1` = round(EN, 1)) |>
        relocate(controls, .after = last_col()), n = Inf)

alpha_fin_scan <- seq(0.02, 0.05, by = 0.0001)
control_scan <- tibble(
  alpha_fin = alpha_fin_scan,
  T1E = map_dbl(alpha_fin_scan, ~cp_boundary_oc(P0, .x)$power)
)
ok <- control_scan |> filter(T1E <= ALPHA)
cat("\nCONCLUSION FOR THE 'UNCALIBRATABLE' CLAIM: ")
if (nrow(ok) > 0) {
  best_alpha <- max(ok$alpha_fin)
  best <- cp_boundary_oc(P1, best_alpha)
  best_t1e <- ok$T1E[which.max(ok$alpha_fin)]
  cat("the CP design can be made conservative by tightening the final ",
      "boundary.\nThe largest controlling alpha_fin on the 0.0001 grid is ",
      sprintf("%.4f", best_alpha), ", with T1E = ", pct(best_t1e),
      ", power = ", pct(best$power), ", and E[N] = ",
      sprintf("%.1f", best$exp_n), ".\n",
      "Because the attainable T1E jumps across discrete critical counts, ",
      "this establishes conservative control rather than exact calibration.\n",
      "The manuscript's claim must be narrowed to: no CP_L on the\n",
      "prespecified grid controls T1E holding CP_U, futility and the\n",
      "final z-test level fixed.\n", sep = "")
} else {
  cat("no alpha_fin on this grid controls T1E; the claim survives\n",
      "the joint calibration comparison and should say so explicitly.\n", sep = "")
}

# ============================================================================
# [3] No-futility sensitivity
# ============================================================================
rule("[3] NO-FUTILITY SENSITIVITY (Referee 1)")

cat("Futility stopping removes null trials before the final test and can\n",
    "therefore depress T1E. Disabling it quantifies the contribution of this\n",
    "component while changing the complete design.\n\n",
    sep = "")

nofut <- map_dfr(CP_L_GRID, function(pl) {
  w <- cp_oc_exact(p_true = P0, cp_lower = pl, cp_upper = CP_UPPER,
                   cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                   n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                   ap = AP, bp = BP, cp_mode = PRIMARY_CP_MODE)
  n <- cp_oc_exact(p_true = P0, cp_lower = pl, cp_upper = CP_UPPER,
                   cp_futility = -Inf, n_init = N_INIT, n_int = N_INT,
                   n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                   ap = AP, bp = BP, cp_mode = PRIMARY_CP_MODE)
  tibble(CP_L = pl, with_fut = w$power, no_fut = n$power)
})
print(nofut |> transmute(CP_L,
                         `T1E with futility` = pct(with_fut),
                         `T1E no futility`   = pct(no_fut),
                         `masked (pp)` = sprintf("%+.2f",
                                                 100 * (no_fut - with_fut))), n = Inf)

bay_fut <- bind_rows(
  bayes_oc_exact(p_true = P0, gamma_fin = GAMMA_CAL, gamma_int = GAMMA_INT,
                 pp_fut = 0.05, pp_upper = PP_UPPER, n_init = N_INIT,
                 n_int = N_INT, n_max = N_MAX, p0 = P0, ap = AP, bp = BP) |>
    mutate(futility = "on"),
  bayes_oc_exact(p_true = P0, gamma_fin = GAMMA_CAL, gamma_int = GAMMA_INT,
                 pp_fut = -Inf, pp_upper = PP_UPPER, n_init = N_INIT,
                 n_int = N_INT, n_max = N_MAX, p0 = P0, ap = AP, bp = BP) |>
    mutate(futility = "off"))
cat("\nBayesian PP design at gamma_final =", GAMMA_CAL, ":\n")
print(bay_fut |> transmute(futility, T1E = pct(power), `E[N]` = round(exp_n, 1)))
cat("With no PP futility boundary, observations formerly classified as futile\n",
    "enter the SSR branch; the resulting increase in E[N] is part of this\n",
    "sensitivity design and should be stated explicitly.\n", sep = "")

# ============================================================================
# [4] Frozen adaptation mapping x final rule
# ============================================================================
rule("[4] ADAPTATION MAPPING x FINAL RULE")

cat("Rows = what drives the sample-size decision.\n",
    "Columns = what decides success at the end.\n",
    "All evaluated at p = p0 (Type I error), base setting.\n\n", sep = "")

cell_A <- cp_oc_exact(p_true = P0, cp_lower = 0.50, cp_upper = CP_UPPER,
                      cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                      n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                      ap = AP, bp = BP,
                      cp_mode = PRIMARY_CP_MODE)$power

cell_B <- cp_oc_post_final(p_true = P0, cp_lower = 0.50, gamma_fin = GAMMA_CAL,
                           cp_upper = CP_UPPER, cp_futility = CP_FUT,
                           n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
                           p0 = P0, p1 = P1, alpha = ALPHA,
                           ap = AP, bp = BP,
                           cp_mode = PRIMARY_CP_MODE)$power

cell_C <- bayes_oc_exact(p_true = P0, gamma_fin = GAMMA_CAL,
                         gamma_int = GAMMA_INT, pp_fut = 0.05,
                         pp_upper = PP_UPPER, n_init = N_INIT, n_int = N_INT,
                         n_max = N_MAX, p0 = P0, ap = AP, bp = BP)$power

# Freeze the same PP mapping/interim decisions and replace only the posterior
# final boundary with the level-alpha exact-binomial boundary.
cell_C_exact <- bayes_oc_exact_final(
  p_true = P0, gamma_fin = GAMMA_CAL,
  gamma_int = GAMMA_INT, pp_fut = 0.05, pp_upper = PP_UPPER,
  n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
  p0 = P0, alpha = ALPHA, ap = AP, bp = BP
)$power

cell_D <- bayes_oc_ztest_final(p_true = P0, gamma_int = GAMMA_INT,
                               pp_fut = 0.05, pp_upper = PP_UPPER,
                               n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
                               p0 = P0, alpha = ALPHA, ap = AP, bp = BP,
                               gamma_fin = GAMMA_CAL,
                               pp_target = "posterior")$power

# Same cell with interim efficacy stopping disabled, so the comparison
# reflects the final rule alone rather than the Bayesian interim stop.
cell_D_noeff <- bayes_oc_ztest_final(p_true = P0, gamma_int = Inf,
                                     pp_fut = 0.05, pp_upper = PP_UPPER,
                                     n_init = N_INIT, n_int = N_INT,
                                     n_max = N_MAX, p0 = P0, alpha = ALPHA,
                                     ap = AP, bp = BP, gamma_fin = GAMMA_CAL,
                                     pp_target = "posterior")$power

cell_C_noeff <- bayes_oc_exact(p_true = P0, gamma_fin = GAMMA_CAL,
                               gamma_int = Inf, pp_fut = 0.05,
                               pp_upper = PP_UPPER, n_init = N_INIT,
                               n_int = N_INT, n_max = N_MAX,
                               p0 = P0, ap = AP, bp = BP)$power

twoby2 <- tribble(
  ~`Adaptation mapping`, ~`z-test final`, ~`posterior final`,
  "CP-driven",   pct(cell_A),            pct(cell_B),
  "PP-driven",   pct(cell_D),            pct(cell_C)
)
print(twoby2)
cat("\nPP mapping with posterior final:", pct(cell_C),
    "| same mapping with exact-binomial final:", pct(cell_C_exact), "\n")

cat("\nWith interim efficacy stopping DISABLED (isolates the final rule):\n")
print(tribble(
  ~`Adaptation mapping`, ~`z-test final`, ~`posterior final`,
  "CP-driven",   pct(cell_A),            pct(cell_B),
  "PP-driven",   pct(cell_D_noeff),      pct(cell_C_noeff)
))

cat("\nMain effects at p0 (no-efficacy-stop version):\n")
cat("  Final rule (z - posterior), averaged over mappings:  ",
    sprintf("%+.2f pp\n",
            100 * (((cell_A + cell_D_noeff) - (cell_B + cell_C_noeff)) / 2)))
cat("  Mapping (CP - PP), averaged over final rules:         ",
    sprintf("%+.2f pp\n",
            100 * (((cell_A + cell_B) - (cell_D_noeff + cell_C_noeff)) / 2)))
cat("\nThese are descriptive contrasts between frozen mappings and final rules.\n",
    "They do not estimate adaptation versus no adaptation; [4c] supplies ",
    "that reference.\n", sep = "")

cell_D_coherent_noeff <- bayes_oc_ztest_final(
  p_true = P0, gamma_int = Inf, pp_fut = 0.05, pp_upper = PP_UPPER,
  n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
  p0 = P0, alpha = ALPHA, ap = AP, bp = BP,
  pp_target = "ztest")$power
cat("\nSeparate coherent design (PP retargeted to the z-test boundary, ",
    "no efficacy stop): ", pct(cell_D_coherent_noeff),
    ". This is not used in the frozen-mapping contrasts because retargeting ",
    "changes the SSR mapping.\n", sep = "")

# ============================================================================
# [M] MECHANISM: is it discreteness, or is it the normal approximation?
# ============================================================================
rule("[M] MECHANISM CHECK: DISCRETENESS vs NORMAL-APPROXIMATION ERROR")

cat("Discreteness alone cannot make a test anticonservative: the exact\n",
    "binomial test is at least as discrete as the z-test and is conservative\n",
    "by construction. If the continuity-corrected z behaves like the exact\n",
    "test, the operative mechanism is APPROXIMATION ERROR, and the paper's\n",
    "causal language must say so.\n\n", sep = "")

cc_tbl <- crit_count_table(c(N_INT, N_INIT, 100L, 120L, 150L, N_MAX),
                           p0 = P0, alpha = ALPHA, ap = AP, bp = BP,
                           gamma_fin = GAMMA_CAL)

cat("Critical counts by rule (rejection iff x_total >= crit):\n")
print(cc_tbl |>
        select(n, rule, crit) |>
        pivot_wider(names_from = rule, values_from = crit), n = Inf)

cat("\nAttained fixed-n Type I error by rule:\n")
print(cc_tbl |>
        mutate(t1e = pct(t1e)) |>
        select(n, rule, t1e) |>
        pivot_wider(names_from = rule, values_from = t1e), n = Inf)

ident_at_init <- cc_tbl |> filter(n == N_INIT)
n_distinct_crit <- length(unique(ident_at_init$crit))
cat("\nAt N_init =", N_INIT, "the four rules give",
    n_distinct_crit, "distinct critical count(s):",
    paste(sort(unique(ident_at_init$crit)), collapse = ", "), "\n")
if (n_distinct_crit == 1L) {
  cat("=> At the planned sample size these are THE SAME TEST. Any apparent\n",
      "   'agreement' between them at N_init is an identity, not evidence of\n",
      "   similar base-design operating points. They can differ through the\n",
      "   extended-N part of the sample space.\n", sep = "")
}

# How much null mass actually reaches an extended N?
mass_extended <- cp_oc_exact(p_true = P0, cp_lower = 0.50,
                             cp_upper = CP_UPPER, cp_futility = CP_FUT,
                             n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
                             p0 = P0, p1 = P1, alpha = ALPHA,
                             ap = AP, bp = BP, cp_mode = PRIMARY_CP_MODE)
cat("\nNull mass entering the SSR (promising) branch at CP_L = 0.50:",
    pct(mass_extended$pr_prom), "\n")
cat("Null mass stopped for futility:", pct(mass_extended$pr_fut), "\n")

cat("\nContinuity-corrected z under the full adaptive design (CP_L = 0.50):\n")
cc_adaptive <- bind_rows(
  cp_oc_exact(p_true = P0, cp_lower = 0.50, cp_upper = CP_UPPER,
              cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
              n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA, ap = AP, bp = BP,
              cp_mode = PRIMARY_CP_MODE, correct = FALSE) |>
    mutate(rule = "z (uncorrected)"),
  cp_oc_exact(p_true = P0, cp_lower = 0.50, cp_upper = CP_UPPER,
              cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
              n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA, ap = AP, bp = BP,
              cp_mode = PRIMARY_CP_MODE, correct = TRUE) |>
    mutate(rule = "z (continuity-corrected)"))
print(cc_adaptive |> transmute(rule, T1E = pct(power),
                               `Pr(prom)` = pct(pr_prom)))

# ============================================================================
# [Z] ZONE OCCUPANCY AND THE CP PROJECTION CONVENTION
# ============================================================================
rule("[Z] Pr(PROMISING | H0) UNDER BOTH CP PROJECTION CONVENTIONS")

cat("If almost no null mass reaches the SSR branch, the adaptation cannot\n",
    "move Type I error much, and a small SSR effect is a property of the\n",
    "projection convention rather than a general finding.\n\n", sep = "")

zones <- map_dfr(CP_L_GRID, function(pl) {
  e <- cp_oc_exact(p_true = P0, cp_lower = pl, cp_upper = CP_UPPER,
                   cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                   n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                   ap = AP, bp = BP, cp_mode = "estimated")
  d <- cp_oc_exact(p_true = P0, cp_lower = pl, cp_upper = CP_UPPER,
                   cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                   n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                   ap = AP, bp = BP, cp_mode = "design")
  tibble(CP_L = pl,
         prom_est = e$pr_prom, t1e_est = e$power,
         prom_des = d$pr_prom, t1e_des = d$power)
})

fixed_z <- fixed_n_t1e(N_INIT, P0, ALPHA)
print(zones |> transmute(
  CP_L,
  `Pr(prom) est`  = pct(prom_est), `T1E est` = pct(t1e_est),
  `est vs fixed`  = sprintf("%+.2f", 100 * (t1e_est - fixed_z)),
  `Pr(prom) des`  = pct(prom_des), `T1E des` = pct(t1e_des),
  `des vs fixed`  = sprintf("%+.2f", 100 * (t1e_des - fixed_z))), n = Inf)

cat("\nFixed-N z-test baseline (no interim, no SSR):", pct(fixed_z), "\n")
cat("Estimated-effect convention: SSR changes T1E by ",
    sprintf("%+.2f to %+.2f pp\n",
            100 * min(zones$t1e_est - fixed_z),
            100 * max(zones$t1e_est - fixed_z)))
cat("Design-alternative convention: SSR changes T1E by ",
    sprintf("%+.2f to %+.2f pp\n",
            100 * min(zones$t1e_des - fixed_z),
            100 * max(zones$t1e_des - fixed_z)))
cat("\nIf these bracket zero with opposite signs, the SSR effect is\n",
    "convention-dependent and the manuscript must condition on it.\n", sep = "")

# ============================================================================
# [4c] FACTORIAL WITH THE NO-ADAPTATION REFERENCE ROW
# ============================================================================
rule("[4c] COMPARISON INCLUDING A NO-ADAPTATION ROW")

cat("The 2x2 in [4] has row levels {CP-triggered, PP-triggered}. That\n",
    "identifies which STATISTIC triggers adaptation, not whether adapting\n",
    "matters. Adding a fixed-N row supplies the missing reference level.\n\n",
    sep = "")

no_ssr_z    <- fixed_n_oc(P0, N_INIT, P0, ALPHA, "z",         AP, BP, GAMMA_CAL)$power
no_ssr_post <- fixed_n_oc(P0, N_INIT, P0, ALPHA, "posterior", AP, BP, GAMMA_CAL)$power

print(tribble(
  ~`Adaptation`,        ~`z final`,          ~`posterior final`,
  "none (fixed N)",     pct(no_ssr_z),       pct(no_ssr_post),
  "CP-triggered",       pct(cell_A),         pct(cell_B),
  "PP-triggered",       pct(cell_D_noeff),   pct(cell_C_noeff)
))

cat("\nEffect of ADAPTING, within each column (vs the fixed-N row):\n")
cat(sprintf("  CP-triggered, z final        : %+.2f pp\n", 100 * (cell_A - no_ssr_z)))
cat(sprintf("  PP-triggered, z final        : %+.2f pp\n", 100 * (cell_D_noeff - no_ssr_z)))
cat(sprintf("  CP-triggered, posterior final: %+.2f pp\n", 100 * (cell_B - no_ssr_post)))
cat(sprintf("  PP-triggered, posterior final: %+.2f pp\n", 100 * (cell_C_noeff - no_ssr_post)))

cat("\nEffect of the FINAL RULE at fixed N (no adaptation involved): ",
    sprintf("%+.2f pp\n", 100 * (no_ssr_z - no_ssr_post)))
cat("Interaction (2x2 adaptive cells only): ",
    sprintf("%+.2f pp\n",
            100 * ((cell_A - cell_B) - (cell_D_noeff - cell_C_noeff))))

# ============================================================================
# [4b] Exact / Clopper-Pearson final-analysis comparator (Referee 2, minor 2)
# ============================================================================
rule("[4b] EXACT BINOMIAL (CLOPPER-PEARSON) FINAL ANALYSIS")

cat("The submitted Discussion asserted that an exact final test 'would push\n",
    "T1E uniformly below alpha'. This section evaluates that assertion.\n",
    "Adaptation is held fixed (CP zone logic, z-boundary projection) so the\n",
    "only change from Section [0]/[1] is the final decision rule.\n\n", sep = "")

cat("Fixed-N reference, no interim and no SSR, at N_init =", N_INIT, ":\n")
cat("  z-test      :", pct(fixed_n_t1e(N_INIT, P0, ALPHA)), "\n")
cat("  exact test  :", pct(fixed_n_t1e_exact(N_INIT, P0, ALPHA)), "\n\n")

exact_final <- map_dfr(CP_L_GRID, function(pl) {
  t0 <- cp_oc_exact_final(p_true = P0, cp_lower = pl, cp_upper = CP_UPPER,
                          cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                          n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                          ap = AP, bp = BP, cp_mode = PRIMARY_CP_MODE)
  t1 <- cp_oc_exact_final(p_true = P1, cp_lower = pl, cp_upper = CP_UPPER,
                          cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                          n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                          ap = AP, bp = BP, cp_mode = PRIMARY_CP_MODE)
  z0 <- cp_oc_exact(p_true = P0, cp_lower = pl, cp_upper = CP_UPPER,
                    cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                    n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                    ap = AP, bp = BP, cp_mode = PRIMARY_CP_MODE)
  z1 <- cp_oc_exact(p_true = P1, cp_lower = pl, cp_upper = CP_UPPER,
                    cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                    n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                    ap = AP, bp = BP, cp_mode = PRIMARY_CP_MODE)
  tibble(CP_L = pl,
         T1E_exact = t0$power, Power_exact = t1$power,
         T1E_z     = z0$power, Power_z     = z1$power)
})

print(exact_final |> transmute(
  CP_L,
  `T1E exact-final` = pct(T1E_exact),
  `T1E z-final`     = pct(T1E_z),
  `Power exact`     = pct(Power_exact),
  `Power z`         = pct(Power_z),
  `power cost (pp)` = sprintf("%+.2f", 100 * (Power_exact - Power_z))), n = Inf)

cat("\nMax T1E with exact final:", pct(max(exact_final$T1E_exact)),
    "| nominal:", pct(ALPHA), "\n")
if (max(exact_final$T1E_exact) <= ALPHA) {
  cat("The assertion HOLDS on this grid: exact final controls T1E at every CP_L.\n")
} else {
  cat("The assertion FAILS: exact final exceeds alpha at some CP_L.\n",
      "The submitted Discussion's claim must be corrected, not merely softened.\n",
      sep = "")
}

lf_exact_final <- least_favorable_null(
  cp_oc_exact_final, p_null = P0,
  cp_lower = 0.50, cp_upper = CP_UPPER, cp_futility = CP_FUT,
  n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
  p0 = P0, p1 = P1, alpha = ALPHA, ap = AP, bp = BP,
  cp_mode = PRIMARY_CP_MODE)
cat("\nsup T1E over p <= p0 (exact final, CP_L = 0.5):",
    pct(lf_exact_final$t1e_worst),
    "at p =", sprintf("%.4f", lf_exact_final$p_worst), "\n")

# ============================================================================
# [5] Combination-test comparator
# ============================================================================
rule("[5] COMBINATION TEST WITH EXACT BINOMIAL STAGEWISE p-VALUES")

cat("Weights fixed at design stage: w1 = sqrt(n_int/N_init) = ",
    sprintf("%.4f", sqrt(N_INT / N_INIT)),
    ", w2 = ", sprintf("%.4f", sqrt(1 - N_INT / N_INIT)), "\n", sep = "")
cat("Conditional on stage-1 data, the selected second-stage size is fixed\n",
    "and its exact p-value is super-uniform. With design-fixed combination\n",
    "weights and futility paths treated as nonrejections, T1E is controlled.\n\n",
    sep = "")

comb <- map_dfr(CP_L_GRID, function(pl) {
  e0 <- comb_oc_exact(p_true = P0, cp_lower = pl, cp_upper = CP_UPPER,
                      cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                      n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                      pval = "exact", cp_mode = PRIMARY_CP_MODE)
  e1 <- comb_oc_exact(p_true = P1, cp_lower = pl, cp_upper = CP_UPPER,
                      cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                      n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                      pval = "exact", cp_mode = PRIMARY_CP_MODE)
  m0 <- comb_oc_exact(p_true = P0, cp_lower = pl, cp_upper = CP_UPPER,
                      cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                      n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                      pval = "midp", cp_mode = PRIMARY_CP_MODE)
  m1 <- comb_oc_exact(p_true = P1, cp_lower = pl, cp_upper = CP_UPPER,
                      cp_futility = CP_FUT, n_init = N_INIT, n_int = N_INT,
                      n_max = N_MAX, p0 = P0, p1 = P1, alpha = ALPHA,
                      pval = "midp", cp_mode = PRIMARY_CP_MODE)
  tibble(CP_L = pl,
         T1E_exact = e0$power, Power_exact = e1$power,
         T1E_midp  = m0$power, Power_midp  = m1$power)
})

print(comb |> transmute(CP_L,
                        `T1E exact`   = pct(T1E_exact),
                        `Power exact` = pct(Power_exact),
                        `T1E mid-p`   = pct(T1E_midp),
                        `Power mid-p` = pct(Power_midp)), n = Inf)

cat("\nMax T1E  exact:", pct(max(comb$T1E_exact)),
    " | mid-p:", pct(max(comb$T1E_midp)),
    " | nominal:", pct(ALPHA), "\n")
cat("Power cost of exactness at CP_L = 0.5: ",
    sprintf("%+.2f pp vs the z-test design\n",
            100 * (comb$Power_exact[comb$CP_L == 0.50] -
                   cp_oc_exact(p_true = P1, cp_lower = 0.50,
                               cp_upper = CP_UPPER, cp_futility = CP_FUT,
                               n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
                               p0 = P0, p1 = P1, alpha = ALPHA,
                               ap = AP, bp = BP,
                               cp_mode = PRIMARY_CP_MODE)$power)))

# ============================================================================
# [6] Least-favourable null
# ============================================================================
rule("[6] LEAST-FAVOURABLE NULL OVER p <= p0 (Referee 2)")

report_lf <- function(label, lf) {
  cat(label, "\n")
  cat("  T1E at p0                  :", pct(lf$t1e_at_p0), "\n")
  cat("  maximum on evaluated grid :", pct(lf$t1e_worst),
      "at p =", sprintf("%.4f", lf$p_worst), "\n")
  cat("  boundary is grid maximum   :", lf$at_boundary, "\n")
  if (!lf$at_boundary) {
    cat("  NOTE: the grid maximum is interior; refine locally before reporting.\n")
  }
  cat("\n")
}

lf_cp <- least_favorable_null(
  cp_oc_exact, p_null = P0,
  cp_lower = 0.50, cp_upper = CP_UPPER, cp_futility = CP_FUT,
  n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
  p0 = P0, p1 = P1, alpha = ALPHA, ap = AP, bp = BP,
  cp_mode = PRIMARY_CP_MODE)
report_lf("CP SSR + z-test final (CP_L = 0.5):", lf_cp)

lf_bayes <- least_favorable_null(
  bayes_oc_exact, p_null = P0,
  gamma_fin = GAMMA_CAL, gamma_int = GAMMA_INT, pp_fut = 0.05,
  pp_upper = PP_UPPER, n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
  p0 = P0, ap = AP, bp = BP)
report_lf("Bayesian PP + posterior final (gamma_final = 0.955):", lf_bayes)

lf_comb <- least_favorable_null(
  comb_oc_exact, p_null = P0,
  cp_lower = 0.50, cp_upper = CP_UPPER, cp_futility = CP_FUT,
  n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
  p0 = P0, p1 = P1, alpha = ALPHA, pval = "exact",
  cp_mode = PRIMARY_CP_MODE)
report_lf("Combination test, exact binomial (CP_L = 0.5):", lf_comb)

lf_midp <- least_favorable_null(
  comb_oc_exact, p_null = P0,
  cp_lower = 0.50, cp_upper = CP_UPPER, cp_futility = CP_FUT,
  n_init = N_INIT, n_int = N_INT, n_max = N_MAX,
  p0 = P0, p1 = P1, alpha = ALPHA, pval = "midp",
  cp_mode = PRIMARY_CP_MODE)
report_lf("Combination test, mid-p (CP_L = 0.5):", lf_midp)

rule("DONE. Numbers above feed Sections 4.x and the response letter.")
