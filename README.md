# cp-ssr-binary-trials

Reproducibility code for the manuscript:

> **Conditional Power Promising Zone Sample Size Re-estimation Inflates
> Type I Error in Single-Arm Binary Trials**

This repository contains the simulation code, exact-enumeration operating
characteristic (OC) calculations, and the Quarto source for the manuscript
and all of its tables and figures.

## Summary of the paper

Mehta and Pocock's promising-zone sample size re-estimation (SSR) was
designed for two-arm continuous endpoints, where the Cui-Hung-Wang weighted
test preserves the type I error rate without alpha adjustment. The same
construction is sometimes carried over to single-arm binary trials with a
plain z-test final analysis. We show by exact enumeration that this carry-
over does not preserve nominal alpha. In our representative oncology base
setting, no value of the promising-zone lower threshold ($\mathrm{CP}_L$)
on a fine grid achieves $\mathrm{T1E} \leq \alpha$. Across a wider
generalizability grid of five $(p_0, p_1, n_\text{int}, N_\text{init},
N_\text{max})$ configurations, the picture is mixed: some configurations
have $\mathrm{CP}_L$ values that control T1E and some do not, because the
discrete z-test critical count's rounding bias changes sign across
configurations. Tuning $\mathrm{CP}_L$ on a fixed grid is therefore
unreliable: a value that controls T1E in one design may not in another.
A Bayesian predictive-probability SSR with a posterior-probability final
rule controls type I error in the same setting, and the symmetric
comparison (CP-driven SSR with a Bayesian posterior final analysis) shows
that swapping the final rule alone is sufficient to rescue calibration,
isolating the final-decision rule rather than the SSR rule itself as the
source of the inflation.

## Repository layout

```
cp-ssr-binary-trials/
├── README.md
├── LICENSE
├── CP_PromisingZone_SSR_TypeIError.qmd   # manuscript source (sources R/analysis.R)
├── references.bib                   # bibliography
└── R/
    └── analysis.R                   # simulator + exact OC functions
```

`R/analysis.R` is self-contained and can be sourced from R for ad-hoc
exploration without rendering the manuscript.

## Software requirements

- R >= 4.4 (developed with R 4.6.0)
- Quarto >= 1.4
- A LaTeX distribution (TinyTeX is sufficient) with `lualatex` for PDF output

CRAN packages:

- tidyverse
- gt (>= 1.0.0)
- scales
- patchwork
- glue
- latex2exp

Install all of them with:

```r
install.packages(c(
  "tidyverse", "gt", "scales", "patchwork", "glue", "latex2exp"
))
```

If you do not yet have a LaTeX engine:

```r
install.packages("tinytex")
tinytex::install_tinytex()
```

## How to render the manuscript

From the repository root:

```bash
quarto render CP_PromisingZone_SSR_TypeIError.qmd --to pdf
quarto render CP_PromisingZone_SSR_TypeIError.qmd --to html
```

The first render performs the full exact-enumeration sweep over the
parameter grid and may take several minutes on a laptop. Subsequent
renders reuse cached chunks unless the parameters change.

## How to use the analysis code on its own

```r
source("R/analysis.R")

# Exact Type I error of a CP promising-zone SSR design with a z-test final.
# cp_oc_exact() enumerates the joint (x_int, x_rem) binomial distribution
# and returns operating characteristics deterministically (no Monte Carlo).
res <- cp_oc_exact(
  p_true   = 0.20,             # evaluate under the null (p_true = p0)
  cp_lower = 0.30,
  cp_upper = 0.80,
  cp_futility = 0.10,
  n_init   = 30L, n_int = 20L, n_max = 60L,
  p0 = 0.20, p1 = 0.40,
  alpha    = 0.05,
  ap = 1, bp = 1
)
res$power   # under p_true = p0 this is the Type I error
```

`bayes_oc_exact()` returns the same shape for the Bayesian PP design,
and `cp_oc_post_final()` does so for the symmetric variant
(CP-driven SSR with a Bayesian posterior final analysis). The Monte
Carlo simulators (`run_cp_scenario`, `run_bayes_scenario`) are useful
for extensions where exact enumeration is awkward (e.g., multiple
interim looks); they take an explicit `n_sim` argument.

## Citation

If you use this code, please cite the manuscript (citation block to be
added on acceptance).

## License

Code is released under the MIT License (see `LICENSE`).
