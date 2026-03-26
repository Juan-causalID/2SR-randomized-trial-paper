# Two-Stage Randomized Trial Design for Antimicrobial Strategies: Simulation Code

Simulation code for the manuscript: *"Advantages of a Two-Stage Randomized Trial Design to Evaluate Antimicrobial Treatment Strategies: a Simulation Study"*

## Overview

This repository contains all code to reproduce the agent-based model (ABM) simulations and figures reported in the manuscript. The ABM simulates a hospital ward with two competing bacterial strains (drug-susceptible and drug-resistant) under two-stage randomization, estimating direct, indirect, total, and overall causal effects on mortality following the Hudgens--Halloran framework.

## Requirements

- R >= 4.4.1
- R packages: `ABM` (>= 0.4.3), `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `scales`

Install packages:
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "scales"))
# ABM package:
install.packages("ABM")
```

## Code

All scripts are in `code/`:

| Script | Description | Replicates | Runtime |
|--------|-------------|-----------|---------|
| `run_sim_extract_100.R` | Main simulation: produces all causal effect estimates, infection/mortality outcomes, risk by infection type, and stratified direct effects (Tables 2--4) | 100 | ~10 min |
| `make_efigure_DE.R` | Time-varying direct effect simulation: produces Figure S3 and Table S3 data | 100 | ~5 min |
| `generate_figures.R` | Figure generation: produces Figures 2--7 | 100 | ~10 min |
| `run_SA.R` | Sensitivity analysis runner (26 scenarios, Tables S1--S2) | 5 | variable |
| `sensitivity_analysis.R` | Sensitivity analysis parameter definitions | -- | -- |

## Figures

Pre-generated figures are in `figures/`:

- `fig3_cumulative_deaths.png` -- Cumulative deaths by allocation strategy (Figure 3)
- `fig4_cumulative_incidence.png` -- Cumulative infections by strategy (Figure 4)
- Additional exploratory figures (fig2, fig5--fig7)

## Reproducing Results

```bash
# 1. Main simulation (Tables 2, 3, 4 + inline results)
Rscript code/run_sim_extract_100.R

# 2. Time-varying DE figure + Table S3
Rscript code/make_efigure_DE.R

# 3. Generate all figures
Rscript code/generate_figures.R

# 4. Sensitivity analyses (Tables S1, S2)
Rscript code/run_SA.R
```

## Notes

- Results will vary slightly across runs due to stochastic simulation. With 100 replicates, estimates are stable to ~0.1--0.3 percentage points.
- Sensitivity analyses use 5 replicates per scenario (supplementary, not main results).
