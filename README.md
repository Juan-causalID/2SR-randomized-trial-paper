# Two-Stage Randomized Trial Design for Antimicrobial Strategies: Simulation Code

Simulation code for the manuscript *"Advantages of a Two-Stage Randomized Trial Design to Evaluate Antimicrobial Treatment Strategies: a Simulation Study."*

## Overview

A two-stage (cluster + individual) randomized design, following Hudgens and Halloran, applied to antibiotic stewardship in a hospital ward. An agent-based model (ABM) simulates a ward with two competing bacterial strains (drug-susceptible and drug-resistant). Wards (clusters) are randomized to an allocation strategy that sets the drug mix, and patients within a ward are randomized to a drug. The design estimates the direct, indirect, total, and overall causal effects on mortality and detects the population-level (spillover) consequences of shifting the drug mix.

- Strategy alpha_1 = 90/10 (90% Drug A, 10% Drug B), the reference.
- Strategy alpha_0 = 50/50 (50% Drug A, 50% Drug B), the intervention.

Drug A and Drug B are model constructs, not specific drugs: both clear the susceptible strain, and only Drug B clears the resistant strain.

Every table and figure in the manuscript is reproduced by a single, self-contained script.

## Requirements

- R >= 4.4.1
- R packages: `ABM` (>= 0.4.3), `tidyverse`, `patchwork`, `scales`

```r
install.packages(c("ABM", "tidyverse", "patchwork", "scales"))
```

## Reproducing the results

```bash
Rscript code/2SR_simulation_consolidated.R
```

The script defines the ABM once, runs it, and derives all outputs from that single run so the tables and figures are mutually consistent. By default it:

- runs the main model (6 wards x 100 replicate trials) and prints **Table 2** (causal estimands with 95% simulation intervals), **Table 3** (per-arm infection and mortality totals), **Table S3** (time-varying direct effect by time window), and the whole-trial stratified direct effects;
- writes **Figures 1-3** (cumulative deaths; cumulative incidence; time-varying direct effect) to the working directory;
- caches the main run to `main_beta1.5_data.rds`. Delete this file to re-simulate from scratch.

Two optional sweeps are gated by flags near the top of the script (`CONFIG` block):

- `RUN_SENSITIVITY <- TRUE` -> one-way sensitivity analysis (**Tables S1/S2**);
- `RUN_EXTENDED <- TRUE` -> imperfect-Drug-B efficacy sweep (**Table S4**).

Runtime: the main run is roughly 20-30 minutes on a laptop; each optional sweep adds a further ~10-20 minutes. Results vary slightly across runs because the simulation is stochastic; with 100 replicates the estimates are stable to within their reported simulation intervals.

## Notes

- Model states, transition rates, and parameter values (with sources) are documented inline in the script header and body.
- No external files or sourcing are required; the single script is fully self-contained.
