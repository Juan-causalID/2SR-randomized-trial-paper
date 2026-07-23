#!/usr/bin/env Rscript
################################################################################
# 2SR TWO-STAGE RANDOMIZED TRIAL -- CONSOLIDATED SIMULATION CODE
#
# Self-contained, single-file reproduction of every simulation result in the
# manuscript. No external sourcing, no stray parameters: the agent-based model
# is defined once, run once, and all tables and figures are derived from that
# single run so they are mutually consistent.
#
# ---------------------------------------------------------------------------
# THE STUDY
#   A two-stage (cluster + individual) randomized design (Hudgens & Halloran)
#   applied to antibiotic stewardship in a hospital ward. Clusters (wards) are
#   randomized to an allocation STRATEGY that sets the probability a treated
#   patient receives Drug A vs Drug B; within a cluster, individuals are
#   randomized to a drug according to that probability.
#     * Strategy alpha_1 = 90/10  (90% Drug A, 10% Drug B)   -- "reference"
#     * Strategy alpha_0 = 50/50  (50% Drug A, 50% Drug B)   -- "intervention"
#   Drug A and Drug B are MODEL CONSTRUCTS, not real drugs. In the model:
#     * both drugs clear the susceptible strain,
#     * only Drug B clears the resistant strain (Drug A does not cover it).
#   The design is meant to detect the population-level (spillover) consequences
#   of shifting the drug mix: treating the resistant strain more (more Drug B)
#   suppresses it, and susceptible-strain infections rise to take its place
#   ("competitive release").
#
# THE MODEL (agent-based, continuous-time, frequency-dependent transmission)
#   States (compartments):
#     X    uncolonized (can acquire carriage)
#     S    colonized by the SUSCEPTIBLE strain (asymptomatic carrier)
#     R    colonized by the RESISTANT   strain (asymptomatic carrier)
#     US   symptomatic susceptible-strain infection, UNtreated
#     UR   symptomatic resistant-strain  infection, UNtreated
#     TS1  symptomatic susceptible-strain infection, treated with Drug A
#     TS2  symptomatic susceptible-strain infection, treated with Drug B
#     TR1  symptomatic resistant-strain  infection, treated with Drug A
#     TR2  symptomatic resistant-strain  infection, treated with Drug B
#     Z    discharged alive (absorbing)
#     D    dead (absorbing)
#
# KEY REGIME CHOICE: beta = 1.5
#   The ward is run NEAR the resistant strain's epidemic threshold. There, the
#   50/50 strategy visibly suppresses R and the full 2SR estimand suite appears,
#   with competitive release showing up in INCIDENCE (new R falls, new S rises).
#   (A much higher beta saturates R and hides this in incidence; a much lower
#   beta produces little transmission at all.)
#
# HOW TO RUN
#   Rscript 2SR_simulation_consolidated.R
#   Edit the CONFIG block below to point OUT_DIR at a writable folder and to
#   toggle the optional supplement analyses. The main run (~a few minutes) is
#   cached to an .rds file; delete it to re-simulate.
################################################################################

suppressPackageStartupMessages({
  library(ABM)          # agent-based modelling engine (provides Simulation, newCounter, ...)
  library(tidyverse)    # data wrangling + ggplot2
  library(patchwork)    # compose multi-panel figures
  library(scales)       # axis formatting
})

## ===========================================================================
## CONFIG  (edit these)
## ===========================================================================
OUT_DIR        <- "."          # where figures / csv / cache are written
N_CLUSTERS     <- 6            # 3 clusters per strategy arm
N_REALIZATIONS <- 100          # independent replicate trials (Monte Carlo sample)
TIMESTEPS      <- 350          # trial length (simulation steps)

RUN_FIGURES     <- TRUE        # Figures: cumulative incidence, cumulative deaths, time-varying DE
RUN_SENSITIVITY <- FALSE       # Tables S1/S2: one-way sensitivity sweep (slow; ~10-20 min)
RUN_EXTENDED    <- FALSE       # Table S4: imperfect Drug B (eff_B sweep)   (slow)

CACHE_RDS <- file.path(OUT_DIR, "main_beta1.5_data.rds")

## ===========================================================================
## 1. THE AGENT-BASED MODEL
## ---------------------------------------------------------------------------
## One function. Defaults ARE the production (main-analysis) parameters, so
## runsim() with no arguments reproduces the main model. `pi` (drug-mix) and,
## for the supplement, `eff_B` are the only knobs the analyses vary.
## ===========================================================================
runsim <- function(
    ## --- allocation strategy knob (the trial's exposure) ---
    pi     = 0.9,     # P(Drug A | treated). 0.9 = 90/10 strategy; 0.5 = 50/50 strategy
    ## --- transmission / epidemic regime ---
    beta   = 1.5,     # transmission rate (frequency-dependent); 1.5 = near R's threshold
    c      = 0.01,    # fitness cost of resistance: R transmits at beta*(1-c)
    lambda = 35,      # admission rate (exponential inter-arrival of new patients)
    ## --- natural history ---
    sigma  = 0.5,     # progression rate: colonized -> symptomatic infection
    rho    = 0.9,     # P(treated | symptomatic)
    gammaS = 0.2,     # baseline clearance rate, susceptible strain
    gammaR = 0.2,     # baseline clearance rate, resistant strain
    tau    = 0.7,     # ADDED clearance from an effective drug
    ## --- discharge / mortality ---
    mu     = 0.095,   # discharge rate, symptomatic patients
    mu_d   = 0.010,   # discharge rate, colonized/uninfected patients ("discharge arrows")
    delta  = 0.025,   # baseline death rate (scaled by a per-state multiplier below)
    ## --- population setup ---
    n      = 1000,    # initial ward census
    S0     = 50,      # initial susceptible-strain carriers  (5%)
    R0     = 50,      # initial resistant-strain  carriers  (5%)
    m      = 0.75,    # fraction of ADMISSIONS arriving uncolonized (X) vs colonized-S
    ## --- Drug B efficacy knob (supplement only; 1.0 = full efficacy = main model) ---
    eff_B  = 1.0,
    timesteps = 350
) {
  ## New patients arrive as a Poisson stream (exponential gaps at rate lambda);
  ## a fraction m arrive uncolonized (X), the rest already carrying the S strain.
  admit.event <- function(time) newEvent(time + rexp(1, lambda), function(time, sim, agent) {
    a <- if (runif(1) < m) newAgent("X") else newAgent("S")
    addAgent(sim, a); schedule(sim, admit.event(time))
  })

  ## Seed the ward: S0 S-carriers, R0 R-carriers, the rest uncolonized.
  sim <- Simulation$new(as.list(c(rep("S", S0), rep("R", R0), rep("X", n - S0 - R0))))
  sim$schedule(admit.event(0))

  ## ---- Loggers ----
  ## Occupancy counters: current census of each state at each timestep.
  for (st in c("X","S","R","US","TS1","TS2","UR","TR1","TR2","Z","D"))
    sim$addLogger(newCounter(st, st))
  ## Incidence counters (from -> to): per-interval FLOW, so summing over time
  ## gives the cumulative number of that transition. Used for "new infections".
  sim$addLogger(newCounter("TS1_inc", "S", "TS1")); sim$addLogger(newCounter("TS2_inc", "S", "TS2"))
  sim$addLogger(newCounter("US_inc",  "S", "US"))
  sim$addLogger(newCounter("TR1_inc", "R", "TR1")); sim$addLogger(newCounter("TR2_inc", "R", "TR2"))
  sim$addLogger(newCounter("UR_inc",  "R", "UR"))
  ## Death counters by originating state (for deaths attributable to each drug x strain).
  sim$addLogger(newCounter("D_TS1", "TS1", "D")); sim$addLogger(newCounter("D_TS2", "TS2", "D"))
  sim$addLogger(newCounter("D_US",  "US",  "D"))
  sim$addLogger(newCounter("D_TR1", "TR1", "D")); sim$addLogger(newCounter("D_TR2", "TR2", "D"))
  sim$addLogger(newCounter("D_UR",  "UR",  "D"))

  mx <- newRandomMixing(); sim$addContact(mx)   # well-mixed ward (mass-action contacts)

  ## ---- Progression to symptomatic infection + treatment assignment ----
  ## A carrier progresses at rate sigma; of those, a fraction rho is treated and,
  ## if treated, gets Drug A with prob pi (-> T*1) or Drug B with prob 1-pi (-> T*2);
  ## the untreated fraction 1-rho goes to U*.
  sim$addTransition("S" -> "TS1", sigma * rho * pi)
  sim$addTransition("S" -> "TS2", sigma * rho * (1 - pi))
  sim$addTransition("S" -> "US",  sigma * (1 - rho))
  sim$addTransition("R" -> "TR1", sigma * rho * pi)
  sim$addTransition("R" -> "TR2", sigma * rho * (1 - pi))
  sim$addTransition("R" -> "UR",  sigma * (1 - rho))

  ## ---- Recovery (return to uncolonized X) ----
  ## An effective drug adds tau to the baseline clearance. THE KEY ASYMMETRY:
  ## Drug A helps against S (TS1 gets +tau) but NOT against R (TR1 gets only gammaR);
  ## Drug B helps against both (TR2 gets +tau*eff_B). eff_B=1 => full Drug B efficacy.
  sim$addTransition("TS1" -> "X", gammaS + tau)
  sim$addTransition("TS2" -> "X", gammaS + tau)
  sim$addTransition("US"  -> "X", gammaS)
  sim$addTransition("TR1" -> "X", gammaR)
  sim$addTransition("TR2" -> "X", gammaR + tau * eff_B)
  sim$addTransition("UR"  -> "X", gammaR)

  ## ---- Discharge alive ----
  ## Symptomatic patients leave at rate mu; colonized/uninfected patients leave at
  ## rate mu_d (the "discharge arrows" -- without them the ward grows unboundedly).
  sim$addTransition("TS1" -> "Z", mu); sim$addTransition("TS2" -> "Z", mu); sim$addTransition("US" -> "Z", mu)
  sim$addTransition("TR1" -> "Z", mu); sim$addTransition("TR2" -> "Z", mu); sim$addTransition("UR" -> "Z", mu)
  sim$addTransition("X"   -> "Z", mu_d); sim$addTransition("S" -> "Z", mu_d); sim$addTransition("R" -> "Z", mu_d)

  ## ---- Death (rate = delta * multiplier) ----
  ## Multipliers encode the clinical ordering:
  ##   concordant (right drug) best; untreated worse; discordant R (wrong drug) worst.
  ##   TR2 multiplier = 2.5 - 1.5*eff_B  =>  eff_B=1 -> 1.0 (Drug B fully works);
  ##                                          eff_B=0 -> 2.5 (Drug B == Drug A on R).
  sim$addTransition("TS1" -> "D", delta * 1)                  # concordant S
  sim$addTransition("TS2" -> "D", delta * 1)                  # concordant S
  sim$addTransition("US"  -> "D", delta * 2)                  # untreated S
  sim$addTransition("TR1" -> "D", delta * 2.5)               # discordant R (Drug A, no coverage)
  sim$addTransition("TR2" -> "D", delta * (2.5 - 1.5 * eff_B)) # concordant R (Drug B)
  sim$addTransition("UR"  -> "D", delta * 3)                  # untreated R

  ## ---- Transmission (any carrier/infected of a strain colonizes an X) ----
  ## Frequency-dependent at rate beta; the resistant strain pays the fitness cost (1-c).
  sim$addTransition("S"   + "X" -> "S"   + "S" ~ mx, beta)
  sim$addTransition("TS1" + "X" -> "TS1" + "S" ~ mx, beta)
  sim$addTransition("TS2" + "X" -> "TS2" + "S" ~ mx, beta)
  sim$addTransition("US"  + "X" -> "US"  + "S" ~ mx, beta)
  sim$addTransition("R"   + "X" -> "R"   + "R" ~ mx, beta * (1 - c))
  sim$addTransition("TR1" + "X" -> "TR1" + "R" ~ mx, beta * (1 - c))
  sim$addTransition("TR2" + "X" -> "TR2" + "R" ~ mx, beta * (1 - c))
  sim$addTransition("UR"  + "X" -> "UR"  + "R" ~ mx, beta * (1 - c))

  sim$run(0:timesteps)   # returns a data.frame: one row per timestep, one column per logger (+ `times`)
}

## Robust wrapper: a rare RNG glitch in the engine can return a short/failed run;
## retry until we get a complete run (correct number of rows).
safe_runsim <- function(...) {
  for (k in 1:12) {
    o <- tryCatch(runsim(...), error = function(e) NULL)
    if (!is.null(o) && nrow(o) == TIMESTEPS + 1) return(o)
  }
  stop("runsim failed 12 times")
}

## ===========================================================================
## 2. REDUCERS  (collapse one raw run into small, analysis-ready summaries)
## ---------------------------------------------------------------------------
## We never keep the full 351-row frames for 600 runs (that caused memory
## blow-ups). Instead each run is reduced immediately to three tibbles.
## ===========================================================================
WIN_BREAKS <- seq(0, TIMESTEPS, 50)                       # 50-step windows for time-varying DE
WIN_MIDS   <- (WIN_BREAKS[-1] + WIN_BREAKS[-length(WIN_BREAKS)]) / 2

## (a) whole-trial totals: cumulative assignments/deaths by drug and by strain,
##     plus end-of-trial census (for the overall-mortality denominator).
reduce_final <- function(x, real, alloc) {
  L <- x[nrow(x), ]
  tibble(
    real = real, alloc = alloc,
    admitted = L$X + L$S + L$R + L$US + L$TS1 + L$TS2 + L$UR + L$TR1 + L$TR2 + L$Z + L$D,
    total_D  = L$D,
    assign_A = sum(x$TS1_inc + x$TR1_inc), death_A = sum(x$D_TS1 + x$D_TR1),  # by drug
    assign_B = sum(x$TS2_inc + x$TR2_inc), death_B = sum(x$D_TS2 + x$D_TR2),
    assign_S = sum(x$TS1_inc + x$TS2_inc + x$US_inc), death_S = sum(x$D_TS1 + x$D_TS2 + x$D_US),  # by strain
    assign_R = sum(x$TR1_inc + x$TR2_inc + x$UR_inc), death_R = sum(x$D_TR1 + x$D_TR2 + x$D_UR),
    new_S = sum(x$TS1_inc + x$TS2_inc + x$US_inc), new_R = sum(x$TR1_inc + x$TR2_inc + x$UR_inc)
  )
}
## (b) per-timestep cumulative series (only needed for the figures).
reduce_series <- function(x, real, alloc) tibble(
  real = real, alloc = alloc, times = x$times,
  new_S = x$TS1_inc + x$TS2_inc + x$US_inc, new_R = x$TR1_inc + x$TR2_inc + x$UR_inc,
  dS = x$D_TS1 + x$D_TS2 + x$D_US, dR = x$D_TR1 + x$D_TR2 + x$D_UR, D = x$D
)
## (c) windowed drug x strain counts (for the time-varying / stratified DE).
reduce_windows <- function(x, real, alloc) {
  wi <- findInterval(x$times, WIN_BREAKS, rightmost.closed = TRUE)
  tibble(real = real, alloc = alloc, w = WIN_MIDS[pmin(pmax(wi, 1), length(WIN_MIDS))],
         nAS = x$TS1_inc, dAS = x$D_TS1, nBS = x$TS2_inc, dBS = x$D_TS2,
         nAR = x$TR1_inc, dAR = x$D_TR1, nBR = x$TR2_inc, dBR = x$D_TR2) %>%
    group_by(real, alloc, w) %>% summarise(across(everything(), sum), .groups = "drop")
}

## ===========================================================================
## 3. RUN THE TRIAL (cached)
## ---------------------------------------------------------------------------
## 6 clusters: clusters 1-3 use the 90/10 strategy, clusters 4-6 use 50/50.
## For each cluster we simulate N_REALIZATIONS replicate wards; a "realization"
## of the whole trial pools the 3 clusters within an arm.
## ===========================================================================
if (file.exists(CACHE_RDS)) {
  cat("Loading cached run:", CACHE_RDS, "(delete this file to re-simulate)\n")
  D <- readRDS(CACHE_RDS); FIN <- D$FIN; TS <- D$TS; WIN <- D$WIN
} else {
  cat(sprintf("Simulating %d clusters x %d realizations = %d runs (beta=1.5)\n",
              N_CLUSTERS, N_REALIZATIONS, N_CLUSTERS * N_REALIZATIONS))
  t0 <- Sys.time()
  FINl <- vector("list", N_CLUSTERS * N_REALIZATIONS); TSl <- FINl; WINl <- FINl; idx <- 1
  for (i in 1:N_CLUSTERS) {
    alloc <- if (i <= N_CLUSTERS / 2) "90/10" else "50/50"
    pi_i  <- if (i <= N_CLUSTERS / 2) 0.9     else 0.5
    for (j in 1:N_REALIZATIONS) {
      x <- safe_runsim(pi = pi_i, timesteps = TIMESTEPS)   # all other params = production defaults
      FINl[[idx]] <- reduce_final(x, j, alloc)
      TSl[[idx]]  <- reduce_series(x, j, alloc)
      WINl[[idx]] <- reduce_windows(x, j, alloc)
      idx <- idx + 1
    }
    cat(sprintf("  cluster %d done (%s), %.1f min elapsed\n",
                i, alloc, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    gc(verbose = FALSE)   # reclaim memory once per cluster (not per run)
  }
  FIN <- bind_rows(FINl); TS <- bind_rows(TSl); WIN <- bind_rows(WINl)
  saveRDS(list(FIN = FIN, TS = TS, WIN = WIN), CACHE_RDS)
  cat("Done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min; cached to", CACHE_RDS, "\n")
}

## ===========================================================================
## 4. CAUSAL ESTIMANDS  (Tables 2 and 3)
## ---------------------------------------------------------------------------
## Risk = deaths / assignments among the relevant treated group. For each
## realization we pool the 3 clusters in an arm (mean of cluster risks), then
## take differences BETWEEN arms per realization, and summarize across
## realizations as mean + 95% simulation interval (2.5/97.5 percentiles).
##
## Sign conventions (as reported in the paper):
##   DE(alpha) = risk(Drug A, alpha) - risk(Drug B, alpha)   > 0: A worse (A misses R)
##   IE(drug)  = risk(drug, alpha_1=90/10) - risk(drug, alpha_0=50/50)   spillover of the mix
##   OE        = risk_all(alpha_0=50/50) - risk_all(alpha_1=90/10)   < 0: 50/50 lowers mortality
## ===========================================================================
estimands <- function(fin) {
  arm <- fin %>%
    mutate(oe = total_D / admitted,
           riskA = death_A / assign_A, riskB = death_B / assign_B,
           riskS = death_S / assign_S, riskR = death_R / assign_R) %>%
    group_by(real, alloc) %>%                       # pool the 3 clusters in each arm
    summarise(oe = mean(oe), riskA = mean(riskA), riskB = mean(riskB),
              riskS = mean(riskS), riskR = mean(riskR),
              new_S = sum(new_S), new_R = sum(new_R),
              death_S = sum(death_S), death_R = sum(death_R), total_D = sum(total_D),
              .groups = "drop")
  wide <- function(col) arm %>% select(real, alloc, all_of(col)) %>%
    pivot_wider(names_from = alloc, values_from = all_of(col))
  a1A <- wide("riskA")$`90/10`; a0A <- wide("riskA")$`50/50`   # Drug A risk at alpha_1 / alpha_0
  a1B <- wide("riskB")$`90/10`; a0B <- wide("riskB")$`50/50`   # Drug B risk at alpha_1 / alpha_0
  o1  <- wide("oe")$`90/10`;    o0  <- wide("oe")$`50/50`      # overall mortality at alpha_1 / alpha_0
  vecs <- list(
    `DE(alpha0=50/50)` = a0A - a0B,   # direct effect under 50/50
    `DE(alpha1=90/10)` = a1A - a1B,   # direct effect under 90/10
    `IE(A)`            = a1A - a0A,   # indirect (spillover) effect on Drug A patients
    `IE(B)`            = a1B - a0B,   # indirect (spillover) effect on Drug B patients
    `TE1`              = a1A - a0B,   # total effect: (A,90/10) vs (B,50/50)
    `TE2`              = a1B - a0A,   # total effect: (B,90/10) vs (A,50/50)
    `OE (alpha0-alpha1)` = o0 - o1)   # overall effect
  tab <- imap_dfr(vecs, ~tibble(estimand = .y,
    pp = mean(.x, na.rm = TRUE) * 100,
    lo = quantile(.x, .025, na.rm = TRUE) * 100,
    hi = quantile(.x, .975, na.rm = TRUE) * 100))
  arm_totals <- arm %>% group_by(alloc) %>%
    summarise(across(c(new_S, new_R, death_S, death_R, total_D), ~round(mean(.))),
              riskA = round(mean(riskA) * 100, 1), riskB = round(mean(riskB) * 100, 1),
              .groups = "drop")
  list(table = tab, arms = arm_totals)
}

E <- estimands(FIN)
cat("\n===== TABLE 2: causal estimands (percentage points, 95% SI) =====\n")
E$table %>% mutate(across(c(pp, lo, hi), ~sprintf("%+.2f", .))) %>%
  transmute(estimand, `estimate [95% SI]` = sprintf("%s [%s, %s]", pp, lo, hi)) %>%
  as.data.frame() %>% print(right = FALSE)
cat("\n===== TABLE 3: per-arm totals (mean over realizations) =====\n")
print(as.data.frame(E$arms), right = FALSE)

## ===========================================================================
## 5. TIME-VARYING & STRATIFIED DIRECT EFFECT  (Table S3 + main-text stratified DE)
## ---------------------------------------------------------------------------
## Within each 50-step window, DE and its strain-specific parts show how the
## direct effect shrinks as the resistant strain washes out of the ward.
## (DE_S / DE_R computed only where each cell has >=30 treated, to avoid noise.)
## ===========================================================================
win_de <- WIN %>%
  mutate(nA = nAS + nAR, dA = dAS + dAR, nB = nBS + nBR, dB = dBS + dBR,
         DE   = (dA / nA - dB / nB) * 100,
         propR = (nAR + nBR) / (nAS + nBS + nAR + nBR) * 100,
         DE_S = ifelse(nAS >= 30 & nBS >= 30, (dAS / nAS - dBS / nBS) * 100, NA_real_),
         DE_R = ifelse(nAR >= 30 & nBR >= 30, (dAR / nAR - dBR / nBR) * 100, NA_real_))
s3 <- win_de %>% group_by(alloc, w) %>%
  summarise(DE = mean(DE, na.rm = TRUE), DE_lo = quantile(DE, .025, na.rm = TRUE), DE_hi = quantile(DE, .975, na.rm = TRUE),
            DE_R = mean(DE_R, na.rm = TRUE), DE_S = mean(DE_S, na.rm = TRUE),
            propR = mean(propR, na.rm = TRUE), .groups = "drop") %>% arrange(alloc, w)
cat("\n===== TABLE S3: marginal DE / DE_R / DE_S / %R by time window =====\n")
print(as.data.frame(s3 %>% mutate(across(where(is.numeric), ~round(., 2)))), right = FALSE)

## Whole-trial stratified DE (main-text "stratified direct effect" numbers).
strat <- WIN %>% group_by(real, alloc) %>%
  summarise(across(c(nAS, dAS, nBS, dBS, nAR, dAR, nBR, dBR), sum), .groups = "drop") %>%
  mutate(DE_S = (dAS / nAS - dBS / nBS) * 100, DE_R = (dAR / nAR - dBR / nBR) * 100) %>%
  group_by(alloc) %>%
  summarise(DE_S = mean(DE_S), DE_R = mean(DE_R),
            DE_R_lo = quantile(DE_R, .025), DE_R_hi = quantile(DE_R, .975), .groups = "drop")
cat("\n===== Whole-trial stratified DE (percentage points) =====\n")
print(as.data.frame(strat %>% mutate(across(where(is.numeric), ~round(., 2)))), right = FALSE)

## ===========================================================================
## 6. FIGURES  (all fed from the single run above -> figures match the tables)
## ===========================================================================
if (RUN_FIGURES) {
  strategy_colors  <- c("90/10" = "#E41A1C", "50/50" = "#377EB8")
  infection_colors <- c("Susceptible" = "#FF7F00", "Resistant" = "#984EA3")
  theme_paper <- theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 11),
          legend.position = "bottom", legend.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 12), axis.title = element_text(size = 10))

  ## Cumulative series per realization (pool 3 clusters, then cumsum over time).
  agg <- TS %>% rename(strategy = alloc) %>%
    group_by(real, strategy, times) %>%
    summarise(new_S = sum(new_S), new_R = sum(new_R), dR = sum(dR), dS = sum(dS), D = sum(D), .groups = "drop") %>%
    arrange(real, strategy, times) %>% group_by(real, strategy) %>%
    mutate(cum_S = cumsum(new_S), cum_R = cumsum(new_R), cum_dR = cumsum(dR), cum_dS = cumsum(dS)) %>% ungroup()

  ## --- Figure: cumulative NEW INFECTIONS by strain and strategy (competitive release) ---
  inc <- agg %>% select(real, strategy, times, Susceptible = cum_S, Resistant = cum_R) %>%
    pivot_longer(c(Susceptible, Resistant), names_to = "Infection", values_to = "cum") %>%
    group_by(strategy, Infection, times) %>%
    summarise(m = mean(cum), lo = quantile(cum, .025), hi = quantile(cum, .975), .groups = "drop")
  fig_inc <- ggplot(inc, aes(times, m, color = Infection, fill = Infection)) +
    geom_line(linewidth = 1.1) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2, color = NA) +
    facet_wrap(~strategy, labeller = labeller(strategy = c("50/50" = "50/50 Strategy", "90/10" = "90/10 Strategy"))) +
    scale_color_manual("Infection", values = infection_colors) +
    scale_fill_manual("Infection", values = infection_colors) +
    labs(title = "Cumulative New Infections by Allocation Strategy",
         subtitle = "Mean across realizations; band = 95% simulation interval",
         x = "Timestep", y = "Cumulative new infections") + theme_paper
  ggsave(file.path(OUT_DIR, "fig_cumulative_incidence.png"), fig_inc, width = 8, height = 4.2, dpi = 300)

  ## --- Figure: cumulative DEATHS (R-strain / S-strain / total) by strategy ---
  sumr <- function(col) agg %>% group_by(strategy, times) %>%
    summarise(m = mean(.data[[col]]), lo = quantile(.data[[col]], .025), hi = quantile(.data[[col]], .975), .groups = "drop")
  sR <- sumr("cum_dR"); sS <- sumr("cum_dS"); sD <- sumr("D"); ymax <- max(sR$hi, sS$hi)
  mkpanel <- function(df, faint, ttl, ylab, ylim = NULL) {
    p <- ggplot() +
      geom_line(data = agg, aes(times, .data[[faint]], group = interaction(real, strategy), color = strategy),
                alpha = .07, linewidth = .25) +
      geom_line(data = df, aes(times, m, color = strategy), linewidth = 1.1) +
      geom_ribbon(data = df, aes(times, ymin = lo, ymax = hi, fill = strategy), alpha = .2, color = NA) +
      scale_color_manual("Strategy", values = strategy_colors) +
      scale_fill_manual("Strategy", values = strategy_colors) +
      labs(title = ttl, x = NULL, y = ylab) + theme_paper + theme(legend.position = "none")
    if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
    p
  }
  pA <- mkpanel(sR, "cum_dR", "A. Resistant-Strain Deaths", "Cumulative deaths", c(0, ymax))
  pB <- mkpanel(sS, "cum_dS", "B. Susceptible-Strain Deaths", NULL, c(0, ymax))
  pC <- mkpanel(sD, "D", "C. Total Deaths", "Cumulative deaths") + labs(x = "Timestep") + theme(legend.position = "bottom")
  fig_deaths <- (pA | pB) / pC + plot_layout(heights = c(1, .85)) +
    plot_annotation(title = "Cumulative Deaths by Allocation Strategy",
                    subtitle = "Faint = individual runs; bold = mean; band = 95% SI",
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
  ggsave(file.path(OUT_DIR, "fig_cumulative_deaths.png"), fig_deaths, width = 8, height = 7, dpi = 300)

  ## --- eFigure: time-varying marginal DE vs the falling %R (dual axis) ---
  sf <- max(s3$propR) / max(s3$DE_hi)
  efig <- ggplot(s3, aes(w)) +
    geom_line(aes(y = propR / sf), color = "#B07AA1", linewidth = .9) +
    geom_point(aes(y = propR / sf), color = "#B07AA1", size = 2.2, shape = 17) +
    geom_ribbon(aes(ymin = DE_lo, ymax = DE_hi), alpha = .12, fill = "#1F77B4") +
    geom_line(aes(y = DE), color = "#1F77B4", linewidth = .9) +
    geom_point(aes(y = DE), color = "#1F77B4", size = 2.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = .3) +
    facet_wrap(~alloc, ncol = 2) +
    scale_y_continuous(name = "Marginal direct effect (percentage points)",
                       sec.axis = sec_axis(~ . * sf, name = "Resistant infections (%)")) +
    labs(x = "Time (simulation steps)") + theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold"),
          axis.title.y.left = element_text(color = "#1F77B4", face = "bold"),
          axis.title.y.right = element_text(color = "#B07AA1", face = "bold"))
  ggsave(file.path(OUT_DIR, "fig_DE_over_time.png"), efig, width = 7.5, height = 3.8, dpi = 300)
  cat("\nFigures written to", OUT_DIR, "(cumulative incidence / deaths / DE-over-time)\n")
}

## ===========================================================================
## 7. SENSITIVITY SWEEP  (Tables S1 / S2)   -- optional, set RUN_SENSITIVITY <- TRUE
## ---------------------------------------------------------------------------
## Re-run the trial while varying one parameter at a time; report the 7
## estimands per scenario. Uses the SAME model and the SAME estimator as above.
## ===========================================================================
if (RUN_SENSITIVITY) {
  SIMS_SENS <- 20
  ## Run one scenario (both arms) and return its estimand table.
  run_scenario <- function(beta, tau, mu_d, pi_hi, pi_lo, sims = SIMS_SENS) {
    fin <- list(); k <- 1
    for (a in c("90/10", "50/50")) {
      pi_a <- if (a == "90/10") pi_hi else pi_lo
      for (r in 1:(3 * sims)) {   # 3 clusters x sims replicates per arm
        x <- safe_runsim(pi = pi_a, beta = beta, tau = tau, mu_d = mu_d, timesteps = TIMESTEPS)
        fin[[k]] <- reduce_final(x, ((r - 1) %/% 3) + 1, a); k <- k + 1
      }
    }
    estimands(bind_rows(fin))$table
  }
  scen <- tribble(
    ~label,                              ~beta, ~tau, ~mu_d, ~pi_hi, ~pi_lo,
    "Baseline (main, beta=1.5)",           1.5,  0.7,  0.010,  0.9,    0.5,
    "Transmission beta=1.0 (lower)",       1.0,  0.7,  0.010,  0.9,    0.5,
    "Transmission beta=2.0 (higher)",      2.0,  0.7,  0.010,  0.9,    0.5,
    "Treatment tau=0.5 (lower)",           1.5,  0.5,  0.010,  0.9,    0.5,
    "Treatment tau=0.9 (higher)",          1.5,  0.9,  0.010,  0.9,    0.5,
    "Length of stay mu_d=0.005 (longer)",  1.5,  0.7,  0.005,  0.9,    0.5,
    "Length of stay mu_d=0.020 (shorter)", 1.5,  0.7,  0.020,  0.9,    0.5,
    "Allocation 90/10 vs 30/70 (wide)",    1.5,  0.7,  0.010,  0.9,    0.3,
    "Allocation 60/40 vs 50/50 (narrow)",  1.5,  0.7,  0.010,  0.6,    0.5)
  cat("\n===== TABLES S1/S2: one-way sensitivity =====\n")
  sens <- pmap_dfr(scen, function(label, beta, tau, mu_d, pi_hi, pi_lo) {
    cat("  scenario:", label, "\n")
    run_scenario(beta, tau, mu_d, pi_hi, pi_lo) %>% mutate(scenario = label, .before = 1)
  })
  print(as.data.frame(sens %>% mutate(across(c(pp, lo, hi), ~round(., 2)))), right = FALSE)
  write_csv(sens, file.path(OUT_DIR, "sensitivity_beta1.5.csv"))
}

## ===========================================================================
## 8. IMPERFECT DRUG B  (Table S4)   -- optional, set RUN_EXTENDED <- TRUE
## ---------------------------------------------------------------------------
## Sweep Drug B's efficacy eff_B from 1.0 (main) down to 0.25; at eff_B=0 Drug B
## would equal Drug A on the resistant strain. Shows how the overall effect and
## the death redistribution attenuate as Drug B's advantage shrinks.
## ===========================================================================
if (RUN_EXTENDED) {
  SIMS_EXT <- 100
  run_effB <- function(eff_B, sims = SIMS_EXT) {
    fin <- list(); k <- 1
    for (a in c("90/10", "50/50")) {
      pi_a <- if (a == "90/10") 0.9 else 0.5
      for (r in 1:(3 * sims)) {
        x <- safe_runsim(pi = pi_a, eff_B = eff_B, timesteps = TIMESTEPS)
        fin[[k]] <- reduce_final(x, ((r - 1) %/% 3) + 1, a); k <- k + 1
      }
    }
    est <- estimands(bind_rows(fin))
    oe <- est$table %>% filter(estimand == "OE (alpha0-alpha1)")
    est$arms %>% select(alloc, death_S, death_R) %>%
      pivot_wider(names_from = alloc, values_from = c(death_S, death_R)) %>%
      mutate(eff_B = eff_B, OE_pp = oe$pp, OE_lo = oe$lo, OE_hi = oe$hi, .before = 1)
  }
  cat("\n===== TABLE S4: imperfect Drug B (eff_B sweep) =====\n")
  ext <- map_dfr(c(1.0, 0.75, 0.50, 0.25), function(e) { cat("  eff_B =", e, "\n"); run_effB(e) })
  print(as.data.frame(ext %>% mutate(across(where(is.numeric), ~round(., 2)))), right = FALSE)
  write_csv(ext, file.path(OUT_DIR, "table_s4_beta1.5.csv"))
}

cat("\n=== DONE ===\n")
