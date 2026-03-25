#!/usr/bin/env Rscript
# Reduced simulation run: extract all key estimands for paper
# Skips all plotting, uses same parameters as main code

suppressPackageStartupMessages({
  library(ABM)
  library(tidyverse)
})

cat("=== ABM Simulation: Extracting Estimands for Paper ===\n\n")

nclus <- 6
sims <- 100

runsim <- function(
    gammaR = 0.2, gammaS = 0.2, beta = 0.4, c = 0.01, n = 1000,
    mu = 0.095, delta = 0.025, lambda = 10, m = 0.75,
    S0 = 50, R0 = 50, sigma = 0.5, rho = 0.9, pi = 0.9,
    tau = 0.7, timesteps = 100
) {
  admit.event <- function(time) {
    newEvent(time + rexp(1, lambda), function(time, sim, agent) {
      if (runif(1) < m) {
        a <- newAgent("X")
      } else {
        a <- newAgent("S")
      }
      addAgent(sim, a)
      schedule(sim, admit.event(time))
    })
  }

  sim <- Simulation$new(as.list(c(
    rep("S", S0), rep("R", R0), rep("X", n - S0 - R0)
  )))

  sim$schedule(admit.event(0))

  # Loggers
  sim$addLogger(newCounter("X", "X"))
  sim$addLogger(newCounter("S", "S"))
  sim$addLogger(newCounter("R", "R"))
  sim$addLogger(newCounter("US", "US"))
  sim$addLogger(newCounter("TS1", "TS1"))
  sim$addLogger(newCounter("TS2", "TS2"))
  sim$addLogger(newCounter("TS1_inc", "S", "TS1"))
  sim$addLogger(newCounter("TS2_inc", "S", "TS2"))
  sim$addLogger(newCounter("UR", "UR"))
  sim$addLogger(newCounter("TR1", "TR1"))
  sim$addLogger(newCounter("TR2", "TR2"))
  sim$addLogger(newCounter("TR1_inc", "R", "TR1"))
  sim$addLogger(newCounter("TR2_inc", "R", "TR2"))
  sim$addLogger(newCounter("Z", "Z"))
  sim$addLogger(newCounter("D", "D"))
  sim$addLogger(newCounter("D_TS1", "TS1", "D"))
  sim$addLogger(newCounter("D_TS2", "TS2", "D"))
  sim$addLogger(newCounter("D_TR1", "TR1", "D"))
  sim$addLogger(newCounter("D_TR2", "TR2", "D"))
  sim$addLogger(newCounter("US_inc", "S", "US"))
  sim$addLogger(newCounter("UR_inc", "R", "UR"))
  sim$addLogger(newCounter("D_US", "US", "D"))
  sim$addLogger(newCounter("D_UR", "UR", "D"))
  sim$addLogger(newCounter("D", "D"))

  mx <- newRandomMixing()
  sim$addContact(mx)

  # Symptom onset & treatment

  sim$addTransition("S" -> "TS1", sigma * rho * pi)
  sim$addTransition("S" -> "TS2", sigma * rho * (1 - pi))
  sim$addTransition("S" -> "US", sigma * (1 - rho))
  sim$addTransition("R" -> "TR1", sigma * rho * pi)
  sim$addTransition("R" -> "TR2", sigma * rho * (1 - pi))
  sim$addTransition("R" -> "UR", sigma * (1 - rho))

  # Recovery
  sim$addTransition("TS1" -> "X", gammaS + tau)
  sim$addTransition("TS2" -> "X", gammaS + tau)
  sim$addTransition("US" -> "X", gammaS)
  sim$addTransition("TR1" -> "X", gammaR)
  sim$addTransition("TR2" -> "X", gammaR + tau)
  sim$addTransition("UR" -> "X", gammaR)

  # Discharge alive
  sim$addTransition("TS1" -> "Z", mu)
  sim$addTransition("TS2" -> "Z", mu)
  sim$addTransition("US" -> "Z", mu)
  sim$addTransition("TR1" -> "Z", mu)
  sim$addTransition("TR2" -> "Z", mu)
  sim$addTransition("UR" -> "Z", mu)

  # Death
  sim$addTransition("TS1" -> "D", delta * 1)     # concordant S: best
  sim$addTransition("TS2" -> "D", delta * 1)     # concordant S: best
  sim$addTransition("US" -> "D", delta * 2)      # untreated S: no drug benefit
  sim$addTransition("TR1" -> "D", delta * 2.5)   # discordant R: wrong drug + R
  sim$addTransition("TR2" -> "D", delta * 1.5)   # concordant R: effective but R harder
  sim$addTransition("UR" -> "D", delta * 3)      # untreated R: worst

  # Transmission
  sim$addTransition("S" + "X" -> "S" + "S" ~ mx, beta)
  sim$addTransition("TS1" + "X" -> "TS1" + "S" ~ mx, beta)
  sim$addTransition("TS2" + "X" -> "TS2" + "S" ~ mx, beta)
  sim$addTransition("US" + "X" -> "US" + "S" ~ mx, beta)
  sim$addTransition("R" + "X" -> "R" + "R" ~ mx, beta * (1 - c))
  sim$addTransition("TR1" + "X" -> "TR1" + "R" ~ mx, beta * (1 - c))
  sim$addTransition("TR2" + "X" -> "TR2" + "R" ~ mx, beta * (1 - c))
  sim$addTransition("UR" + "X" -> "UR" + "R" ~ mx, beta * (1 - c))

  x <- sim$run(0:timesteps)
  return(x)
}

# ===== RUN SIMULATIONS =====
cat("Running", nclus, "clusters x", sims, "simulations x 350 timesteps...\n")
t0 <- Sys.time()

res <- list()
for (i in 1:nclus) {
  if (i <= floor(nclus / 2)) {
    pi <- 0.9
  } else {
    pi <- 0.5
  }
  for (j in 1:sims) {
    x <- runsim(pi = pi, timesteps = 350)
    x$clus <- i
    x$sim <- j
    res[[(i - 1) * sims + j]] <- x
  }
  cat("  Cluster", i, "done (pi =", pi, ")\n")
}

t1 <- Sys.time()
cat("Simulation complete in", round(difftime(t1, t0, units = "secs"), 1), "seconds\n\n")

res <- bind_rows(res, .id = "id")

# ===== EXTRACT ALL KEY NUMBERS =====
cat("========================================\n")
cat("EXTRACTING KEY ESTIMANDS\n")
cat("========================================\n\n")

# --- 1) Overall Risk by Treatment (Drug A vs B) ---
rd_cluster_df <- res %>%
  group_by(sim, clus) %>%
  summarise(
    assign_A = sum(TS1_inc + TR1_inc),
    assign_B = sum(TS2_inc + TR2_inc),
    death_A = sum(D_TS1 + D_TR1),
    death_B = sum(D_TS2 + D_TR2),
    .groups = "drop"
  ) %>%
  mutate(
    risk_A = death_A / assign_A,
    risk_B = death_B / assign_B,
    rd_cluster = risk_A - risk_B,
    alloc_grp = if_else(clus <= floor(nclus / 2), "90/10", "50/50")
  )

risk_by_group <- rd_cluster_df %>%
  group_by(sim, alloc_grp) %>%
  summarise(
    total_deaths = sum(death_A + death_B),
    total_assign = sum(assign_A + assign_B),
    risk = total_deaths / total_assign,
    .groups = "drop"
  )

group_summary_risk <- risk_by_group %>%
  group_by(alloc_grp) %>%
  summarise(
    mean_risk = mean(risk),
    se_risk = sd(risk) / sqrt(n()),
    .groups = "drop"
  )

cat("--- Group Risk Summary (treated patients) ---\n")
print(group_summary_risk)

rd_overall <- with(group_summary_risk,
  mean_risk[alloc_grp == "90/10"] - mean_risk[alloc_grp == "50/50"]
)
cat("Overall RD (90/10 - 50/50):", round(rd_overall, 4), "\n\n")

# --- 2) Overall Deaths (all patients, not just treated) ---
end_deaths <- res %>%
  filter(times == max(times)) %>%
  group_by(sim, clus) %>%
  summarise(deaths = sum(D), .groups = "drop") %>%
  mutate(alloc_grp = if_else(clus <= floor(nclus / 2), "90/10", "50/50"))

group_deaths <- end_deaths %>%
  group_by(sim, alloc_grp) %>%
  summarise(total_deaths = sum(deaths), .groups = "drop")

group_summary_deaths <- group_deaths %>%
  group_by(alloc_grp) %>%
  summarise(
    mean_deaths = mean(total_deaths),
    se_deaths = sd(total_deaths) / sqrt(n()),
    .groups = "drop"
  )

cat("--- Total Deaths by Strategy (all patients) ---\n")
print(group_summary_deaths)

rd_overall_deaths <- with(group_summary_deaths,
  mean_deaths[alloc_grp == "90/10"] - mean_deaths[alloc_grp == "50/50"]
)
cat("Difference in mean total deaths (90/10 - 50/50):", round(rd_overall_deaths, 2), "\n\n")

# --- 3) Spillover Effect: Drug B treated ---
spill_df <- res %>%
  group_by(sim, clus) %>%
  summarise(
    assign_B = sum(TS2_inc + TR2_inc),
    death_B = sum(D_TS2 + D_TR2),
    .groups = "drop"
  ) %>%
  filter(assign_B > 0) %>%
  mutate(
    risk_B = death_B / assign_B,
    alloc_grp = if_else(clus <= floor(nclus / 2), "90/10", "50/50")
  )

spill_summary <- spill_df %>%
  group_by(alloc_grp) %>%
  summarise(
    mean_riskB = mean(risk_B),
    se_riskB = sd(risk_B) / sqrt(n()),
    .groups = "drop"
  )

cat("--- Spillover: Risk among B-treated ---\n")
print(spill_summary)

spill_rd <- with(spill_summary,
  mean_riskB[alloc_grp == "90/10"] - mean_riskB[alloc_grp == "50/50"]
)
cat("Spillover RD (B-treated, 90/10 vs 50/50):", round(spill_rd, 4), "\n\n")

# --- 4) Spillover Effect: Drug A treated ---
spill_A_df <- res %>%
  group_by(sim, clus) %>%
  summarise(
    assign_A = sum(TS1_inc + TR1_inc),
    death_A = sum(D_TS1 + D_TR1),
    .groups = "drop"
  ) %>%
  filter(assign_A > 0) %>%
  mutate(
    risk_A = death_A / assign_A,
    alloc_grp = if_else(clus <= floor(nclus / 2), "90/10", "50/50")
  )

spill_summary_A <- spill_A_df %>%
  group_by(alloc_grp) %>%
  summarise(
    mean_riskA = mean(risk_A),
    se_riskA = sd(risk_A) / sqrt(n()),
    .groups = "drop"
  )

cat("--- Spillover: Risk among A-treated ---\n")
print(spill_summary_A)

spill_rd_A <- with(spill_summary_A,
  mean_riskA[alloc_grp == "90/10"] - mean_riskA[alloc_grp == "50/50"]
)
cat("Spillover RD (A-treated, 90/10 vs 50/50):", round(spill_rd_A, 4), "\n\n")

# --- 5) Risk by Infection Type (S vs R) ---
sus_df <- res %>%
  group_by(sim, clus) %>%
  summarise(
    assign_S = sum(TS1_inc + TS2_inc + US_inc),
    death_S = sum(D_TS1 + D_TS2 + D_US),
    .groups = "drop"
  ) %>%
  filter(assign_S > 0) %>%
  mutate(
    risk_S = death_S / assign_S,
    alloc_grp = if_else(clus <= floor(nclus / 2), "90/10", "50/50")
  )

sus_summary <- sus_df %>%
  group_by(alloc_grp) %>%
  summarise(
    mean_riskS = mean(risk_S),
    se_riskS = sd(risk_S) / sqrt(n()),
    .groups = "drop"
  )

res_df <- res %>%
  group_by(sim, clus) %>%
  summarise(
    assign_R = sum(TR1_inc + TR2_inc + UR_inc),
    death_R = sum(D_TR1 + D_TR2 + D_UR),
    .groups = "drop"
  ) %>%
  filter(assign_R > 0) %>%
  mutate(
    risk_R = death_R / assign_R,
    alloc_grp = if_else(clus <= floor(nclus / 2), "90/10", "50/50")
  )

res_summary <- res_df %>%
  group_by(alloc_grp) %>%
  summarise(
    mean_riskR = mean(risk_R),
    se_riskR = sd(risk_R) / sqrt(n()),
    .groups = "drop"
  )

# Simulation intervals for S
rdS_sim <- sus_df %>%
  group_by(sim, alloc_grp) %>%
  summarise(mean_r = mean(risk_S), .groups = "drop") %>%
  pivot_wider(names_from = alloc_grp, values_from = mean_r) %>%
  mutate(rd_S = `90/10` - `50/50`) %>%
  pull(rd_S)

ci_S <- quantile(rdS_sim, c(0.025, 0.975))

# Simulation intervals for R
rdR_sim <- res_df %>%
  group_by(sim, alloc_grp) %>%
  summarise(mean_r = mean(risk_R), .groups = "drop") %>%
  pivot_wider(names_from = alloc_grp, values_from = mean_r) %>%
  mutate(rd_R = `90/10` - `50/50`) %>%
  pull(rd_R)

ci_R <- quantile(rdR_sim, c(0.025, 0.975))

rd_S <- mean(rdS_sim)
rd_R <- mean(rdR_sim)

cat("--- Risk by Infection Type ---\n")
cat("Susceptible infections:\n")
print(sus_summary)
cat("Resistant infections:\n")
print(res_summary)
cat("RD_susceptible (90/10 - 50/50):", round(rd_S, 4),
    " 95% SI [", round(ci_S[1], 4), ",", round(ci_S[2], 4), "]\n")
cat("RD_resistant   (90/10 - 50/50):", round(rd_R, 4),
    " 95% SI [", round(ci_R[1], 4), ",", round(ci_R[2], 4), "]\n\n")

# --- 6) Stratified Spillover (Drug x Infection Type) ---
spill_strat <- res %>%
  group_by(sim, clus) %>%
  summarise(
    n_B_S = sum(TS2_inc), d_B_S = sum(D_TS2),
    n_B_R = sum(TR2_inc), d_B_R = sum(D_TR2),
    n_A_S = sum(TS1_inc), d_A_S = sum(D_TS1),
    n_A_R = sum(TR1_inc), d_A_R = sum(D_TR1),
    .groups = "drop"
  ) %>%
  filter(n_B_S + n_B_R + n_A_S + n_A_R > 0) %>%
  mutate(
    alloc = if_else(clus <= floor(nclus / 2), "90/10", "50/50"),
    risk_B_S = if_else(n_B_S > 0, d_B_S / n_B_S, NA_real_),
    risk_B_R = if_else(n_B_R > 0, d_B_R / n_B_R, NA_real_),
    risk_A_S = if_else(n_A_S > 0, d_A_S / n_A_S, NA_real_),
    risk_A_R = if_else(n_A_R > 0, d_A_R / n_A_R, NA_real_)
  )

spill_summary_strat <- spill_strat %>%
  group_by(alloc) %>%
  summarise(
    mean_B_S = mean(risk_B_S, na.rm = TRUE),
    se_B_S = sd(risk_B_S, na.rm = TRUE) / sqrt(sum(!is.na(risk_B_S))),
    mean_B_R = mean(risk_B_R, na.rm = TRUE),
    se_B_R = sd(risk_B_R, na.rm = TRUE) / sqrt(sum(!is.na(risk_B_R))),
    mean_A_S = mean(risk_A_S, na.rm = TRUE),
    se_A_S = sd(risk_A_S, na.rm = TRUE) / sqrt(sum(!is.na(risk_A_S))),
    mean_A_R = mean(risk_A_R, na.rm = TRUE),
    se_A_R = sd(risk_A_R, na.rm = TRUE) / sqrt(sum(!is.na(risk_A_R))),
    .groups = "drop"
  )

cat("--- Stratified Spillover (Drug x Infection Type) ---\n")
print(spill_summary_strat)

with(spill_summary_strat, {
  rd_B_S <- mean_B_S[alloc == "90/10"] - mean_B_S[alloc == "50/50"]
  rd_B_R <- mean_B_R[alloc == "90/10"] - mean_B_R[alloc == "50/50"]
  rd_A_S <- mean_A_S[alloc == "90/10"] - mean_A_S[alloc == "50/50"]
  rd_A_R <- mean_A_R[alloc == "90/10"] - mean_A_R[alloc == "50/50"]
  cat("Spillover RD among B-treated, Susceptible:", round(rd_B_S, 4), "\n")
  cat("Spillover RD among B-treated, Resistant:  ", round(rd_B_R, 4), "\n")
  cat("Spillover RD among A-treated, Susceptible:", round(rd_A_S, 4), "\n")
  cat("Spillover RD among A-treated, Resistant:  ", round(rd_A_R, 4), "\n\n")
})

# --- 7) Cumulative Risk by Drug x Infection x Strategy ---
final_counts <- res %>%
  group_by(sim, clus) %>%
  summarise(
    tot_SA = sum(TS1_inc), dead_SA = sum(D_TS1),
    tot_SB = sum(TS2_inc), dead_SB = sum(D_TS2),
    tot_RA = sum(TR1_inc), dead_RA = sum(D_TR1),
    tot_RB = sum(TR2_inc), dead_RB = sum(D_TR2),
    .groups = "drop"
  ) %>%
  mutate(alloc = if_else(clus <= floor(nclus / 2), "90/10", "50/50"))

cumrisk_summary <- final_counts %>%
  group_by(alloc) %>%
  summarise(
    risk_SA = mean(dead_SA / tot_SA, na.rm = TRUE),
    se_SA = sd(dead_SA / tot_SA, na.rm = TRUE) / sqrt(sum(!is.na(tot_SA))),
    risk_SB = mean(dead_SB / tot_SB, na.rm = TRUE),
    se_SB = sd(dead_SB / tot_SB, na.rm = TRUE) / sqrt(sum(!is.na(tot_SB))),
    risk_RA = mean(dead_RA / tot_RA, na.rm = TRUE),
    se_RA = sd(dead_RA / tot_RA, na.rm = TRUE) / sqrt(sum(!is.na(tot_RA))),
    risk_RB = mean(dead_RB / tot_RB, na.rm = TRUE),
    se_RB = sd(dead_RB / tot_RB, na.rm = TRUE) / sqrt(sum(!is.na(tot_RB))),
    .groups = "drop"
  )

cat("--- Cumulative Death Risk (Drug x Infection x Strategy) ---\n")
print(cumrisk_summary)

# --- 8) New Infections & Deaths by Type ---
clus_allocs <- tibble(
  clus = 1:nclus,
  alloc_grp = if_else(clus <= floor(nclus / 2), "90/10", "50/50")
)
res2 <- res %>% left_join(clus_allocs, by = "clus")

new_inf_summary <- res2 %>%
  group_by(sim, clus) %>%
  summarise(
    tot_new_S = sum(TS1_inc + TS2_inc + US_inc),
    tot_new_R = sum(TR1_inc + TR2_inc + UR_inc),
    .groups = "drop"
  ) %>%
  left_join(clus_allocs, by = "clus") %>%
  group_by(sim, alloc_grp) %>%
  summarise(
    strategy_new_S = sum(tot_new_S),
    strategy_new_R = sum(tot_new_R),
    .groups = "drop"
  ) %>%
  group_by(alloc_grp) %>%
  summarise(
    mean_new_S = mean(strategy_new_S),
    se_new_S = sd(strategy_new_S) / sqrt(n()),
    mean_new_R = mean(strategy_new_R),
    se_new_R = sd(strategy_new_R) / sqrt(n()),
    .groups = "drop"
  )

cat("--- New Infections by Type & Strategy ---\n")
print(new_inf_summary)

deaths_by_type <- res2 %>%
  group_by(sim, clus) %>%
  summarise(
    deaths_S = sum(D_TS1 + D_TS2 + D_US),
    deaths_R = sum(D_TR1 + D_TR2 + D_UR),
    .groups = "drop"
  ) %>%
  left_join(clus_allocs, by = "clus") %>%
  group_by(sim, alloc_grp) %>%
  summarise(
    strat_deaths_S = sum(deaths_S),
    strat_deaths_R = sum(deaths_R),
    .groups = "drop"
  ) %>%
  group_by(alloc_grp) %>%
  summarise(
    mean_deaths_S = mean(strat_deaths_S),
    se_deaths_S = sd(strat_deaths_S) / sqrt(n()),
    mean_deaths_R = mean(strat_deaths_R),
    se_deaths_R = sd(strat_deaths_R) / sqrt(n()),
    .groups = "drop"
  )

cat("--- Deaths by Infection Type & Strategy ---\n")
print(deaths_by_type)

# --- 9) Baseline totals ---
baseline_totals <- res2 %>%
  filter(times == 0) %>%
  group_by(sim, alloc_grp) %>%
  summarise(
    total_agents = sum(X + S + R),
    total_S_inf = sum(S),
    total_R_inf = sum(R),
    prev_S = total_S_inf / total_agents,
    prev_R = total_R_inf / total_agents,
    .groups = "drop"
  ) %>%
  group_by(alloc_grp) %>%
  summarise(
    mean_agents = mean(total_agents),
    mean_S_inf = mean(total_S_inf),
    mean_R_inf = mean(total_R_inf),
    mean_prev_S = mean(prev_S),
    mean_prev_R = mean(prev_R),
    .groups = "drop"
  )

cat("\n--- Baseline Totals ---\n")
print(baseline_totals)

# --- 10) Agent counts per strategy ---
agent_counts_strategy <- res2 %>%
  group_by(sim, clus) %>%
  summarise(
    discharges_cluster = max(Z),
    agents_cluster = 1000 + discharges_cluster,
    deaths_cluster = max(D),
    .groups = "drop"
  ) %>%
  left_join(clus_allocs, by = "clus") %>%
  group_by(sim, alloc_grp) %>%
  summarise(
    total_agents_strategy = sum(agents_cluster),
    total_deaths_strategy = sum(deaths_cluster),
    .groups = "drop"
  ) %>%
  group_by(alloc_grp) %>%
  summarise(
    mean_total_agents = mean(total_agents_strategy),
    sd_total_agents = sd(total_agents_strategy),
    mean_deaths = mean(total_deaths_strategy),
    sd_deaths = sd(total_deaths_strategy),
    .groups = "drop"
  )

cat("\n--- Agent Counts & Deaths per Strategy ---\n")
print(agent_counts_strategy)

# --- 11) Compute formal causal estimands for Table 3 ---
cat("\n========================================\n")
cat("FORMAL CAUSAL ESTIMANDS (Table 3)\n")
cat("========================================\n\n")

# DE at 90/10: E[Y(A, alpha_1)] - E[Y(B, alpha_1)]
de_90 <- with(cumrisk_summary, {
  # Average over S and R for each drug
  # Weight by number of assignments
  list(
    risk_A_90 = (risk_SA[alloc == "90/10"] + risk_RA[alloc == "90/10"]) / 2,
    risk_B_90 = (risk_SB[alloc == "90/10"] + risk_RB[alloc == "90/10"]) / 2
  )
})

# Use the overall risk by drug within each strategy
de_data <- rd_cluster_df %>%
  group_by(sim, alloc_grp) %>%
  summarise(
    mean_riskA = mean(risk_A),
    mean_riskB = mean(risk_B),
    .groups = "drop"
  )

de_90_sim <- de_data %>%
  filter(alloc_grp == "90/10") %>%
  mutate(de = mean_riskA - mean_riskB) %>%
  pull(de)

de_50_sim <- de_data %>%
  filter(alloc_grp == "50/50") %>%
  mutate(de = mean_riskA - mean_riskB) %>%
  pull(de)

cat("DE(alpha_1 = 90/10): ", round(mean(de_90_sim), 4),
    "  95% SI [", round(quantile(de_90_sim, 0.025), 4), ",",
    round(quantile(de_90_sim, 0.975), 4), "]\n")

cat("DE(alpha_0 = 50/50): ", round(mean(de_50_sim), 4),
    "  95% SI [", round(quantile(de_50_sim, 0.025), 4), ",",
    round(quantile(de_50_sim, 0.975), 4), "]\n")

# IE for Drug A: E[Y(A, alpha_1)] - E[Y(A, alpha_0)]
ie_A_sim <- de_data %>%
  select(sim, alloc_grp, mean_riskA) %>%
  pivot_wider(names_from = alloc_grp, values_from = mean_riskA) %>%
  mutate(ie = `90/10` - `50/50`) %>%
  pull(ie)

cat("IE(Drug A):           ", round(mean(ie_A_sim), 4),
    "  95% SI [", round(quantile(ie_A_sim, 0.025), 4), ",",
    round(quantile(ie_A_sim, 0.975), 4), "]\n")

# IE for Drug B: E[Y(B, alpha_1)] - E[Y(B, alpha_0)]
ie_B_sim <- de_data %>%
  select(sim, alloc_grp, mean_riskB) %>%
  pivot_wider(names_from = alloc_grp, values_from = mean_riskB) %>%
  mutate(ie = `90/10` - `50/50`) %>%
  pull(ie)

cat("IE(Drug B):           ", round(mean(ie_B_sim), 4),
    "  95% SI [", round(quantile(ie_B_sim, 0.025), 4), ",",
    round(quantile(ie_B_sim, 0.975), 4), "]\n")

# TE: E[Y(A, alpha_1)] - E[Y(B, alpha_0)]
te_sim <- de_data %>%
  select(sim, alloc_grp, mean_riskA, mean_riskB) %>%
  pivot_wider(names_from = alloc_grp, values_from = c(mean_riskA, mean_riskB)) %>%
  mutate(te = `mean_riskA_90/10` - `mean_riskB_50/50`) %>%
  pull(te)

cat("TE:                   ", round(mean(te_sim), 4),
    "  95% SI [", round(quantile(te_sim, 0.025), 4), ",",
    round(quantile(te_sim, 0.975), 4), "]\n")

# OE: E[Y(alpha_1)] - E[Y(alpha_0)]  (marginal over drug assignment)
oe_sim <- risk_by_group %>%
  select(sim, alloc_grp, risk) %>%
  pivot_wider(names_from = alloc_grp, values_from = risk) %>%
  mutate(oe = `90/10` - `50/50`) %>%
  pull(oe)

cat("OE:                   ", round(mean(oe_sim), 4),
    "  95% SI [", round(quantile(oe_sim, 0.025), 4), ",",
    round(quantile(oe_sim, 0.975), 4), "]\n")

# TE_2: E[Y(B, alpha_1)] - E[Y(A, alpha_0)]
te2_sim <- de_data %>%
  select(sim, alloc_grp, mean_riskA, mean_riskB) %>%
  pivot_wider(names_from = alloc_grp, values_from = c(mean_riskA, mean_riskB)) %>%
  mutate(te2 = `mean_riskB_90/10` - `mean_riskA_50/50`) %>%
  pull(te2)

cat("TE_2:                 ", round(mean(te2_sim), 4),
    "  95% SI [", round(quantile(te2_sim, 0.025), 4), ",",
    round(quantile(te2_sim, 0.975), 4), "]\n")

# Verify decomposition: TE = DE(alpha_1) + IE(B)
cat("\n--- Decomposition Check ---\n")
cat("TE (computed):               ", round(mean(te_sim), 4), "\n")
cat("DE(alpha_1) + IE(B):         ", round(mean(de_90_sim) + mean(ie_B_sim), 4), "\n")
cat("IE(A) + DE(alpha_0):         ", round(mean(ie_A_sim) + mean(de_50_sim), 4), "\n")
cat("TE_2 (computed):             ", round(mean(te2_sim), 4), "\n")
cat("IE(A) - DE(alpha_1):         ", round(mean(ie_A_sim) - mean(de_90_sim), 4), "\n")
cat("IE(B) - DE(alpha_0):         ", round(mean(ie_B_sim) - mean(de_50_sim), 4), "\n")

# --- 12) Stratified DE: DE among S-infected and R-infected separately ---
cat("\n========================================\n")
cat("STRATIFIED DE (by infection type)\n")
cat("========================================\n\n")

# DE_S(alpha): risk(Drug A | S-infected, alpha) - risk(Drug B | S-infected, alpha)
# DE_R(alpha): risk(Drug A | R-infected, alpha) - risk(Drug B | R-infected, alpha)
de_strat_data <- spill_strat %>%
  group_by(sim, alloc) %>%
  summarise(
    mean_risk_A_S = mean(risk_A_S, na.rm = TRUE),
    mean_risk_B_S = mean(risk_B_S, na.rm = TRUE),
    mean_risk_A_R = mean(risk_A_R, na.rm = TRUE),
    mean_risk_B_R = mean(risk_B_R, na.rm = TRUE),
    .groups = "drop"
  )

de_S_90_sim <- de_strat_data %>%
  filter(alloc == "90/10") %>%
  mutate(de = mean_risk_A_S - mean_risk_B_S) %>%
  pull(de)

de_S_50_sim <- de_strat_data %>%
  filter(alloc == "50/50") %>%
  mutate(de = mean_risk_A_S - mean_risk_B_S) %>%
  pull(de)

de_R_90_sim <- de_strat_data %>%
  filter(alloc == "90/10") %>%
  mutate(de = mean_risk_A_R - mean_risk_B_R) %>%
  pull(de)

de_R_50_sim <- de_strat_data %>%
  filter(alloc == "50/50") %>%
  mutate(de = mean_risk_A_R - mean_risk_B_R) %>%
  pull(de)

cat("DE_S(alpha_1 = 90/10): ", round(mean(de_S_90_sim, na.rm=TRUE)*100, 2), " pp",
    "  95% SI [", round(quantile(de_S_90_sim, 0.025, na.rm=TRUE)*100, 2), ",",
    round(quantile(de_S_90_sim, 0.975, na.rm=TRUE)*100, 2), "]\n")
cat("DE_S(alpha_0 = 50/50): ", round(mean(de_S_50_sim, na.rm=TRUE)*100, 2), " pp",
    "  95% SI [", round(quantile(de_S_50_sim, 0.025, na.rm=TRUE)*100, 2), ",",
    round(quantile(de_S_50_sim, 0.975, na.rm=TRUE)*100, 2), "]\n")
cat("DE_R(alpha_1 = 90/10): ", round(mean(de_R_90_sim, na.rm=TRUE)*100, 2), " pp",
    "  95% SI [", round(quantile(de_R_90_sim, 0.025, na.rm=TRUE)*100, 2), ",",
    round(quantile(de_R_90_sim, 0.975, na.rm=TRUE)*100, 2), "]\n")
cat("DE_R(alpha_0 = 50/50): ", round(mean(de_R_50_sim, na.rm=TRUE)*100, 2), " pp",
    "  95% SI [", round(quantile(de_R_50_sim, 0.025, na.rm=TRUE)*100, 2), ",",
    round(quantile(de_R_50_sim, 0.975, na.rm=TRUE)*100, 2), "]\n")

cat("\nInterpretation: DE_S should be small (both drugs cover S).\n")
cat("DE_R should be large (Drug A does not cover R, Drug B does).\n")
cat("The marginal DE is a mixture of DE_S and DE_R, weighted by strain prevalence.\n")

cat("\n=== DONE ===\n")
