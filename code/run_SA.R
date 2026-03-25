#!/usr/bin/env Rscript
# Sensitivity analyses for the paper
suppressPackageStartupMessages({
  library(ABM)
  library(tidyverse)
})

cat("=== SENSITIVITY ANALYSES ===\n\n")

# Copy the runsim function from the main code
runsim <- function(
    gammaR = 0.2, gammaS = 0.2, beta = 0.4, c = 0.02, n = 1000,
    mu = 0.095, delta = 0.025, lambda = 10, m = 0.75,
    S0 = 50, R0 = 50, sigma = 1, rho = 0.9, pi = 0.9,
    tau = 0.7, timesteps = 100
) {
  admit.event <- function(time) {
    newEvent(time + rexp(1, lambda), function(time, sim, agent) {
      if (runif(1) < m) a <- newAgent("X") else a <- newAgent("S")
      addAgent(sim, a)
      schedule(sim, admit.event(time))
    })
  }
  sim <- Simulation$new(as.list(c(rep("S", S0), rep("R", R0), rep("X", n - S0 - R0))))
  sim$schedule(admit.event(0))
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
  sim$addTransition("S" -> "TS1", sigma * rho * pi)
  sim$addTransition("S" -> "TS2", sigma * rho * (1 - pi))
  sim$addTransition("S" -> "US", sigma * (1 - rho))
  sim$addTransition("R" -> "TR1", sigma * rho * pi)
  sim$addTransition("R" -> "TR2", sigma * rho * (1 - pi))
  sim$addTransition("R" -> "UR", sigma * (1 - rho))
  sim$addTransition("TS1" -> "X", gammaS + tau)
  sim$addTransition("TS2" -> "X", gammaS + tau)
  sim$addTransition("US" -> "X", gammaS)
  sim$addTransition("TR1" -> "X", gammaR)
  sim$addTransition("TR2" -> "X", gammaR + tau)
  sim$addTransition("UR" -> "X", gammaR)
  sim$addTransition("TS1" -> "Z", mu)
  sim$addTransition("TS2" -> "Z", mu)
  sim$addTransition("US" -> "Z", mu)
  sim$addTransition("TR1" -> "Z", mu)
  sim$addTransition("TR2" -> "Z", mu)
  sim$addTransition("UR" -> "Z", mu)
  sim$addTransition("TS1" -> "D", delta * 1)
  sim$addTransition("TS2" -> "D", delta * 1)
  sim$addTransition("US" -> "D", delta * 1.5)
  sim$addTransition("TR1" -> "D", delta * 2)
  sim$addTransition("TR2" -> "D", delta * 1)
  sim$addTransition("UR" -> "D", delta * 3)
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

# Helper: run a full scenario and return key estimands
run_scenario <- function(nclus = 6, sims = 10, ...) {
  res <- list()
  for (i in 1:nclus) {
    pi_val <- if (i <= floor(nclus / 2)) 0.9 else 0.5
    for (j in 1:sims) {
      x <- runsim(pi = pi_val, timesteps = 350, ...)
      x$clus <- i
      x$sim <- j
      res[[(i - 1) * sims + j]] <- x
    }
  }
  res <- bind_rows(res, .id = "id")

  # Compute key estimands
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
      alloc_grp = if_else(clus <= floor(nclus / 2), "90/10", "50/50")
    )

  de_data <- rd_cluster_df %>%
    group_by(sim, alloc_grp) %>%
    summarise(mean_riskA = mean(risk_A), mean_riskB = mean(risk_B), .groups = "drop")

  # IE(A)
  ie_A_sim <- de_data %>%
    select(sim, alloc_grp, mean_riskA) %>%
    pivot_wider(names_from = alloc_grp, values_from = mean_riskA) %>%
    mutate(ie = `90/10` - `50/50`) %>%
    pull(ie)

  # OE
  risk_by_group <- rd_cluster_df %>%
    group_by(sim, alloc_grp) %>%
    summarise(risk = sum(death_A + death_B) / sum(assign_A + assign_B), .groups = "drop")

  oe_sim <- risk_by_group %>%
    select(sim, alloc_grp, risk) %>%
    pivot_wider(names_from = alloc_grp, values_from = risk) %>%
    mutate(oe = `50/50` - `90/10`) %>%
    pull(oe)

  # Total deaths
  clus_allocs <- tibble(clus = 1:nclus, alloc_grp = if_else(clus <= floor(nclus / 2), "90/10", "50/50"))
  res2 <- res %>% left_join(clus_allocs, by = "clus")

  deaths_by_type <- res2 %>%
    group_by(sim, clus) %>%
    summarise(deaths_S = sum(D_TS1 + D_TS2 + D_US), deaths_R = sum(D_TR1 + D_TR2 + D_UR), .groups = "drop") %>%
    left_join(clus_allocs, by = "clus") %>%
    group_by(sim, alloc_grp) %>%
    summarise(strat_deaths_S = sum(deaths_S), strat_deaths_R = sum(deaths_R), .groups = "drop") %>%
    group_by(alloc_grp) %>%
    summarise(
      mean_deaths_S = mean(strat_deaths_S),
      mean_deaths_R = mean(strat_deaths_R),
      .groups = "drop"
    )

  list(
    IE_A = round(mean(ie_A_sim) * 100, 2),
    IE_A_si = round(quantile(ie_A_sim, c(0.025, 0.975)) * 100, 2),
    OE = round(mean(oe_sim) * 100, 2),
    OE_si = round(quantile(oe_sim, c(0.025, 0.975)) * 100, 2),
    deaths = deaths_by_type
  )
}

# ============ SA 1: Treatment Effect Magnitude (tau) ============
cat("--- SA 1: Treatment Effect Magnitude (tau) ---\n")
for (tau_val in c(0.3, 0.7, 1.0)) {
  cat("\ntau =", tau_val, "\n")
  t0 <- Sys.time()
  result <- run_scenario(tau = tau_val)
  t1 <- Sys.time()
  cat("  IE(A):", result$IE_A, "pp  SI [", result$IE_A_si[1], ",", result$IE_A_si[2], "]\n")
  cat("  OE:   ", result$OE, "pp  SI [", result$OE_si[1], ",", result$OE_si[2], "]\n")
  cat("  Deaths:\n")
  print(result$deaths)
  cat("  Time:", round(difftime(t1, t0, units = "secs"), 1), "s\n")
}

# ============ SA 2: Transmission Rate (beta) ============
cat("\n--- SA 2: Transmission Rate (beta) ---\n")
for (beta_val in c(0.2, 0.4, 0.6)) {
  cat("\nbeta =", beta_val, "\n")
  t0 <- Sys.time()
  result <- run_scenario(beta = beta_val)
  t1 <- Sys.time()
  cat("  IE(A):", result$IE_A, "pp  SI [", result$IE_A_si[1], ",", result$IE_A_si[2], "]\n")
  cat("  OE:   ", result$OE, "pp  SI [", result$OE_si[1], ",", result$OE_si[2], "]\n")
  cat("  Deaths:\n")
  print(result$deaths)
  cat("  Time:", round(difftime(t1, t0, units = "secs"), 1), "s\n")
}

# ============ SA 3: Cluster Configuration ============
cat("\n--- SA 3: Cluster Configuration ---\n")
for (nclus_val in c(6, 12)) {
  cat("\nnclus =", nclus_val, "\n")
  t0 <- Sys.time()
  result <- run_scenario(nclus = nclus_val, sims = 5)  # reduce sims for 12 clusters
  t1 <- Sys.time()
  cat("  IE(A):", result$IE_A, "pp  SI [", result$IE_A_si[1], ",", result$IE_A_si[2], "]\n")
  cat("  OE:   ", result$OE, "pp  SI [", result$OE_si[1], ",", result$OE_si[2], "]\n")
  cat("  Time:", round(difftime(t1, t0, units = "secs"), 1), "s\n")
}

cat("\n=== SA COMPLETE ===\n")
