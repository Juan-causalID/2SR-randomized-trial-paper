#!/usr/bin/env Rscript
# Generate publication-quality figures for the 2SR ABM paper
# Saves figures as PDF files for LaTeX inclusion

suppressPackageStartupMessages({
  library(ABM)
  library(tidyverse)
  library(patchwork)
  library(scales)
})

cat("=== Generating Figures for Paper ===\n\n")

# ── Consistent theme for all figures ────────────────────────────────────────
theme_paper <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    axis.title = element_text(size = 10)
  )

# Color palettes
strategy_colors <- c("90/10" = "#E41A1C", "50/50" = "#377EB8")
infection_colors <- c("Susceptible" = "#FF7F00", "Resistant" = "#984EA3")
drug_colors <- c("Drug A" = "#E41A1C", "Drug B" = "#377EB8")

# ── Run simulation ──────────────────────────────────────────────────────────
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
  sim$addTransition("TS1" -> "D", delta * 1)     # concordant S
  sim$addTransition("TS2" -> "D", delta * 1)     # concordant S
  sim$addTransition("US" -> "D", delta * 2)      # untreated S
  sim$addTransition("TR1" -> "D", delta * 2.5)   # discordant R
  sim$addTransition("TR2" -> "D", delta * 1.5)   # concordant R
  sim$addTransition("UR" -> "D", delta * 3)      # untreated R
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

cat("Running simulations...\n")
t0 <- Sys.time()
res <- list()
for (i in 1:nclus) {
  pi_val <- if (i <= floor(nclus / 2)) 0.9 else 0.5
  for (j in 1:sims) {
    x <- runsim(pi = pi_val, timesteps = 350)
    x$clus <- i
    x$sim <- j
    res[[(i - 1) * sims + j]] <- x
  }
  cat("  Cluster", i, "done\n")
}
t1 <- Sys.time()
cat("Simulations done in", round(difftime(t1, t0, units = "secs"), 1), "s\n\n")

res <- bind_rows(res, .id = "id")

# Label allocation strategy
res <- res %>%
  mutate(strategy = if_else(clus <= floor(nclus / 2), "90/10", "50/50"))


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Infection Dynamics — Prevalence of S, R, X by Strategy
# ════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 2: Infection dynamics...\n")

# Aggregate all S and R compartments
dynamics_df <- res %>%
  mutate(
    all_R = R + TR1 + TR2 + UR,
    all_S = S + TS1 + TS2 + US,
    all_X = X
  ) %>%
  select(sim, times, strategy, all_R, all_S, all_X) %>%
  group_by(sim, strategy, times) %>%
  summarise(
    Susceptible = sum(all_S),
    Resistant   = sum(all_R),
    Uninfected  = sum(all_X),
    .groups = "drop"
  ) %>%
  pivot_longer(c(Susceptible, Resistant, Uninfected),
               names_to = "Compartment", values_to = "Count")

# Summary
dynamics_summary <- dynamics_df %>%
  group_by(strategy, times, Compartment) %>%
  summarise(
    mean_count = mean(Count),
    se_count   = sd(Count) / sqrt(n()),
    .groups = "drop"
  )

fig2 <- ggplot(dynamics_summary,
       aes(x = times, y = mean_count, color = Compartment, fill = Compartment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean_count - se_count, ymax = mean_count + se_count),
              alpha = 0.15, color = NA) +
  facet_wrap(~ strategy, labeller = labeller(
    strategy = c("50/50" = "50/50 Strategy (\u03B1\u2080)",
                 "90/10" = "90/10 Strategy (\u03B1\u2081)")
  )) +
  scale_color_manual(values = c("Susceptible" = "#FF7F00", "Resistant" = "#984EA3",
                                "Uninfected" = "#4DAF4A")) +
  scale_fill_manual(values = c("Susceptible" = "#FF7F00", "Resistant" = "#984EA3",
                               "Uninfected" = "#4DAF4A")) +
  labs(
    title = "Infection Dynamics by Allocation Strategy",
    subtitle = "Mean across 100 simulation replicates per cluster (3 clusters per arm); ribbon = \u00B1 1 SE",
    x = "Timestep",
    y = "Number of agents (across 3 clusters)",
    color = "Compartment",
    fill  = "Compartment"
  ) +
  theme_paper

ggsave("fig2_dynamics.pdf", fig2, width = 8, height = 4.5)
ggsave("fig2_dynamics.png", fig2, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig2_dynamics.pdf/.png\n")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Cumulative Deaths — Combined panel (Resistant / Susceptible / Total)
# ════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 3: Cumulative deaths...\n")

# Resistant deaths
resistant_deaths_ts <- res %>%
  group_by(sim, strategy, times) %>%
  summarise(deaths_R = sum(D_TR1 + D_TR2 + D_UR), .groups = "drop")

resistant_cum_ts <- resistant_deaths_ts %>%
  arrange(sim, strategy, times) %>%
  group_by(sim, strategy) %>%
  mutate(cum_R = cumsum(deaths_R)) %>%
  ungroup()

resistant_cum_summary <- resistant_cum_ts %>%
  group_by(strategy, times) %>%
  summarise(
    mean_cum = mean(cum_R),
    se_cum   = sd(cum_R) / sqrt(n()),
    .groups  = "drop"
  )

# Susceptible deaths
susceptible_deaths_ts <- res %>%
  group_by(sim, strategy, times) %>%
  summarise(deaths_S = sum(D_TS1 + D_TS2 + D_US), .groups = "drop")

susceptible_cum_ts <- susceptible_deaths_ts %>%
  arrange(sim, strategy, times) %>%
  group_by(sim, strategy) %>%
  mutate(cum_S = cumsum(deaths_S)) %>%
  ungroup()

susceptible_cum_summary <- susceptible_cum_ts %>%
  group_by(strategy, times) %>%
  summarise(
    mean_cum = mean(cum_S),
    se_cum   = sd(cum_S) / sqrt(n()),
    .groups  = "drop"
  )

# Total deaths — D logger is already cumulative (absorbing state), no cumsum needed
deaths_cum_ts <- res %>%
  group_by(sim, strategy, times) %>%
  summarise(cum_d = sum(D), .groups = "drop")

deaths_cum_summary <- deaths_cum_ts %>%
  group_by(strategy, times) %>%
  summarise(
    mean_cum = mean(cum_d),
    se_cum   = sd(cum_d) / sqrt(n()),
    .groups  = "drop"
  )

# Common y-axis for top panels
y_max_top <- max(
  resistant_cum_summary$mean_cum + resistant_cum_summary$se_cum,
  susceptible_cum_summary$mean_cum + susceptible_cum_summary$se_cum
)

# Panel A: Resistant
pA <- ggplot(resistant_cum_summary, aes(x = times, y = mean_cum,
                                         color = strategy, fill = strategy)) +
  geom_line(data = resistant_cum_ts, aes(y = cum_R, group = interaction(sim, strategy)),
            alpha = 0.12, linewidth = 0.3) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = mean_cum - se_cum, ymax = mean_cum + se_cum),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = strategy_colors) +
  scale_fill_manual(values = strategy_colors) +
  scale_y_continuous(limits = c(0, y_max_top), expand = c(0, 5)) +
  labs(title = "A. Resistant-Strain Deaths", x = NULL, y = "Cumulative deaths") +
  theme_paper + theme(legend.position = "none")

# Panel B: Susceptible
pB <- ggplot(susceptible_cum_summary, aes(x = times, y = mean_cum,
                                           color = strategy, fill = strategy)) +
  geom_line(data = susceptible_cum_ts, aes(y = cum_S, group = interaction(sim, strategy)),
            alpha = 0.12, linewidth = 0.3) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = mean_cum - se_cum, ymax = mean_cum + se_cum),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = strategy_colors) +
  scale_fill_manual(values = strategy_colors) +
  scale_y_continuous(limits = c(0, y_max_top), expand = c(0, 5)) +
  labs(title = "B. Susceptible-Strain Deaths", x = NULL, y = NULL) +
  theme_paper + theme(legend.position = "none")

# Panel C: Total
pC <- ggplot(deaths_cum_summary, aes(x = times, y = mean_cum,
                                      color = strategy, fill = strategy)) +
  geom_line(data = deaths_cum_ts, aes(y = cum_d, group = interaction(sim, strategy)),
            alpha = 0.12, linewidth = 0.3) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = mean_cum - se_cum, ymax = mean_cum + se_cum),
              alpha = 0.2, color = NA) +
  scale_color_manual("Strategy", values = strategy_colors) +
  scale_fill_manual("Strategy", values = strategy_colors) +
  labs(title = "C. Total Deaths", x = "Timestep", y = "Cumulative deaths") +
  theme_paper

fig3 <- (pA | pB) / pC +
  plot_layout(heights = c(1, 0.85)) +
  plot_annotation(
    title = "Cumulative Deaths by Allocation Strategy",
    subtitle = "Faint lines = individual simulation runs; bold = mean across replicates; ribbon = \u00B1 1 SE",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )

ggsave("fig3_cumulative_deaths.pdf", fig3, width = 8, height = 7, )
ggsave("fig3_cumulative_deaths.png", fig3, width = 8, height = 7, dpi = 300)
cat("  Saved fig3_cumulative_deaths.pdf/.png\n")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 4: Cumulative Incidence by Infection Type and Strategy
# ════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 4: Cumulative incidence...\n")

inc_S_ts <- res %>%
  group_by(sim, strategy, times) %>%
  summarise(new_S = sum(TS1_inc + TS2_inc + US_inc), .groups = "drop")

inc_R_ts <- res %>%
  group_by(sim, strategy, times) %>%
  summarise(new_R = sum(TR1_inc + TR2_inc + UR_inc), .groups = "drop")

cum_S_ts <- inc_S_ts %>%
  arrange(sim, strategy, times) %>%
  group_by(sim, strategy) %>%
  mutate(cum = cumsum(new_S)) %>%
  ungroup() %>%
  mutate(Infection = "Susceptible")

cum_R_ts <- inc_R_ts %>%
  arrange(sim, strategy, times) %>%
  group_by(sim, strategy) %>%
  mutate(cum = cumsum(new_R)) %>%
  ungroup() %>%
  mutate(Infection = "Resistant")

cum_all <- bind_rows(cum_S_ts %>% select(sim, strategy, times, cum, Infection),
                     cum_R_ts %>% select(sim, strategy, times, cum, Infection))

cum_summary <- cum_all %>%
  group_by(strategy, Infection, times) %>%
  summarise(
    mean_cum = mean(cum),
    se_cum   = sd(cum) / sqrt(n()),
    .groups = "drop"
  )

fig4 <- ggplot() +
  geom_line(data = cum_all,
            aes(x = times, y = cum, group = interaction(sim, Infection),
                color = Infection),
            alpha = 0.1, linewidth = 0.3) +
  geom_line(data = cum_summary,
            aes(x = times, y = mean_cum, color = Infection),
            linewidth = 1.1) +
  geom_ribbon(data = cum_summary,
              aes(x = times, ymin = mean_cum - se_cum, ymax = mean_cum + se_cum,
                  fill = Infection),
              alpha = 0.15, color = NA) +
  facet_wrap(~ strategy, labeller = labeller(
    strategy = c("50/50" = "50/50 Strategy (\u03B1\u2080)",
                 "90/10" = "90/10 Strategy (\u03B1\u2081)")
  )) +
  scale_color_manual(values = infection_colors) +
  scale_fill_manual(values = infection_colors) +
  labs(
    title = "Cumulative New Infections by Allocation Strategy",
    subtitle = "Faint lines = individual simulation runs; bold = mean across replicates; ribbon = \u00B1 1 SE",
    x = "Timestep",
    y = "Cumulative new infections (across 3 clusters)",
    color = "Infection Type",
    fill  = "Infection Type"
  ) +
  theme_paper

ggsave("fig4_cumulative_incidence.pdf", fig4, width = 8, height = 4.5, )
ggsave("fig4_cumulative_incidence.png", fig4, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig4_cumulative_incidence.pdf/.png\n")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 5: Death Risk Over Time by Drug — Direct + Indirect Effects
# ════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 5: Death risk by drug...\n")

ab_ts <- res %>%
  group_by(sim, strategy, times) %>%
  summarise(
    new_A = sum(TS1_inc + TR1_inc),
    d_A   = sum(D_TS1 + D_TR1),
    new_B = sum(TS2_inc + TR2_inc),
    d_B   = sum(D_TS2 + D_TR2),
    .groups = "drop"
  )

ab_cum <- ab_ts %>%
  arrange(sim, strategy, times) %>%
  group_by(sim, strategy) %>%
  mutate(
    cumA   = cumsum(new_A),
    cumdA  = cumsum(d_A),
    riskA  = cumdA / cumA,
    cumB   = cumsum(new_B),
    cumdB  = cumsum(d_B),
    riskB  = cumdB / cumB
  ) %>%
  ungroup()

cum_risk_summary <- ab_cum %>%
  group_by(strategy, times) %>%
  summarise(
    mean_riskA = mean(riskA, na.rm = TRUE),
    se_riskA   = sd(riskA, na.rm = TRUE) / sqrt(n()),
    mean_riskB = mean(riskB, na.rm = TRUE),
    se_riskB   = sd(riskB, na.rm = TRUE) / sqrt(n()),
    .groups    = "drop"
  ) %>%
  pivot_longer(
    cols = c(mean_riskA, mean_riskB, se_riskA, se_riskB),
    names_to = c(".value", "Drug"),
    names_pattern = "(mean_risk|se_risk)(A|B)"
  ) %>%
  mutate(Drug = recode(Drug, A = "Drug A", B = "Drug B"))

fig5 <- ggplot(cum_risk_summary %>% filter(times > 5),
       aes(x = times, y = mean_risk, color = strategy, fill = strategy)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean_risk - se_risk, ymax = mean_risk + se_risk),
              alpha = 0.15, color = NA) +
  facet_wrap(~ Drug) +
  scale_color_manual("Strategy", values = strategy_colors) +
  scale_fill_manual("Strategy", values = strategy_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "Cumulative Mortality Risk Over Time by Drug and Strategy",
    subtitle = "Vertical gap within a panel = indirect effect (IE); horizontal gap across panels = direct effect (DE)",
    x = "Timestep",
    y = "Cumulative death risk (deaths \u00F7 assignments)"
  ) +
  theme_paper

ggsave("fig5_risk_by_drug.pdf", fig5, width = 8, height = 4.5, )
ggsave("fig5_risk_by_drug.png", fig5, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig5_risk_by_drug.pdf/.png\n")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 6: Death Risk Stratified by Drug × Infection Type (2×2 grid)
# ════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 6: Stratified death risk...\n")

strat_df <- res %>%
  group_by(sim, strategy, times) %>%
  summarise(
    assn_A_S = sum(TS1_inc), death_A_S = sum(D_TS1),
    assn_B_S = sum(TS2_inc), death_B_S = sum(D_TS2),
    assn_A_R = sum(TR1_inc), death_A_R = sum(D_TR1),
    assn_B_R = sum(TR2_inc), death_B_R = sum(D_TR2),
    .groups = "drop"
  ) %>%
  arrange(sim, strategy, times) %>%
  group_by(sim, strategy) %>%
  mutate(
    cum_A_S = cumsum(assn_A_S), cum_dA_S = cumsum(death_A_S),
    cum_B_S = cumsum(assn_B_S), cum_dB_S = cumsum(death_B_S),
    cum_A_R = cumsum(assn_A_R), cum_dA_R = cumsum(death_A_R),
    cum_B_R = cumsum(assn_B_R), cum_dB_R = cumsum(death_B_R),
    risk_A_S = cum_dA_S / cum_A_S,
    risk_B_S = cum_dB_S / cum_B_S,
    risk_A_R = cum_dA_R / cum_A_R,
    risk_B_R = cum_dB_R / cum_B_R
  ) %>%
  ungroup()

strat_summary <- strat_df %>%
  group_by(strategy, times) %>%
  summarise(
    mean_A_S = mean(risk_A_S, na.rm = TRUE),
    mean_B_S = mean(risk_B_S, na.rm = TRUE),
    mean_A_R = mean(risk_A_R, na.rm = TRUE),
    mean_B_R = mean(risk_B_R, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = c("Drug", "Infection"),
    names_pattern = "mean_([AB])_([SR])",
    values_to = "risk"
  ) %>%
  mutate(
    Drug = recode(Drug, A = "Drug A", B = "Drug B"),
    Infection = recode(Infection, S = "Susceptible", R = "Resistant")
  )

fig6 <- ggplot(strat_summary %>% filter(times > 5),
       aes(x = times, y = risk, color = strategy)) +
  geom_line(linewidth = 1) +
  facet_grid(Infection ~ Drug, scales = "free_y") +
  scale_color_manual("Strategy", values = strategy_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = "Cumulative Death Risk by Drug, Infection Type, and Strategy",
    subtitle = "Rows = infection type; columns = drug assignment",
    x = "Timestep",
    y = "Cumulative death risk"
  ) +
  theme_paper

ggsave("fig6_stratified_risk.pdf", fig6, width = 7, height = 6, )
ggsave("fig6_stratified_risk.png", fig6, width = 7, height = 6, dpi = 300)
cat("  Saved fig6_stratified_risk.pdf/.png\n")


# ════════════════════════════════════════════════════════════════════════════
# FIGURE 7: Bar chart — End-of-trial Death Risk by Infection & Drug
# ════════════════════════════════════════════════════════════════════════════
cat("Generating Figure 7: End-of-trial death risk bar chart...\n")

final_counts <- res %>%
  group_by(sim, clus) %>%
  summarise(
    tot_SA = sum(TS1_inc), dead_SA = sum(D_TS1),
    tot_SB = sum(TS2_inc), dead_SB = sum(D_TS2),
    tot_RA = sum(TR1_inc), dead_RA = sum(D_TR1),
    tot_RB = sum(TR2_inc), dead_RB = sum(D_TR2),
    .groups = "drop"
  ) %>%
  mutate(strategy = if_else(clus <= floor(nclus / 2), "90/10", "50/50"))

cumrisk_summary <- final_counts %>%
  group_by(strategy) %>%
  summarise(
    risk_SA = mean(dead_SA / tot_SA, na.rm = TRUE),
    se_SA   = sd(dead_SA / tot_SA, na.rm = TRUE) / sqrt(sum(!is.na(tot_SA))),
    risk_SB = mean(dead_SB / tot_SB, na.rm = TRUE),
    se_SB   = sd(dead_SB / tot_SB, na.rm = TRUE) / sqrt(sum(!is.na(tot_SB))),
    risk_RA = mean(dead_RA / tot_RA, na.rm = TRUE),
    se_RA   = sd(dead_RA / tot_RA, na.rm = TRUE) / sqrt(sum(!is.na(tot_RA))),
    risk_RB = mean(dead_RB / tot_RB, na.rm = TRUE),
    se_RB   = sd(dead_RB / tot_RB, na.rm = TRUE) / sqrt(sum(!is.na(tot_RB))),
    .groups = "drop"
  )

bar_df <- cumrisk_summary %>%
  pivot_longer(
    cols = -strategy,
    names_to = c("metric", "inf", "drug"),
    names_pattern = "(risk|se)_(S|R)(A|B)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(
    Infection = recode(inf, S = "Susceptible", R = "Resistant"),
    Drug = recode(drug, A = "Drug A", B = "Drug B")
  )

fig7 <- ggplot(bar_df, aes(x = Infection, y = risk, fill = Drug)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = risk - se, ymax = risk + se),
                position = position_dodge(0.8), width = 0.25) +
  facet_wrap(~ strategy, labeller = labeller(
    strategy = c("50/50" = "50/50 Strategy (\u03B1\u2080)",
                 "90/10" = "90/10 Strategy (\u03B1\u2081)")
  )) +
  scale_fill_manual(values = drug_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "End-of-Trial Cumulative Death Risk by Infection Type and Drug",
    subtitle = "Error bars = \u00B1 1 SE across simulation replicates",
    x = "Infection Type",
    y = "Cumulative death risk",
    fill = "Drug"
  ) +
  theme_paper

ggsave("fig7_bar_risk.pdf", fig7, width = 7, height = 4.5, )
ggsave("fig7_bar_risk.png", fig7, width = 7, height = 4.5, dpi = 300)
cat("  Saved fig7_bar_risk.pdf/.png\n")


cat("\n=== All figures generated ===\n")
cat("Files saved in:", getwd(), "\n")
