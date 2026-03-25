#!/usr/bin/env Rscript
# Generate eFigure: Time-varying IE, TE, TE2, OE

suppressPackageStartupMessages({
  library(ABM)
  library(tidyverse)
})

source_lines <- readLines("run_sim_extract.R")
func_lines <- source_lines[15:113]
eval(parse(text = paste(func_lines, collapse = "\n")))

nclus <- 6; sims <- 10

cat("Running simulation...\n")
set.seed(42)
res <- list()
for (i in 1:nclus) {
  pi <- if (i <= 3) 0.9 else 0.5
  for (j in 1:sims) {
    x <- runsim(pi = pi, timesteps = 350)
    x$clus <- i; x$sim <- j
    x$alloc <- if (i <= 3) "90/10" else "50/50"
    res[[(i - 1) * sims + j]] <- x
  }
  cat("  Cluster", i, "done\n")
}
res <- bind_rows(res)
cat("Simulation done.\n")

# 50-step windows
breaks <- seq(0, 350, 50)
mids <- (breaks[-1] + breaks[-length(breaks)]) / 2
res <- res %>% mutate(
  window_idx = findInterval(times, breaks, rightmost.closed = TRUE),
  window_mid = mids[pmin(window_idx, length(mids))]
) %>% filter(window_idx >= 1, window_idx <= length(mids))

# Per-window aggregation per sim and alloc
window_data <- res %>%
  group_by(sim, alloc, window_mid) %>%
  summarise(
    n_A = sum(TS1_inc + TR1_inc),      # total Drug A treated
    d_A = sum(D_TS1 + D_TR1),          # deaths among Drug A
    n_B = sum(TS2_inc + TR2_inc),      # total Drug B treated
    d_B = sum(D_TS2 + D_TR2),          # deaths among Drug B
    n_total = sum(TS1_inc + TR1_inc + TS2_inc + TR2_inc),
    d_total = sum(D_TS1 + D_TR1 + D_TS2 + D_TR2),
    n_R = sum(TR1_inc + TR2_inc),
    n_all_inf = sum(TS1_inc + TR1_inc + TS2_inc + TR2_inc),
    .groups = "drop"
  ) %>%
  mutate(
    risk_A = d_A / n_A,    # Y(A, alloc)
    risk_B = d_B / n_B,    # Y(B, alloc)
    risk_all = d_total / n_total,  # Y(alloc)
    prop_R = n_R / n_all_inf
  )

# Split by strategy and join on sim x window
d90 <- window_data %>% filter(alloc == "90/10") %>%
  select(sim, window_mid, rA_90 = risk_A, rB_90 = risk_B, rAll_90 = risk_all, pR_90 = prop_R)
d50 <- window_data %>% filter(alloc == "50/50") %>%
  select(sim, window_mid, rA_50 = risk_A, rB_50 = risk_B, rAll_50 = risk_all, pR_50 = prop_R)

wide <- inner_join(d90, d50, by = c("sim", "window_mid"))

# Compute all estimands per sim per window (in percentage points)
effects <- wide %>%
  mutate(
    IE_A = (rA_90 - rA_50) * 100,
    IE_B = (rB_90 - rB_50) * 100,
    TE   = (rA_90 - rB_50) * 100,
    TE2  = (rB_90 - rA_50) * 100,
    OE   = (rAll_50 - rAll_90) * 100,
    prop_R_90 = pR_90 * 100,
    prop_R_50 = pR_50 * 100
  )

# Summary across sims
effect_cols <- c("IE_A", "IE_B", "TE", "TE2", "OE", "prop_R_90", "prop_R_50")

eff_summary <- effects %>%
  group_by(window_mid) %>%
  summarise(
    across(all_of(effect_cols),
           list(mean = ~mean(.x, na.rm = TRUE),
                lo = ~quantile(.x, 0.025, na.rm = TRUE),
                hi = ~quantile(.x, 0.975, na.rm = TRUE))),
    .groups = "drop"
  )

cat("\n=== Summary: Time-Varying Effects ===\n")
print(eff_summary %>%
  select(window_mid,
         IE_A_mean, IE_B_mean, TE_mean, TE2_mean, OE_mean,
         prop_R_90_mean, prop_R_50_mean), n = 20)

# --- Figure: 3-panel ---
# Panel A: IE(A) and IE(B)
# Panel B: OE
# Panel C: Prop R under each strategy

# Prepare long data for IE panel
ie_long <- bind_rows(
  eff_summary %>% transmute(window_mid, effect = "IE(A)",
    mean = IE_A_mean, lo = IE_A_lo, hi = IE_A_hi),
  eff_summary %>% transmute(window_mid, effect = "IE(B)",
    mean = IE_B_mean, lo = IE_B_lo, hi = IE_B_hi)
)

# TE long
te_long <- bind_rows(
  eff_summary %>% transmute(window_mid, effect = "TE",
    mean = TE_mean, lo = TE_lo, hi = TE_hi),
  eff_summary %>% transmute(window_mid, effect = "TE2",
    mean = TE2_mean, lo = TE2_lo, hi = TE2_hi)
)

# OE data
oe_data <- eff_summary %>%
  transmute(window_mid, mean = OE_mean, lo = OE_lo, hi = OE_hi)

# Prop R data
propR_long <- bind_rows(
  eff_summary %>% transmute(window_mid, effect = "90/10",
    mean = prop_R_90_mean, lo = prop_R_90_lo, hi = prop_R_90_hi),
  eff_summary %>% transmute(window_mid, effect = "50/50",
    mean = prop_R_50_mean, lo = prop_R_50_lo, hi = prop_R_50_hi)
)

# Colors
col_ie_a <- "#E15759"   # red
col_ie_b <- "#4E79A7"   # blue
col_te   <- "#E15759"
col_te2  <- "#4E79A7"
col_oe   <- "#59A14F"   # green
col_90   <- "#F28E2B"   # orange
col_50   <- "#76B7B2"   # teal

# Panel A: IE(A) and IE(B)
p_ie <- ggplot(ie_long, aes(x = window_mid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = effect), alpha = 0.12) +
  geom_line(aes(y = mean, color = effect), linewidth = 0.9) +
  geom_point(aes(y = mean, color = effect, shape = effect), size = 2.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c("IE(A)" = col_ie_a, "IE(B)" = col_ie_b),
                     labels = c("IE(A)" = "IE(A)", "IE(B)" = "IE(B)")) +
  scale_fill_manual(values = c("IE(A)" = col_ie_a, "IE(B)" = col_ie_b),
                    labels = c("IE(A)" = "IE(A)", "IE(B)" = "IE(B)")) +
  scale_shape_manual(values = c("IE(A)" = 16, "IE(B)" = 17),
                     labels = c("IE(A)" = "IE(A)", "IE(B)" = "IE(B)")) +
  labs(x = NULL, y = "Effect (pp)", title = "A. Indirect Effects") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.82, 0.85),
        legend.title = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0.8)),
        plot.title = element_text(face = "bold", size = 11))

# Panel B: TE and TE2
p_te <- ggplot(te_long, aes(x = window_mid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = effect), alpha = 0.12) +
  geom_line(aes(y = mean, color = effect), linewidth = 0.9) +
  geom_point(aes(y = mean, color = effect, shape = effect), size = 2.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c("TE" = col_te, "TE2" = col_te2),
                     labels = c("TE" = expression(TE), "TE2" = expression(TE[2]))) +
  scale_fill_manual(values = c("TE" = col_te, "TE2" = col_te2),
                    labels = c("TE" = expression(TE), "TE2" = expression(TE[2]))) +
  scale_shape_manual(values = c("TE" = 16, "TE2" = 17),
                     labels = c("TE" = expression(TE), "TE2" = expression(TE[2]))) +
  labs(x = NULL, y = "Effect (pp)", title = "B. Total Effects") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.82, 0.85),
        legend.title = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0.8)),
        plot.title = element_text(face = "bold", size = 11))

# Panel C: OE
p_oe <- ggplot(oe_data, aes(x = window_mid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, fill = col_oe) +
  geom_line(aes(y = mean), color = col_oe, linewidth = 0.9) +
  geom_point(aes(y = mean), color = col_oe, size = 2.2, shape = 15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  labs(x = NULL, y = "Effect (pp)", title = "C. Overall Effect") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

# Panel D: Proportion R under each strategy
p_pr <- ggplot(propR_long, aes(x = window_mid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = effect), alpha = 0.12) +
  geom_line(aes(y = mean, color = effect), linewidth = 0.9) +
  geom_point(aes(y = mean, color = effect, shape = effect), size = 2.2) +
  scale_color_manual(values = c("90/10" = col_90, "50/50" = col_50)) +
  scale_fill_manual(values = c("90/10" = col_90, "50/50" = col_50)) +
  scale_shape_manual(values = c("90/10" = 16, "50/50" = 17)) +
  labs(x = "Time (simulation steps)", y = "Resistant infections (%)",
       title = "D. Resistant Strain Prevalence") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.82, 0.85),
        legend.title = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0.8)),
        plot.title = element_text(face = "bold", size = 11))

# Combine 4 panels
library(patchwork)
p_all <- (p_ie | p_te) / (p_oe | p_pr) +
  plot_annotation(
    theme = theme(plot.margin = margin(5, 5, 5, 5))
  )

outpath <- "/Users/jug242/Dropbox/Apps/Overleaf/2SR new version/eFigure_all_effects_time.pdf"
ggsave(outpath, p_all, width = 8, height = 6, dpi = 300)
cat("\nSaved figure to:", outpath, "\n")

# --- eTable: Time-Varying IE, TE, OE ---
etable_effects <- effects %>%
  group_by(window_mid) %>%
  summarise(
    IE_A_m = sprintf("%.2f", mean(IE_A)),
    IE_A_si = sprintf("[%.2f, %.2f]", quantile(IE_A, 0.025), quantile(IE_A, 0.975)),
    IE_B_m = sprintf("%.2f", mean(IE_B)),
    IE_B_si = sprintf("[%.2f, %.2f]", quantile(IE_B, 0.025), quantile(IE_B, 0.975)),
    TE_m = sprintf("%.2f", mean(TE)),
    TE_si = sprintf("[%.2f, %.2f]", quantile(TE, 0.025), quantile(TE, 0.975)),
    TE2_m = sprintf("%.2f", mean(TE2)),
    TE2_si = sprintf("[%.2f, %.2f]", quantile(TE2, 0.025), quantile(TE2, 0.975)),
    OE_m = sprintf("%.2f", mean(OE)),
    OE_si = sprintf("[%.2f, %.2f]", quantile(OE, 0.025), quantile(OE, 0.975)),
    pR_90 = sprintf("%.1f", mean(prop_R_90)),
    pR_50 = sprintf("%.1f", mean(prop_R_50)),
    .groups = "drop"
  )

cat("\n=== eTable: Time-Varying Causal Effects ===\n\n")
print(as.data.frame(etable_effects), right = FALSE)
