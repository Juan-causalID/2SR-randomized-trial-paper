#!/usr/bin/env Rscript
# Generate eFigure: Time-varying Direct Effect

suppressPackageStartupMessages({
  library(ABM)
  library(tidyverse)
})

source_lines <- readLines("run_sim_extract.R")
func_lines <- source_lines[15:113]
eval(parse(text = paste(func_lines, collapse = "\n")))

nclus <- 6; sims <- 100

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

# Per-window aggregation
window_data <- res %>%
  group_by(sim, alloc, window_mid) %>%
  summarise(
    n_A_S = sum(TS1_inc), d_A_S = sum(D_TS1),
    n_B_S = sum(TS2_inc), d_B_S = sum(D_TS2),
    n_A_R = sum(TR1_inc), d_A_R = sum(D_TR1),
    n_B_R = sum(TR2_inc), d_B_R = sum(D_TR2),
    .groups = "drop"
  )

MIN_COUNT <- 30

de_by_window <- window_data %>%
  group_by(sim, alloc, window_mid) %>%
  summarise(
    tot_A = sum(n_A_S + n_A_R), tot_d_A = sum(d_A_S + d_A_R),
    tot_B = sum(n_B_S + n_B_R), tot_d_B = sum(d_B_S + d_B_R),
    tot_A_S = sum(n_A_S), tot_d_A_S = sum(d_A_S),
    tot_B_S = sum(n_B_S), tot_d_B_S = sum(d_B_S),
    tot_A_R = sum(n_A_R), tot_d_A_R = sum(d_A_R),
    tot_B_R = sum(n_B_R), tot_d_B_R = sum(d_B_R),
    prop_R = sum(n_A_R + n_B_R) / sum(n_A_S + n_B_S + n_A_R + n_B_R),
    .groups = "drop"
  ) %>%
  mutate(
    risk_A = tot_d_A / tot_A,
    risk_B = tot_d_B / tot_B,
    DE = (risk_A - risk_B) * 100,
    DE_S = if_else(tot_A_S >= MIN_COUNT & tot_B_S >= MIN_COUNT,
                   (tot_d_A_S/tot_A_S - tot_d_B_S/tot_B_S) * 100, NA_real_),
    DE_R = if_else(tot_A_R >= MIN_COUNT & tot_B_R >= MIN_COUNT,
                   (tot_d_A_R/tot_A_R - tot_d_B_R/tot_B_R) * 100, NA_real_),
    prop_R_pct = prop_R * 100
  )

# Summary across sims
de_summary <- de_by_window %>%
  group_by(alloc, window_mid) %>%
  summarise(
    across(c(DE, DE_S, DE_R, prop_R_pct),
           list(mean = ~mean(.x, na.rm = TRUE),
                lo = ~quantile(.x, 0.025, na.rm = TRUE),
                hi = ~quantile(.x, 0.975, na.rm = TRUE),
                n = ~sum(!is.na(.x)))),
    .groups = "drop"
  ) %>%
  mutate(
    DE_S_mean = if_else(DE_S_n >= 5, DE_S_mean, NA_real_),
    DE_S_lo = if_else(DE_S_n >= 5, DE_S_lo, NA_real_),
    DE_S_hi = if_else(DE_S_n >= 5, DE_S_hi, NA_real_),
    DE_R_mean = if_else(DE_R_n >= 5, DE_R_mean, NA_real_),
    DE_R_lo = if_else(DE_R_n >= 5, DE_R_lo, NA_real_),
    DE_R_hi = if_else(DE_R_n >= 5, DE_R_hi, NA_real_)
  )

cat("\n=== Summary ===\n")
print(de_summary %>% select(alloc, window_mid, DE_mean, DE_S_mean, DE_R_mean,
                             DE_R_n, prop_R_pct_mean), n = 20)

# --- Single figure: Marginal DE + Proportion R on dual y-axes ---
# Push R% line well above DE so they only converge at zero
# DE axis: -2 to 26pp; R% axis: 0 to ~75%

scale_factor <- 75 / 28  # smaller factor → R% line sits higher

plot_data <- de_summary %>%
  select(alloc, window_mid,
         DE_mean, DE_lo, DE_hi,
         prop_R_pct_mean, prop_R_pct_lo, prop_R_pct_hi)

p <- ggplot(plot_data, aes(x = window_mid)) +
  # Proportion R (rescaled to DE axis) — offset slightly upward
  geom_ribbon(aes(ymin = prop_R_pct_lo / scale_factor,
                  ymax = prop_R_pct_hi / scale_factor),
              alpha = 0.10, fill = "#B07AA1") +
  geom_line(aes(y = prop_R_pct_mean / scale_factor),
            color = "#B07AA1", linewidth = 0.9, linetype = "solid") +
  geom_point(aes(y = prop_R_pct_mean / scale_factor),
             color = "#B07AA1", size = 2.2, shape = 17) +
  # Marginal DE
  geom_ribbon(aes(ymin = DE_lo, ymax = DE_hi),
              alpha = 0.12, fill = "#1F77B4") +
  geom_line(aes(y = DE_mean), color = "#1F77B4", linewidth = 0.9) +
  geom_point(aes(y = DE_mean), color = "#1F77B4", size = 2.2) +
  # Reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  # Facet by strategy
  facet_wrap(~alloc, ncol = 2) +
  # Dual y-axes
  scale_y_continuous(
    name = "Marginal direct effect (percentage points)",
    limits = c(-2, 26),
    breaks = seq(0, 25, 5),
    sec.axis = sec_axis(~ . * scale_factor,
                        name = "Resistant infections (%)",
                        breaks = seq(0, 60, 20))
  ) +
  labs(x = "Time (simulation steps)") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 11),
    axis.title.y.left = element_text(color = "#1F77B4", face = "bold"),
    axis.text.y.left = element_text(color = "#1F77B4"),
    axis.title.y.right = element_text(color = "#B07AA1", face = "bold"),
    axis.text.y.right = element_text(color = "#B07AA1")
  )

outpath <- "/Users/jug242/Dropbox/Apps/Overleaf/2SR new version/eFigure_DE_time.pdf"
ggsave(outpath, p, width = 7.5, height = 3.8, dpi = 300)
cat("\nSaved to:", outpath, "\n")

# --- eTable: Time-varying DE summary ---
etable <- de_by_window %>%
  group_by(alloc, window_mid) %>%
  summarise(
    DE_mean = sprintf("%.2f", mean(DE, na.rm = TRUE)),
    DE_95SI = sprintf("[%.2f, %.2f]",
                      quantile(DE, 0.025, na.rm = TRUE),
                      quantile(DE, 0.975, na.rm = TRUE)),
    DE_S_mean = sprintf("%.2f", mean(DE_S, na.rm = TRUE)),
    DE_S_95SI = sprintf("[%.2f, %.2f]",
                        quantile(DE_S, 0.025, na.rm = TRUE),
                        quantile(DE_S, 0.975, na.rm = TRUE)),
    DE_R_mean = if_else(sum(!is.na(DE_R)) >= 5,
                        sprintf("%.2f", mean(DE_R, na.rm = TRUE)),
                        "---"),
    DE_R_95SI = if_else(sum(!is.na(DE_R)) >= 5,
                        sprintf("[%.2f, %.2f]",
                                quantile(DE_R, 0.025, na.rm = TRUE),
                                quantile(DE_R, 0.975, na.rm = TRUE)),
                        "---"),
    prop_R = sprintf("%.1f", mean(prop_R_pct, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(Strategy = alloc, `Time window` = window_mid,
         `Marginal DE (pp)` = DE_mean, `95% SI` = DE_95SI,
         `DE_S (pp)` = DE_S_mean, `95% SI (S)` = DE_S_95SI,
         `DE_R (pp)` = DE_R_mean, `95% SI (R)` = DE_R_95SI,
         `% R infections` = prop_R)

cat("\n=== eTable: Time-Varying Direct Effects ===\n\n")
print(as.data.frame(etable), right = FALSE)
