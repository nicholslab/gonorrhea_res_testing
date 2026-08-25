# ============================================================
# Figure 4: Weighted-average cost-neutrality threshold prices
# Weights: MSM=0.53, MSW=0.21, Women=0.26
# Layout: Strategy (rows) × NDA scenario (columns)
# Color: diverging blue=negative / white=0 / red=positive
# Data source: Figure_4.xlsx (Col E = 0% NDA, Col F = 5% NDA)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ---- Hardcoded weighted-average data ----
# All values are population-weighted (MSM=0.53, MSW=0.21, Women=0.26)

df <- tribble(
  ~strategy,       ~nda_label,          ~cip_res, ~cro_res, ~threshold_price,
  
  # ── CIP, 0% NDA ──────────────────────────────────────────
  "CIP","0% NDA resistance", 0.1, 0.000,  0.117,
  "CIP","0% NDA resistance", 0.1, 0.025,  6.924,
  "CIP","0% NDA resistance", 0.1, 0.050, 13.307,
  "CIP","0% NDA resistance", 0.1, 0.075, 20.819,
  "CIP","0% NDA resistance", 0.1, 0.100, 26.700,
  "CIP","0% NDA resistance", 0.1, 0.125, 31.711,
  "CIP","0% NDA resistance", 0.1, 0.150, 39.747,
  "CIP","0% NDA resistance", 0.3, 0.000,  0.403,
  "CIP","0% NDA resistance", 0.3, 0.025,  6.788,
  "CIP","0% NDA resistance", 0.3, 0.050, 11.008,
  "CIP","0% NDA resistance", 0.3, 0.075, 14.242,
  "CIP","0% NDA resistance", 0.3, 0.100, 18.832,
  "CIP","0% NDA resistance", 0.3, 0.125, 25.021,
  "CIP","0% NDA resistance", 0.3, 0.150, 29.504,
  "CIP","0% NDA resistance", 0.5, 0.000, -0.158,
  "CIP","0% NDA resistance", 0.5, 0.025,  2.940,
  "CIP","0% NDA resistance", 0.5, 0.050,  4.347,
  "CIP","0% NDA resistance", 0.5, 0.075, 10.969,
  "CIP","0% NDA resistance", 0.5, 0.100, 13.635,
  "CIP","0% NDA resistance", 0.5, 0.125, 18.114,
  "CIP","0% NDA resistance", 0.5, 0.150, 21.128,
  "CIP","0% NDA resistance", 0.7, 0.000, -0.199,
  "CIP","0% NDA resistance", 0.7, 0.025,  0.239,
  "CIP","0% NDA resistance", 0.7, 0.050,  2.997,
  "CIP","0% NDA resistance", 0.7, 0.075,  5.814,
  "CIP","0% NDA resistance", 0.7, 0.100,  6.841,
  "CIP","0% NDA resistance", 0.7, 0.125,  9.929,
  "CIP","0% NDA resistance", 0.7, 0.150, 11.704,
  
  # ── CIP, 5% NDA ──────────────────────────────────────────
  "CIP","5% NDA resistance", 0.1, 0.000,  1.210,
  "CIP","5% NDA resistance", 0.1, 0.025,  7.292,
  "CIP","5% NDA resistance", 0.1, 0.050, 15.415,
  "CIP","5% NDA resistance", 0.1, 0.075, 19.241,
  "CIP","5% NDA resistance", 0.1, 0.100, 26.560,
  "CIP","5% NDA resistance", 0.1, 0.125, 32.512,
  "CIP","5% NDA resistance", 0.1, 0.150, 38.832,
  "CIP","5% NDA resistance", 0.3, 0.000, -1.339,
  "CIP","5% NDA resistance", 0.3, 0.025,  3.687,
  "CIP","5% NDA resistance", 0.3, 0.050,  8.437,
  "CIP","5% NDA resistance", 0.3, 0.075, 15.110,
  "CIP","5% NDA resistance", 0.3, 0.100, 17.725,
  "CIP","5% NDA resistance", 0.3, 0.125, 26.452,
  "CIP","5% NDA resistance", 0.3, 0.150, 29.973,
  "CIP","5% NDA resistance", 0.5, 0.000,  0.139,
  "CIP","5% NDA resistance", 0.5, 0.025,  1.266,
  "CIP","5% NDA resistance", 0.5, 0.050,  8.487,
  "CIP","5% NDA resistance", 0.5, 0.075, 10.366,
  "CIP","5% NDA resistance", 0.5, 0.100, 13.306,
  "CIP","5% NDA resistance", 0.5, 0.125, 18.863,
  "CIP","5% NDA resistance", 0.5, 0.150, 20.557,
  "CIP","5% NDA resistance", 0.7, 0.000, -2.359,
  "CIP","5% NDA resistance", 0.7, 0.025,  0.521,
  "CIP","5% NDA resistance", 0.7, 0.050,  2.921,
  "CIP","5% NDA resistance", 0.7, 0.075,  5.468,
  "CIP","5% NDA resistance", 0.7, 0.100,  8.158,
  "CIP","5% NDA resistance", 0.7, 0.125,  9.493,
  "CIP","5% NDA resistance", 0.7, 0.150, 10.167,
  
  # ── CIP+CRO, 0% NDA ──────────────────────────────────────
  "CIP+CRO","0% NDA resistance", 0.1, 0.000,  0.877,
  "CIP+CRO","0% NDA resistance", 0.1, 0.025,  7.612,
  "CIP+CRO","0% NDA resistance", 0.1, 0.050, 14.494,
  "CIP+CRO","0% NDA resistance", 0.1, 0.075, 22.697,
  "CIP+CRO","0% NDA resistance", 0.1, 0.100, 30.313,
  "CIP+CRO","0% NDA resistance", 0.1, 0.125, 37.664,
  "CIP+CRO","0% NDA resistance", 0.1, 0.150, 43.754,
  "CIP+CRO","0% NDA resistance", 0.3, 0.000, -1.407,
  "CIP+CRO","0% NDA resistance", 0.3, 0.025,  6.207,
  "CIP+CRO","0% NDA resistance", 0.3, 0.050, 13.297,
  "CIP+CRO","0% NDA resistance", 0.3, 0.075, 20.424,
  "CIP+CRO","0% NDA resistance", 0.3, 0.100, 26.239,
  "CIP+CRO","0% NDA resistance", 0.3, 0.125, 33.691,
  "CIP+CRO","0% NDA resistance", 0.3, 0.150, 41.499,
  "CIP+CRO","0% NDA resistance", 0.5, 0.000, -1.245,
  "CIP+CRO","0% NDA resistance", 0.5, 0.025,  5.790,
  "CIP+CRO","0% NDA resistance", 0.5, 0.050, 12.780,
  "CIP+CRO","0% NDA resistance", 0.5, 0.075, 18.929,
  "CIP+CRO","0% NDA resistance", 0.5, 0.100, 25.618,
  "CIP+CRO","0% NDA resistance", 0.5, 0.125, 33.338,
  "CIP+CRO","0% NDA resistance", 0.5, 0.150, 38.903,
  "CIP+CRO","0% NDA resistance", 0.7, 0.000, -4.901,
  "CIP+CRO","0% NDA resistance", 0.7, 0.025,  4.423,
  "CIP+CRO","0% NDA resistance", 0.7, 0.050,  9.434,
  "CIP+CRO","0% NDA resistance", 0.7, 0.075, 19.146,
  "CIP+CRO","0% NDA resistance", 0.7, 0.100, 22.466,
  "CIP+CRO","0% NDA resistance", 0.7, 0.125, 30.931,
  "CIP+CRO","0% NDA resistance", 0.7, 0.150, 34.573,
  
  # ── CIP+CRO, 5% NDA ──────────────────────────────────────
  "CIP+CRO","5% NDA resistance", 0.1, 0.000, -1.712,
  "CIP+CRO","5% NDA resistance", 0.1, 0.025,  6.914,
  "CIP+CRO","5% NDA resistance", 0.1, 0.050, 15.991,
  "CIP+CRO","5% NDA resistance", 0.1, 0.075, 21.999,
  "CIP+CRO","5% NDA resistance", 0.1, 0.100, 31.240,
  "CIP+CRO","5% NDA resistance", 0.1, 0.125, 36.985,
  "CIP+CRO","5% NDA resistance", 0.1, 0.150, 43.755,
  "CIP+CRO","5% NDA resistance", 0.3, 0.000, -2.331,
  "CIP+CRO","5% NDA resistance", 0.3, 0.025,  4.543,
  "CIP+CRO","5% NDA resistance", 0.3, 0.050, 13.062,
  "CIP+CRO","5% NDA resistance", 0.3, 0.075, 19.326,
  "CIP+CRO","5% NDA resistance", 0.3, 0.100, 25.830,
  "CIP+CRO","5% NDA resistance", 0.3, 0.125, 33.581,
  "CIP+CRO","5% NDA resistance", 0.3, 0.150, 40.509,
  "CIP+CRO","5% NDA resistance", 0.5, 0.000, -1.305,
  "CIP+CRO","5% NDA resistance", 0.5, 0.025,  5.082,
  "CIP+CRO","5% NDA resistance", 0.5, 0.050, 11.006,
  "CIP+CRO","5% NDA resistance", 0.5, 0.075, 18.720,
  "CIP+CRO","5% NDA resistance", 0.5, 0.100, 26.081,
  "CIP+CRO","5% NDA resistance", 0.5, 0.125, 30.304,
  "CIP+CRO","5% NDA resistance", 0.5, 0.150, 35.735,
  "CIP+CRO","5% NDA resistance", 0.7, 0.000, -2.345,
  "CIP+CRO","5% NDA resistance", 0.7, 0.025,  2.315,
  "CIP+CRO","5% NDA resistance", 0.7, 0.050, 10.771,
  "CIP+CRO","5% NDA resistance", 0.7, 0.075, 14.910,
  "CIP+CRO","5% NDA resistance", 0.7, 0.100, 23.588,
  "CIP+CRO","5% NDA resistance", 0.7, 0.125, 27.738,
  "CIP+CRO","5% NDA resistance", 0.7, 0.150, 33.561,
  
  # ── CIP+CRO+NDA, 0% NDA ──────────────────────────────────
  "CIP+CRO+NDA","0% NDA resistance", 0.1, 0.000,  2.264,
  "CIP+CRO+NDA","0% NDA resistance", 0.1, 0.025,  8.143,
  "CIP+CRO+NDA","0% NDA resistance", 0.1, 0.050, 14.593,
  "CIP+CRO+NDA","0% NDA resistance", 0.1, 0.075, 21.381,
  "CIP+CRO+NDA","0% NDA resistance", 0.1, 0.100, 30.048,
  "CIP+CRO+NDA","0% NDA resistance", 0.1, 0.125, 36.967,
  "CIP+CRO+NDA","0% NDA resistance", 0.1, 0.150, 44.895,
  "CIP+CRO+NDA","0% NDA resistance", 0.3, 0.000, -0.402,
  "CIP+CRO+NDA","0% NDA resistance", 0.3, 0.025,  6.791,
  "CIP+CRO+NDA","0% NDA resistance", 0.3, 0.050, 13.267,
  "CIP+CRO+NDA","0% NDA resistance", 0.3, 0.075, 20.459,
  "CIP+CRO+NDA","0% NDA resistance", 0.3, 0.100, 26.766,
  "CIP+CRO+NDA","0% NDA resistance", 0.3, 0.125, 33.701,
  "CIP+CRO+NDA","0% NDA resistance", 0.3, 0.150, 41.396,
  "CIP+CRO+NDA","0% NDA resistance", 0.5, 0.000, -2.039,
  "CIP+CRO+NDA","0% NDA resistance", 0.5, 0.025,  5.425,
  "CIP+CRO+NDA","0% NDA resistance", 0.5, 0.050, 11.325,
  "CIP+CRO+NDA","0% NDA resistance", 0.5, 0.075, 18.168,
  "CIP+CRO+NDA","0% NDA resistance", 0.5, 0.100, 25.865,
  "CIP+CRO+NDA","0% NDA resistance", 0.5, 0.125, 32.558,
  "CIP+CRO+NDA","0% NDA resistance", 0.5, 0.150, 37.485,
  "CIP+CRO+NDA","0% NDA resistance", 0.7, 0.000, -2.987,
  "CIP+CRO+NDA","0% NDA resistance", 0.7, 0.025,  4.368,
  "CIP+CRO+NDA","0% NDA resistance", 0.7, 0.050, 10.688,
  "CIP+CRO+NDA","0% NDA resistance", 0.7, 0.075, 15.881,
  "CIP+CRO+NDA","0% NDA resistance", 0.7, 0.100, 22.174,
  "CIP+CRO+NDA","0% NDA resistance", 0.7, 0.125, 30.337,
  "CIP+CRO+NDA","0% NDA resistance", 0.7, 0.150, 34.742,
  
  # ── CIP+CRO+NDA, 5% NDA ──────────────────────────────────
  "CIP+CRO+NDA","5% NDA resistance", 0.1, 0.000, -0.720,
  "CIP+CRO+NDA","5% NDA resistance", 0.1, 0.025,  7.902,
  "CIP+CRO+NDA","5% NDA resistance", 0.1, 0.050, 13.641,
  "CIP+CRO+NDA","5% NDA resistance", 0.1, 0.075, 21.724,
  "CIP+CRO+NDA","5% NDA resistance", 0.1, 0.100, 28.774,
  "CIP+CRO+NDA","5% NDA resistance", 0.1, 0.125, 36.724,
  "CIP+CRO+NDA","5% NDA resistance", 0.1, 0.150, 43.650,
  "CIP+CRO+NDA","5% NDA resistance", 0.3, 0.000, -0.100,
  "CIP+CRO+NDA","5% NDA resistance", 0.3, 0.025,  7.204,
  "CIP+CRO+NDA","5% NDA resistance", 0.3, 0.050, 13.959,
  "CIP+CRO+NDA","5% NDA resistance", 0.3, 0.075, 22.035,
  "CIP+CRO+NDA","5% NDA resistance", 0.3, 0.100, 27.667,
  "CIP+CRO+NDA","5% NDA resistance", 0.3, 0.125, 33.233,
  "CIP+CRO+NDA","5% NDA resistance", 0.3, 0.150, 40.763,
  "CIP+CRO+NDA","5% NDA resistance", 0.5, 0.000, -2.589,
  "CIP+CRO+NDA","5% NDA resistance", 0.5, 0.025,  8.019,
  "CIP+CRO+NDA","5% NDA resistance", 0.5, 0.050, 12.207,
  "CIP+CRO+NDA","5% NDA resistance", 0.5, 0.075, 19.244,
  "CIP+CRO+NDA","5% NDA resistance", 0.5, 0.100, 25.589,
  "CIP+CRO+NDA","5% NDA resistance", 0.5, 0.125, 31.903,
  "CIP+CRO+NDA","5% NDA resistance", 0.5, 0.150, 38.773,
  "CIP+CRO+NDA","5% NDA resistance", 0.7, 0.000, -3.431,
  "CIP+CRO+NDA","5% NDA resistance", 0.7, 0.025,  2.897,
  "CIP+CRO+NDA","5% NDA resistance", 0.7, 0.050, 10.869,
  "CIP+CRO+NDA","5% NDA resistance", 0.7, 0.075, 15.270,
  "CIP+CRO+NDA","5% NDA resistance", 0.7, 0.100, 24.004,
  "CIP+CRO+NDA","5% NDA resistance", 0.7, 0.125, 29.951,
  "CIP+CRO+NDA","5% NDA resistance", 0.7, 0.150, 36.686
) %>%
  mutate(
    Strategy  = factor(strategy,  levels = c("CIP", "CIP+CRO", "CIP+CRO+NDA")),
    NDA_label = factor(nda_label, levels = c("0% NDA resistance", "5% NDA resistance")),
    label     = sprintf("%.0f", threshold_price)
  )

# ---- Axis breaks ----
x_breaks <- c(0.1, 0.3, 0.5, 0.7)
y_breaks <- c(0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15)

# ---- Shared diverging color scale ----
# Fixed to match original figure (lim_abs = 45)
lim_abs <- 45

# ---- Plot ----
p <- ggplot(df, aes(x = cip_res, y = cro_res, fill = threshold_price)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = label,
        color  = abs(threshold_price) < lim_abs * 0.3),
    size = 3, fontface = "bold"
  ) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(
    low      = "#0072B2",   # Okabe-Ito blue  (negative = testing adds cost)
    mid      = "white",
    high     = "#E69F00",   # Okabe-Ito orange (positive = testing saves money)
    midpoint = 0,
    limits   = c(-lim_abs, lim_abs),
    oob      = scales::squish,
    name     = "Threshold\nPrice (USD)",
    labels   = scales::dollar_format(accuracy = 1)
  ) +
  scale_x_continuous(
    "CIP Resistance Prevalence",
    breaks = x_breaks, labels = scales::percent_format(accuracy = 1)
  ) +
  scale_y_continuous(
    "CRO Resistance Prevalence",
    breaks = y_breaks, labels = scales::percent_format(accuracy = 0.1),
    expand = c(0, 0)
  ) +
  facet_grid(Strategy ~ NDA_label) +
  labs(
    title    = "Figure 4. Cost-Neutrality Threshold Price for Reflex NG Resistance Testing",
    subtitle = "Population-weighted average (MSM 53%, MSW 21%, Women 26%)  |  Orange = testing saves money; Blue = testing adds cost",
    caption  = paste0(
      "Threshold price: maximum per-test cost at which resistance-guided therapy is cost-neutral vs. empiric ceftriaxone.\n",
      "Negative values indicate scenarios where resistance testing increases expected costs (testing not cost-neutral at any price).\n",
      "Results based on 1,000 PSA draws \u00d7 10,000 individuals per scenario."
    )
  ) 
  theme_minimal(base_size = 12) +
  theme(
    panel.grid     = element_blank(),
    axis.text      = element_text(color = "black", size = 9),
    axis.text.x    = element_text(angle = 45, hjust = 1),
    plot.title     = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle  = element_text(hjust = 0.5, size = 10, color = "grey30"),
    plot.caption   = element_text(hjust = 0, size = 8, color = "grey50", lineheight = 1.3),
    strip.text     = element_text(face = "bold", size = 11),
    legend.position     = "right",
    legend.title        = element_text(size = 10),
    legend.key.height   = unit(1.5, "cm")
  )

ggsave("Figure4_weighted_threshold_heatmap.png", p, width = 11, height = 9, dpi = 300)
cat("Saved: Figure4_weighted_threshold_heatmap.png\n")
cat(sprintf("Color scale: -$%.0f to +$%.0f (diverging, centered at 0)\n", lim_abs, lim_abs))
cat("Weights applied: MSM=0.53, MSW=0.21, Women=0.26\n")