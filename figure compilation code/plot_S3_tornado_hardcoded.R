# ============================================================
# S3 Tornado diagram — fully hardcoded from updated DSA results
# Requires: tidyverse, ggplot2, patchwork, cowplot
# ============================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

BASE_CIP    <- 6.801
BASE_TRIPLE <- 11.912

dsa_df <- tribble(
  ~panel,        ~param,            ~label,                        ~label_unit, ~base_val, ~low_val, ~high_val, ~low_thresh, ~high_thresh,
  "CIP",         "p_partner_treat", "Partner notification rate",   "%",          0.30,      0.20,     0.50,       5.767,      7.573,
  "CIP+CRO+NDA", "p_partner_treat", "Partner notification rate",   "%",          0.30,      0.20,     0.50,      10.566,     13.697,
  "CIP",         "cost_nda_ndb",    "Cost of NDA/NDB",             "$",          50,        25,       500,        5.965,      15.399,
  "CIP+CRO+NDA", "cost_nda_ndb",    "Cost of NDA/NDB",             "$",          50,        25,       500,       12.002,      9.674,
  "CIP",         "p_return_sympt",  "P(return | symptoms)",        "%",          0.90,      0.70,     0.95,       5.712,      6.574,
  "CIP+CRO+NDA", "p_return_sympt",  "P(return | symptoms)",        "%",          0.90,      0.70,     0.95,      10.333,     12.299,
  "CIP",         "sens_res_test",   "Resistance test sensitivity", "%",          0.99,      0.92,     1.00,      -2.031,      7.787,
  "CIP+CRO+NDA", "sens_res_test",   "Resistance test sensitivity", "%",          0.99,      0.92,     1.00,       2.242,     12.900,
  "CIP",         "spec_res_test",   "Resistance test specificity", "%",          0.98,      0.96,     1.00,       6.543,      6.661,
  "CIP+CRO+NDA", "spec_res_test",   "Resistance test specificity", "%",          0.98,      0.96,     1.00,      11.160,     12.400
)

# ---- y-axis labels: "Parameter name (low; high)" ----
dsa_df <- dsa_df %>%
  mutate(
    y_label = case_when(
      label_unit == "$" ~
        sprintf("%s ($%.0f; $%.0f)", label, low_val, high_val),
      label_unit == "%" ~
        sprintf("%s (%.0f%%; %.0f%%)", label, low_val * 100, high_val * 100),
      TRUE ~
        sprintf("%s (%.2f; %.2f)", label, low_val, high_val)
    )
  )

# ---- Tornado geometry ----
tornado_df <- dsa_df %>%
  mutate(
    bar_left    = pmin(low_thresh, high_thresh),
    bar_right   = pmax(low_thresh, high_thresh),
    total_swing = bar_right - bar_left,
    low_is_left = low_thresh <= high_thresh
  )

# ---- Colours ----
COL_HIGH <- "#4682B4"   # steelblue  — higher parameter value
COL_LOW  <- "#C06B5E"   # muted coral — lower parameter value
COL_BAR  <- "#7AAFC8"   # bar fill

# ---- Plot function ----
make_tornado <- function(data, panel_name, base_val) {
  
  param_order <- data %>%
    arrange(total_swing) %>%
    pull(y_label)
  
  data <- data %>%
    mutate(y_label = factor(y_label, levels = param_order))
  
  # Build annotation df — colour passed OUTSIDE aes() to avoid check_subclass error
  ann <- bind_rows(
    data %>% transmute(
      y_label, x = bar_left,
      txt   = sprintf("$%.1f", bar_left),
      hjust = 1.12,
      col   = ifelse(low_is_left, COL_LOW, COL_HIGH)
    ),
    data %>% transmute(
      y_label, x = bar_right,
      txt   = sprintf("$%.1f", bar_right),
      hjust = -0.12,
      col   = ifelse(low_is_left, COL_HIGH, COL_LOW)
    )
  )
  
  ggplot(data, aes(y = y_label)) +
    geom_segment(
      aes(x = bar_left, xend = bar_right, yend = y_label),
      linewidth = 9, color = COL_BAR, alpha = 0.78
    ) +
    geom_vline(
      xintercept = base_val,
      color = "black", linetype = "dashed", linewidth = 0.8
    ) +
    geom_text(
      data    = ann,
      mapping = aes(x = x, y = y_label, label = txt, hjust = hjust),
      color   = ann$col,   # <-- outside aes(), avoids check_subclass error
      size    = 3.3,
      fontface = "bold"
    ) +
    scale_x_continuous(
      labels = scales::dollar,
      expand = expansion(mult = 0.18)
    ) +
    labs(
      title    = sprintf("One-Way Sensitivity Analysis: %s Panel", panel_name),
      subtitle = sprintf(
        "Reference scenario: MSM, CIP=50%%, CRO=5%%, NDA=0%%  |  Base case threshold: $%.1f",
        base_val
      ),
      x = "Threshold Price (USD)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y        = element_text(size = 10.5),
      plot.subtitle      = element_text(size = 9, color = "grey40"),
      plot.title         = element_text(size = 12)
    )
}

# ---- Make panels ----
p_cip <- make_tornado(
  tornado_df %>% filter(panel == "CIP"),
  "CIP-only", BASE_CIP
)

p_triple <- make_tornado(
  tornado_df %>% filter(panel == "CIP+CRO+NDA"),
  "CIP+CRO+NDA", BASE_TRIPLE
)

# ---- Shared colour legend ----
legend_df <- tibble(
  x   = c(1, 2),
  y   = c(1, 1),
  lbl = c("Higher parameter value", "Lower parameter value")
)

p_leg <- ggplot(legend_df, aes(x = x, y = y, color = lbl)) +
  geom_point(size = 5, shape = 15) +
  scale_color_manual(
    values = c("Higher parameter value" = COL_HIGH,
               "Lower parameter value"  = COL_LOW),
    name = NULL
  ) +
  theme_void() +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text      = element_text(size = 10)
  )

shared_legend <- cowplot::get_legend(p_leg)

# ---- Combine ----
p_stacked <- p_cip / p_triple +
  plot_annotation(
    title = "Figure S3. Deterministic Sensitivity Analysis \u2014 Tornado Diagrams",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

p_final <- cowplot::plot_grid(
  p_stacked, shared_legend,
  ncol = 1, rel_heights = c(1, 0.04)
)

ggsave("gonorrhea_S3_tornado_diagram.png", p_final,
       width = 11, height = 10, dpi = 300)
cat("Saved: gonorrhea_S3_tornado_diagram.png\n")