# ============================================================
# Figure 2 & Figure 3 — Updated with NDA=0% output
# Reads: gonorrhea_summary_table_for_manuscript_fast_fp_treated_with_categories.csv
#
# Figure 2: Per-person cost by group & ceftriaxone resistance, NO testing
#           (replicates original Figure 2 — bar chart by group)
# Figure 3: Per-person cost by group & ceftriaxone resistance, comparing
#           No testing vs Ciprofloxacin vs Ciprofloxacin+Ceftriaxone testing strategies
#           (replicates original Figure 3 — A/B/C panel)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

df <- read_csv("gonorrhea_summary_table_for_manuscript_fast_fp_treated_with_categories.csv",
               show_col_types = FALSE)

# ── Shared constants ──────────────────────────────────────────────────────────
# Fix ciprofloxacin resistance at 50% (reference scenario used in original figures)
CIP_REF <- 0.5

# Ceftriaxone resistance levels to show (0%, 5%, 15%)
CRO_SHOW <- c(0, 0.05, 0.15)

# Cost category colours — match original figure exactly
cat_colors <- c(
  "Initial visit"              = "#B8A9C9",   # dusty lavender
  "Initial diagnostic"         = "#6BAED6",   # steel blue
  "Ceftriaxone"                = "#D9795A",   # terracotta
  "Ciprofloxacin"              = "#E0A96D",   # warm sand
  "Novel Drug A"               = "#C994B0",   # mauve
  "Follow-up visit"            = "#74A9C8",   # muted sky blue
  "Follow-up diagnostic"       = "#74B49B",   # sage green
  "Follow-up meds"             = "#9DC88D",   # soft green
  "Partner treatment"          = "#C97A9A",   # dusty rose
  "Complications"              = "#A8B86A"    # muted olive
)

cat_order <- names(cat_colors)

# ── Helper: pivot WITHOUT-test cost categories for one row ────────────────────
pivot_without <- function(df_in) {
  df_in %>%
    select(Group, Strategy, cip_res, cro_res,
           "Initial visit"        = cost_without_test_ng_init_visit_mean,
           "Initial diagnostic"   = cost_without_test_ng_init_diag_mean,
           "Ceftriaxone"          = cost_without_test_ng_cro_mean,
           "Ciprofloxacin"        = cost_without_test_ng_cip_mean,
           "Novel Drug A"         = cost_without_test_ng_nda_mean,
           "Follow-up visit"      = cost_without_test_ng_follow_visit_mean,
           "Follow-up diagnostic" = cost_without_test_ng_follow_diag_mean,
           "Follow-up meds"       = cost_without_test_ng_follow_drugs_mean,
           "Partner treatment"    = cost_without_test_ng_partner_total_mean,
           "Complications"        = cost_without_test_ng_complications_mean
    ) %>%
    pivot_longer(cols = all_of(cat_order),
                 names_to = "category", values_to = "cost") %>%
    mutate(category = factor(category, levels = cat_order),
           testing = "No resistance testing")
}

# ── Helper: pivot WITH-test cost categories ───────────────────────────────────
pivot_with <- function(df_in) {
  df_in %>%
    select(Group, Strategy, cip_res, cro_res,
           "Initial visit"        = cost_with_test_ng_init_visit_mean,
           "Initial diagnostic"   = cost_with_test_ng_init_diag_mean,
           "Ceftriaxone"          = cost_with_test_ng_cro_mean,
           "Ciprofloxacin"        = cost_with_test_ng_cip_mean,
           "Novel Drug A"         = cost_with_test_ng_nda_mean,
           "Follow-up visit"      = cost_with_test_ng_follow_visit_mean,
           "Follow-up diagnostic" = cost_with_test_ng_follow_diag_mean,
           "Follow-up meds"       = cost_with_test_ng_follow_drugs_mean,
           "Partner treatment"    = cost_with_test_ng_partner_total_mean,
           "Complications"        = cost_with_test_ng_complications_mean
    ) %>%
    pivot_longer(cols = all_of(cat_order),
                 names_to = "category", values_to = "cost") %>%
    mutate(category = factor(category, levels = cat_order),
           testing = Strategy)
}

# ── Filter to reference ciprofloxacin resistance & selected ceftriaxone levels ─
df_ref <- df %>%
  filter(cip_res == CIP_REF, cro_res %in% CRO_SHOW)

# ── x-axis labels ─────────────────────────────────────────────────────────────
cro_labels <- c("0"    = "0% Ceftriaxone resistance",
                "0.05" = "5% Ceftriaxone resistance",
                "0.15" = "15% Ceftriaxone resistance")

# ── Recode Group labels (WSM → Women) ─────────────────────────────────────────
df_ref <- df_ref %>%
  mutate(Group = recode(Group, "WSM" = "Women"))

# =============================================================================
# FIGURE 2
# Per-person cost (NG+ individuals) by group & ceftriaxone resistance, NO testing
# Bars: one per (Group × ceftriaxone level), faceted by Group
# =============================================================================

fig2_data <- df_ref %>%
  # No testing = "without test" columns; strategy doesn't matter for no-test,
  # so just take ciprofloxacin rows to avoid tripling bars
  filter(Strategy == "CIP") %>%
  pivot_without() %>%
  mutate(
    cro_label = factor(cro_labels[as.character(cro_res)],
                       levels = cro_labels),
    x_label = paste0(Group, "\n", cro_labels[as.character(cro_res)])
  )

fig2 <- ggplot(fig2_data,
               aes(x = factor(cro_res), y = cost, fill = category)) +
  geom_col(width = 0.7, position = "stack") +
  scale_fill_manual(values = cat_colors, name = "Cost category") +
  scale_x_discrete(labels = c("0" = "0%", "0.05" = "5%", "0.15" = "15%")) +
  scale_y_continuous("Cost per person (NG+ individuals)", labels = dollar) +
  facet_wrap(~ Group, nrow = 1) +
  labs(
    x       = "Ceftriaxone Resistance Prevalence",
    title   = "Figure 2. Average per-person cost of gonorrhea care by population group\nand ceftriaxone resistance, disaggregated by cost category",
    subtitle = paste0("Ciprofloxacin resistance fixed at ", percent(CIP_REF),
                      " | No resistance testing | Novel Drug A resistance = 0%")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle      = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position    = "right"
  )

ggsave("Figure2_cost_by_group_NDA0.png", fig2, width = 11, height = 6, dpi = 300)
cat("Saved: Figure2_cost_by_group_NDA0.png\n")

# =============================================================================
# FIGURE 3
# Per-person cost by group & ceftriaxone resistance, comparing testing strategies
# Panel A = MSM, B = MSW, C = Women
# Within each panel: 3 facets (0% / 5% / 15% ceftriaxone) × 3 bars
#                    (No testing / Ciprofloxacin / Ciprofloxacin+Ceftriaxone)
# =============================================================================

# No-testing arm (same across strategies — use CIP row to get "without" cols)
no_test_data <- df_ref %>%
  filter(Strategy == "CIP") %>%
  pivot_without() %>%
  mutate(testing = "No resistance testing")

# With-testing arms (ciprofloxacin and ciprofloxacin+ceftriaxone strategies)
with_test_data <- df_ref %>%
  filter(Strategy %in% c("CIP", "CIP+CRO")) %>%
  pivot_with() %>%
  mutate(testing = recode(Strategy,
                          "CIP"     = "Ciprofloxacin testing",
                          "CIP+CRO" = "Ciprofloxacin+Ceftriaxone testing"))

fig3_data <- bind_rows(no_test_data, with_test_data) %>%
  mutate(
    testing = factor(testing,
                     levels = c("No resistance testing",
                                "Ciprofloxacin testing",
                                "Ciprofloxacin+Ceftriaxone testing")),
    cro_label = factor(cro_labels[as.character(cro_res)],
                       levels = cro_labels),
    Group = factor(Group, levels = c("MSM", "MSW", "Women"))
  )

make_fig3_panel <- function(grp, panel_letter) {
  ggplot(filter(fig3_data, Group == grp),
         aes(x = testing, y = cost, fill = category)) +
    geom_col(width = 0.7, position = "stack") +
    scale_fill_manual(values = cat_colors, name = "Cost category") +
    scale_y_continuous("Cost per person (NG+ individuals)", labels = dollar) +
    scale_x_discrete(labels = c(
      "No resistance testing"             = "No resistance\ntesting",
      "Ciprofloxacin testing"            = "Ciprofloxacin\ntesting",
      "Ciprofloxacin+Ceftriaxone testing" = "Ciprofloxacin+\nCeftriaxone\ntesting"
    )) +
    facet_wrap(~ cro_label, nrow = 1) +
    labs(
      x        = NULL,
      subtitle = paste0(panel_letter, "  —  ", grp)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      strip.text         = element_text(face = "bold", size = 10),
      plot.subtitle      = element_text(face = "bold", size = 11),
      axis.text.x        = element_text(size = 8),
      legend.position    = "none"
    )
}

p_a <- make_fig3_panel("MSM", "A")
p_b <- make_fig3_panel("MSW", "B")
p_c <- make_fig3_panel("Women", "C")

# Shared legend
legend_plot <- ggplot(fig3_data,
                      aes(x = testing, y = cost, fill = category)) +
  geom_col() +
  scale_fill_manual(values = cat_colors, name = "Cost category") +
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))
shared_legend <- cowplot::get_legend(legend_plot)

# Stack panels with shared legend
library(cowplot)
fig3_panels <- plot_grid(p_a, p_b, p_c,
                         ncol = 1, align = "v", axis = "lr")

fig3 <- plot_grid(
  fig3_panels, shared_legend,
  ncol = 2, rel_widths = c(1, 0.22)
)

title_row <- ggdraw() +
  draw_label(
    paste0("Figure 3. Per-person cost of gonorrhea care across population groups\n",
           "and increasing ceftriaxone resistance, by testing strategy\n",
           "(Ciprofloxacin resistance fixed at ", percent(CIP_REF),
           ", Novel Drug A resistance = 0%)"),
    fontface = "bold", size = 11, hjust = 0.5
  )

fig3_final <- plot_grid(title_row, fig3, ncol = 1, rel_heights = c(0.1, 1))

ggsave("Figure3_cost_by_strategy_NDA0.png", fig3_final,
       width = 11, height = 12, dpi = 300)
cat("Saved: Figure3_cost_by_strategy_NDA0.png\n")

cat("\nDone! Both figures saved.\n")
cat("Install cowplot if needed: install.packages('cowplot')\n")