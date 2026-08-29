# ============================================================
# Updated Figures 2, 3, & 4
# Uses new CSV outputs directly
# Ready to run: no external file dependencies
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(cowplot)
})

# ============================================================
# LOAD DATA
# ============================================================

# Read the summary table (contains cost breakdowns for Figures 2 & 3)
df_summary <- read_csv("gonorrhea_summary_table_v3.csv",
                       show_col_types = FALSE)

# Read the thresholds (for Figure 4)
# Note: The uploaded file should contain both 0% and 5% NDA resistance scenarios
df_thresholds <- read_csv("gonorrhea_thresholds_v3.csv",
                          show_col_types = FALSE)

# ============================================================
# FIGURE 2 & 3: COST BREAKDOWN DATA
# ============================================================

# Fix ciprofloxacin resistance at 50% (reference scenario used in original figures)
CIP_REF <- 0.5

# Ceftriaxone resistance levels to show (0%, 5%, 15%)
CRO_SHOW <- c(0, 0.05, 0.15)

# NDA resistance to show in Figures 2 & 3.
# The model now emits BOTH 0% and 5% NDA in a single run. Figures 2 and 3 show
# one NDA scenario at a time, so this MUST be filtered — without it every
# (Group, Strategy, cro_res) combination matches two rows and geom_col()
# silently stacks them, drawing every bar at roughly twice its true height.
NDA_REF <- 0

# Cost category colours — match original figure exactly
#
# NOTE (v4): "Partner treatment" is now partner care EXCLUDING complications.
# In the model, partner complication cost is added to BOTH c_partner_total and
# c_complications (it is only added to `cost` once, so scenario totals were
# always right). Stacking the two raw columns therefore double-counted partner
# complications and inflated every bar — by ~$209 of a $529 total for MSW.
# The pivots below subtract the overlap so the stack sums to the reported total.
cat_colors <- c(
  "Initial visit"              = "#B8A9C9",   # dusty lavender
  "Initial diagnostic"         = "#6BAED6",   # steel blue
  "Ceftriaxone"                = "#D9795A",   # terracotta
  "Ciprofloxacin"              = "#E0A96D",   # warm sand
  "Novel Drug A"               = "#C994B0",   # mauve
  "Follow-up visit"            = "#74A9C8",   # muted sky blue
  "Follow-up diagnostic"       = "#74B49B",   # sage green
  "Follow-up meds"             = "#9DC88D",   # soft green
  "Partner care (excl. complications)" = "#C97A9A",   # dusty rose
  "Complications (index + partner)"    = "#A8B86A"    # muted olive
)

cat_order <- names(cat_colors)

# ── Prepare summary table: recode group & strategy names ───────────────────────
df_costs <- df_summary %>%
  mutate(
    Group = recode(group, "MSM" = "MSM", "MSW" = "MSW", "WSM" = "Women"),
    Strategy = toupper(strategy)
  ) %>%
  filter(cip_res == CIP_REF, cro_res %in% CRO_SHOW, nda_res == NDA_REF)

# Guard: exactly one row per Group x Strategy x cro_res. If this ever fires,
# a filter is missing and the stacked bars would be silently inflated.
.dup <- df_costs %>% count(Group, Strategy, cro_res) %>% filter(n > 1)
if (nrow(.dup) > 0) {
  print(.dup)
  stop("df_costs has duplicate rows per Group/Strategy/cro_res - bars would be double-counted.")
}
cat(sprintf("Figures 2-3: %d rows, NDA resistance = %s\n",
            nrow(df_costs), scales::percent(NDA_REF)))

# ── Helper: pivot WITHOUT-test cost categories ────────────────────────────────
# The overlap (partner complications, counted in both partner_total and
# complications) is recovered as (sum of raw components - reported total) and
# subtracted from partner care, so the stack reconciles to the total exactly.
pivot_without <- function(df_in) {
  df_in %>%
    mutate(
      .overlap = cost_without_test_ng_init_visit_mean +
                 cost_without_test_ng_init_diag_mean +
                 cost_without_test_ng_cro_mean +
                 cost_without_test_ng_cip_mean +
                 cost_without_test_ng_nda_mean +
                 cost_without_test_ng_follow_visit_mean +
                 cost_without_test_ng_follow_diag_mean +
                 cost_without_test_ng_follow_drugs_mean +
                 cost_without_test_ng_partner_total_mean +
                 cost_without_test_ng_complications_mean -
                 cost_without_test_mean_ng,
      .partner_net = cost_without_test_ng_partner_total_mean - .overlap
    ) %>%
    select(Group, Strategy, cip_res, cro_res,
           "Initial visit"        = cost_without_test_ng_init_visit_mean,
           "Initial diagnostic"   = cost_without_test_ng_init_diag_mean,
           "Ceftriaxone"          = cost_without_test_ng_cro_mean,
           "Ciprofloxacin"        = cost_without_test_ng_cip_mean,
           "Novel Drug A"         = cost_without_test_ng_nda_mean,
           "Follow-up visit"      = cost_without_test_ng_follow_visit_mean,
           "Follow-up diagnostic" = cost_without_test_ng_follow_diag_mean,
           "Follow-up meds"       = cost_without_test_ng_follow_drugs_mean,
           "Partner care (excl. complications)" = .partner_net,
           "Complications (index + partner)"    = cost_without_test_ng_complications_mean
    ) %>%
    pivot_longer(cols = all_of(cat_order),
                 names_to = "category", values_to = "cost") %>%
    mutate(category = factor(category, levels = cat_order),
           testing = "No resistance testing")
}

# ── Helper: pivot WITH-test cost categories ───────────────────────────────────
pivot_with <- function(df_in) {
  df_in %>%
    mutate(
      .overlap = cost_with_test_ng_init_visit_mean +
                 cost_with_test_ng_init_diag_mean +
                 cost_with_test_ng_cro_mean +
                 cost_with_test_ng_cip_mean +
                 cost_with_test_ng_nda_mean +
                 cost_with_test_ng_follow_visit_mean +
                 cost_with_test_ng_follow_diag_mean +
                 cost_with_test_ng_follow_drugs_mean +
                 cost_with_test_ng_partner_total_mean +
                 cost_with_test_ng_complications_mean -
                 cost_with_test_mean_ng,
      .partner_net = cost_with_test_ng_partner_total_mean - .overlap
    ) %>%
    select(Group, Strategy, cip_res, cro_res,
           "Initial visit"        = cost_with_test_ng_init_visit_mean,
           "Initial diagnostic"   = cost_with_test_ng_init_diag_mean,
           "Ceftriaxone"          = cost_with_test_ng_cro_mean,
           "Ciprofloxacin"        = cost_with_test_ng_cip_mean,
           "Novel Drug A"         = cost_with_test_ng_nda_mean,
           "Follow-up visit"      = cost_with_test_ng_follow_visit_mean,
           "Follow-up diagnostic" = cost_with_test_ng_follow_diag_mean,
           "Follow-up meds"       = cost_with_test_ng_follow_drugs_mean,
           "Partner care (excl. complications)" = .partner_net,
           "Complications (index + partner)"    = cost_with_test_ng_complications_mean
    ) %>%
    pivot_longer(cols = all_of(cat_order),
                 names_to = "category", values_to = "cost") %>%
    mutate(category = factor(category, levels = cat_order),
           testing = Strategy)
}

# Guard: stacked categories must reconcile to the reported scenario total.
# This protects BOTH Figure 2 and Figure 3 - they share pivot_without().
.recon <- df_costs %>%
  pivot_without() %>%
  group_by(Group, Strategy, cro_res) %>%
  summarise(stacked = sum(cost), .groups = "drop") %>%
  left_join(df_costs %>% select(Group, Strategy, cro_res,
                                reported = cost_without_test_mean_ng),
            by = c("Group", "Strategy", "cro_res")) %>%
  mutate(err = abs(stacked - reported))
if (max(.recon$err) > 0.01) {
  print(arrange(.recon, desc(err)))
  stop("Stacked cost categories do not sum to the scenario total.")
}
cat(sprintf("Figures 2-3 stack reconciles to total (max error $%.4f)\n", max(.recon$err)))

# Suffix and label driven by NDA_REF so filenames can never disagree with
# their contents. Previously both were hard-coded to "NDA0"/"0%", so setting
# NDA_REF <- 0.05 wrote 5% data into a file named NDA0 with a subtitle saying 0%.
NDA_SUFFIX <- sprintf("NDA%g", NDA_REF * 100)
NDA_TEXT   <- sprintf("Novel Drug A resistance = %s", percent(NDA_REF))

# ── x-axis labels ─────────────────────────────────────────────────────────────
cro_labels <- c("0"    = "0% Ceftriaxone resistance",
                "0.05" = "5% Ceftriaxone resistance",
                "0.15" = "15% Ceftriaxone resistance")

# =============================================================================
# FIGURE 2
# Per-person cost (NG+ individuals) by group & ceftriaxone resistance, NO testing
# =============================================================================

fig2_data <- df_costs %>%
  # No testing = "without test" columns; strategy doesn't matter for no-test,
  # so just take CIP rows to avoid tripling bars
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
                      " | No resistance testing | ", NDA_TEXT)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.x    = unit(0.5, "lines"),
    strip.text         = element_text(face = "bold"),
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle      = element_text(hjust = 0.5, size = 9, color = "grey40"),
    plot.margin        = margin(t = 10, r = 10, b = 80, l = 10, unit = "pt"),
    legend.position    = "right"
  )

.f2 <- sprintf("Figure2_cost_by_group_%s.png", NDA_SUFFIX)
ggsave(.f2, fig2, width = 8, height = 7, dpi = 300)
cat(sprintf("Saved: %s\n", .f2))

# =============================================================================
# FIGURE 3
# Per-person cost by group & ceftriaxone resistance, comparing testing strategies
# Panel A = MSM, B = MSW, C = Women
# =============================================================================

# No-testing arm (same across strategies — use CIP row)
no_test_data <- df_costs %>%
  filter(Strategy == "CIP") %>%
  pivot_without() %>%
  mutate(testing = "No resistance testing")

# With-testing arms (ciprofloxacin and ciprofloxacin+ceftriaxone strategies)
with_test_data <- df_costs %>%
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
shared_legend <- get_legend(legend_plot)

# Stack panels with shared legend
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
           ", ", NDA_TEXT, ")"),
    fontface = "bold", size = 11, hjust = 0.5
  )

fig3_final <- plot_grid(title_row, fig3, ncol = 1, rel_heights = c(0.1, 1))

.f3 <- sprintf("Figure3_cost_by_strategy_%s.png", NDA_SUFFIX)
ggsave(.f3, fig3_final, width = 11, height = 12, dpi = 300)
cat(sprintf("Saved: %s\n", .f3))

# =============================================================================
# FIGURE 4
# Weighted-average cost-neutrality threshold prices
# Weights: MSM=0.53, MSW=0.21, Women=0.26
# Layout: Strategy (rows) × NDA scenario (columns)
# =============================================================================

# (v4) Figure 4 now reads gonorrhea_thresholds_weighted_v3.csv, produced by the
# model itself. That file weights the per-draw thresholds and THEN takes
# quantiles, so it carries a valid 95% UI. Re-weighting point estimates here
# gave the identical central value (the mean is linear) but no usable interval,
# and there is no correct way to recover one from group-level bounds.
df_weighted <- read_csv("gonorrhea_thresholds_weighted_v3.csv",
                        show_col_types = FALSE)

fig4_data <- df_weighted %>%
  filter(nda_res %in% c(0, 0.05)) %>%
  mutate(
    Strategy  = factor(toupper(strategy),
                       levels = c("CIP", "CIP+CRO", "CIP+CRO+NDA")),
    nda_label = ifelse(nda_res == 0, "0% NDA resistance", "5% NDA resistance"),
    NDA_label = factor(nda_label,
                       levels = c("0% NDA resistance", "5% NDA resistance")),
    label     = sprintf("%.0f", threshold_price)
  )

# Guard: one tile per cip x cro x nda x strategy.
.dup4 <- fig4_data %>% count(cip_res, cro_res, nda_res, Strategy) %>% filter(n > 1)
if (nrow(.dup4) > 0) {
  print(.dup4)
  stop("fig4_data has duplicate tiles - the heatmap would overplot.")
}
cat(sprintf("Figure 4: %d tiles from the model's weighted-threshold file\n",
            nrow(fig4_data)))

# Reference values quoted in the Results text (50% CIP, 5% CRO).
fig4_data %>%
  filter(cip_res == 0.5, cro_res == 0.05) %>%
  transmute(nda_label, Strategy,
            text = sprintf("$%.0f (95%% UI $%.0f-$%.0f)",
                           threshold_price, threshold_lower, threshold_upper)) %>%
  arrange(nda_label, Strategy) %>%
  print(n = Inf)

# Axis breaks
x_breaks <- c(0.1, 0.3, 0.5, 0.7)
y_breaks <- c(0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15)

# Diverging colour scale.
# (v4) Derived from the data rather than hard-coded. Under v3 the weighted
# thresholds reach $47.6, which the previous fixed limit of 45 would have
# clipped — the highest-resistance tiles would all have rendered the same
# saturated colour despite differing by several dollars.
lim_abs <- ceiling(max(abs(fig4_data$threshold_price), na.rm = TRUE) / 5) * 5
cat(sprintf("Figure 4 colour limit: +/- $%d (data max |value| = $%.1f)\n",
            lim_abs, max(abs(fig4_data$threshold_price), na.rm = TRUE)))

fig4 <- ggplot(fig4_data, aes(x = cip_res, y = cro_res, fill = threshold_price)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = label,
        color  = abs(threshold_price) < lim_abs * 0.3),
    size = 3, fontface = "bold"
  ) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(
    low      = "#E69F00",   # Okabe-Ito orange/red (negative = testing adds cost)
    mid      = "white",
    high     = "#0072B2",   # Okabe-Ito blue (positive = testing saves money)
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
    subtitle = "Population-weighted average (MSM 53%, MSW 21%, Women 26%) | Blue = testing saves money; Red = testing adds cost",
    caption  = paste0(
      "Threshold price: maximum per-test cost at which resistance-guided therapy is cost-neutral vs. empiric ceftriaxone.\n",
      "Negative values indicate scenarios where resistance testing increases expected costs (testing not cost-neutral at any price).\n",
      "Results based on 1,000 PSA draws × 10,000 individuals per scenario."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
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

ggsave("Figure4_weighted_threshold_heatmap.png", fig4,
       width = 11, height = 9, dpi = 300)
cat("✓ Saved: Figure4_weighted_threshold_heatmap.png\n")
cat(sprintf("\nColor scale: -$%.0f to +$%.0f (diverging, centered at 0)\n", lim_abs, lim_abs))
cat("Weights applied: MSM=0.53, MSW=0.21, Women=0.26\n")
cat("\n✓ All figures generated successfully!\n")
