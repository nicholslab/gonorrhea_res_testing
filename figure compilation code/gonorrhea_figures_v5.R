# ============================================================
# Figures 2, 3 & 4 (v5) — as submitted to PLOS Medicine (R2)
# Inputs (same folder): gonorrhea_summary_table_v3.csv,
#                       gonorrhea_thresholds_weighted_v3.csv
# Outputs: Figure 2.tiff, Figure 3.tiff, Figure 4.tiff (300 dpi)
# Figure titles and legends are in the manuscript, not in the images.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ============================================================
# FIGURES 2 & 3: COST BREAKDOWN
# ============================================================

CIP_REF  <- 0.5               # ciprofloxacin resistance (reference scenario)
CRO_SHOW <- c(0, 0.05, 0.15)  # ceftriaxone resistance levels shown
NDA_REF  <- 0                 # novel drug A resistance. The model emits 0% and 5%
                              # in one run; without this filter every bar would be
                              # stacked twice.

# NOTE on partner complications: the model adds each partner complication cost
# to the total cost once, but records it in BOTH c_partner_total and
# c_complications. Summing all raw cost columns therefore double-counts it.
# Here the overlap (= sum of raw components - reported total) is split out as
# its own category, so partner treatment and index complications exclude it and
# every stack sums exactly to the reported total.
cat_colors <- c(
  "Initial visit"           = "#B8A9C9",
  "Initial diagnostic"      = "#6BAED6",
  "Ceftriaxone"             = "#D9795A",
  "Ciprofloxacin"           = "#E0A96D",
  "Novel drug A"            = "#C994B0",
  "Follow-up visit"         = "#74A9C8",
  "Follow-up diagnostic"    = "#74B49B",
  "Follow-up medication"    = "#9DC88D",
  "Partner treatment"       = "#C97A9A",
  "Complications (index)"   = "#A8B86A",
  "Complications (partner)" = "#5E7A3A"
)

df_costs <- read_csv("gonorrhea_summary_table_v3.csv", show_col_types = FALSE) %>%
  mutate(Group    = recode(group, "MSM" = "MSM", "MSW" = "MSW", "WSM" = "Women"),
         Strategy = toupper(strategy)) %>%
  filter(cip_res == CIP_REF, cro_res %in% CRO_SHOW, nda_res == NDA_REF)

.dup <- df_costs %>% count(Group, Strategy, cro_res) %>% filter(n > 1)
if (nrow(.dup) > 0) { print(.dup); stop("Duplicate rows - bars would be double-counted.") }

# arm = "without" or "with" (resistance testing)
pivot_costs <- function(d, arm) {
  p   <- function(x) d[[sprintf("cost_%s_test_ng_%s_mean", arm, x)]]
  tot <- d[[sprintf("cost_%s_test_mean_ng", arm)]]
  ov  <- p("init_visit") + p("init_diag") + p("cro") + p("cip") + p("nda") +
         p("follow_visit") + p("follow_diag") + p("follow_drugs") +
         p("partner_total") + p("complications") - tot
  tibble(Group = d$Group, Strategy = d$Strategy, cro_res = d$cro_res, total = tot,
    `Initial visit`           = p("init_visit"),
    `Initial diagnostic`      = p("init_diag"),
    Ceftriaxone               = p("cro"),
    Ciprofloxacin             = p("cip"),
    `Novel drug A`            = p("nda"),
    `Follow-up visit`         = p("follow_visit"),
    `Follow-up diagnostic`    = p("follow_diag"),
    `Follow-up medication`    = p("follow_drugs"),
    `Partner treatment`       = p("partner_total") - ov,
    `Complications (index)`   = p("complications") - ov,
    `Complications (partner)` = ov) %>%
    pivot_longer(all_of(names(cat_colors)), names_to = "category", values_to = "cost") %>%
    mutate(category = factor(category, levels = names(cat_colors)))
}

theme_fig <- theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
y_cost <- scale_y_continuous("Cost per person diagnosed with gonorrhea (2025 USD)",
                             labels = dollar, expand = expansion(mult = c(0, .03)))

# ---- Figure 2: no resistance testing ----
fig2_data <- df_costs %>% filter(Strategy == "CIP") %>% pivot_costs("without") %>%
  mutate(cro = factor(paste0(cro_res * 100, "%"), levels = c("0%", "5%", "15%")))

fig2 <- ggplot(fig2_data, aes(cro, cost, fill = category)) +
  geom_col(width = 0.6, position = position_stack(reverse = TRUE)) +
  facet_wrap(~ Group, nrow = 1) +
  scale_fill_manual("Cost category", values = cat_colors) +
  y_cost + labs(x = "Prevalence of ceftriaxone resistance") + theme_fig
ggsave("Figure 2.tiff", fig2, width = 10, height = 5.5, dpi = 300, compression = "lzw", bg = "white")

# ---- Figure 3: by testing strategy ----
fig3_data <- bind_rows(
  df_costs %>% filter(Strategy == "CIP") %>% pivot_costs("without") %>%
    mutate(arm = "No testing"),
  df_costs %>% filter(Strategy %in% c("CIP", "CIP+CRO")) %>% pivot_costs("with") %>%
    mutate(arm = recode(Strategy, "CIP" = "CIP testing", "CIP+CRO" = "CIP+CRO testing"))) %>%
  mutate(arm   = factor(arm, levels = c("No testing", "CIP testing", "CIP+CRO testing")),
         Group = factor(recode(Group, MSM = "A. MSM", MSW = "B. MSW", Women = "C. Women")),
         cro   = factor(paste0(cro_res * 100, "% ceftriaxone resistance"),
                        levels = paste0(c(0, 5, 15), "% ceftriaxone resistance")))

fig3 <- ggplot(fig3_data, aes(arm, cost, fill = category)) +
  geom_col(width = 0.6, position = position_stack(reverse = TRUE)) +
  facet_grid(Group ~ cro, scales = "free_y") +
  scale_fill_manual("Cost category", values = cat_colors) +
  y_cost + labs(x = NULL) + theme_fig +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9))
ggsave("Figure 3.tiff", fig3, width = 11, height = 9, dpi = 300, compression = "lzw", bg = "white")

# Guard: every stack must sum to the reported total
recon <- bind_rows(fig2_data %>% mutate(arm = "none"), fig3_data) %>%
  group_by(Group, cro_res, arm) %>%
  summarise(err = abs(sum(cost) - first(total)), .groups = "drop")
if (max(recon$err) > 0.01) stop("Stacked categories do not sum to the reported total.")
cat(sprintf("Figures 2-3 saved; stacks reconcile to totals (max error $%.4f)\n", max(recon$err)))

# ============================================================
# FIGURE 4: population-weighted cost-neutral threshold prices
# Weights (applied in the model): MSM 0.53, MSW 0.21, women 0.26
# Uses the model's weighted-threshold file, which weights per draw and then
# takes quantiles, so 95% UIs are valid (values provided in S1 Data).
# ============================================================

fig4_data <- read_csv("gonorrhea_thresholds_weighted_v3.csv", show_col_types = FALSE) %>%
  filter(nda_res %in% c(0, 0.05)) %>%
  mutate(Strategy  = factor(toupper(strategy), levels = c("CIP", "CIP+CRO", "CIP+CRO+NDA")),
         NDA_label = factor(ifelse(nda_res == 0, "0% NDA resistance", "5% NDA resistance")),
         label     = sprintf("%.0f", threshold_price))

.dup4 <- fig4_data %>% count(cip_res, cro_res, nda_res, Strategy) %>% filter(n > 1)
if (nrow(.dup4) > 0) { print(.dup4); stop("Duplicate tiles in Figure 4.") }

# Values quoted in the Results text (50% CIP, 5% CRO)
fig4_data %>% filter(cip_res == 0.5, cro_res == 0.05) %>%
  transmute(NDA_label, Strategy,
            text = sprintf("$%.0f (95%% UI [$%.0f, $%.0f])",
                           threshold_price, threshold_lower, threshold_upper)) %>%
  arrange(NDA_label, Strategy) %>% print(n = Inf)

# Colour limit derived from the data so no tiles are clipped
lim_abs <- ceiling(max(abs(fig4_data$threshold_price), na.rm = TRUE) / 5) * 5

fig4 <- ggplot(fig4_data, aes(cip_res, cro_res, fill = threshold_price)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label, color = abs(threshold_price) < lim_abs * 0.3),
            size = 3, fontface = "bold") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
  scale_fill_gradient2(low = "#E69F00", mid = "white", high = "#0072B2", midpoint = 0,
                       limits = c(-lim_abs, lim_abs), oob = squish,
                       name = "Threshold\nPrice (USD)", labels = dollar_format(accuracy = 1)) +
  scale_x_continuous("CIP Resistance Prevalence", breaks = c(0.1, 0.3, 0.5, 0.7),
                     labels = percent_format(accuracy = 1)) +
  scale_y_continuous("CRO Resistance Prevalence",
                     breaks = c(0, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15),
                     labels = percent_format(accuracy = 0.1), expand = c(0, 0)) +
  facet_grid(Strategy ~ NDA_label) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = 11),
        legend.key.height = unit(1.5, "cm"))
ggsave("Figure 4.tiff", fig4, width = 11, height = 9, dpi = 300, compression = "lzw", bg = "white")
cat("Figure 4 saved\n")
