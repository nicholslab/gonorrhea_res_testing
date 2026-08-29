# =============================================================================
# Figure S1 and S2 - tornado diagrams
#
# Standalone: reads gonorrhea_S3_dsa_results_v3.csv, so the model does NOT
# need to be re-run. Plot styling is identical to the model's version.
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

IN_FILE  <- "gonorrhea_S3_dsa_results_v3.csv"
OUT_FILE <- "gonorrhea_S3_tornado_diagram_v4.png"

# Parameters excluded from the figure (inert under an MSM reference case).
DROP_PARAMS <- c("p_pid_asympt", "pid_partner_trunc")

if (!file.exists(IN_FILE)) {
  cat("Working directory:", getwd(), "\n")
  print(list.files(pattern = "\\.csv$"))
  stop("Cannot find ", IN_FILE, " - setwd() to the folder holding the model output.")
}

dsa_df <- read_csv(IN_FILE, show_col_types = FALSE)

# ---- Verify the dropped parameters really are inert before removing them ----
inert_check <- dsa_df %>%
  filter(param %in% DROP_PARAMS) %>%
  mutate(swing = abs(high_thresh - low_thresh))

if (nrow(inert_check) == 0) {
  message("Note: none of the excluded parameters were present in ", IN_FILE)
} else {
  cat("Checking excluded parameters are inert:\n")
  inert_check %>%
    transmute(panel, param, low_thresh = round(low_thresh, 4),
              high_thresh = round(high_thresh, 4), swing = round(swing, 4)) %>%
    print(n = Inf)

  if (max(inert_check$swing) > 0.005) {
    stop("An excluded parameter moved the threshold by $",
         sprintf("%.4f", max(inert_check$swing)),
         ". It is NOT inert and must not be dropped - investigate before plotting.")
  }
  cat("Confirmed inert (max swing $",
      sprintf("%.4f", max(inert_check$swing)), ").\n\n", sep = "")
}

# ---- Build plotting frame ----
tornado_df <- dsa_df %>%
  filter(!param %in% DROP_PARAMS) %>%
  mutate(
    y_label = case_when(
      label_unit == "$" ~ sprintf("%s ($%.0f; $%.0f)", label, low_val, high_val),
      label_unit == "%" ~ sprintf("%s (%.1f%%; %.1f%%)", label, low_val * 100, high_val * 100),
      TRUE              ~ sprintf("%s (%.2f; %.2f)", label, low_val, high_val)
    ),
    bar_left    = pmin(low_thresh, high_thresh),
    bar_right   = pmax(low_thresh, high_thresh),
    total_swing = bar_right - bar_left,
    left_is_low = low_thresh < high_thresh
  )

cat(sprintf("Plotting %d parameters per panel (excluded %d).\n",
            n_distinct(tornado_df$param), length(DROP_PARAMS)))

COL_RIGHT <- "#4682B4"
COL_LEFT  <- "#B05C57"
COL_BAR   <- "#6FA3C7"

make_tornado <- function(data, panel_name, base_thresh_val) {
  param_order <- data %>% arrange(total_swing) %>% pull(y_label)
  data <- data %>% mutate(y_label = factor(y_label, levels = param_order))

  ann_left <- data %>%
    mutate(x = bar_left, txt = sprintf("$%.1f", bar_left), hjust = 1.15,
           col = ifelse(left_is_low, COL_LEFT, COL_RIGHT))
  ann_right <- data %>%
    mutate(x = bar_right, txt = sprintf("$%.1f", bar_right), hjust = -0.15,
           col = ifelse(left_is_low, COL_RIGHT, COL_LEFT))
  ann <- bind_rows(ann_left, ann_right)

  ggplot(data, aes(y = y_label)) +
    geom_segment(aes(x = bar_left, xend = bar_right, yend = y_label),
                 linewidth = 9, color = COL_BAR, alpha = 0.75) +
    geom_vline(xintercept = base_thresh_val,
               color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_text(data = ann, aes(x = x, y = y_label, label = txt, hjust = hjust, color = col),
              size = 3.2, fontface = "bold", inherit.aes = FALSE) +
    scale_color_identity() +
    scale_x_continuous(labels = scales::dollar, expand = expansion(mult = 0.15)) +
    labs(
      title = sprintf("One-Way Sensitivity Analysis: %s Panel", panel_name),
      subtitle = sprintf(
        "Reference: MSM, CIP=50%%, CRO=5%%, NDA=0%%\nBase case threshold: $%.1f  |  Blue = higher value; Red = lower value",
        base_thresh_val),
      x = "Threshold Price (USD)", y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 10),
          plot.subtitle = element_text(size = 9, color = "grey40"))
}

# base_thresh is constant within panel; take the first value
base_cip    <- tornado_df %>% filter(panel == "CIP")         %>% slice(1) %>% pull(base_thresh)
base_triple <- tornado_df %>% filter(panel == "CIP+CRO+NDA") %>% slice(1) %>% pull(base_thresh)

p_cip    <- make_tornado(filter(tornado_df, panel == "CIP"),
                         "CIP-only", base_cip)
p_triple <- make_tornado(filter(tornado_df, panel == "CIP+CRO+NDA"),
                         "CIP+CRO+NDA", base_triple)

p_combined <- p_cip / p_triple +
  plot_annotation(
    title = "Figure S3. Deterministic Sensitivity Analysis - Tornado Diagrams",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

# Six parameters instead of eight, so the panel is a little shorter.
ok <- tryCatch({
  ggsave(OUT_FILE, p_combined, width = 11, height = 9, dpi = 300); TRUE
}, error = function(e) { message("FAILED to write ", OUT_FILE, "\n  ",
                                 conditionMessage(e)); FALSE })
if (ok && file.exists(OUT_FILE)) {
  cat(sprintf("Saved: %s (%.0f KB)\n",
              normalizePath(OUT_FILE), file.size(OUT_FILE) / 1024))
}

# ---- Values for the supplement text ----
cat("\nParameter ranges, ordered by influence:\n")
tornado_df %>%
  arrange(panel, desc(total_swing)) %>%
  transmute(panel, label,
            range = sprintf("$%.1f to $%.1f", low_thresh, high_thresh),
            swing = sprintf("$%.1f", total_swing)) %>%
  print(n = Inf)
