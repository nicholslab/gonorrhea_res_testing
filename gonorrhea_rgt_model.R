# =============================================================================
# Economic value of gonorrhea resistance testing in the United States
# Individual-level stochastic simulation model
#
# This script reproduces the per-episode cost estimates, cost-neutral threshold
# prices (Figures 2-4, Table 2 dimension A), and the one-way deterministic
# sensitivity analysis (Supplement S3) reported in the manuscript.
#
# The system-level and drug-development value dimensions (Table 2 B, C, D) are
# derived from these outputs together with published transmission-model
# trajectories; see the README and Supplement S2 for those calculations.
#
# Reproducibility: the global seed is set below (set.seed(123)); the sensitivity
# analysis uses its own per-row seeds so that low/base/high estimates differ only
# in the parameter being varied. Runtime is on the order of a few hours on a
# desktop machine (1,000 PSA draws x 20,000 individuals per scenario).
#
# Author: Brooke E. Nichols
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

set.seed(123)

# ------------------------------
# controls
# ------------------------------
n_psa      <- 1000     # PSA parameter draws
n_indiv    <- 20000    # individuals per PSA per (group, resistance, strategy)
nda_res_fixed <- 0.00  # fixed NDA resistance for heatmap (tunable)

# resistance grids (CIP x CRO) for heatmaps
cip_res_vals <- c(0.10, 0.30, 0.50, 0.70)
cro_res_vals <- c(0.00, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15)

# groups and NG prevalence among those tested
groups <- c("MSM","MSW","WSM")
ng_prev_by_group <- c(MSM = 0.0454, MSW = 0.0244, WSM = 0.0103)

# strategies (which resistance results are available/used programmatically)
resistance_panels <- list(
  `CIP`         = c("cip"),
  `CIP+CRO`     = c("cip","cro"),
  `CIP+CRO+NDA` = c("cip","cro","nda")
)

# Resistance test performance
sens_res_test <- list(cip = 0.99, cro = 0.99, nda = 0.99)
spec_res_test <- list(cip = 0.98, cro = 0.98, nda = 0.98)

# ------------------------------
# Draw PSA parameters
# ------------------------------
draw_psa <- function() {
  # NG test sensitivity (specificity ~1 unless noted)
  sens_ng_men   <- rbeta(1, 99.64, 6.36)       # MSM + MSW
  sens_ng_women <- rbeta(1, 11.89, 0.50)       # WSM
  spec_ng_men   <- 0.999
  spec_ng_women <- 0.999
  
  # Return if symptomatic
  p_return_sympt    <- rbeta(1, 54.45, 6.05)     # mean ~0.9
  
  # Partner treatment probability (engagement)
  p_partner_treat   <- rbeta(1, 10.47, 24.43)    # mean ~0.30
  
  # Probability symptomatic by group (asymptomatic = 1 - symptomatic)
  p_sympt <- c(
    MSW = rbeta(1, 5.7, 1.6),
    MSM = rbeta(1, 8.0, 3.8),
    WSM = rbeta(1, 9.2, 13.6)
  )
  
  # Complication probabilities
  # Women: PID risk depends on symptom status
  p_pid_sympt      <- runif(1, 0.00092, 0.0055)
  p_pid_asympt     <- runif(1, 0.037,  0.16)
  p_cpp_given_pid  <- runif(1, 0.23,   0.29)
  p_ep_given_pid   <- runif(1, 0.049,  0.098)
  p_ti_given_pid   <- runif(1, 0.12,   0.23)
  
  # Men (if untreated)
  p_epi_untreated  <- runif(1, 0.0012, 0.14)
  p_dgi_untreated  <- runif(1, 0.0075, 0.013)
  
  # Complication costs
  c_pid <- runif(1, 1081, 3090)
  c_cpp <- runif(1,  788, 2245)
  c_ep  <- runif(1, 3013, 8689)
  c_ti  <- runif(1, 4108, 11710)
  c_epi <- runif(1,  302,  877)
  c_dgi <- runif(1, 1308, 5264)
  
  # Base service costs (fixed)
  visit_init          <- 167.10
  diag_men            <- 133.49
  diag_women          <- 435.39
  visit_followup      <- 125.18
  diag_followup_men   <- 154.83
  diag_followup_women <- 495.64
  
  # Drugs
  cost_cip      <- runif(1, 0.21, 0.96)
  cost_cro      <- 1.11                 # ceftriaxone
  cost_cefixime <- runif(1, 14.26, 31.09)  # retained for scenario use; not on the base-case pathway
  cost_nda      <- 50
  cost_ndb      <- 50
  
  list(
    sens_ng_men = sens_ng_men,
    sens_ng_women = sens_ng_women,
    spec_ng_men = spec_ng_men,
    spec_ng_women = spec_ng_women,
    p_return_sympt = p_return_sympt,
    p_partner_treat = p_partner_treat,
    p_sympt = p_sympt,
    
    # Complications
    p_pid_sympt = p_pid_sympt,
    p_pid_asympt = p_pid_asympt,
    p_cpp_given_pid = p_cpp_given_pid,
    p_ep_given_pid  = p_ep_given_pid,
    p_ti_given_pid  = p_ti_given_pid,
    p_epi_untreated = p_epi_untreated,
    p_dgi_untreated = p_dgi_untreated,
    
    # Costs
    visit_init = visit_init,
    diag_men = diag_men,
    diag_women = diag_women,
    visit_followup = visit_followup,
    diag_followup_men = diag_followup_men,
    diag_followup_women = diag_followup_women,
    c_pid = c_pid, c_cpp = c_cpp, c_ep = c_ep, c_ti = c_ti, c_epi = c_epi, c_dgi = c_dgi,
    cost_cip = cost_cip, cost_cro = cost_cro, cost_cefixime = cost_cefixime,
    cost_nda = cost_nda, cost_ndb = cost_ndb
  )
}

# ----------------------------------
# Drug selection helpers
# ----------------------------------
select_first_drug <- function(panel, cip_test_pos, cro_test_pos, nda_test_pos) {
  # No resistance testing -> default ceftriaxone
  if (length(panel) == 0) return("cro")
  # Prefer CIP, then CRO, then NDA if each appears susceptible (test negative)
  if ("cip" %in% panel && !cip_test_pos) return("cip")
  if ("cro" %in% panel && !cro_test_pos) return("cro")
  if ("nda" %in% panel && !nda_test_pos) return("nda")
  # If all tested appear resistant, panel-specific fallbacks
  if (identical(sort(panel), "cip")) return("cro")
  if (identical(sort(panel), c("cip","cro"))) return("nda")
  if (identical(sort(panel), c("cip","cro","nda"))) return("ndb")
  "cro"
}

select_next_drug <- function(prev_drug, cip_r, cro_r, nda_r) {
  # Salvage cascade is clinical and resistance-driven only.
  # panel is NOT used here — what was tested doesn't change
  # what drugs are available for salvage.
  # CIP fail -> CRO (if susceptible) -> NDA (if susceptible) -> NDB
  # CRO fail -> NDA (if susceptible) -> NDB
  # NDA fail -> NDB (last resort; assumed always susceptible)
  if (prev_drug == "cip") {
    if (!cro_r) return("cro")
    if (!nda_r) return("nda")
    return("ndb")
  }
  if (prev_drug == "cro") {
    if (!nda_r) return("nda")
    return("ndb")
  }
  if (prev_drug == "nda") {
    return("ndb")
  }
  return(NA_character_)
}

drug_cost <- function(drug, psa) {
  switch(drug,
         cip = psa$cost_cip,
         cro = psa$cost_cro,
         nda = psa$cost_nda,
         ndb = psa$cost_ndb,
         0)
}

drug_success <- function(drug, cip_r, cro_r, nda_r) {
  switch(drug,
         cip = !cip_r,
         cro = !cro_r,
         nda = !nda_r,
         ndb = TRUE,
         FALSE)
}

# ------------------------------
# Microsim for one PSA draw, group, resistance, and strategy
# Returns mean cost (overall) + per-NG+ cost & categories
# ------------------------------
simulate_once <- function(grp, cip_res, cro_res, nda_res, panel, psa, n = 5000, test_cost = 0) {
  is_woman <- grp == "WSM"
  is_male  <- !is_woman
  ng_prev  <- ng_prev_by_group[[grp]]
  
  # Index symptoms
  sympt <- rbinom(n, 1, psa$p_sympt[[grp]]) == 1
  
  # Test perf by sex
  sens_ng <- if (is_male) psa$sens_ng_men else psa$sens_ng_women
  spec_ng <- if (is_male) psa$spec_ng_men else psa$spec_ng_women
  
  # True infection + test outcome
  true_ng     <- rbinom(n, 1, ng_prev) == 1
  test_ng_pos <- rbinom(n, 1, ifelse(true_ng, sens_ng, 1 - spec_ng)) == 1
  
  # True resistance (index & partner share)
  cip_r <- rbinom(n, 1, cip_res) == 1
  cro_r <- rbinom(n, 1, cro_res) == 1
  nda_r <- rbinom(n, 1, nda_res) == 1
  
  # Observed resistance test results (if tested)
  cip_test_pos <- rbinom(n, 1, ifelse(cip_r, sens_res_test$cip, 1 - spec_res_test$cip)) == 1
  cro_test_pos <- rbinom(n, 1, ifelse(cro_r, sens_res_test$cro, 1 - spec_res_test$cro)) == 1
  nda_test_pos <- rbinom(n, 1, ifelse(nda_r, sens_res_test$nda, 1 - spec_res_test$nda)) == 1
  
  # Cost-category trackers (per individual)
  # Index, initial
  c_init_visit <- rep(psa$visit_init, n)
  diag_cost    <- if (is_woman) psa$diag_women else psa$diag_men
  c_init_diag  <- rep(diag_cost, n)
  
  # Index, follow-up
  c_follow_visit <- numeric(n)
  c_follow_diag  <- numeric(n)
  c_follow_drugs <- numeric(n)   # any second+ line drugs (incl. NDB)
  
  # Index drug categories (first-line CIP, CRO, NDA)
  c_drug_cip <- numeric(n)
  c_drug_cro <- numeric(n)
  c_drug_nda <- numeric(n)
  
  # Partner total (all partner-related costs)
  c_partner_total <- numeric(n)
  
  # Complications (index + partner)
  c_complications <- numeric(n)
  
  # Total cost (index + partner)
  cost <- c_init_visit + c_init_diag
  
  # Partner group mapping (determines partner symptom probability)
  partner_group <- if (grp == "MSM") {
    "MSM"   # male partner
  } else if (grp == "MSW") {
    "WSM"   # female partner
  } else {   # grp == "WSM"
    "MSW"   # male partner
  }
  partner_is_male <- partner_group %in% c("MSM","MSW")
  partner_sympt_p <- psa$p_sympt[[partner_group]]
  
  index_final_treated   <- rep(FALSE, n)
  partner_final_treated <- rep(FALSE, n)
  partner_engaged       <- rbinom(n, 1, psa$p_partner_treat) == 1
  
  # Resistance test cost applies to any NG test-positive (true or false)
  cost <- cost + ifelse(test_ng_pos, test_cost, 0)
  test_pos_share <- mean(test_ng_pos)
  
  # Treat if test positive (true or false positive)
  treat_now <- test_ng_pos
  
  for (i in seq_len(n)) {
    if (treat_now[i]) {
      # ---------- INDEX: FIRST-LINE ----------
      drug1 <- select_first_drug(panel, cip_test_pos[i], cro_test_pos[i], nda_test_pos[i])
      d1    <- drug_cost(drug1, psa)
      cost[i] <- cost[i] + d1
      
      # First-line drug categories for index
      if (drug1 == "cip") {
        c_drug_cip[i] <- c_drug_cip[i] + d1
      } else if (drug1 == "cro") {
        c_drug_cro[i] <- c_drug_cro[i] + d1
      } else if (drug1 == "nda") {
        c_drug_nda[i] <- c_drug_nda[i] + d1
      } else {
        # e.g., ndb as first line (rare)
        c_follow_drugs[i] <- c_follow_drugs[i] + d1
      }
      
      if (true_ng[i]) {
        # Only true infections can "fail" and return
        success1 <- drug_success(drug1, cip_r[i], cro_r[i], nda_r[i])
        
        if (!success1 && sympt[i] && (rbinom(1,1,psa$p_return_sympt)==1)) {
          # ---------- INDEX: FOLLOW-UP ----------
          cost[i] <- cost[i] + psa$visit_followup
          c_follow_visit[i] <- c_follow_visit[i] + psa$visit_followup
          
          diag_fu <- if (is_woman) psa$diag_followup_women else psa$diag_followup_men
          cost[i] <- cost[i] + diag_fu
          c_follow_diag[i] <- c_follow_diag[i] + diag_fu
          
          drug2 <- select_next_drug(drug1, cip_r[i], cro_r[i], nda_r[i])
          if (!is.na(drug2)) {
            d2 <- drug_cost(drug2, psa)
            cost[i] <- cost[i] + d2
            c_follow_drugs[i] <- c_follow_drugs[i] + d2
            success2 <- drug_success(drug2, cip_r[i], cro_r[i], nda_r[i])
            index_final_treated[i] <- success2
          } else {
            index_final_treated[i] <- success1
          }
        } else {
          index_final_treated[i] <- success1
        }
      } else {
        # False positive: treated but cannot fail (no infection)
        index_final_treated[i] <- TRUE
      }
      
      # ---------- PARTNER PATHWAY ----------
      if (partner_engaged[i]) {
        partner_sympt <- rbinom(1, 1, partner_sympt_p) == 1
        
        # Same first-line choice as index for partner
        partner_drug1 <- drug1
        pd1 <- drug_cost(partner_drug1, psa)
        cost[i]            <- cost[i] + pd1
        c_partner_total[i] <- c_partner_total[i] + pd1
        
        if (true_ng[i]) {
          # Assume partner infected when index truly infected
          partner_success1 <- drug_success(partner_drug1, cip_r[i], cro_r[i], nda_r[i])
          if (!partner_success1 && partner_sympt) {
            cost[i]            <- cost[i] + psa$visit_followup
            c_partner_total[i] <- c_partner_total[i] + psa$visit_followup
            
            diag_fu_p <- if (partner_is_male) psa$diag_followup_men else psa$diag_followup_women
            cost[i]            <- cost[i] + diag_fu_p
            c_partner_total[i] <- c_partner_total[i] + diag_fu_p
            
            partner_drug2 <- select_next_drug(partner_drug1, cip_r[i], cro_r[i], nda_r[i])
            if (!is.na(partner_drug2)) {
              pd2 <- drug_cost(partner_drug2, psa)
              cost[i]            <- cost[i] + pd2
              c_partner_total[i] <- c_partner_total[i] + pd2
              partner_final_treated[i] <- drug_success(partner_drug2, cip_r[i], cro_r[i], nda_r[i])
            } else {
              partner_final_treated[i] <- partner_success1
            }
          } else {
            partner_final_treated[i] <- partner_success1
          }
        } else {
          # Index false positive: partner treated unnecessarily but no infection / complications
          partner_final_treated[i] <- TRUE
        }
      } else {
        partner_final_treated[i] <- FALSE
      }
    } # end treat_now
    
    # ---------- COMPLICATION COSTS ----------
    if (true_ng[i]) {
      # Index complications
      if (is_woman) {
        pid_prob <- if (sympt[i]) psa$p_pid_sympt else psa$p_pid_asympt
        if (rbinom(1,1,pid_prob)==1) {
          comp_cost <- psa$c_pid +
            psa$p_cpp_given_pid * psa$c_cpp +
            psa$p_ep_given_pid  * psa$c_ep  +
            psa$p_ti_given_pid  * psa$c_ti
          cost[i]           <- cost[i] + comp_cost
          c_complications[i]<- c_complications[i] + comp_cost
        }
      } else {
        if (!index_final_treated[i]) {
          comp_cost <- (psa$p_epi_untreated * psa$c_epi) +
            (psa$p_dgi_untreated * psa$c_dgi)
          cost[i]           <- cost[i] + comp_cost
          c_complications[i]<- c_complications[i] + comp_cost
        }
      }
      
      # Partner complications (only if index true infection)
      if (partner_is_male) {
        if (!partner_final_treated[i]) {
          comp_cost_p <- (psa$p_epi_untreated * psa$c_epi) +
            (psa$p_dgi_untreated * psa$c_dgi)
          cost[i]            <- cost[i] + comp_cost_p
          c_partner_total[i] <- c_partner_total[i] + comp_cost_p
          c_complications[i] <- c_complications[i] + comp_cost_p
        }
      } else {
        # female partner: PID risk by partner symptom status
        partner_sympt_here <- rbinom(1,1,partner_sympt_p)==1
        pid_prob_partner <- if (partner_sympt_here) psa$p_pid_sympt else psa$p_pid_asympt
        if (rbinom(1,1,pid_prob_partner)==1) {
          comp_cost_p <- psa$c_pid +
            psa$p_cpp_given_pid * psa$c_cpp +
            psa$p_ep_given_pid  * psa$c_ep  +
            psa$p_ti_given_pid  * psa$c_ti
          cost[i]            <- cost[i] + comp_cost_p
          c_partner_total[i] <- c_partner_total[i] + comp_cost_p
          c_complications[i] <- c_complications[i] + comp_cost_p
        }
      }
    }
  }
  
  # Per-person means among those WITH gonorrhea
  if (any(true_ng)) {
    idx <- which(true_ng)
    mean_cost_ng            <- mean(cost[idx])
    mean_init_visit_ng      <- mean(c_init_visit[idx])
    mean_init_diag_ng       <- mean(c_init_diag[idx])
    mean_drug_cip_ng        <- mean(c_drug_cip[idx])
    mean_drug_cro_ng        <- mean(c_drug_cro[idx])
    mean_drug_nda_ng        <- mean(c_drug_nda[idx])
    mean_follow_visit_ng    <- mean(c_follow_visit[idx])
    mean_follow_diag_ng     <- mean(c_follow_diag[idx])
    mean_follow_drugs_ng    <- mean(c_follow_drugs[idx])
    mean_partner_total_ng   <- mean(c_partner_total[idx])
    mean_complications_ng   <- mean(c_complications[idx])
  } else {
    mean_cost_ng          <- NA_real_
    mean_init_visit_ng    <- NA_real_
    mean_init_diag_ng     <- NA_real_
    mean_drug_cip_ng      <- NA_real_
    mean_drug_cro_ng      <- NA_real_
    mean_drug_nda_ng      <- NA_real_
    mean_follow_visit_ng  <- NA_real_
    mean_follow_diag_ng   <- NA_real_
    mean_follow_drugs_ng  <- NA_real_
    mean_partner_total_ng <- NA_real_
    mean_complications_ng <- NA_real_
  }
  
  list(
    mean_cost = mean(cost),
    mean_cost_ng          = mean_cost_ng,
    mean_init_visit_ng    = mean_init_visit_ng,
    mean_init_diag_ng     = mean_init_diag_ng,
    mean_drug_cip_ng      = mean_drug_cip_ng,
    mean_drug_cro_ng      = mean_drug_cro_ng,
    mean_drug_nda_ng      = mean_drug_nda_ng,
    mean_follow_visit_ng  = mean_follow_visit_ng,
    mean_follow_diag_ng   = mean_follow_diag_ng,
    mean_follow_drugs_ng  = mean_follow_drugs_ng,
    mean_partner_total_ng = mean_partner_total_ng,
    mean_complications_ng = mean_complications_ng,
    test_pos_share        = test_pos_share
  )
}

# ------------------------------
# Evaluate a scenario: summaries + FAST THRESHOLD
# ------------------------------
evaluate_scenario <- function(grp, cip_res, cro_res, nda_res, panel) {
  cost_with0      <- numeric(n_psa)
  cost_without    <- numeric(n_psa)
  testpos_share_w <- numeric(n_psa)
  
  # Per-NG+ total cost
  cost_with0_ng   <- numeric(n_psa)
  cost_without_ng <- numeric(n_psa)
  
  # Per-NG+ cost categories (with / without resistance testing)
  ng_init_visit_with      <- numeric(n_psa)
  ng_init_visit_without   <- numeric(n_psa)
  ng_init_diag_with       <- numeric(n_psa)
  ng_init_diag_without    <- numeric(n_psa)
  ng_cip_with             <- numeric(n_psa)
  ng_cip_without          <- numeric(n_psa)
  ng_cro_with             <- numeric(n_psa)
  ng_cro_without          <- numeric(n_psa)
  ng_nda_with             <- numeric(n_psa)
  ng_nda_without          <- numeric(n_psa)
  ng_follow_visit_with    <- numeric(n_psa)
  ng_follow_visit_without <- numeric(n_psa)
  ng_follow_diag_with     <- numeric(n_psa)
  ng_follow_diag_without  <- numeric(n_psa)
  ng_follow_drugs_with    <- numeric(n_psa)
  ng_follow_drugs_without <- numeric(n_psa)
  ng_partner_with         <- numeric(n_psa)
  ng_partner_without      <- numeric(n_psa)
  ng_comp_with            <- numeric(n_psa)
  ng_comp_without         <- numeric(n_psa)
  
  for (k in seq_len(n_psa)) {
    psa <- draw_psa()
    
    # No resistance testing
    s0 <- simulate_once(grp, cip_res, cro_res, nda_res,
                        panel = character(0), psa = psa,
                        n = n_indiv, test_cost = 0)
    cost_without[k]    <- s0$mean_cost
    cost_without_ng[k] <- s0$mean_cost_ng
    
    ng_init_visit_without[k]   <- s0$mean_init_visit_ng
    ng_init_diag_without[k]    <- s0$mean_init_diag_ng
    ng_cip_without[k]          <- s0$mean_drug_cip_ng
    ng_cro_without[k]          <- s0$mean_drug_cro_ng
    ng_nda_without[k]          <- s0$mean_drug_nda_ng
    ng_follow_visit_without[k] <- s0$mean_follow_visit_ng
    ng_follow_diag_without[k]  <- s0$mean_follow_diag_ng
    ng_follow_drugs_without[k] <- s0$mean_follow_drugs_ng
    ng_partner_without[k]      <- s0$mean_partner_total_ng
    ng_comp_without[k]         <- s0$mean_complications_ng
    
    # With resistance testing (test cost = 0) to get baseline mean and test-positive share
    s1 <- simulate_once(grp, cip_res, cro_res, nda_res,
                        panel = panel, psa = psa,
                        n = n_indiv, test_cost = 0)
    cost_with0[k]      <- s1$mean_cost
    cost_with0_ng[k]   <- s1$mean_cost_ng
    testpos_share_w[k] <- s1$test_pos_share
    
    ng_init_visit_with[k]    <- s1$mean_init_visit_ng
    ng_init_diag_with[k]     <- s1$mean_init_diag_ng
    ng_cip_with[k]           <- s1$mean_drug_cip_ng
    ng_cro_with[k]           <- s1$mean_drug_cro_ng
    ng_nda_with[k]           <- s1$mean_drug_nda_ng
    ng_follow_visit_with[k]  <- s1$mean_follow_visit_ng
    ng_follow_diag_with[k]   <- s1$mean_follow_diag_ng
    ng_follow_drugs_with[k]  <- s1$mean_follow_drugs_ng
    ng_partner_with[k]       <- s1$mean_partner_total_ng
    ng_comp_with[k]          <- s1$mean_complications_ng
  }
  
  out <- tibble(
    group = grp,
    panel = paste(panel, collapse = "+") %>% { if (. == "") "NoTest" else . },
    cip_res = cip_res,
    cro_res = cro_res,
    nda_res = nda_res,
    
    # Whole-cohort mean costs (everyone tested)
    cost_with_test_mean     = mean(cost_with0),
    cost_with_test_lower    = quantile(cost_with0, 0.025),
    cost_with_test_upper    = quantile(cost_with0, 0.975),
    cost_without_test_mean  = mean(cost_without),
    cost_without_test_lower = quantile(cost_without, 0.025),
    cost_without_test_upper = quantile(cost_without, 0.975),
    
    # Per-person cost among those WITH gonorrhea (total)
    cost_with_test_mean_ng     = mean(cost_with0_ng,   na.rm = TRUE),
    cost_with_test_lower_ng    = quantile(cost_with0_ng, 0.025, na.rm = TRUE),
    cost_with_test_upper_ng    = quantile(cost_with0_ng, 0.975, na.rm = TRUE),
    cost_without_test_mean_ng  = mean(cost_without_ng,  na.rm = TRUE),
    cost_without_test_lower_ng = quantile(cost_without_ng, 0.025, na.rm = TRUE),
    cost_without_test_upper_ng = quantile(cost_without_ng, 0.975, na.rm = TRUE),
    
    # Per-NG+ cost categories WITH resistance testing
    cost_with_test_ng_init_visit_mean      = mean(ng_init_visit_with,    na.rm = TRUE),
    cost_with_test_ng_init_visit_lower     = quantile(ng_init_visit_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_init_visit_upper     = quantile(ng_init_visit_with, 0.975, na.rm = TRUE),
    
    cost_with_test_ng_init_diag_mean       = mean(ng_init_diag_with,     na.rm = TRUE),
    cost_with_test_ng_init_diag_lower      = quantile(ng_init_diag_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_init_diag_upper      = quantile(ng_init_diag_with, 0.975, na.rm = TRUE),
    
    cost_with_test_ng_cip_mean             = mean(ng_cip_with,           na.rm = TRUE),
    cost_with_test_ng_cip_lower            = quantile(ng_cip_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_cip_upper            = quantile(ng_cip_with, 0.975, na.rm = TRUE),
    
    cost_with_test_ng_cro_mean             = mean(ng_cro_with,           na.rm = TRUE),
    cost_with_test_ng_cro_lower            = quantile(ng_cro_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_cro_upper            = quantile(ng_cro_with, 0.975, na.rm = TRUE),
    
    cost_with_test_ng_nda_mean             = mean(ng_nda_with,           na.rm = TRUE),
    cost_with_test_ng_nda_lower            = quantile(ng_nda_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_nda_upper            = quantile(ng_nda_with, 0.975, na.rm = TRUE),
    
    cost_with_test_ng_follow_visit_mean    = mean(ng_follow_visit_with,  na.rm = TRUE),
    cost_with_test_ng_follow_visit_lower   = quantile(ng_follow_visit_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_follow_visit_upper   = quantile(ng_follow_visit_with, 0.975, na.rm = TRUE),
    
    cost_with_test_ng_follow_diag_mean     = mean(ng_follow_diag_with,   na.rm = TRUE),
    cost_with_test_ng_follow_diag_lower    = quantile(ng_follow_diag_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_follow_diag_upper    = quantile(ng_follow_diag_with, 0.975, na.rm = TRUE),
    
    cost_with_test_ng_follow_drugs_mean    = mean(ng_follow_drugs_with,  na.rm = TRUE),
    cost_with_test_ng_follow_drugs_lower   = quantile(ng_follow_drugs_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_follow_drugs_upper   = quantile(ng_follow_drugs_with, 0.975, na.rm = TRUE),
    
    cost_with_test_ng_partner_total_mean   = mean(ng_partner_with,       na.rm = TRUE),
    cost_with_test_ng_partner_total_lower  = quantile(ng_partner_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_partner_total_upper  = quantile(ng_partner_with, 0.975, na.rm = TRUE),
    
    cost_with_test_ng_complications_mean   = mean(ng_comp_with,          na.rm = TRUE),
    cost_with_test_ng_complications_lower  = quantile(ng_comp_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_complications_upper  = quantile(ng_comp_with, 0.975, na.rm = TRUE),
    
    # Per-NG+ cost categories WITHOUT resistance testing
    cost_without_test_ng_init_visit_mean      = mean(ng_init_visit_without,    na.rm = TRUE),
    cost_without_test_ng_init_visit_lower     = quantile(ng_init_visit_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_init_visit_upper     = quantile(ng_init_visit_without, 0.975, na.rm = TRUE),
    
    cost_without_test_ng_init_diag_mean       = mean(ng_init_diag_without,     na.rm = TRUE),
    cost_without_test_ng_init_diag_lower      = quantile(ng_init_diag_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_init_diag_upper      = quantile(ng_init_diag_without, 0.975, na.rm = TRUE),
    
    cost_without_test_ng_cip_mean             = mean(ng_cip_without,           na.rm = TRUE),
    cost_without_test_ng_cip_lower            = quantile(ng_cip_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_cip_upper            = quantile(ng_cip_without, 0.975, na.rm = TRUE),
    
    cost_without_test_ng_cro_mean             = mean(ng_cro_without,           na.rm = TRUE),
    cost_without_test_ng_cro_lower            = quantile(ng_cro_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_cro_upper            = quantile(ng_cro_without, 0.975, na.rm = TRUE),
    
    cost_without_test_ng_nda_mean             = mean(ng_nda_without,           na.rm = TRUE),
    cost_without_test_ng_nda_lower            = quantile(ng_nda_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_nda_upper            = quantile(ng_nda_without, 0.975, na.rm = TRUE),
    
    cost_without_test_ng_follow_visit_mean    = mean(ng_follow_visit_without,  na.rm = TRUE),
    cost_without_test_ng_follow_visit_lower   = quantile(ng_follow_visit_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_follow_visit_upper   = quantile(ng_follow_visit_without, 0.975, na.rm = TRUE),
    
    cost_without_test_ng_follow_diag_mean     = mean(ng_follow_diag_without,   na.rm = TRUE),
    cost_without_test_ng_follow_diag_lower    = quantile(ng_follow_diag_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_follow_diag_upper    = quantile(ng_follow_diag_without, 0.975, na.rm = TRUE),
    
    cost_without_test_ng_follow_drugs_mean    = mean(ng_follow_drugs_without,  na.rm = TRUE),
    cost_without_test_ng_follow_drugs_lower   = quantile(ng_follow_drugs_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_follow_drugs_upper   = quantile(ng_follow_drugs_without, 0.975, na.rm = TRUE),
    
    cost_without_test_ng_partner_total_mean   = mean(ng_partner_without,       na.rm = TRUE),
    cost_without_test_ng_partner_total_lower  = quantile(ng_partner_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_partner_total_upper  = quantile(ng_partner_without, 0.975, na.rm = TRUE),
    
    cost_without_test_ng_complications_mean   = mean(ng_comp_without,          na.rm = TRUE),
    cost_without_test_ng_complications_lower  = quantile(ng_comp_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_complications_upper  = quantile(ng_comp_without, 0.975, na.rm = TRUE)
  )
  
  # FAST THRESHOLD: per resistance test among NG test-positive specimens
  # Compute per-PSA-draw threshold to get proper uncertainty intervals
  mean_share <- mean(testpos_share_w)
  
  # Per-draw threshold (same formula applied to each PSA draw)
  thresh_draws <- ifelse(
    testpos_share_w > 0,
    (cost_without - cost_with0) / testpos_share_w,
    NA_real_
  )
  
  threshold_cost  <- if (mean_share > 0) mean(thresh_draws, na.rm = TRUE) else NA_real_
  threshold_lower <- quantile(thresh_draws, 0.025, na.rm = TRUE)
  threshold_upper <- quantile(thresh_draws, 0.975, na.rm = TRUE)
  
  list(
    summary = out,
    threshold = tibble(
      group = grp,
      panel = paste(panel, collapse = "+") %>% { if (. == "") "NoTest" else . },
      cip_res = cip_res,
      cro_res = cro_res,
      nda_res = nda_res,
      threshold_price = threshold_cost,
      threshold_lower = threshold_lower,
      threshold_upper = threshold_upper
    )
  )
}

# ------------------------------
# Run across groups × resistance grids × strategies (with progress bar)
# ------------------------------
all_summaries  <- list()
all_thresholds <- list()

cat("Running PSA across groups, resistance grids, and strategies…\n")
pb <- txtProgressBar(min = 0,
                     max = length(groups)*length(cip_res_vals)*length(cro_res_vals)*length(resistance_panels),
                     style = 3)
step <- 0

for (grp in groups) {
  for (cipr in cip_res_vals) {
    for (cror in cro_res_vals) {
      for (panel_name in names(resistance_panels)) {
        panel <- resistance_panels[[panel_name]]
        res <- evaluate_scenario(grp, cipr, cror, nda_res_fixed, panel)
        all_summaries[[length(all_summaries)+1]]  <- res$summary %>% mutate(strategy = panel_name)
        all_thresholds[[length(all_thresholds)+1]]<- res$threshold %>% mutate(strategy = panel_name)
        step <- step + 1
        setTxtProgressBar(pb, step)
      }
    }
  }
}
close(pb)

summary_df   <- bind_rows(all_summaries)
threshold_df <- bind_rows(all_thresholds)

# Save raw outputs (including cost categories)
readr::write_csv(summary_df,   "gonorrhea_costs_summary.csv")
readr::write_csv(threshold_df, "gonorrhea_thresholds.csv")

# ------------------------------
# Clean, manuscript-ready summary table (subset of key fields)
# ------------------------------
nice_summary <- summary_df %>%
  mutate(
    Group = factor(group, levels = c("MSM","MSW","WSM")),
    Strategy = factor(strategy, levels = c("CIP","CIP+CRO","CIP+CRO+NDA","NoTest"))
  ) %>%
  select(
    Group, Strategy, cip_res, cro_res, nda_res,
    # Whole-cohort cost per person tested
    cost_with_test_mean,  cost_with_test_lower,  cost_with_test_upper,
    cost_without_test_mean, cost_without_test_lower, cost_without_test_upper,
    
    # Total cost per NG+ person
    cost_with_test_mean_ng,  cost_with_test_lower_ng,  cost_with_test_upper_ng,
    cost_without_test_mean_ng, cost_without_test_lower_ng, cost_without_test_upper_ng,
    
    # Cost categories per NG+ person (with resistance testing)
    cost_with_test_ng_init_visit_mean,      cost_with_test_ng_init_visit_lower,      cost_with_test_ng_init_visit_upper,
    cost_with_test_ng_init_diag_mean,       cost_with_test_ng_init_diag_lower,       cost_with_test_ng_init_diag_upper,
    cost_with_test_ng_cip_mean,             cost_with_test_ng_cip_lower,             cost_with_test_ng_cip_upper,
    cost_with_test_ng_cro_mean,             cost_with_test_ng_cro_lower,             cost_with_test_ng_cro_upper,
    cost_with_test_ng_nda_mean,             cost_with_test_ng_nda_lower,             cost_with_test_ng_nda_upper,
    cost_with_test_ng_follow_visit_mean,    cost_with_test_ng_follow_visit_lower,    cost_with_test_ng_follow_visit_upper,
    cost_with_test_ng_follow_diag_mean,     cost_with_test_ng_follow_diag_lower,     cost_with_test_ng_follow_diag_upper,
    cost_with_test_ng_follow_drugs_mean,    cost_with_test_ng_follow_drugs_lower,    cost_with_test_ng_follow_drugs_upper,
    cost_with_test_ng_partner_total_mean,   cost_with_test_ng_partner_total_lower,   cost_with_test_ng_partner_total_upper,
    cost_with_test_ng_complications_mean,   cost_with_test_ng_complications_lower,   cost_with_test_ng_complications_upper,
    
    # Cost categories per NG+ person (without resistance testing)
    cost_without_test_ng_init_visit_mean,      cost_without_test_ng_init_visit_lower,      cost_without_test_ng_init_visit_upper,
    cost_without_test_ng_init_diag_mean,       cost_without_test_ng_init_diag_lower,       cost_without_test_ng_init_diag_upper,
    cost_without_test_ng_cip_mean,             cost_without_test_ng_cip_lower,             cost_without_test_ng_cip_upper,
    cost_without_test_ng_cro_mean,             cost_without_test_ng_cro_lower,             cost_without_test_ng_cro_upper,
    cost_without_test_ng_nda_mean,             cost_without_test_ng_nda_lower,             cost_without_test_ng_nda_upper,
    cost_without_test_ng_follow_visit_mean,    cost_without_test_ng_follow_visit_lower,    cost_without_test_ng_follow_visit_upper,
    cost_without_test_ng_follow_diag_mean,     cost_without_test_ng_follow_diag_lower,     cost_without_test_ng_follow_diag_upper,
    cost_without_test_ng_follow_drugs_mean,    cost_without_test_ng_follow_drugs_lower,    cost_without_test_ng_follow_drugs_upper,
    cost_without_test_ng_partner_total_mean,   cost_without_test_ng_partner_total_lower,   cost_without_test_ng_partner_total_upper,
    cost_without_test_ng_complications_mean,   cost_without_test_ng_complications_lower,   cost_without_test_ng_complications_upper
  ) %>%
  arrange(Group, Strategy, cip_res, cro_res)

readr::write_csv(nice_summary, "gonorrhea_summary_table.csv")

# ------------------------------
# 3×3 heatmap: threshold price vs CIP × CRO, by Strategy × Group
# ------------------------------
heat_df <- threshold_df %>%
  filter(strategy %in% c("CIP","CIP+CRO","CIP+CRO+NDA")) %>%
  mutate(
    Group = factor(group, levels = c("MSM","MSW","WSM")),
    Strategy = factor(strategy, levels = c("CIP","CIP+CRO","CIP+CRO+NDA")),
    label = ifelse(is.finite(threshold_price),
                   sprintf("$%.1f", threshold_price),
                   "")
  )

p_heat <- ggplot(heat_df, aes(x = cip_res, y = cro_res, fill = threshold_price)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), color = "white", fontface = "bold", size = 3) +
  scale_x_continuous("Ciprofloxacin resistance", breaks = cip_res_vals, labels = scales::percent) +
  scale_y_continuous("Ceftriaxone resistance", breaks = cro_res_vals, labels = scales::percent) +
  labs(fill = "Threshold price ($)",
       title = "Cost-neutral threshold prices for resistance testing",
       subtitle = paste0("Panels: 3 strategies × 3 groups (NDA resistance fixed at ",
                         scales::percent(nda_res_fixed), ")")) +
  facet_grid(Strategy ~ Group) +
  theme_minimal(base_size = 12)

ggsave("threshold_heatmap.png", p_heat, width = 10, height = 8, dpi = 300)

cat("\nDone. Files written:\n",
    "- gonorrhea_costs_summary.csv\n",
    "- gonorrhea_thresholds.csv\n",
    "- gonorrhea_summary_table.csv\n",
    "- threshold_heatmap.png\n")

# ============================================================
# Supplement S3: one-way deterministic sensitivity analysis
# Reference scenario: MSM, CIP = 50%, CRO = 5%, NDA = 0%
# Outcome: cost-neutral threshold price (CIP-only and CIP+CRO+NDA panels)
#
# Each parameter is varied one at a time between its low and high value.
# The base, low, and high runs for a given parameter share a PSA seed so that
# the difference between them is driven only by the parameter being varied;
# the base-case reference line is computed once per panel (seed = 42).
# ============================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

# ---- Reference scenario parameters ----
REF_GRP     <- "MSM"
REF_CIP_RES <- 0.50
REF_CRO_RES <- 0.05
REF_NDA_RES <- 0.00
N_DSA_PSA   <- 1000
N_DSA_INDIV <- 20000

# ---- draw_psa_override: set one parameter to a fixed value ----
draw_psa_override <- function(override_name = NULL, override_val = NULL) {
  psa <- draw_psa()
  if (!is.null(override_name)) {
    if (override_name == "p_partner_treat") {
      psa$p_partner_treat <- override_val
    } else if (override_name == "cost_nda_ndb") {
      psa$cost_nda <- override_val
      psa$cost_ndb <- override_val
      if ("cost_nda_ndb" %in% names(psa)) psa$cost_nda_ndb <- override_val
    } else if (override_name == "p_return_sympt") {
      psa$p_return_sympt <- override_val
    }
    # sens/spec overridden via globals in evaluate_scenario_dsa()
  }
  psa
}

# ---- evaluate_scenario_dsa: run PSA and return threshold ----
# seed_offset: set before the k-loop so base/low/high are comparable
evaluate_scenario_dsa <- function(grp, cip_res, cro_res, nda_res, panel,
                                  override_name = NULL, override_val = NULL,
                                  n_psa_dsa = N_DSA_PSA, n_indiv_dsa = N_DSA_INDIV,
                                  row_seed = 42) {
  cost_with0      <- numeric(n_psa_dsa)
  cost_without    <- numeric(n_psa_dsa)
  testpos_share_w <- numeric(n_psa_dsa)
  
  orig_sens <- sens_res_test
  orig_spec <- spec_res_test
  if (!is.null(override_name) && override_name == "sens_res_test")
    sens_res_test <<- list(cip = override_val, cro = override_val, nda = override_val)
  if (!is.null(override_name) && override_name == "spec_res_test")
    spec_res_test <<- list(cip = override_val, cro = override_val, nda = override_val)
  
  set.seed(row_seed)   # same seed for base/low/high within a parameter row
  for (k in seq_len(n_psa_dsa)) {
    psa <- draw_psa_override(override_name, override_val)
    
    s0 <- simulate_once(grp, cip_res, cro_res, nda_res,
                        panel = character(0), psa = psa,
                        n = n_indiv_dsa, test_cost = 0)
    cost_without[k] <- s0$mean_cost
    
    s1 <- simulate_once(grp, cip_res, cro_res, nda_res,
                        panel = panel, psa = psa,
                        n = n_indiv_dsa, test_cost = 0)
    cost_with0[k]      <- s1$mean_cost
    testpos_share_w[k] <- s1$test_pos_share
  }
  
  sens_res_test <<- orig_sens
  spec_res_test <<- orig_spec
  
  mean_share <- mean(testpos_share_w)
  if (mean_share > 0) (mean(cost_without) - mean(cost_with0)) / mean_share else NA_real_
}

# ---- Panel definitions ----
panels_dsa <- list(
  "CIP"         = c("cip"),
  "CIP+CRO+NDA" = c("cip", "cro", "nda")
)

# ---- Base-case reference threshold, one per panel (seed = 42) ----
cat("Computing base case thresholds (seed=42)...\n")
fixed_base_thresh <- list()
for (pname in names(panels_dsa)) {
  fixed_base_thresh[[pname]] <- evaluate_scenario_dsa(
    REF_GRP, REF_CIP_RES, REF_CRO_RES, REF_NDA_RES,
    panels_dsa[[pname]], NULL, NULL, row_seed = 42
  )
  cat(sprintf("  %s: $%.2f\n", pname, fixed_base_thresh[[pname]]))
}

# ---- DSA parameters ----
dsa_params <- tribble(
  ~param,            ~label,                       ~label_unit, ~base, ~low,  ~high,
  "p_partner_treat", "Partner notification rate",  "%",          0.30,  0.20,  0.50,
  "cost_nda_ndb",    "Cost of NDA/NDB",            "$",         50,    25,    500,
  "p_return_sympt",  "P(return | symptoms)",       "%",          0.90,  0.70,  0.95,
  "sens_res_test",   "Resistance test sensitivity","%" ,         0.99,  0.92,  1.00,
  "spec_res_test",   "Resistance test specificity","%",          0.98,  0.96,  1.00
)

# ---- One-way sweeps: shared seed per parameter row ----
# base/low/high all call set.seed(row_seed) at the top of their PSA loop,
# so the stochastic draws are identical and the only difference between the
# three runs is the value of the parameter being varied.
dsa_results <- list()
cat("\nRunning DSA...\n")

for (i in seq_len(nrow(dsa_params))) {
  row      <- dsa_params[i, ]
  row_seed <- 100 + i   # unique per parameter, shared across base/low/high
  cat(sprintf("  [%d/%d] %s\n", i, nrow(dsa_params), row$label))
  
  for (pname in names(panels_dsa)) {
    panel <- panels_dsa[[pname]]
    
    base_thresh <- evaluate_scenario_dsa(
      REF_GRP, REF_CIP_RES, REF_CRO_RES, REF_NDA_RES,
      panel, NULL, NULL, row_seed = row_seed
    )
    low_thresh <- evaluate_scenario_dsa(
      REF_GRP, REF_CIP_RES, REF_CRO_RES, REF_NDA_RES,
      panel, row$param, row$low, row_seed = row_seed
    )
    high_thresh <- evaluate_scenario_dsa(
      REF_GRP, REF_CIP_RES, REF_CRO_RES, REF_NDA_RES,
      panel, row$param, row$high, row_seed = row_seed
    )
    
    if (!is.na(low_thresh) && !is.na(high_thresh) && low_thresh > high_thresh)
      cat(sprintf("    [!] %s | %s: low=$%.2f > high=$%.2f — check economic direction\n",
                  pname, row$param, low_thresh, high_thresh))
    
    dsa_results[[length(dsa_results) + 1]] <- tibble(
      panel       = pname,
      param       = row$param,
      label       = row$label,
      label_unit  = row$label_unit,
      base_val    = row$base,
      low_val     = row$low,
      high_val    = row$high,
      base_thresh = fixed_base_thresh[[pname]],  # global fixed base
      row_base_thresh = base_thresh,             # per-row base (for QC)
      low_thresh  = low_thresh,
      high_thresh = high_thresh
    )
  }
}

dsa_df <- bind_rows(dsa_results)
readr::write_csv(dsa_df, "gonorrhea_S3_dsa_results.csv")

# ---- Build tornado plot data ----
tornado_df <- dsa_df %>%
  mutate(
    # Format y-axis label with low and high values in parentheses
    # e.g. "Partner notification rate (20%; 50%)"
    y_label = case_when(
      label_unit == "$" ~ sprintf("%s ($%.0f; $%.0f)", label, low_val, high_val),
      label_unit == "%" ~ sprintf("%s (%.0f%%; %.0f%%)", label, low_val * 100, high_val * 100),
      TRUE              ~ sprintf("%s (%s%.2f; %s%.2f)", label, label_unit, low_val, label_unit, high_val)
    ),
    bar_left    = pmin(low_thresh, high_thresh),
    bar_right   = pmax(low_thresh, high_thresh),
    total_swing = bar_right - bar_left,
    # Which endpoint is left vs right of dashed line?
    left_is_low  = low_thresh < high_thresh   # TRUE => low value produced left bar end
  )

# ---- Tornado plot function ----
# Two-colour annotation: value left of dashed = coral, right of dashed = steelblue
# Y-axis label includes (low; high) parenthetical

COL_RIGHT <- "#4682B4"   # steelblue — value to the right (higher threshold)
COL_LEFT  <- "#B05C57"   # muted coral — value to the left (lower threshold)
COL_BAR   <- "#6FA3C7"   # bar fill

make_tornado <- function(data, panel_name, base_thresh_val) {
  
  param_order <- data %>%
    arrange(total_swing) %>%
    pull(y_label)
  
  data <- data %>%
    mutate(y_label = factor(y_label, levels = param_order))
  
  # Two annotation rows: one for bar_left, one for bar_right
  ann_left <- data %>%
    mutate(
      x     = bar_left,
      txt   = sprintf("$%.1f", bar_left),
      hjust = 1.15,
      # colour: if left end is the low-value end, use COL_LEFT; else COL_RIGHT
      col   = ifelse(left_is_low, COL_LEFT, COL_RIGHT)
    )
  ann_right <- data %>%
    mutate(
      x     = bar_right,
      txt   = sprintf("$%.1f", bar_right),
      hjust = -0.15,
      col   = ifelse(left_is_low, COL_RIGHT, COL_LEFT)
    )
  ann <- bind_rows(ann_left, ann_right)
  
  ggplot(data, aes(y = y_label)) +
    geom_segment(
      aes(x = bar_left, xend = bar_right, yend = y_label),
      linewidth = 9, color = COL_BAR, alpha = 0.75
    ) +
    geom_vline(
      xintercept = base_thresh_val,
      color = "black", linetype = "dashed", linewidth = 0.8
    ) +
    geom_text(
      data = ann,
      aes(x = x, label = txt, hjust = hjust, color = col),
      size = 3.2, fontface = "bold", inherit.aes = FALSE,
      mapping = aes(y = y_label)
    ) +
    scale_color_identity() +
    scale_x_continuous(labels = scales::dollar, expand = expansion(mult = 0.15)) +
    labs(
      title    = sprintf("One-Way Sensitivity Analysis: %s Panel", panel_name),
      subtitle = sprintf(
        "Reference scenario: MSM, CIP=50%%, CRO=5%%, NDA=0%%\nBase case threshold: $%.1f  |  Blue = higher value; Red = lower value",
        base_thresh_val
      ),
      x = "Threshold Price (USD)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y  = element_text(size = 10.5),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
}

p_cip <- make_tornado(
  tornado_df %>% filter(panel == "CIP"),
  "CIP-only",
  fixed_base_thresh[["CIP"]]
)

p_triple <- make_tornado(
  tornado_df %>% filter(panel == "CIP+CRO+NDA"),
  "CIP+CRO+NDA",
  fixed_base_thresh[["CIP+CRO+NDA"]]
)

p_combined <- p_cip / p_triple +
  plot_annotation(
    title = "Figure S3. Deterministic Sensitivity Analysis — Tornado Diagrams",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave("gonorrhea_S3_tornado_diagram.png", p_combined, width = 11, height = 10, dpi = 300)
cat("\nSaved: gonorrhea_S3_tornado_diagram.png\n")

# ---- QC: compare fixed base vs per-row base ----
cat("\n--- QC: base threshold consistency ---\n")
dsa_df %>%
  select(panel, label, base_thresh, row_base_thresh) %>%
  mutate(diff = round(row_base_thresh - base_thresh, 3)) %>%
  print(n = Inf)

# ---- Summary table ----
cat("\n--- S3 Summary Table ---\n")
dsa_df %>%
  select(panel, label, low_val, base_val, high_val, low_thresh, base_thresh, high_thresh) %>%
  mutate(across(ends_with("thresh"), ~round(.x, 1))) %>%
  print(n = Inf)