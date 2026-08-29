# =============================================================================
# Economic value of gonorrhea resistance testing in the United States
# Individual-level stochastic simulation model
#
# VERSION 2 — revisions following co-author review:
#
#   (A) PARTNER INFECTION IS NOW PROBABILISTIC.
#       Previously partners were assumed infected whenever the index case was
#       truly infected. Partners are now infected with a per-partnership
#       transmission/concordance probability (~0.455), following fitted
#       estimates for MSM in Tuite et al. 2017 and Reichert et al. 2023 (with wide range).
#
#   (B) PID IS UPDATED TO BE DURATION-DRIVEN RATHER THAN SYMPTOM-DRIVEN.
#       Li et al. 2022 (Lancet Reg Health Am) parameterise PID risk via a
#       continuous-time hazard applied to duration of infection; their
#       symptomatic / asymptomatic split is a proxy for duration, not an
#       independent risk factor. Dividing their p(PID) by matching duration
#       gives a constant hazard of ~0.147/yr across all three strata:
#           symptomatic      0.0025 / 0.017 = 0.147
#           asympt 15-24y    0.075  / 0.53  = 0.147
#           asympt 25-39y    0.091  / 0.65  = 0.147
#       We therefore model lambda directly and let duration carry the
#       treatment effect. Successful treatment does not undo PID risk already
#       accrued before presentation (correct, and unchanged from v1), but
#       failed treatment now accrues additional duration and hence additional
#       PID risk. This is the channel through which resistance-guided therapy
#       averts PID; it was absent in v1.
#
#   (C) AGE STRATA BLENDED.
#       v1 used the 25-39y asymptomatic values only. We now blend the two
#       age strata (see PID_W_15_24 below).
#
#   (D) PARTNER PATHWAY CORRECTIONS.
#       - Partner symptom status is drawn ONCE per individual (v1 drew it
#         twice, independently, for return-to-care and for PID).
#       - Partner PID is now gated on the partner actually being infected
#         (v1 applied it whenever the index was infected).
#       - Partner PID is now averted-in-part by successful partner treatment
#         (v1 applied PID regardless of whether the partner was ever
#         notified or treated, so notification averted nothing).
#
# Reproducibility: global seed set below (set.seed(123)); the sensitivity
# analysis uses per-row seeds so low/base/high differ only in the parameter
# being varied. 
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


# resistance grids (CIP x CRO) for heatmaps
cip_res_vals <- c(0.10, 0.30, 0.50, 0.70)
cro_res_vals <- c(0.00, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15)

nda_res_vals <- c(0.00, 0.05)
nda_res_fixed <- 0.05   # retained only for backward compatibility; unused below

PID_AS_EXPECTED_COST <- TRUE

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

# =============================================================================
# (C) AGE BLENDING FOR WOMEN'S PID PARAMETERS
# -----------------------------------------------------------------------------
# Li et al. 2022 report PID probability and infection duration separately for
# ages 15-24 and 25-39. Reported gonorrhoea in women skews young, so we blend
# with more weight on the younger stratum. Adjust PID_W_15_24 to re-weight.
# Note the implied hazard is invariant to this weight (both strata give 0.147),
# so the blend affects the reference duration and probability but not lambda.
# =============================================================================
PID_W_15_24 <- 0.55    # weight on ages 15-24; (1 - this) on ages 25-39

# Li et al. Table 1 — asymptomatic PID probability, by age stratum (mean, UI)
.p_asym_1524 <- c(mean = 0.075, lo = 0.031, hi = 0.13)
.p_asym_2539 <- c(mean = 0.091, lo = 0.037, hi = 0.16)
# Li et al. Table 1 — duration of asymptomatic infection, years (mean, UI)
.d_asym_1524 <- c(mean = 0.53,  lo = 0.39,  hi = 0.66)
.d_asym_2539 <- c(mean = 0.65,  lo = 0.47,  hi = 0.83)

.blend <- function(a, b, w = PID_W_15_24) w * a + (1 - w) * b

P_PID_ASYMPT_MEAN <- .blend(.p_asym_1524["mean"], .p_asym_2539["mean"])  # 0.0822
P_PID_ASYMPT_LO   <- .blend(.p_asym_1524["lo"],   .p_asym_2539["lo"])    # 0.0337
P_PID_ASYMPT_HI   <- .blend(.p_asym_1524["hi"],   .p_asym_2539["hi"])    # 0.1435

D_ASYMPT_REF      <- unname(.blend(.d_asym_1524["mean"], .d_asym_2539["mean"]))  # 0.584

# Duration of symptomatic infection in women (Li et al.; not age-stratified)
D_SYMPT_REF       <- 0.017

# Additional duration accrued when first-line treatment FAILS:
#   - symptomatic index who returns: ~1-3 weeks before re-presenting
#   - anyone who does not return: infection persists; by the memoryless
#     property the expected further duration is the asymptomatic reference
D_FAIL_RETURN     <- 0.038          # ~2 weeks
D_FAIL_PERSIST    <- D_ASYMPT_REF   # clock effectively restarts

# (D) Fraction of her pre-treatment duration a successfully-treated female
# partner still accrues. Partner notification reaches her earlier than she
# would otherwise present. This is an assumption and is varied.

PID_PARTNER_TRUNC <- 0.5

# (A) Per-partnership probability that the partner is infected given the index
# is infected. Tuite 2017 / Reichert 2023 fit ~0.44-0.47 for MSM. Applied to
# all groups in the base case; see manuscript limitations. Beta(20,24) has
# mean 0.455, 95% ~ 0.31-0.60.
PARTNER_INF_A <- 20
PARTNER_INF_B <- 24

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

  # Partner treatment probability (engagement / notification)
  p_partner_treat   <- rbeta(1, 10.47, 24.43)    # mean ~0.30

  # (A) NEW: per-partnership probability partner is infected | index infected
  p_partner_infected <- rbeta(1, PARTNER_INF_A, PARTNER_INF_B)   # mean ~0.455

  # Probability symptomatic by group (asymptomatic = 1 - symptomatic)
  p_sympt <- c(
    MSW = rbeta(1, 5.7, 1.6),
    MSM = rbeta(1, 8.0, 3.8),
    WSM = rbeta(1, 9.2, 13.6)
  )

  # ---------------------------------------------------------------------------
  # (B) PID hazard, women.
  # Draw the asymptomatic PID probability on the age-blended UI range, then
  # back out the hazard against the FIXED blended reference duration:
  #     lambda = -log(1 - p) / d_ref
  # By construction, a cured asymptomatic woman then has PID probability
  # exactly equal to the drawn p, so the cured arm reproduces Li et al.
  # Uncertainty is carried entirely by p (as in v1); duration is a structural
  # anchor and not  an independently varied quantity, preserving the
  # p-duration correlation implied by Li et al.'s Markov fit.
  # ---------------------------------------------------------------------------
  p_pid_asympt_ref <- runif(1, P_PID_ASYMPT_LO, P_PID_ASYMPT_HI)
  lambda_pid       <- -log(1 - p_pid_asympt_ref) / D_ASYMPT_REF

  # PID sequelae (unchanged)
  p_cpp_given_pid  <- runif(1, 0.23,   0.29)
  p_ep_given_pid   <- runif(1, 0.049,  0.098)
  p_ti_given_pid   <- runif(1, 0.12,   0.23)

  # Men (if untreated) — unchanged
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
  cost_cro      <- 1.11                    # ceftriaxone
  cost_cefixime <- runif(1, 14.26, 31.09)  # retained for scenario use
  cost_nda      <- 50
  cost_ndb      <- 50

  list(
    sens_ng_men = sens_ng_men,
    sens_ng_women = sens_ng_women,
    spec_ng_men = spec_ng_men,
    spec_ng_women = spec_ng_women,
    p_return_sympt = p_return_sympt,
    p_partner_treat = p_partner_treat,
    p_partner_infected = p_partner_infected,     # (A)
    p_sympt = p_sympt,

    # PID hazard machinery (B)
    lambda_pid        = lambda_pid,
    p_pid_asympt_ref  = p_pid_asympt_ref,
    d_asympt          = D_ASYMPT_REF,
    d_sympt           = D_SYMPT_REF,
    d_fail_return     = D_FAIL_RETURN,
    d_fail_persist    = D_FAIL_PERSIST,
    pid_partner_trunc = PID_PARTNER_TRUNC,       # (D)

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
    c_pid = c_pid, c_cpp = c_cpp, c_ep = c_ep, c_ti = c_ti,
    c_epi = c_epi, c_dgi = c_dgi,
    cost_cip = cost_cip, cost_cro = cost_cro, cost_cefixime = cost_cefixime,
    cost_nda = cost_nda, cost_ndb = cost_ndb
  )
}

# ---- (B) PID probability from accrued duration ----
pid_prob_from_duration <- function(d, psa) 1 - exp(-psa$lambda_pid * d)

# ---- Expected cost of a PID episode including sequelae ----
pid_episode_cost <- function(psa) {
  psa$c_pid +
    psa$p_cpp_given_pid * psa$c_cpp +
    psa$p_ep_given_pid  * psa$c_ep  +
    psa$p_ti_given_pid  * psa$c_ti
}

male_comp_cost <- function(psa) {
  (psa$p_epi_untreated * psa$c_epi) + (psa$p_dgi_untreated * psa$c_dgi)
}

# ----------------------------------
# Drug selection helpers (unchanged)
# ----------------------------------
select_first_drug <- function(panel, cip_test_pos, cro_test_pos, nda_test_pos) {
  if (length(panel) == 0) return("cro")
  if ("cip" %in% panel && !cip_test_pos) return("cip")
  if ("cro" %in% panel && !cro_test_pos) return("cro")
  if ("nda" %in% panel && !nda_test_pos) return("nda")
  if (identical(sort(panel), "cip")) return("cro")
  if (identical(sort(panel), c("cip","cro"))) return("nda")
  if (identical(sort(panel), c("cip","cro","nda"))) return("ndb")
  "cro"
}

select_next_drug <- function(prev_drug, cip_r, cro_r, nda_r) {
  if (prev_drug == "cip") {
    if (!cro_r) return("cro")
    if (!nda_r) return("nda")
    return("ndb")
  }
  if (prev_drug == "cro") {
    if (!nda_r) return("nda")
    return("ndb")
  }
  if (prev_drug == "nda") return("ndb")
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

# =============================================================================
# (E) COMMON RANDOM NUMBERS  [NEW IN V3]
# -----------------------------------------------------------------------------
# In v2 the with-test and without-test arms each drew their own cohort inside
# simulate_once. The two arms therefore compared different sets of 20,000
# people, so everything unrelated to the testing decision (who is infected,
# who is symptomatic, who gets PID) failed to cancel in the difference
#
#     (cost_without - cost_with) / test_pos_share
#
# Because test_pos_share is ~0.01 for WSM, that residual noise was multiplied
# by ~95, producing threshold intervals of +/- $200 that were essentially
# independent of the point estimate (a Monte Carlo artefact, not uncertainty).
#
# draw_cohort() now generates every individual-level random variable once per
# PSA iteration; both arms are run through the identical cohort. Two details
# matter for the coupling to be exact:
#
#   1. Quantities whose probability depends on the arm (PID, which depends on
#      accrued duration, and return-to-care) are drawn as uniforms and compared
#      against the probability, rather than as rbinom() calls. This gives a
#      monotone coupling: a woman who develops PID under the shorter (treated)
#      duration also develops it under the longer (failed) duration, so only
#      the genuinely marginal cases flip between arms.
#
#   2. All draws are vectors indexed by i, never conditional rbinom(1, ...)
#      calls inside the loop. Conditional draws consume the RNG stream at
#      arm-dependent rates and would desynchronise the two arms even if they
#      started from the same seed.
#
# Parameter uncertainty is untouched: draw_psa() is still called once per PSA
# iteration and both arms share it, exactly as before.
# =============================================================================
draw_cohort <- function(grp, cip_res, cro_res, nda_res, psa, n) {
  is_woman <- grp == "WSM"
  is_male  <- !is_woman
  ng_prev  <- ng_prev_by_group[[grp]]

  sens_ng <- if (is_male) psa$sens_ng_men else psa$sens_ng_women
  spec_ng <- if (is_male) psa$spec_ng_men else psa$spec_ng_women

  partner_group   <- if (grp == "MSM") "MSM" else if (grp == "MSW") "WSM" else "MSW"
  partner_sympt_p <- psa$p_sympt[[partner_group]]

  true_ng <- rbinom(n, 1, ng_prev) == 1

  list(
    # ---- index characteristics ----
    sympt       = rbinom(n, 1, psa$p_sympt[[grp]]) == 1,
    true_ng     = true_ng,
    test_ng_pos = rbinom(n, 1, ifelse(true_ng, sens_ng, 1 - spec_ng)) == 1,

    # ---- true resistance profile of the infecting strain ----
    cip_r = rbinom(n, 1, cip_res) == 1,
    cro_r = rbinom(n, 1, cro_res) == 1,
    nda_r = rbinom(n, 1, nda_res) == 1,

    # ---- uniforms for resistance-test readout (sens/spec vary in the DSA) ----
    u_cip_test = runif(n),
    u_cro_test = runif(n),
    u_nda_test = runif(n),

    # ---- partner characteristics ----
    partner_engaged      = rbinom(n, 1, psa$p_partner_treat) == 1,
    partner_transmit     = rbinom(n, 1, psa$p_partner_infected) == 1,
    partner_sympt_status = rbinom(n, 1, partner_sympt_p) == 1,

    # ---- uniforms for arm-dependent events ----
    u_return      = runif(n),
    u_pid_index   = runif(n),
    u_pid_partner = runif(n)
  )
}

# ------------------------------
# Microsim for one PSA draw, group, resistance, and strategy
# ------------------------------
# cohort: a list from draw_cohort(). If NULL, one is drawn internally (this
# reproduces the v2 behaviour and is only kept so the function can still be
# called standalone; the analysis code always passes a shared cohort).
simulate_once <- function(grp, cip_res, cro_res, nda_res, panel, psa,
                          n = 5000, test_cost = 0, cohort = NULL) {
  is_woman <- grp == "WSM"
  is_male  <- !is_woman
  ng_prev  <- ng_prev_by_group[[grp]]

  if (is.null(cohort)) cohort <- draw_cohort(grp, cip_res, cro_res, nda_res, psa, n)

  # ---- unpack the shared cohort (identical across arms) ----
  sympt       <- cohort$sympt
  true_ng     <- cohort$true_ng
  test_ng_pos <- cohort$test_ng_pos

  cip_r <- cohort$cip_r
  cro_r <- cohort$cro_r
  nda_r <- cohort$nda_r

  cip_test_pos <- cohort$u_cip_test < ifelse(cip_r, sens_res_test$cip, 1 - spec_res_test$cip)
  cro_test_pos <- cohort$u_cro_test < ifelse(cro_r, sens_res_test$cro, 1 - spec_res_test$cro)
  nda_test_pos <- cohort$u_nda_test < ifelse(nda_r, sens_res_test$nda, 1 - spec_res_test$nda)

  u_return      <- cohort$u_return
  u_pid_index   <- cohort$u_pid_index
  u_pid_partner <- cohort$u_pid_partner

  c_init_visit <- rep(psa$visit_init, n)
  diag_cost    <- if (is_woman) psa$diag_women else psa$diag_men
  c_init_diag  <- rep(diag_cost, n)

  c_follow_visit <- numeric(n)
  c_follow_diag  <- numeric(n)
  c_follow_drugs <- numeric(n)

  c_drug_cip <- numeric(n)
  c_drug_cro <- numeric(n)
  c_drug_nda <- numeric(n)

  c_partner_total <- numeric(n)
  c_complications <- numeric(n)

  cost <- c_init_visit + c_init_diag

  partner_group <- if (grp == "MSM") "MSM" else if (grp == "MSW") "WSM" else "MSW"
  partner_is_male <- partner_group %in% c("MSM","MSW")
  partner_sympt_p <- psa$p_sympt[[partner_group]]

  index_final_treated   <- rep(FALSE, n)
  partner_final_treated <- rep(FALSE, n)
  partner_engaged       <- cohort$partner_engaged

  # (A) partner infected only if index truly infected AND transmission
  #     occurred across the partnership. (E) both draws come from the shared
  #     cohort, so the same partnerships are infected in both arms.
  partner_infected <- true_ng & cohort$partner_transmit

  # (D) partner symptom status drawn ONCE and reused for both
  #     return-to-care and PID duration.
  partner_sympt_status <- cohort$partner_sympt_status

  cost <- cost + ifelse(test_ng_pos, test_cost, 0)
  test_pos_share <- mean(test_ng_pos)

  treat_now <- test_ng_pos

  for (i in seq_len(n)) {
    if (treat_now[i]) {
      # ---------- INDEX: FIRST-LINE ----------
      drug1 <- select_first_drug(panel, cip_test_pos[i], cro_test_pos[i], nda_test_pos[i])
      d1    <- drug_cost(drug1, psa)
      cost[i] <- cost[i] + d1

      if (drug1 == "cip") {
        c_drug_cip[i] <- c_drug_cip[i] + d1
      } else if (drug1 == "cro") {
        c_drug_cro[i] <- c_drug_cro[i] + d1
      } else if (drug1 == "nda") {
        c_drug_nda[i] <- c_drug_nda[i] + d1
      } else {
        c_follow_drugs[i] <- c_follow_drugs[i] + d1
      }

      if (true_ng[i]) {
        success1 <- drug_success(drug1, cip_r[i], cro_r[i], nda_r[i])

        # (E) uniform comparison, not rbinom: the same symptomatic person
        #     returns to care in both arms.
        if (!success1 && sympt[i] && (u_return[i] < psa$p_return_sympt)) {
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
            index_final_treated[i] <- drug_success(drug2, cip_r[i], cro_r[i], nda_r[i])
          } else {
            index_final_treated[i] <- success1
          }
          # mark that this person did return (affects PID duration below)
          attr(index_final_treated, "returned") <- NULL  # (no-op; see idx_returned)
          idx_returned_flag <- TRUE
        } else {
          index_final_treated[i] <- success1
          idx_returned_flag <- FALSE
        }
      } else {
        index_final_treated[i] <- TRUE
        idx_returned_flag <- FALSE
      }

      # ---------- PARTNER PATHWAY ----------
      if (partner_engaged[i]) {
        partner_sympt <- partner_sympt_status[i]     # (D) reuse single draw

        partner_drug1 <- drug1
        pd1 <- drug_cost(partner_drug1, psa)
        cost[i]            <- cost[i] + pd1
        c_partner_total[i] <- c_partner_total[i] + pd1

        # (A) partner can only fail treatment if actually infected
        if (partner_infected[i]) {
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
          # Partner uninfected (index false positive, or no transmission):
          # treated unnecessarily, drug cost already counted, no sequelae.
          partner_final_treated[i] <- TRUE
        }
      } else {
        partner_final_treated[i] <- FALSE
      }
    } else {
      idx_returned_flag <- FALSE
    } # end treat_now

    # ---------- COMPLICATION COSTS ----------
    if (true_ng[i]) {

      # ===== INDEX =====
      if (is_woman) {
        # (B) Duration accrued before presentation. Symptom status is the
        #     proxy for how long she was infected before being detected.
        d_total <- if (sympt[i]) psa$d_sympt else psa$d_asympt

        # (B) additional duration if treatment ultimately failed. This is the
        #     channel through which resistance-guided therapy averts PID.
        if (!index_final_treated[i]) {
          d_total <- d_total + if (idx_returned_flag) psa$d_fail_return else psa$d_fail_persist
        }

        # (G) Cost PID as an expectation, matching how male complications are
        #     already handled. The monotone uniform coupling is kept as the
        #     fallback so PID_AS_EXPECTED_COST = FALSE reproduces v2.
        pid_p <- pid_prob_from_duration(d_total, psa)
        if (PID_AS_EXPECTED_COST) {
          comp_cost <- pid_p * pid_episode_cost(psa)
          cost[i]            <- cost[i] + comp_cost
          c_complications[i] <- c_complications[i] + comp_cost
        } else if (u_pid_index[i] < pid_p) {
          comp_cost <- pid_episode_cost(psa)
          cost[i]            <- cost[i] + comp_cost
          c_complications[i] <- c_complications[i] + comp_cost
        }
      } else {
        if (!index_final_treated[i]) {
          comp_cost <- male_comp_cost(psa)
          cost[i]            <- cost[i] + comp_cost
          c_complications[i] <- c_complications[i] + comp_cost
        }
      }

      # ===== PARTNER =====
      # (A)(D) Only if the partner is actually infected.
      if (partner_infected[i]) {
        if (partner_is_male) {
          if (!partner_final_treated[i]) {
            comp_cost_p <- male_comp_cost(psa)
            cost[i]            <- cost[i] + comp_cost_p
            c_partner_total[i] <- c_partner_total[i] + comp_cost_p
            c_complications[i] <- c_complications[i] + comp_cost_p
          }
        } else {
          # Female partner: duration-based PID, using the single symptom draw.
          d_p <- if (partner_sympt_status[i]) psa$d_sympt else psa$d_asympt

          if (partner_engaged[i] && partner_final_treated[i]) {
            # (D) Notification + successful treatment reaches her earlier than
            #     she would otherwise have presented -> truncated duration.
            d_p <- d_p * psa$pid_partner_trunc
          } else {
            # Never notified, or treatment failed -> infection persists.
            d_p <- d_p + psa$d_fail_persist
          }

          # (G) same expected-value treatment for the female partner.
          pid_p_p <- pid_prob_from_duration(d_p, psa)
          if (PID_AS_EXPECTED_COST) {
            comp_cost_p <- pid_p_p * pid_episode_cost(psa)
            cost[i]            <- cost[i] + comp_cost_p
            c_partner_total[i] <- c_partner_total[i] + comp_cost_p
            c_complications[i] <- c_complications[i] + comp_cost_p
          } else if (u_pid_partner[i] < pid_p_p) {
            comp_cost_p <- pid_episode_cost(psa)
            cost[i]            <- cost[i] + comp_cost_p
            c_partner_total[i] <- c_partner_total[i] + comp_cost_p
            c_complications[i] <- c_complications[i] + comp_cost_p
          }
        }
      }
    }
  }

  # Per-person means among those WITH gonorrhea
  if (any(true_ng)) {
    idx <- which(true_ng)
    mean_cost_ng          <- mean(cost[idx])
    mean_init_visit_ng    <- mean(c_init_visit[idx])
    mean_init_diag_ng     <- mean(c_init_diag[idx])
    mean_drug_cip_ng      <- mean(c_drug_cip[idx])
    mean_drug_cro_ng      <- mean(c_drug_cro[idx])
    mean_drug_nda_ng      <- mean(c_drug_nda[idx])
    mean_follow_visit_ng  <- mean(c_follow_visit[idx])
    mean_follow_diag_ng   <- mean(c_follow_diag[idx])
    mean_follow_drugs_ng  <- mean(c_follow_drugs[idx])
    mean_partner_total_ng <- mean(c_partner_total[idx])
    mean_complications_ng <- mean(c_complications[idx])
  } else {
    mean_cost_ng <- mean_init_visit_ng <- mean_init_diag_ng <-
      mean_drug_cip_ng <- mean_drug_cro_ng <- mean_drug_nda_ng <-
      mean_follow_visit_ng <- mean_follow_diag_ng <- mean_follow_drugs_ng <-
      mean_partner_total_ng <- mean_complications_ng <- NA_real_
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
# Evaluate a scenario: summaries + threshold
# ------------------------------
evaluate_scenario <- function(grp, cip_res, cro_res, nda_res, panel) {
  cost_with0      <- numeric(n_psa)
  cost_without    <- numeric(n_psa)
  testpos_share_w <- numeric(n_psa)

  cost_with0_ng   <- numeric(n_psa)
  cost_without_ng <- numeric(n_psa)

  ng_init_visit_with      <- numeric(n_psa); ng_init_visit_without   <- numeric(n_psa)
  ng_init_diag_with       <- numeric(n_psa); ng_init_diag_without    <- numeric(n_psa)
  ng_cip_with             <- numeric(n_psa); ng_cip_without          <- numeric(n_psa)
  ng_cro_with             <- numeric(n_psa); ng_cro_without          <- numeric(n_psa)
  ng_nda_with             <- numeric(n_psa); ng_nda_without          <- numeric(n_psa)
  ng_follow_visit_with    <- numeric(n_psa); ng_follow_visit_without <- numeric(n_psa)
  ng_follow_diag_with     <- numeric(n_psa); ng_follow_diag_without  <- numeric(n_psa)
  ng_follow_drugs_with    <- numeric(n_psa); ng_follow_drugs_without <- numeric(n_psa)
  ng_partner_with         <- numeric(n_psa); ng_partner_without      <- numeric(n_psa)
  ng_comp_with            <- numeric(n_psa); ng_comp_without         <- numeric(n_psa)

  for (k in seq_len(n_psa)) {
    psa <- draw_psa()

    # (E) ONE cohort per PSA iteration, shared by both arms.
    cohort <- draw_cohort(grp, cip_res, cro_res, nda_res, psa, n_indiv)

    s0 <- simulate_once(grp, cip_res, cro_res, nda_res,
                        panel = character(0), psa = psa,
                        n = n_indiv, test_cost = 0, cohort = cohort)
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

    s1 <- simulate_once(grp, cip_res, cro_res, nda_res,
                        panel = panel, psa = psa,
                        n = n_indiv, test_cost = 0, cohort = cohort)
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
    cip_res = cip_res, cro_res = cro_res, nda_res = nda_res,

    cost_with_test_mean     = mean(cost_with0),
    cost_with_test_lower    = quantile(cost_with0, 0.025),
    cost_with_test_upper    = quantile(cost_with0, 0.975),
    cost_without_test_mean  = mean(cost_without),
    cost_without_test_lower = quantile(cost_without, 0.025),
    cost_without_test_upper = quantile(cost_without, 0.975),

    cost_with_test_mean_ng     = mean(cost_with0_ng,   na.rm = TRUE),
    cost_with_test_lower_ng    = quantile(cost_with0_ng, 0.025, na.rm = TRUE),
    cost_with_test_upper_ng    = quantile(cost_with0_ng, 0.975, na.rm = TRUE),
    cost_without_test_mean_ng  = mean(cost_without_ng,  na.rm = TRUE),
    cost_without_test_lower_ng = quantile(cost_without_ng, 0.025, na.rm = TRUE),
    cost_without_test_upper_ng = quantile(cost_without_ng, 0.975, na.rm = TRUE),

    cost_with_test_ng_init_visit_mean    = mean(ng_init_visit_with,    na.rm = TRUE),
    cost_with_test_ng_init_visit_lower   = quantile(ng_init_visit_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_init_visit_upper   = quantile(ng_init_visit_with, 0.975, na.rm = TRUE),
    cost_with_test_ng_init_diag_mean     = mean(ng_init_diag_with,     na.rm = TRUE),
    cost_with_test_ng_init_diag_lower    = quantile(ng_init_diag_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_init_diag_upper    = quantile(ng_init_diag_with, 0.975, na.rm = TRUE),
    cost_with_test_ng_cip_mean           = mean(ng_cip_with,           na.rm = TRUE),
    cost_with_test_ng_cip_lower          = quantile(ng_cip_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_cip_upper          = quantile(ng_cip_with, 0.975, na.rm = TRUE),
    cost_with_test_ng_cro_mean           = mean(ng_cro_with,           na.rm = TRUE),
    cost_with_test_ng_cro_lower          = quantile(ng_cro_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_cro_upper          = quantile(ng_cro_with, 0.975, na.rm = TRUE),
    cost_with_test_ng_nda_mean           = mean(ng_nda_with,           na.rm = TRUE),
    cost_with_test_ng_nda_lower          = quantile(ng_nda_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_nda_upper          = quantile(ng_nda_with, 0.975, na.rm = TRUE),
    cost_with_test_ng_follow_visit_mean  = mean(ng_follow_visit_with,  na.rm = TRUE),
    cost_with_test_ng_follow_visit_lower = quantile(ng_follow_visit_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_follow_visit_upper = quantile(ng_follow_visit_with, 0.975, na.rm = TRUE),
    cost_with_test_ng_follow_diag_mean   = mean(ng_follow_diag_with,   na.rm = TRUE),
    cost_with_test_ng_follow_diag_lower  = quantile(ng_follow_diag_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_follow_diag_upper  = quantile(ng_follow_diag_with, 0.975, na.rm = TRUE),
    cost_with_test_ng_follow_drugs_mean  = mean(ng_follow_drugs_with,  na.rm = TRUE),
    cost_with_test_ng_follow_drugs_lower = quantile(ng_follow_drugs_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_follow_drugs_upper = quantile(ng_follow_drugs_with, 0.975, na.rm = TRUE),
    cost_with_test_ng_partner_total_mean  = mean(ng_partner_with,      na.rm = TRUE),
    cost_with_test_ng_partner_total_lower = quantile(ng_partner_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_partner_total_upper = quantile(ng_partner_with, 0.975, na.rm = TRUE),
    cost_with_test_ng_complications_mean  = mean(ng_comp_with,         na.rm = TRUE),
    cost_with_test_ng_complications_lower = quantile(ng_comp_with, 0.025, na.rm = TRUE),
    cost_with_test_ng_complications_upper = quantile(ng_comp_with, 0.975, na.rm = TRUE),

    cost_without_test_ng_init_visit_mean    = mean(ng_init_visit_without,    na.rm = TRUE),
    cost_without_test_ng_init_visit_lower   = quantile(ng_init_visit_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_init_visit_upper   = quantile(ng_init_visit_without, 0.975, na.rm = TRUE),
    cost_without_test_ng_init_diag_mean     = mean(ng_init_diag_without,     na.rm = TRUE),
    cost_without_test_ng_init_diag_lower    = quantile(ng_init_diag_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_init_diag_upper    = quantile(ng_init_diag_without, 0.975, na.rm = TRUE),
    cost_without_test_ng_cip_mean           = mean(ng_cip_without,           na.rm = TRUE),
    cost_without_test_ng_cip_lower          = quantile(ng_cip_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_cip_upper          = quantile(ng_cip_without, 0.975, na.rm = TRUE),
    cost_without_test_ng_cro_mean           = mean(ng_cro_without,           na.rm = TRUE),
    cost_without_test_ng_cro_lower          = quantile(ng_cro_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_cro_upper          = quantile(ng_cro_without, 0.975, na.rm = TRUE),
    cost_without_test_ng_nda_mean           = mean(ng_nda_without,           na.rm = TRUE),
    cost_without_test_ng_nda_lower          = quantile(ng_nda_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_nda_upper          = quantile(ng_nda_without, 0.975, na.rm = TRUE),
    cost_without_test_ng_follow_visit_mean  = mean(ng_follow_visit_without,  na.rm = TRUE),
    cost_without_test_ng_follow_visit_lower = quantile(ng_follow_visit_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_follow_visit_upper = quantile(ng_follow_visit_without, 0.975, na.rm = TRUE),
    cost_without_test_ng_follow_diag_mean   = mean(ng_follow_diag_without,   na.rm = TRUE),
    cost_without_test_ng_follow_diag_lower  = quantile(ng_follow_diag_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_follow_diag_upper  = quantile(ng_follow_diag_without, 0.975, na.rm = TRUE),
    cost_without_test_ng_follow_drugs_mean  = mean(ng_follow_drugs_without,  na.rm = TRUE),
    cost_without_test_ng_follow_drugs_lower = quantile(ng_follow_drugs_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_follow_drugs_upper = quantile(ng_follow_drugs_without, 0.975, na.rm = TRUE),
    cost_without_test_ng_partner_total_mean  = mean(ng_partner_without,      na.rm = TRUE),
    cost_without_test_ng_partner_total_lower = quantile(ng_partner_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_partner_total_upper = quantile(ng_partner_without, 0.975, na.rm = TRUE),
    cost_without_test_ng_complications_mean  = mean(ng_comp_without,         na.rm = TRUE),
    cost_without_test_ng_complications_lower = quantile(ng_comp_without, 0.025, na.rm = TRUE),
    cost_without_test_ng_complications_upper = quantile(ng_comp_without, 0.975, na.rm = TRUE)
  )

  mean_share <- mean(testpos_share_w)
  thresh_draws <- ifelse(testpos_share_w > 0,
                         (cost_without - cost_with0) / testpos_share_w,
                         NA_real_)

  threshold_cost  <- if (mean_share > 0) mean(thresh_draws, na.rm = TRUE) else NA_real_
  threshold_lower <- quantile(thresh_draws, 0.025, na.rm = TRUE)
  threshold_upper <- quantile(thresh_draws, 0.975, na.rm = TRUE)

  list(
    summary = out,
    threshold = tibble(
      group = grp,
      panel = paste(panel, collapse = "+") %>% { if (. == "") "NoTest" else . },
      cip_res = cip_res, cro_res = cro_res, nda_res = nda_res,
      threshold_price = threshold_cost,
      threshold_lower = threshold_lower,
      threshold_upper = threshold_upper
    ),
    # (E) the raw per-draw thresholds. Population-weighted intervals must be
    # built by weighting these draws and THEN taking quantiles; weighting the
    # already-summarised bounds sums three independent noise terms linearly
    # and overstates the interval.
    thresh_draws = thresh_draws
  )
}

# =============================================================================
# QC: verify the PID reparameterisation reproduces Li et al. in the cured arm
# =============================================================================
cat("\n=== QC: PID hazard reparameterisation ===\n")
cat(sprintf("Age blend: %.0f%% ages 15-24 / %.0f%% ages 25-39\n",
            PID_W_15_24*100, (1-PID_W_15_24)*100))
cat(sprintf("Blended asymptomatic p(PID): %.4f  [%.4f - %.4f]\n",
            P_PID_ASYMPT_MEAN, P_PID_ASYMPT_LO, P_PID_ASYMPT_HI))
cat(sprintf("Blended reference duration : %.3f years\n", D_ASYMPT_REF))
.lam_check <- -log(1 - P_PID_ASYMPT_MEAN) / D_ASYMPT_REF
cat(sprintf("Implied hazard lambda      : %.4f /yr  (Li et al. implies ~0.147)\n", .lam_check))
cat(sprintf("  -> cured asymptomatic PID: %.4f  (target %.4f)\n",
            1 - exp(-.lam_check * D_ASYMPT_REF), P_PID_ASYMPT_MEAN))
cat(sprintf("  -> cured symptomatic  PID: %.5f (Li et al. reports 0.0025)\n",
            1 - exp(-.lam_check * D_SYMPT_REF)))
cat(sprintf("  -> FAILED, no return     : %.4f  (was %.4f in v1)\n",
            1 - exp(-.lam_check * (D_ASYMPT_REF + D_FAIL_PERSIST)), P_PID_ASYMPT_MEAN))
cat(sprintf("Partner infection prob     : %.3f  (Beta(%d,%d))\n\n",
            PARTNER_INF_A/(PARTNER_INF_A+PARTNER_INF_B), PARTNER_INF_A, PARTNER_INF_B))

# ------------------------------
# Run across groups x resistance grids x strategies
# ------------------------------
all_summaries  <- list()
all_thresholds <- list()
all_draws      <- list()   # (E) per-draw thresholds, keyed by scenario

cat("Running PSA across groups, resistance grids, and strategies...\n")
pb <- txtProgressBar(min = 0,
                     max = length(groups)*length(cip_res_vals)*length(cro_res_vals)*
                           length(nda_res_vals)*length(resistance_panels),
                     style = 3)
step <- 0

for (grp in groups) {
  for (cipr in cip_res_vals) {
    for (cror in cro_res_vals) {
      for (ndar in nda_res_vals) {
        for (panel_name in names(resistance_panels)) {
          panel <- resistance_panels[[panel_name]]
          res <- evaluate_scenario(grp, cipr, cror, ndar, panel)
          all_summaries[[length(all_summaries)+1]]   <- res$summary   %>% mutate(strategy = panel_name)
          all_thresholds[[length(all_thresholds)+1]] <- res$threshold %>% mutate(strategy = panel_name)
          all_draws[[paste(grp, cipr, cror, ndar, panel_name, sep = "|")]] <- res$thresh_draws
          step <- step + 1
          setTxtProgressBar(pb, step)
        }
      }
    }
  }
}
close(pb)

summary_df   <- bind_rows(all_summaries)
threshold_df <- bind_rows(all_thresholds)

readr::write_csv(summary_df,   "gonorrhea_costs_summary_v3.csv")
readr::write_csv(threshold_df, "gonorrhea_thresholds_v3.csv")

# =============================================================================
# (E) POPULATION-WEIGHTED THRESHOLDS
# -----------------------------------------------------------------------------
# Weight the per-draw thresholds, then take quantiles of the weighted
# distribution. Because the three groups' draws are independent, their noise
# adds in quadrature here rather than linearly, which is both correct and
# narrower than weighting the summarised bounds.
# =============================================================================
POP_WEIGHTS <- c(MSM = 0.53, MSW = 0.21, WSM = 0.26)

weighted_thresholds <- list()

for (cipr in cip_res_vals) {
  for (cror in cro_res_vals) {
    for (ndar in nda_res_vals) {
      for (panel_name in names(resistance_panels)) {

        keys <- sprintf("%s|%s|%s|%s|%s", names(POP_WEIGHTS), cipr, cror, ndar, panel_name)
        if (!all(keys %in% names(all_draws))) next

        wd <- Reduce(`+`, Map(function(g, k) POP_WEIGHTS[[g]] * all_draws[[k]],
                              names(POP_WEIGHTS), keys))

        weighted_thresholds[[length(weighted_thresholds) + 1]] <- tibble(
          cip_res = cipr, cro_res = cror, nda_res = ndar,
          strategy = panel_name,
          threshold_price = mean(wd, na.rm = TRUE),
          threshold_lower = unname(quantile(wd, 0.025, na.rm = TRUE)),
          threshold_upper = unname(quantile(wd, 0.975, na.rm = TRUE))
        )
      }
    }
  }
}

weighted_df <- bind_rows(weighted_thresholds)
readr::write_csv(weighted_df, "gonorrhea_thresholds_weighted_v3.csv")

cat("\n=== Population-weighted thresholds (MSM 53% / MSW 21% / Women 26%) ===\n")
weighted_df %>%
  filter(cip_res == 0.50, cro_res == 0.05) %>%
  mutate(across(starts_with("threshold"), ~round(.x, 1))) %>%
  print()

# ---- Diagnostic: is the interval still dominated by Monte Carlo noise? ----
# Under parameter uncertainty the half-width should scale with the point
# estimate. A flat half-width across scenarios means residual simulation
# noise; if that shows up here, raise n_indiv.
cat("\n=== Noise diagnostic: correlation of half-width with point estimate ===\n")
cat("   (near +1 = parameter-driven, good; near 0 = still noise-dominated)\n")
threshold_df %>%
  mutate(halfwidth = (threshold_upper - threshold_lower) / 2) %>%
  group_by(group, strategy) %>%
  summarise(cor_point_width = round(cor(threshold_price, halfwidth), 2),
            median_halfwidth = round(median(halfwidth), 1),
            .groups = "drop") %>%
  print(n = Inf)

# ------------------------------
# Manuscript-ready summary table
# ------------------------------
nice_summary <- summary_df %>%
  mutate(
    Group = factor(group, levels = c("MSM","MSW","WSM")),
    Strategy = factor(strategy, levels = c("CIP","CIP+CRO","CIP+CRO+NDA","NoTest"))
  ) %>%
  arrange(Group, Strategy, cip_res, cro_res)

readr::write_csv(nice_summary, "gonorrhea_summary_table_v3.csv")

# ------------------------------
# Heatmap
# ------------------------------
heat_df <- threshold_df %>%
  filter(strategy %in% c("CIP","CIP+CRO","CIP+CRO+NDA")) %>%
  mutate(
    Group = factor(group, levels = c("MSM","MSW","WSM")),
    Strategy = factor(strategy, levels = c("CIP","CIP+CRO","CIP+CRO+NDA")),
    NDA = factor(ifelse(nda_res == 0, "0% NDA resistance", "5% NDA resistance"),
                 levels = c("0% NDA resistance", "5% NDA resistance")),
    label = ifelse(is.finite(threshold_price), sprintf("$%.1f", threshold_price), "")
  )

# One panel file per NDA scenario (both are now produced in a single run).
for (nda_lvl in levels(heat_df$NDA)) {
  p_heat <- ggplot(filter(heat_df, NDA == nda_lvl),
                   aes(x = cip_res, y = cro_res, fill = threshold_price)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), color = "white", fontface = "bold", size = 3) +
    scale_x_continuous("Ciprofloxacin resistance", breaks = cip_res_vals, labels = scales::percent) +
    scale_y_continuous("Ceftriaxone resistance", breaks = cro_res_vals, labels = scales::percent) +
    labs(fill = "Threshold price ($)",
         title = "Cost-neutral threshold prices for resistance testing",
         subtitle = paste0(nda_lvl,
                           "; duration-based PID, probabilistic partner infection")) +
    facet_grid(Strategy ~ Group) +
    theme_minimal(base_size = 12)

  fname <- sprintf("threshold_heatmap_v3_NDA%s.png",
                   ifelse(grepl("^0", nda_lvl), "0", "5"))
  ggsave(fname, p_heat, width = 10, height = 8, dpi = 300)
  cat(sprintf("Saved: %s\n", fname))
}

cat("\nDone. Files written:\n",
    "- gonorrhea_costs_summary_v2.csv\n",
    "- gonorrhea_thresholds_v2.csv\n",
    "- gonorrhea_summary_table_v2.csv\n",
    "- threshold_heatmap_v2.png\n")

# =============================================================================
# Supplement S3: one-way deterministic sensitivity analysis
# Reference scenario: MSM, CIP = 50%, CRO = 5%, NDA = 0%
# =============================================================================

library(ggplot2)
library(patchwork)

REF_GRP     <- "MSM"
REF_CIP_RES <- 0.50
REF_CRO_RES <- 0.05
REF_NDA_RES <- 0.00
N_DSA_PSA   <- 1000
N_DSA_INDIV <- 20000

draw_psa_override <- function(override_name = NULL, override_val = NULL) {
  psa <- draw_psa()
  if (!is.null(override_name)) {
    if (override_name == "p_partner_treat") {
      psa$p_partner_treat <- override_val
    } else if (override_name == "cost_nda_ndb") {
      psa$cost_nda <- override_val
      psa$cost_ndb <- override_val
    } else if (override_name == "p_return_sympt") {
      psa$p_return_sympt <- override_val
    } else if (override_name == "p_partner_infected") {
      # (A) NEW
      psa$p_partner_infected <- override_val
    } else if (override_name == "pid_partner_trunc") {
      # (D) NEW
      psa$pid_partner_trunc <- override_val
    } else if (override_name == "p_pid_asympt") {
      # (B) NEW: override the PID probability and recompute the hazard
      psa$p_pid_asympt_ref <- override_val
      psa$lambda_pid <- -log(1 - override_val) / psa$d_asympt
    }
  }
  psa
}

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

  set.seed(row_seed)
  for (k in seq_len(n_psa_dsa)) {
    psa <- draw_psa_override(override_name, override_val)

    # (E) shared cohort, as in evaluate_scenario()
    cohort <- draw_cohort(grp, cip_res, cro_res, nda_res, psa, n_indiv_dsa)

    s0 <- simulate_once(grp, cip_res, cro_res, nda_res,
                        panel = character(0), psa = psa,
                        n = n_indiv_dsa, test_cost = 0, cohort = cohort)
    cost_without[k] <- s0$mean_cost

    s1 <- simulate_once(grp, cip_res, cro_res, nda_res,
                        panel = panel, psa = psa,
                        n = n_indiv_dsa, test_cost = 0, cohort = cohort)
    cost_with0[k]      <- s1$mean_cost
    testpos_share_w[k] <- s1$test_pos_share
  }

  sens_res_test <<- orig_sens
  spec_res_test <<- orig_spec

  # (H) v2 returned a RATIO OF MEANS here while evaluate_scenario() returned a
  # MEAN OF RATIOS. The two are not equal, so the DSA base case never quite
  # matched the corresponding PSA cell and the tornado midline sat slightly off
  # the reported threshold. Matched to evaluate_scenario() below.
  thresh_draws <- ifelse(testpos_share_w > 0,
                         (cost_without - cost_with0) / testpos_share_w,
                         NA_real_)
  if (all(is.na(thresh_draws))) NA_real_ else mean(thresh_draws, na.rm = TRUE)
}

panels_dsa <- list(
  "CIP"         = c("cip"),
  "CIP+CRO+NDA" = c("cip", "cro", "nda")
)

cat("\nComputing base case thresholds (seed=42)...\n")
fixed_base_thresh <- list()
for (pname in names(panels_dsa)) {
  fixed_base_thresh[[pname]] <- evaluate_scenario_dsa(
    REF_GRP, REF_CIP_RES, REF_CRO_RES, REF_NDA_RES,
    panels_dsa[[pname]], NULL, NULL, row_seed = 42
  )
  cat(sprintf("  %s: $%.2f\n", pname, fixed_base_thresh[[pname]]))
}

# ---- DSA parameters (three new rows added) ----
dsa_params <- tribble(
  ~param,               ~label,                          ~label_unit, ~base,  ~low,   ~high,
  "p_partner_treat",    "Partner notification rate",     "%",          0.30,   0.20,   0.50,
  "p_partner_infected", "P(partner infected | index)",   "%",          0.455,  0.30,   0.65,
  "pid_partner_trunc",  "Partner PID duration truncation","%",         0.50,   0.25,   0.75,
  "p_pid_asympt",       "P(PID | asymptomatic, cured)",  "%",          0.0822, 0.0337, 0.1435,
  "cost_nda_ndb",       "Cost of NDA/NDB",               "$",         50,     25,     500,
  "p_return_sympt",     "P(return | symptoms)",          "%",          0.90,   0.70,   0.95,
  "sens_res_test",      "Resistance test sensitivity",   "%",          0.99,   0.92,   1.00,
  "spec_res_test",      "Resistance test specificity",   "%",          0.98,   0.96,   1.00
)

dsa_results <- list()
cat("\nRunning DSA...\n")

for (i in seq_len(nrow(dsa_params))) {
  row      <- dsa_params[i, ]
  row_seed <- 100 + i
  cat(sprintf("  [%d/%d] %s\n", i, nrow(dsa_params), row$label))

  for (pname in names(panels_dsa)) {
    panel <- panels_dsa[[pname]]

    base_thresh <- evaluate_scenario_dsa(
      REF_GRP, REF_CIP_RES, REF_CRO_RES, REF_NDA_RES,
      panel, NULL, NULL, row_seed = row_seed)
    low_thresh <- evaluate_scenario_dsa(
      REF_GRP, REF_CIP_RES, REF_CRO_RES, REF_NDA_RES,
      panel, row$param, row$low, row_seed = row_seed)
    high_thresh <- evaluate_scenario_dsa(
      REF_GRP, REF_CIP_RES, REF_CRO_RES, REF_NDA_RES,
      panel, row$param, row$high, row_seed = row_seed)

    if (!is.na(low_thresh) && !is.na(high_thresh) && low_thresh > high_thresh)
      cat(sprintf("    [!] %s | %s: low=$%.2f > high=$%.2f - check direction\n",
                  pname, row$param, low_thresh, high_thresh))

    dsa_results[[length(dsa_results) + 1]] <- tibble(
      panel = pname, param = row$param, label = row$label,
      label_unit = row$label_unit,
      base_val = row$base, low_val = row$low, high_val = row$high,
      base_thresh = fixed_base_thresh[[pname]],
      row_base_thresh = base_thresh,
      low_thresh = low_thresh, high_thresh = high_thresh
    )
  }
}

dsa_df <- bind_rows(dsa_results)
readr::write_csv(dsa_df, "gonorrhea_S3_dsa_results_v3.csv")

# ---- Tornado plots ----
tornado_df <- dsa_df %>%
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

p_cip <- make_tornado(tornado_df %>% filter(panel == "CIP"),
                      "CIP-only", fixed_base_thresh[["CIP"]])
p_triple <- make_tornado(tornado_df %>% filter(panel == "CIP+CRO+NDA"),
                         "CIP+CRO+NDA", fixed_base_thresh[["CIP+CRO+NDA"]])

p_combined <- p_cip / p_triple +
  plot_annotation(
    title = "Figure S3. Deterministic Sensitivity Analysis - Tornado Diagrams",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave("gonorrhea_S3_tornado_diagram_v3.png", p_combined, width = 11, height = 11, dpi = 300)
cat("\nSaved: gonorrhea_S3_tornado_diagram_v2.png\n")

cat("\n--- QC: base threshold consistency ---\n")
dsa_df %>%
  select(panel, label, base_thresh, row_base_thresh) %>%
  mutate(diff = round(row_base_thresh - base_thresh, 3)) %>%
  print(n = Inf)

cat("\n--- S3 Summary Table ---\n")
dsa_df %>%
  select(panel, label, low_val, base_val, high_val, low_thresh, base_thresh, high_thresh) %>%
  mutate(across(ends_with("thresh"), ~round(.x, 1))) %>%
  print(n = Inf)
