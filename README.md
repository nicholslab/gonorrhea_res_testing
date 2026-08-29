# gonorrhea_res_testing
Economic value of gonorrhea resistance testing in the United States

Code accompanying the manuscript *"Economic value of gonorrhea resistance testing
in the United States: a multi-perspective analysis of when testing pays for
itself"* ([Nichols et al.).] (https://www.medrxiv.org/content/10.64898/2026.04.07.26350302v1) 

The repository contains an individual-level stochastic simulation model of
gonorrhea diagnosis and treatment that estimates, for a range of resistance
scenarios and testing strategies, the per-episode cost of care and the
cost-neutral price at which resistance-guided therapy (RGT) breaks even against
standard empiric care.

> **This release (v2.0.0) accompanies the revised manuscript submitted in response to peer review.** The model has changed substantively since v1.0.0 and the numerical results differ from those in the original submission. See [Version history](#version-history) below. The code and results underlying the original submission remain permanently available at the `v1.0.0` tag and its Zenodo DOI.

## What changed in v2

Five changes to the model, made in response to reviewer and co-author comments. Each is documented in the header of the model script.

**A. Partner infection is now probabilistic.** v1 assumed a partner was infected whenever the index case was truly infected. Partners are now infected with a per-partnership transmission/concordance probability drawn from Beta(20, 24) (mean 0.455, 95% ~0.31–0.60), following fitted estimates for MSM in Tuite et al. 2017 and Reichert et al. 2023.

**B. PID risk is duration-driven rather than symptom-driven.** Li et al. 2022 parameterise PID risk via a continuous-time hazard applied to duration of infection; their symptomatic/asymptomatic split proxies duration rather than acting as an independent risk factor. Dividing their p(PID) by matching duration gives a constant hazard of ~0.147/yr across all three strata, so the model now applies that hazard directly and lets duration carry the treatment effect. Failed treatment accrues additional duration and therefore additional PID risk — the channel through which RGT averts PID, which was absent in v1. The script prints a QC block at run time verifying that the reparameterisation reproduces the Li et al. probabilities in the cured arm.

**C. Age strata blended.** v1 used the 25–39y asymptomatic values only; the two strata are now blended with 55% weight on ages 15–24 (`PID_W_15_24`). The implied hazard is invariant to this weight.

**D. Partner pathway corrections.** Partner symptom status is drawn once per individual rather than twice independently; partner PID is gated on the partner actually being infected; and partner PID is partly averted by successful partner treatment, so notification now averts something.

**E. Common random numbers across arms.** The with-test and without-test arms previously drew separate cohorts, so everything unrelated to the testing decision failed to cancel in the difference. Because that difference is divided by the test-positive share (~0.01 for women), the residual noise was amplified roughly 95-fold, producing threshold intervals of ±$200 that were a Monte Carlo artefact rather than uncertainty. `draw_cohort()` now generates every individual-level random variable once per PSA iteration and runs both arms through the identical cohort, using uniform draws compared against arm-dependent probabilities to keep the coupling monotone. Parameter uncertainty is untouched.

Additionally, **population-weighted thresholds are now computed from the joint per-draw samples** (weighted, then quantiled) rather than by weighting the summarised bounds, so the three groups' independent noise adds in quadrature. This addresses the reviewer's point about naive summation of interval bounds.

## Contents

| File / folder | Description |
| --- | --- |
| `gonorrhea_model_v2_29 Aug.R` | Full analysis script: PSA over resistance scenarios and testing panels, cost-neutral threshold calculation, population-weighted thresholds, cost-category disaggregation, and the one-way deterministic sensitivity analysis (Supplement S3). |
| `figure compilation code/gonorrhea_figures_v4.R` | Compiles the manuscript figures from the saved model output. |
| `figure compilation code/make_tornado_v4.R` | Standalone tornado-diagram script; reads the saved DSA results, so the model does not need to be re-run to regenerate the figure. |
| `results/` | Model output files and final figures used in the manuscript. |
| `LICENSE` | MIT license. |
| `README.md` | This file. |

## What the scripts produce

`gonorrhea_model_v2_29 Aug.R` writes the following to the working directory (copied into `results/` in this repository):

| Output | Content |
| --- | --- |
| `gonorrhea_costs_summary_v3.csv` | Full per-scenario cost summary (whole-cohort and per-gonorrhea-positive means with 95% uncertainty intervals), disaggregated by cost category. |
| `gonorrhea_thresholds_v3.csv` | Cost-neutral threshold price per resistance test, by group, resistance scenario, and testing panel (dimension A in Table 2). |
| `gonorrhea_thresholds_weighted_v3.csv` | Population-weighted thresholds (MSM 53% / MSW 21% / women 26%), computed from the joint per-draw samples. |
| `gonorrhea_summary_table_v3.csv` | Trimmed, manuscript-oriented subset of the summary. |
| `gonorrhea_S3_dsa_results_v3.csv` | One-way deterministic sensitivity analysis over eight parameters (Supplement S3). |
| `threshold_heatmap_v3_NDA0.png`, `threshold_heatmap_v3_NDA5.png` | Cost-neutral threshold prices across ciprofloxacin × ceftriaxone resistance, by strategy and population group; one file per novel-drug-A resistance scenario (both produced in a single run). |
| `gonorrhea_S3_tornado_diagram_v3.png` | Tornado diagrams for the DSA. |


The per-episode cost figures (Figures 2 and 3) and the dimension-A rows of Table 2 are read directly from these outputs. The system-level (dimensions B and C) and drug-development (dimension D) values in Table 2 are computed from these per-episode outputs combined with published transmission-model trajectories and national test volumes, as described in the Methods and Supplement S2.


## Requirements

- R (>= 4.2 recommended)
- `tidyverse` (dplyr, ggplot2, tidyr, readr, tibble, purrr)
- `patchwork` (tornado figures)
- `scales` (installed with the tidyverse)

```r
install.packages(c("tidyverse", "patchwork"))
```

## Reproducing the results

From the repository root:

```bash
Rscript "gonorrhea_model_v2_29 Aug.R"
```

or, from an interactive R session in the repository root:

```r
source("gonorrhea_model_v2_29 Aug.R")
```

(The quotes are required — the filename contains a space.)

To regenerate only the tornado figure from saved DSA output, set the working directory to `results/` and run:

```r
source("../figure compilation code/make_tornado_v4.R")
```

**Seeds.** The global seed is set at the top of the model script (`set.seed(123)`), so the main PSA is fully reproducible. The deterministic sensitivity analysis uses per-row seeds so that the low, base, and high estimates for each parameter differ only in the parameter being varied, not in the underlying stochastic draws.

**Runtime.** The full run is 1,000 PSA draws × 20,000 simulated individuals per group × resistance scenario × testing panel (3 groups × 4 ciprofloxacin levels × 7 ceftriaxone levels × 2 novel-drug levels × 3 panels), and takes on the order of several hours on a typical desktop. To check the pipeline quickly, reduce `n_psa` and `n_indiv` near the top of the script; threshold estimates will be noisier but the structure of the results is preserved.

**Diagnostics.** The script prints two checks at run time: the PID hazard QC block described above, and a noise diagnostic reporting the correlation between each scenario's interval half-width and its point estimate. Values near +1 indicate parameter-driven uncertainty; values near 0 indicate residual simulation noise and are a signal to raise `n_indiv`.

## Model structure (brief)

Each simulated individual is assigned an infection status, symptom status, and (if infected) an independent resistance profile to ciprofloxacin, ceftriaxone, and a novel drug A. Individuals are tested for gonorrhea; those testing positive receive a resistance panel (CIP, CIP+CRO, or CIP+CRO+NDA) that guides first-line drug choice. Treatment failure among symptomatic, returning individuals triggers a follow-up visit, repeat diagnostics, and salvage therapy. Partners are infected probabilistically, notified and treated with the modelled notification probability, and their own outcomes are tracked and costed. Sex-specific complications (PID and sequelae for women; epididymitis and DGI for men) are accrued as a function of infection duration. Full equations are given in Supplement S1.

Parameters are drawn from published sources (Table 1 of the manuscript); the model is parameterized rather than statistically fitted or calibrated to a target dataset.

## Version history

| Release | DOI | Corresponds to |
| --- | --- | --- |
| `v2.0.0` | [10.5281/zenodo.22166671] (https://doi.org/10.5281/zenodo.22166671) | Revised manuscript (response to peer review), August 2026 |
| `v1.0.0` | [10.5281/zenodo.22093355](https://doi.org/10.5281/zenodo.22093355) | Original submitted manuscript |

Each Zenodo release is a permanent, immutable snapshot. Results reported in a given version of the manuscript are reproducible from the correspondingly tagged release.

## Citation

If you use this code, please cite the accompanying manuscript and the archived release corresponding to the version you used.


