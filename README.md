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

## Contents

| File | Description |
|------|-------------|
| `R/gonorrhea_rgt_model.R` | Full analysis script: PSA over resistance scenarios and testing panels, cost-neutral threshold calculation, cost-category disaggregation, and the one-way deterministic sensitivity analysis (Supplement S3). |
| `README.md` | This file. |
| `LICENSE` | MIT license. |

## What the script produces

Running the script writes the following files to the working directory:

- `gonorrhea_costs_summary.csv` — full per-scenario cost summary (whole-cohort and per-gonorrhea-positive means with 95% uncertainty intervals), disaggregated by cost category.
- `gonorrhea_thresholds.csv` — cost-neutral threshold price per resistance test, by group, resistance scenario, and testing panel (dimension A in Table 2).
- `gonorrhea_summary_table.csv` — a trimmed, manuscript-oriented subset of the summary.
- `threshold_heatmap.png` — cost-neutral threshold prices across ciprofloxacin x ceftriaxone resistance, by strategy and population group (Figure 4).
- `gonorrhea_S3_dsa_results.csv` and `gonorrhea_S3_tornado_diagram.png` — one-way deterministic sensitivity analysis (Supplement S3).

The per-episode cost figures (Figures 2 and 3) and the dimension-A rows of Table 2
are read directly from these outputs. The system-level (dimensions B and C) and
drug-development (dimension D) values in Table 2 are computed from these per-episode
outputs combined with published transmission-model trajectories and national test
volumes, as described in the manuscript Methods and Supplement S2.

## Requirements

- R (>= 4.2 recommended)
- The [`tidyverse`](https://www.tidyverse.org/) meta-package (`dplyr`, `ggplot2`, `tidyr`, `readr`, `tibble`, `purrr`)
- [`patchwork`](https://patchwork.data-imaginist.com/) (for the S3 tornado figure)
- `scales` (installed with the tidyverse)

Install the dependencies with:

```r
install.packages(c("tidyverse", "patchwork"))
```

## Reproducing the results

```bash
Rscript R/gonorrhea_rgt_model.R
```

or, from an interactive R session in the repository root:

```r
source("R/gonorrhea_rgt_model.R")
```

The global random seed is set at the top of the script (`set.seed(123)`), so the
main PSA is fully reproducible. The deterministic sensitivity analysis sets its own
per-parameter seeds so that the low, base, and high estimates for each parameter
differ only in the parameter being varied (not in the underlying stochastic draws).

**Runtime.** The full run is 1,000 PSA draws x 20,000 simulated individuals per
group x resistance scenario x testing panel, and takes on the order of a few hours
on a typical desktop. To check the pipeline quickly, reduce `n_psa` and `n_indiv`
near the top of the script; threshold estimates will be noisier but the structure
of the results is preserved.

## Model structure (brief)

Each simulated individual is assigned an infection status, symptom status, and (if
infected) an independent resistance profile to ciprofloxacin, ceftriaxone, and a
novel drug A. Individuals are tested for gonorrhea; those testing positive receive
a resistance panel (single-, dual-, or triple-target) that guides first-line drug
choice. Treatment failure among symptomatic, returning individuals triggers a
follow-up visit, repeat diagnostics, and salvage therapy. Partner treatment and
sex-specific complications (PID and sequelae for women; epididymitis and DGI for
men) are tracked and costed. Full equations are given in Supplement S1.

Parameters are drawn from published sources (Table 1 of the manuscript); the model
is parameterized rather than statistically fitted or calibrated to a target
dataset.

## Citation

If you use this code, please cite the accompanying manuscript. A permanent
archived release of this repository is available via Zenodo (DOI to be added on
publication).

