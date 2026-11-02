# AdoPA: Adolescent Physical Activity — Hierarchical Ordinal Model

Replication code for Strain T, Flaxman SR, Ricardo LIC, et al, and the Country Data Author Group. National, regional, and global trends in adolescent physical activity from 2010 to 2025: a pooled analysis of 389 surveys with 3·2 million participants from 140 countries. Lancet Glob Health 2026; published online 4 November 2026. https://doi.org/10.1016/j.langlo.2026.104099 

## Overview

This repository contains the statistical model and analysis script for estimating trends in
adolescent physical activity, using a hierarchical ordinal Bayesian model, with inflated tails,
 fit with Stan.

## Required R packages

```r
install.packages(c("rstan", "data.table", "dplyr", "tidyr", "here", "broom.mixed"))
```

## Data

- `data/surveys.csv` — survey-level physical activity data
- `data/covariates.csv` — country-level covariates and region mappings

## Running the analysis

```bash
Rscript fit.r
```

This script:
1. Loads and prepares the data
2. Compiles the Stan model `model.stan` 
3. Fits separate models for boys and girls (4 chains × 500 post-warmup iterations, thinned by 2)
4. Extracts predictions for all in-sample countries and out-of-sample countries
5. Writes results to `results/`

_Note that the runtime is 10+ hours!_

## Output

All output is written to `results/` (created at runtime):

| File | Contents |
|------|----------|
| `predictions_<model>_<timestamp>.csv` | Posterior means and 95% CIs for all countries × ages × years × sex |
| `mcmc_draws_<model>_<timestamp>_all_countries.csv` | Full per-draw probability matrix (all 8 physical activity categories) |
| `ppc_<model>_<timestamp>.csv` | Posterior predictive check data |

## Model

The Stan model (`model.stan`) is a hierarchical ordinal regression with:
- 8 ordered categories (0–7 days/week meeting physical activity recommendations)
- Regional random slopes for age and year trends
- Country-level random intercepts and year slopes
- Survey-level random intercepts
- Out-of-sample predictions for countries not in the training data

