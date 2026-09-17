functions {
  real ordered_logistic_hierarchical_lpmf(int y, real mu, vector thres) {
    return ordered_logistic_lpmf(y | mu, thres);
  }
}

data {
  int<lower=1> N;
  array[N] int Y;
  vector<lower=0>[N] weights;

  // Design matrices
  int<lower=1> K;  // number of population-level effects (age_centered, year_scaled)
  matrix[N, K] X;  // population-level design matrix

  // Hierarchical structure
  int<lower=1> nthres;  // number of thresholds (7 thresholds for 8 categories: Day 0,1,2,...,7)
  int<lower=1> nregions;
  int<lower=1> ncountries;
  array[N] int<lower=1, upper=ncountries> country;
  array[ncountries] int<lower=1, upper=nregions> country_to_region;

  // Region random slopes design matrix (age_centered, year_scaled)
  int<lower=1> K_region;  // should be 2 (age_centered, year_scaled)
  matrix[N, K_region] X_region;  // design matrix for region random effects
  array[N] int<lower=1, upper=nregions> region;  // region index for each observation

  // Survey random effects with national grouping
  int<lower=1> nsurvey_groups;  // number of unique survey:national combinations
  array[N] int<lower=1, upper=nsurvey_groups> survey_group;  // survey group index

  // Trend prediction grid data
  int<lower=0> N_trends;
  matrix[N_trends, K] X_trends;
  matrix[N_trends, K_region] X_region_trends;
  array[N_trends] int<lower=1, upper=nregions> region_trends;
  array[N_trends] int<lower=1, upper=ncountries> country_trends;
  array[N_trends] int<lower=1, upper=nsurvey_groups> survey_group_trends;

  // Prediction data for new countries
  int<lower=0> N_new;
  matrix[N_new, K] X_new;
  matrix[N_new, K_region] X_region_new;
  array[N_new] int<lower=1, upper=nregions> region_new;
  array[N_new] int<lower=1, upper=nsurvey_groups> survey_group_new;
  int<lower=0, upper=1> sample_new_levels;

  // Number of distinct new countries (iso3) represented in the N_new grid and a mapping from each grid row to its new-country id.
  // These should be constructed in R so that all rows for a given new iso3 share the same new_country_id, matching brms' per-level behavior.
  int<lower=0> N_new_countries;
  array[N_new] int<lower=1> new_country_id;
  array[N_new_countries] int<lower=1, upper=nregions> new_country_region;

  int prior_only;
}

parameters {
  // Population-level effects
  vector[K] b;  // age_centered, year_scaled coefficients

  // Regional threshold baselines
  array[nregions] ordered[nthres] region_thresholds;

  // Country extreme category inflation parameters: positive values make the category more likely
  vector[ncountries] country_inflate_cat0;
  vector[ncountries] country_inflate_cat7;

  // Hierarchical variance for country inflation effects
  real<lower=0> sigma_country_inflate_cat0;
  real<lower=0> sigma_country_inflate_cat7;

  // Region random slopes, i.e.: (-1 + age_centered + year_scaled | adopa_region)
  vector<lower=0>[K_region] sigma_region_slopes;
  cholesky_factor_corr[K_region] L_region_slopes;
  matrix[nregions, K_region] z_region_slopes;

  // Country random effects, i.e.: (intercept + year slope | iso3)
  vector<lower=0>[2] sigma_country_re;
  cholesky_factor_corr[2] L_country_re;
  matrix[ncountries, 2] z_country_re;

  // Survey random intercepts grouped by national, i.e.: (1 | gr(survey, by = national))
  real<lower=0> sigma_survey;
  vector[nsurvey_groups] z_survey;
}

transformed parameters {
  // Final country thresholds
  array[ncountries] vector[nthres] country_thresholds;

  // Region random slopes
  matrix[nregions, K_region] region_slopes;

  // Country random intercepts and year slopes
  vector[ncountries] country_intercepts;
  vector[ncountries] country_year_slopes;
  matrix[ncountries, 2] country_re;

  // Survey random intercepts
  vector[nsurvey_groups] survey_intercepts;

  real lprior = 0;

  // Region random slopes
  region_slopes = z_region_slopes * diag_pre_multiply(sigma_region_slopes, L_region_slopes)';

  // Country random intercepts and year slopes
  country_re = (diag_pre_multiply(sigma_country_re, L_country_re) * z_country_re')';
  country_intercepts = country_re[, 1];
  country_year_slopes = country_re[, 2];

  // Survey random intercepts
  survey_intercepts = z_survey * sigma_survey;

  // Country thresholds (same construction as slopes model)
  for (c in 1:ncountries) {
    int r = country_to_region[c];

    // Step 1: Calculate new extreme thresholds
    real new_first_threshold = region_thresholds[r][1] + country_inflate_cat0[c];
    real new_last_threshold = region_thresholds[r][nthres] - country_inflate_cat7[c];

    // Step 2: Calculate new range
    real new_range = new_last_threshold - new_first_threshold;

    // Step 3: Ensure minimum range for identifiability
    real min_total_range = 1.0;
    real min_spacing = 0.1;

    if (new_range < min_total_range) {
      real expansion_needed = (min_total_range - new_range) * 0.5;
      new_first_threshold -= expansion_needed;
      new_last_threshold += expansion_needed;
      new_range = min_total_range;
    }

    // Step 4: Set extreme thresholds
    country_thresholds[c][1] = new_first_threshold;
    country_thresholds[c][nthres] = new_last_threshold;

    // Step 5: Interpolate middle thresholds proportionally
    if (nthres > 2) {
      real regional_range = region_thresholds[r][nthres] - region_thresholds[r][1];

      for (t in 2:(nthres - 1)) {
        real regional_prop = (region_thresholds[r][t] - region_thresholds[r][1]) / regional_range;
        real proposed_threshold = new_first_threshold + regional_prop * new_range;
        country_thresholds[c][t] = fmax(proposed_threshold,
                                        country_thresholds[c][t - 1] + min_spacing);
      }

      for (t in 2:nthres) {
        country_thresholds[c][t] = fmax(country_thresholds[c][t],
                                        country_thresholds[c][t - 1] + min_spacing);
      }
    }
  }

  // Priors
  for (r in 1:nregions) {
    lprior += student_t_lpdf(region_thresholds[r] | 3, 0, 2.5);
  }

  lprior += normal_lpdf(country_inflate_cat0 | 0, sigma_country_inflate_cat0);
  lprior += normal_lpdf(country_inflate_cat7 | 0, sigma_country_inflate_cat7);

  lprior += student_t_lpdf(sigma_country_inflate_cat0 | 3, 0, 1);
  lprior += student_t_lpdf(sigma_country_inflate_cat7 | 3, 0, 1);
  lprior += student_t_lpdf(sigma_region_slopes | 3, 0, 2.5);
  lprior += student_t_lpdf(sigma_country_re | 3, 0, 2.5);
  lprior += student_t_lpdf(sigma_survey | 3, 0, 2.5);
  lprior += lkj_corr_cholesky_lpdf(L_region_slopes | 1);
  lprior += lkj_corr_cholesky_lpdf(L_country_re | 1);

  for (r in 1:nregions) {
    lprior += std_normal_lpdf(z_region_slopes[r, ]);
  }
  lprior += std_normal_lpdf(to_vector(z_country_re));
  lprior += std_normal_lpdf(z_survey);
}

model {
  if (!prior_only) {
    for (n in 1:N) {
      real mu = dot_product(X[n, ], b) +
                dot_product(X_region[n, ], region_slopes[region[n], ]) +
                country_intercepts[country[n]] +
                country_year_slopes[country[n]] * X[n, 2] +
                survey_intercepts[survey_group[n]];

      target += weights[n] *
                ordered_logistic_hierarchical_lpmf(Y[n] | mu, country_thresholds[country[n]]);
    }
  }

  target += lprior;
  target += std_normal_lpdf(b);
}

generated quantities {
  // Output for analysis
  array[ncountries] vector[nthres] b_country_thresholds;
  array[nregions] vector[nthres] b_region_thresholds;
  vector[ncountries] b_country_inflate_cat0;
  vector[ncountries] b_country_inflate_cat7;
  matrix[nregions, K_region] b_region_slopes;
  vector[ncountries] b_country_intercepts;
  vector[ncountries] b_country_year_slopes;
  vector[nsurvey_groups] b_survey_intercepts;

  // Trend predictions
  array[N_trends] vector[8] trend_category_probs;
  vector[N_trends] trend_expected_days;

  // New country predictions
  array[N_new] vector[8] new_category_probs;
  vector[N_new] y_pred_new;

  // Copy effects
  for (c in 1:ncountries) {
    b_country_thresholds[c] = country_thresholds[c];
    b_country_inflate_cat0[c] = country_inflate_cat0[c];
    b_country_inflate_cat7[c] = country_inflate_cat7[c];
    b_country_intercepts[c] = country_intercepts[c];
    b_country_year_slopes[c] = country_year_slopes[c];
  }

  for (r in 1:nregions) {
    b_region_thresholds[r] = region_thresholds[r];
    for (k in 1:K_region) {
      b_region_slopes[r, k] = region_slopes[r, k];
    }
  }

  for (sg in 1:nsurvey_groups) {
    b_survey_intercepts[sg] = survey_intercepts[sg];
  }

  // Trend predictions
  if (N_trends > 0) {
    for (n in 1:N_trends) {
      real mu_trend = dot_product(X_trends[n, ], b) +
                      dot_product(X_region_trends[n, ], region_slopes[region_trends[n], ]) +
                      country_intercepts[country_trends[n]] +
                      country_year_slopes[country_trends[n]] * X_trends[n, 2];

      vector[nthres] thresholds = country_thresholds[country_trends[n]];
      vector[8] probs;

      probs[1] = inv_logit(thresholds[1] - mu_trend);
      for (k in 2:7) {
        probs[k] = inv_logit(thresholds[k] - mu_trend) -
                   inv_logit(thresholds[k - 1] - mu_trend);
      }
      probs[8] = 1 - inv_logit(thresholds[nthres] - mu_trend);

      for (k in 1:8) {
        probs[k] = fmax(probs[k], 1e-10);
      }
      probs = probs / sum(probs);

      trend_category_probs[n] = probs;

      trend_expected_days[n] = 0.0;
      for (k in 1:8) {
        trend_expected_days[n] += (k - 1) * probs[k];
      }
    }
  }

  // New country predictions with per iso3 random effects
  if (N_new > 0) {
    // Per-new-country random effects, drawn once per posterior draw
    array[N_new_countries] vector[nthres] new_country_thresholds;
    vector[N_new_countries] new_country_int;
    vector[N_new_countries] new_country_year;

    if (sample_new_levels == 1 && N_new_countries > 0) {
      for (c_new in 1:N_new_countries) {
        int r = new_country_region[c_new];

        // Sample new-country threshold inflations
        real new_inflate_cat0 = normal_rng(0, sigma_country_inflate_cat0);
        real new_inflate_cat7 = normal_rng(0, sigma_country_inflate_cat7);

        // Build thresholds for this new country (corrected algorithm)
        real new_first = region_thresholds[r][1] + new_inflate_cat0;
        real new_last = region_thresholds[r][nthres] - new_inflate_cat7;
        real new_range = fmax(new_last - new_first, 0.1);

        new_country_thresholds[c_new][1] = new_first;
        new_country_thresholds[c_new][nthres] = new_first + new_range;

        if (nthres > 2) {
          real regional_range = region_thresholds[r][nthres] - region_thresholds[r][1];
          for (t in 2:(nthres - 1)) {
            real regional_prop =
              (region_thresholds[r][t] - region_thresholds[r][1]) / regional_range;
            new_country_thresholds[c_new][t] =
              new_first + regional_prop * new_range;
          }
        }

        // Sample country intercept and year slope
        {
          vector[2] new_country_effect = multi_normal_cholesky_rng(
            rep_vector(0, 2),
            diag_pre_multiply(sigma_country_re, L_country_re)
          );
          new_country_int[c_new] = new_country_effect[1];
          new_country_year[c_new] = new_country_effect[2];
        }
      }
    }

    for (n in 1:N_new) {
      vector[nthres] thresholds_n;
      real int_n;
      real year_n;
      real mu_new;
      vector[8] probs;

      if (sample_new_levels == 1 && N_new_countries > 0) {
        int cid = new_country_id[n];
        int r = new_country_region[cid];

        thresholds_n = new_country_thresholds[cid];
        int_n = new_country_int[cid];
        year_n = new_country_year[cid];

        mu_new = dot_product(X_new[n, ], b) +
                 dot_product(X_region_new[n, ], region_slopes[r, ]) +
                 int_n +
                 year_n * X_new[n, 2];
      } else {
        int r = region_new[n];
        thresholds_n = region_thresholds[r];
        int_n = 0;
        year_n = 0;

        mu_new = dot_product(X_new[n, ], b) +
                 dot_product(X_region_new[n, ], region_slopes[r, ]);
      }

      // Category probabilities
      probs[1] = inv_logit(thresholds_n[1] - mu_new);
      for (k in 2:7) {
        probs[k] = inv_logit(thresholds_n[k] - mu_new) -
                   inv_logit(thresholds_n[k - 1] - mu_new);
      }
      probs[8] = 1 - inv_logit(thresholds_n[nthres] - mu_new);

      for (k in 1:8) {
        probs[k] = fmax(probs[k], 1e-10);
      }
      probs = probs / sum(probs);

      new_category_probs[n] = probs;

      // Categorical sample for compatibility
      y_pred_new[n] = categorical_rng(probs);
    }
  }
}

