################################################################################
# fit.r — AdoPA model fitting and prediction export
#
# Fits a hierarchical ordinal Stan model for adolescent physical activity data.
# Saves prediction summaries, PPC data, diagnostics, and an MCMC draws CSV
# to the results/ directory.
#
# Usage:
#   Rscript fit.r
#
# Output:
#   results/predictions_<model>_<timestamp>.csv
#   results/mcmc_draws_<model>_<timestamp>_all_countries.csv
#   results/ppc_<model>_<timestamp>.csv
################################################################################

options(stringsAsFactors = FALSE)
set.seed(42)

library(rstan)
library(data.table)
library(dplyr)
library(tidyr)
options(mc.cores = parallel::detectCores())

################################################################################
# MODEL CONFIGURATION
################################################################################

config <- list(
  name       = "model",
  stan_file  = "model.stan",
  chains     = 4,
  iter       = 400,
  warmup     = 200,
  thin       = 1,
  control    = list(adapt_delta = 0.9, max_treedepth = 15),
  seed       = 42
)
config$timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
prefix <- paste0(config$name, "_", config$timestamp)
cat("Model:", config$name, "\n")
cat("Stan file:", config$stan_file, "\n")
cat("Timestamp:", config$timestamp, "\n")

cat("Compiling Stan model...\n")
stan_model_obj <- stan_model(file = config$stan_file, model_name = config$name)
cat("Compilation complete.\n")

################################################################################
# CSV UTILITIES
################################################################################

generate_csv_from_pred_data <- function(pred_data, model_config, timestamp = NULL,
                                        output_dir = "results") {
  if (is.null(timestamp)) timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  csv_filename <- file.path(output_dir,
    paste0("predictions_", model_config, "_", timestamp, ".csv"))

  if (file.exists(csv_filename)) {
    cat("CSV already exists:", csv_filename, "\n")
    return(csv_filename)
  }
  cat("Generating CSV:", csv_filename, "\n")
  if (nrow(pred_data) == 0) {
    cat("Warning: pred_data is empty, cannot generate CSV\n")
    return(NULL)
  }

  tryCatch({
    csv_data <- data.frame(
      iso3              = pred_data$iso3,
      country_name      = pred_data$whoname,
      region            = pred_data$adopa_region,
      year              = pred_data$year,
      age               = pred_data$age,
      sex               = pred_data$sex,
      day0_prob_mean    = pred_data$day0_prob,
      day1_prob_mean    = pred_data$day1_prob,
      day2_prob_mean    = pred_data$day2_prob,
      day3_prob_mean    = pred_data$day3_prob,
      day4_prob_mean    = pred_data$day4_prob,
      day5_prob_mean    = pred_data$day5_prob,
      day6_prob_mean    = pred_data$day6_prob,
      day7_prob_mean    = pred_data$day7_prob,
      day1to3_prob_mean = pred_data$day1to3_prob,
      day4to6_prob_mean = pred_data$day4to6_prob,
      mean_days         = pred_data$pred_mean,
      mean_days_lower   = pred_data$pred_q025,
      mean_days_upper   = pred_data$pred_q975,
      day0_prob_lower   = pred_data$day0_prob_q025,
      day0_prob_upper   = pred_data$day0_prob_q975,
      day7_prob_lower   = pred_data$day7_prob_q025,
      day7_prob_upper   = pred_data$day7_prob_q975,
      stringsAsFactors  = FALSE
    )
    csv_data <- csv_data[order(csv_data$iso3, csv_data$year, csv_data$age, csv_data$sex), ]
    write.csv(csv_data, csv_filename, row.names = FALSE)
    cat("Wrote predictions CSV:", csv_filename, "(", nrow(csv_data), "rows)\n")
    return(csv_filename)
  }, error = function(e) {
    cat("Error generating predictions CSV:", e$message, "\n")
    return(NULL)
  })
}

################################################################################
# MCMC DRAWS HELPER
################################################################################

make_prob_df <- function(prob_array, grid, sex_value) {
  # prob_array: [n_draws, N_points, 8]
  base_df <- data.frame(
    iso3 = grid$iso3, year = grid$year, age = grid$age, sex = sex_value,
    stringsAsFactors = FALSE
  )
  prob_cols <- list()
  for (day in 0:7) {
    mat <- t(prob_array[, , day + 1])   # [N_points, n_draws]
    colnames(mat) <- paste0("days", day, "_prob_", seq_len(ncol(mat)))
    prob_cols[[length(prob_cols) + 1]] <- as.data.frame(mat)
  }
  weights   <- 0:7
  mean_mat  <- matrix(0, nrow = dim(prob_array)[2], ncol = dim(prob_array)[1])
  for (k in seq_len(8)) {
    mean_mat <- mean_mat + weights[k] * t(prob_array[, , k])
  }
  colnames(mean_mat) <- paste0("mean_days_", seq_len(ncol(mean_mat)))
  cbind(base_df, do.call(cbind, prob_cols), as.data.frame(mean_mat))
}

################################################################################
# DATA LOADING AND PREPARATION
################################################################################

data_file <- "data/surveys.csv"
df <- fread(data_file)
cat("Data loaded from:", data_file, "\n")

df <- df %>% mutate(year = floor((beginyear + endyear) / 2))
covars <- fread("data/covariates.csv")
df <- df %>% select(-whoname)
df <- df %>% left_join(covars, by = c("iso3", "adopa_region", "year", "low_school_attendance"))

df_clean <- df %>%
  filter(exclPE == 0, year >= 2008) %>%
  mutate(
    trials       = ess,
    count_day0   = ess * day0,
    count_day1   = ess * day1,
    count_day2   = ess * day2,
    count_day3   = ess * day3,
    count_day4   = ess * day4,
    count_day5   = ess * day5,
    count_day6   = ess * day6,
    count_day7   = ess * day7,
    year_scaled       = (year - 2011) / 10,
    age_centered      = age - 14,
    adminlevel_factor = factor(adminlevel),
    sample_frame_factor = factor(sample_frame),
    adopa_region      = as.factor(adopa_region),
    iso3              = as.factor(iso3),
    national          = adminlevel == "N",
    regional          = adminlevel == "R"
  ) %>%
  filter(!is.na(count_day0), !is.na(count_day1), !is.na(count_day2), !is.na(count_day3),
         !is.na(count_day4), !is.na(count_day5), !is.na(count_day6), !is.na(count_day7))

cat("Running on ALL countries. Total countries:", length(unique(df_clean$iso3)), "\n")

create_aggregated_ordinal_data <- function(df_subset) {
  core_vars <- c("age_centered", "year_scaled", "adopa_region", "iso3",
                 "survey", "sdi",
                 "oosrate", "probitobesity", "low_school_attendance", "veh_per_cap")
  long_data <- df_subset %>%
    select(any_of(core_vars), starts_with("count_day")) %>%
    mutate(obs_id = row_number()) %>%
    pivot_longer(cols = starts_with("count_day"),
                 names_to  = "days_activity_chr",
                 values_to = "count") %>%
    mutate(
      days_activity_numeric = as.numeric(gsub("count_day", "", days_activity_chr)),
      response = factor(days_activity_numeric, levels = 0:7, ordered = TRUE),
      ageCat   = cut(age_centered + 14, c(11, 13, 16, 18), right = FALSE)
    ) %>%
    select(obs_id, any_of(core_vars), response, days_activity_numeric, count, ageCat)
  return(long_data)
}

################################################################################
# STAN DATA PREPARATION
################################################################################

prepare_custom_stan_data <- function(data, sex_label, config) {
  cat("Preparing Stan data for", sex_label, "...\n")

  countries <- levels(data$iso3)
  regions   <- levels(data$adopa_region)

  country_region_map <- data %>%
    select(iso3, adopa_region) %>%
    distinct() %>%
    arrange(iso3) %>%
    mutate(
      country_idx = as.numeric(iso3),
      region_idx  = as.numeric(adopa_region)
    )

  all_vars <- c("age_centered", "year_scaled")
  X_combined     <- as.matrix(data[, all_vars, drop = FALSE])

  stan_data <- list(
    N       = nrow(data),
    Y       = as.numeric(data$response),
    weights = data$count,
    K       = length(all_vars),
    X       = X_combined,
    nthres  = 7L,
    nregions   = length(regions),
    ncountries = length(countries),
    country    = as.numeric(data$iso3),
    country_to_region = country_region_map$region_idx,
    prior_only = 0L
  )

  # Survey groupings (one group per unique survey ID)
  survey_groups <- data %>%
    select(survey) %>%
    distinct() %>%
    mutate(survey_group_idx = row_number())

  data_with_groups <- data %>%
    left_join(survey_groups, by = "survey")

  X_region_matrix <- cbind(data$age_centered, data$year_scaled)

  stan_data <- c(stan_data, list(
    K_region       = 2L,
    X_region       = X_region_matrix,
    region         = as.numeric(data$adopa_region),
    nsurvey_groups = nrow(survey_groups),
    survey_group   = data_with_groups$survey_group_idx,
    N_new          = 0L,
    X_new          = matrix(0, nrow = 0, ncol = length(all_vars)),
    X_region_new   = matrix(0, nrow = 0, ncol = 2),
    region_new     = integer(0),
    survey_group_new  = integer(0),
    sample_new_levels = 0L
  ))

  # Trend prediction grid: all ages (11-17), years 2008-2025, all countries
  years <- seq(2008, 2025, by = 1)
  ages  <- seq(11, 17,   by = 1)
  trend_grid <- expand.grid(
    iso3 = countries, age = ages, year = years,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      age_centered = age - 14,
      year_scaled  = (year - 2011) / 10
    )
  cat("Trend prediction grid:", nrow(trend_grid), "points (",
      length(countries), "countries x", length(years), "years x", length(ages), "ages)\n")

  country_indices      <- match(trend_grid$iso3, countries)
  region_indices       <- stan_data$country_to_region[country_indices]
  survey_group_indices <- rep(1L, nrow(trend_grid))

  X_trends_combined <- as.matrix(trend_grid[, all_vars, drop = FALSE])

  stan_data <- c(stan_data, list(
    N_trends             = nrow(trend_grid),
    X_trends             = X_trends_combined,
    X_region_trends      = cbind(trend_grid$age_centered, trend_grid$year_scaled),
    region_trends        = region_indices,
    country_trends       = country_indices,
    survey_group_trends  = survey_group_indices
  ))

  trend_info <- list(trend_grid = trend_grid, countries = countries,
                     years = years, ages = ages)

  # New-country prediction grid
  new_info <- NULL
  covar_countries   <- sort(unique(as.character(covars$iso3)))
  training_countries <- as.character(countries)
  new_countries     <- setdiff(covar_countries, training_countries)

  if (length(new_countries) > 0) {
    covar_region_map <- covars %>%
      filter(iso3 %in% new_countries) %>%
      select(iso3, adopa_region) %>%
      distinct()

    new_grid <- expand.grid(
      iso3 = new_countries, age = ages, year = years,
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        age_centered = age - 14,
        year_scaled  = (year - 2011) / 10,
        adopa_region = covar_region_map$adopa_region[match(iso3, covar_region_map$iso3)],
        region_idx   = as.numeric(factor(adopa_region,
                         levels = levels(data$adopa_region)))
      )

    new_grid <- new_grid[!is.na(new_grid$region_idx) & new_grid$region_idx >= 1, , drop = FALSE]

    if (nrow(new_grid) > 0) {
      stan_data$N_new           <- nrow(new_grid)
      stan_data$X_new           <- as.matrix(new_grid[, all_vars, drop = FALSE])
      stan_data$X_region_new    <- cbind(new_grid$age_centered, new_grid$year_scaled)
      stan_data$region_new      <- as.integer(new_grid$region_idx)
      stan_data$survey_group_new <- rep(1L, nrow(new_grid))
      stan_data$sample_new_levels <- 1L

      # Map each new-country row to a new-country id (per iso3), for shared RE per iso3 per draw
      new_iso_levels     <- sort(unique(new_grid$iso3))
      new_country_factor <- factor(new_grid$iso3, levels = new_iso_levels)
      stan_data$N_new_countries <- length(new_iso_levels)
      stan_data$new_country_id  <- as.integer(new_country_factor)

      new_country_region <- tapply(
        new_grid$region_idx, new_country_factor, function(x) x[1]
      )
      stan_data$new_country_region <- as.integer(new_country_region)

      new_info <- list(new_grid = new_grid)
    }
  }

  return(list(
    stan_data     = stan_data,
    country_names = countries,
    region_names  = regions,
    survey_groups = survey_groups,
    trend_info    = trend_info,
    new_info      = new_info
  ))
}

################################################################################
# MODEL FITTING (rstan)
################################################################################

fit_model <- function(data, sex_label, config) {
  cat("Fitting", sex_label, "Stan model...\n")
  stan_prep <- prepare_custom_stan_data(data, sex_label, config)

  cat("Sampling with", config$chains, "chains x", config$iter, "iter (thin =",
      config$thin, ", warmup =", config$warmup, ")...\n")

  fit <- sampling(
    stan_model_obj,
    data    = stan_prep$stan_data,
    chains  = config$chains,
    iter    = config$iter,
    warmup  = config$warmup,
    thin    = config$thin,
    control = config$control,
    seed    = config$seed,
    verbose = FALSE
  )

  cat(sex_label, "model fitted successfully\n")
  return(list(
    fit           = fit,
    stan_data     = stan_prep$stan_data,
    country_names = stan_prep$country_names,
    region_names  = stan_prep$region_names,
    survey_groups = stan_prep$survey_groups,
    trend_info    = stan_prep$trend_info,
    new_info      = stan_prep$new_info
  ))
}

################################################################################
# PREPARE DATA AND FIT
################################################################################

df_boys  <- df_clean %>% filter(sex == 1)
df_girls <- df_clean %>% filter(sex == 2)

boys_data  <- create_aggregated_ordinal_data(df_boys)
girls_data <- create_aggregated_ordinal_data(df_girls)

start_time <- Sys.time()

boys_results  <- fit_model(boys_data,  "Boys",  config)
girls_results <- fit_model(girls_data, "Girls", config)

elapsed_time <- Sys.time() - start_time
cat("Total elapsed time:", round(elapsed_time, 2), attr(elapsed_time, "units"), "\n")

if (!dir.exists("results")) dir.create("results")

################################################################################
# PREDICTION EXTRACTION
################################################################################

extract_stan_trend_predictions <- function(stan_results, sex_value, sex_label) {
  cat("Extracting trend predictions for", sex_label, "...\n")

  tryCatch({
    # extract() returns [n_draws, N_trends] and [n_draws, N_trends, 8]
    trend_expected <- rstan::extract(stan_results$fit, "trend_expected_days")$trend_expected_days
    if (is.null(trend_expected) || ncol(trend_expected) == 0)
      stop("No trend_expected_days in fit")

    trend_probs <- rstan::extract(stan_results$fit, "trend_category_probs")$trend_category_probs

    n_trends     <- ncol(trend_expected)
    n_categories <- 8L

    # Reconstruct trend grid
    countries <- stan_results$country_names
    years     <- seq(2008, 2025, by = 1)
    ages      <- seq(11, 17, by = 1)
    trend_grid <- expand.grid(
      iso3 = countries, age = ages, year = years,
      stringsAsFactors = FALSE
    ) %>%
      mutate(age_centered = age - 14, year_scaled = (year - 2011) / 10)

    cat("Extracted", n_trends, "trend points for", sex_label, "\n")

    covars_local   <- fread("data/covariates.csv")
    country_meta   <- covars_local %>% select(iso3, whoname, adopa_region) %>% distinct()

    pred_summary <- trend_grid %>%
      mutate(
        sex       = sex_value,
        pred_mean = apply(trend_expected, 2, mean, na.rm = TRUE),
        pred_q025 = apply(trend_expected, 2, quantile, 0.025, na.rm = TRUE),
        pred_q975 = apply(trend_expected, 2, quantile, 0.975, na.rm = TRUE),
        day0_prob = apply(trend_probs[,, 1], 2, mean, na.rm = TRUE),
        day1_prob = apply(trend_probs[,, 2], 2, mean, na.rm = TRUE),
        day2_prob = apply(trend_probs[,, 3], 2, mean, na.rm = TRUE),
        day3_prob = apply(trend_probs[,, 4], 2, mean, na.rm = TRUE),
        day4_prob = apply(trend_probs[,, 5], 2, mean, na.rm = TRUE),
        day5_prob = apply(trend_probs[,, 6], 2, mean, na.rm = TRUE),
        day6_prob = apply(trend_probs[,, 7], 2, mean, na.rm = TRUE),
        day7_prob = apply(trend_probs[,, 8], 2, mean, na.rm = TRUE),
        day1to3_prob = day1_prob + day2_prob + day3_prob,
        day4to6_prob = day4_prob + day5_prob + day6_prob,
        day0_prob_q025 = apply(trend_probs[,,1], 2, quantile, 0.025, na.rm = TRUE),
        day0_prob_q975 = apply(trend_probs[,,1], 2, quantile, 0.975, na.rm = TRUE),
        day1_prob_q025 = apply(trend_probs[,,2], 2, quantile, 0.025, na.rm = TRUE),
        day1_prob_q975 = apply(trend_probs[,,2], 2, quantile, 0.975, na.rm = TRUE),
        day2_prob_q025 = apply(trend_probs[,,3], 2, quantile, 0.025, na.rm = TRUE),
        day2_prob_q975 = apply(trend_probs[,,3], 2, quantile, 0.975, na.rm = TRUE),
        day3_prob_q025 = apply(trend_probs[,,4], 2, quantile, 0.025, na.rm = TRUE),
        day3_prob_q975 = apply(trend_probs[,,4], 2, quantile, 0.975, na.rm = TRUE),
        day4_prob_q025 = apply(trend_probs[,,5], 2, quantile, 0.025, na.rm = TRUE),
        day4_prob_q975 = apply(trend_probs[,,5], 2, quantile, 0.975, na.rm = TRUE),
        day5_prob_q025 = apply(trend_probs[,,6], 2, quantile, 0.025, na.rm = TRUE),
        day5_prob_q975 = apply(trend_probs[,,6], 2, quantile, 0.975, na.rm = TRUE),
        day6_prob_q025 = apply(trend_probs[,,7], 2, quantile, 0.025, na.rm = TRUE),
        day6_prob_q975 = apply(trend_probs[,,7], 2, quantile, 0.975, na.rm = TRUE),
        day7_prob_q025 = apply(trend_probs[,,8], 2, quantile, 0.025, na.rm = TRUE),
        day7_prob_q975 = apply(trend_probs[,,8], 2, quantile, 0.975, na.rm = TRUE)
      ) %>%
      left_join(country_meta, by = "iso3") %>%
      mutate(whoname = ifelse(is.na(whoname) | whoname == "", iso3, whoname))

    cat("Trend predictions extracted:", nrow(pred_summary), "rows for", sex_label, "\n")
    return(pred_summary)

  }, error = function(e) {
    stop("Stan trend extraction failed for ", sex_label, ": ", conditionMessage(e))
  })
}

extract_stan_new_predictions <- function(stan_results, sex_value, sex_label) {
  if (is.null(stan_results$stan_data$N_new) || stan_results$stan_data$N_new <= 0)
    return(data.frame())

  cat("Extracting new-country predictions for", sex_label, "...\n")
  N_new        <- stan_results$stan_data$N_new
  n_categories <- 8L

  tryCatch({
    # model.stan outputs new_category_probs
    new_probs <- rstan::extract(stan_results$fit, "new_category_probs")$new_category_probs
    # shape: [n_draws, N_new, 8]  — already correctly shaped with RStan extract()

    if (is.null(new_probs) || dim(new_probs)[2] != N_new)
      stop("new_category_probs dimensions don't match N_new")

    # Compute expected days from per-draw probabilities
    new_expected_array <- matrix(0, nrow = dim(new_probs)[1], ncol = N_new)
    for (k in seq_len(n_categories)) {
      new_expected_array <- new_expected_array + (k - 1) * new_probs[,, k]
    }
    pred_mean <- apply(new_expected_array, 2, mean,             na.rm = TRUE)
    pred_q025 <- apply(new_expected_array, 2, quantile, 0.025,  na.rm = TRUE)
    pred_q975 <- apply(new_expected_array, 2, quantile, 0.975,  na.rm = TRUE)

    day_prob_mean <- function(k) apply(new_probs[,, k], 2, mean,            na.rm = TRUE)
    day_prob_q025 <- function(k) apply(new_probs[,, k], 2, quantile, 0.025, na.rm = TRUE)
    day_prob_q975 <- function(k) apply(new_probs[,, k], 2, quantile, 0.975, na.rm = TRUE)

    if (is.null(stan_results$new_info) || is.null(stan_results$new_info$new_grid))
      stop("new_info$new_grid not available")

    new_grid     <- stan_results$new_info$new_grid
    covars_local <- fread("data/covariates.csv")
    meta <- covars_local %>% select(iso3, whoname, adopa_region) %>% distinct()

    out <- new_grid %>%
      select(iso3, age, year, age_centered, year_scaled) %>%
      mutate(sex = sex_value) %>%
      mutate(
        pred_mean = pred_mean,
        pred_q025 = pred_q025,
        pred_q975 = pred_q975,
        day0_prob = day_prob_mean(1), day1_prob = day_prob_mean(2),
        day2_prob = day_prob_mean(3), day3_prob = day_prob_mean(4),
        day4_prob = day_prob_mean(5), day5_prob = day_prob_mean(6),
        day6_prob = day_prob_mean(7), day7_prob = day_prob_mean(8),
        day1to3_prob = day1_prob + day2_prob + day3_prob,
        day4to6_prob = day4_prob + day5_prob + day6_prob,
        day0_prob_q025 = day_prob_q025(1), day0_prob_q975 = day_prob_q975(1),
        day1_prob_q025 = day_prob_q025(2), day1_prob_q975 = day_prob_q975(2),
        day2_prob_q025 = day_prob_q025(3), day2_prob_q975 = day_prob_q975(3),
        day3_prob_q025 = day_prob_q025(4), day3_prob_q975 = day_prob_q975(4),
        day4_prob_q025 = day_prob_q025(5), day4_prob_q975 = day_prob_q975(5),
        day5_prob_q025 = day_prob_q025(6), day5_prob_q975 = day_prob_q975(6),
        day6_prob_q025 = day_prob_q025(7), day6_prob_q975 = day_prob_q975(7),
        day7_prob_q025 = day_prob_q025(8), day7_prob_q975 = day_prob_q975(8)
      ) %>%
      left_join(meta, by = "iso3") %>%
      mutate(whoname = ifelse(is.na(whoname) | whoname == "", iso3, whoname))

    cat("New-country predictions extracted:", nrow(out), "rows for", sex_label, "\n")
    return(out)

  }, error = function(e) {
    cat("Error extracting new-country predictions for", sex_label, ":", e$message, "\n")
    return(data.frame())
  })
}

# Generate predictions for both sexes
boys_pred  <- extract_stan_trend_predictions(boys_results,  1, "Boys")
girls_pred <- extract_stan_trend_predictions(girls_results, 2, "Girls")
pred_data <- rbind(boys_pred, girls_pred)
cat("Total in-sample predictions:", nrow(pred_data), "\n")

# Enrich with country metadata
covars_df   <- covars %>% as.data.frame()
covar_meta  <- covars_df %>% select(iso3, whoname, adopa_region) %>% distinct()

pred_data <- pred_data %>%
  mutate(
    adopa_region = as.character(adopa_region),
    whoname = ifelse(is.na(whoname) | whoname == "",
      covar_meta$whoname[match(iso3, covar_meta$iso3)], whoname),
    whoname = ifelse(is.na(whoname) | whoname == "", iso3, whoname),
    adopa_region = ifelse(is.na(adopa_region) | adopa_region == "",
      covar_meta$adopa_region[match(iso3, covar_meta$iso3)], adopa_region)
  )
pred_data$adopa_region <- factor(
  pred_data$adopa_region,
  levels = sort(unique(na.omit(covar_meta$adopa_region)))
)

# Add new-country predictions
existing_iso <- unique(as.character(df_clean$iso3))
covar_iso    <- unique(covar_meta$iso3)
missing_iso  <- setdiff(covar_iso, existing_iso)
missing_iso  <- missing_iso[!is.na(missing_iso) & nzchar(missing_iso)]

if (length(missing_iso) > 0 &&
    (!is.null(boys_results$stan_data$N_new) && boys_results$stan_data$N_new > 0)) {
  boys_new  <- extract_stan_new_predictions(boys_results,  1, "Boys")
  girls_new <- extract_stan_new_predictions(girls_results, 2, "Girls")
  synthetic_predictions <- rbind(boys_new, girls_new)
  if (nrow(synthetic_predictions) > 0) {
    pred_data <- rbind(pred_data, synthetic_predictions)
    cat("Added", nrow(synthetic_predictions), "new-country rows for",
        length(unique(synthetic_predictions$iso3)), "new countries\n")
  }
} else {
  cat("All covariate countries represented in survey data; no new-country rows added\n")
}

cat("Total predictions:", nrow(pred_data), "\n")

################################################################################
# POSTERIOR PREDICTIVE CHECKS
################################################################################

generate_real_ppc <- function(stan_results, ordinal_data, sex_value, sex_label) {
  cat("Generating PPC for", sex_label, "...\n")

  covars_local <- fread("data/covariates.csv")
  country_names_meta <- covars_local %>% select(iso3, whoname) %>% distinct()

  # names(fit) gives element-level names with indices, e.g. "b[1]", "country_thresholds[1,1]"
  available_vars         <- names(stan_results$fit)
  has_country_year_slopes <- any(grepl("^b_country_year_slopes\\[|^country_year_slopes\\[",
                                        available_vars))

  # Build set of base parameter names to extract as a draws matrix
  pars_to_extract <- c("b")
  thresh_param <- if (any(grepl("^country_thresholds\\[", available_vars))) "country_thresholds"
                  else "b_country_thresholds"
  int_param    <- if (any(grepl("^country_intercepts\\[", available_vars))) "country_intercepts"
                  else "b_country_intercepts"
  rslope_param <- if (any(grepl("^region_slopes\\[", available_vars))) "region_slopes"
                  else "b_region_slopes"

  pars_to_extract <- c(pars_to_extract, thresh_param, int_param, rslope_param)
  if (has_country_year_slopes) {
    slope_param <- if (any(grepl("^country_year_slopes\\[", available_vars))) "country_year_slopes"
                   else "b_country_year_slopes"
    pars_to_extract <- c(pars_to_extract, slope_param)
  }
  # Only extract params actually present in fit
  pars_to_extract <- intersect(pars_to_extract, stan_results$fit@model_pars)

  # RStan as.matrix returns named matrix with columns like "b[1]", "b_country_thresholds[1,1]"
  cat("Extracting", length(pars_to_extract), "parameter groups for PPC\n")
  draws   <- as.matrix(stan_results$fit, pars = pars_to_extract)
  n_draws <- nrow(draws)

  get_draw_vector <- function(column_name) {
    if (length(column_name) == 0 || is.na(column_name) || !nzchar(column_name))
      return(rep(0, n_draws))
    if (column_name %in% colnames(draws))
      return(draws[, column_name])
    rep(0, n_draws)
  }

  get_threshold_matrix <- function(column_names) {
    available <- column_names %in% colnames(draws)
    if (!all(available)) return(NULL)
    as.matrix(draws[, column_names, drop = FALSE])
  }

  sample_multinom_counts <- function(size_vec, prob_mat) {
    n_draws_local <- length(size_vec)
    n_cat         <- ncol(prob_mat)
    if (n_draws_local == 0L || n_cat == 0L) return(matrix(0, n_draws_local, n_cat))

    size_vec <- round(size_vec)
    size_vec[!is.finite(size_vec) | size_vec < 0] <- 0
    counts         <- matrix(0, n_draws_local, n_cat)
    remaining      <- as.integer(size_vec)
    prob_remaining <- prob_mat
    prob_remaining[!is.finite(prob_remaining)] <- 0
    prob_remaining <- pmax(prob_remaining, 0)

    for (k in seq_len(n_cat - 1)) {
      current_prob <- ifelse(remaining > 0, prob_remaining[, k], 0)
      current_prob[!is.finite(current_prob)] <- 0
      current_prob <- pmin(pmax(current_prob, 0), 1)
      draws_k      <- stats::rbinom(n_draws_local, size = remaining, prob = current_prob)
      draws_k[!is.finite(draws_k)] <- 0
      counts[, k] <- draws_k
      remaining   <- pmax(remaining - draws_k, 0)
      if (k < n_cat - 1) {
        next_cols <- prob_remaining[, (k + 1):n_cat, drop = FALSE]
        stay      <- pmax(1 - current_prob, 0)
        stay[!is.finite(stay)] <- 0
        non_zero  <- stay > 0 & remaining > 0
        if (any(non_zero)) {
          next_cols[non_zero, ] <- next_cols[non_zero, , drop = FALSE] / stay[non_zero]
        }
        next_cols[!non_zero, ]       <- 0
        next_cols[!is.finite(next_cols)] <- 0
        prob_remaining[, (k + 1):n_cat] <- next_cols
      }
      prob_remaining[, k] <- 0
    }
    counts[, n_cat] <- remaining
    counts
  }

  # Standard model: uses b[1], b[2]
  base_b_age  <- get_draw_vector("b[1]")
  base_b_year <- get_draw_vector("b[2]")

  if (!"age" %in% names(ordinal_data))
    ordinal_data <- ordinal_data %>% mutate(age = age_centered + 14)

  surveys_data <- ordinal_data %>%
    filter(age >= 13, age <= 15) %>%
    group_by(iso3, survey) %>%
    summarise(
      total_count  = sum(count),
      age_centered = mean(age_centered),
      year_scaled  = mean(year_scaled),
      adopa_region = first(adopa_region),
      region_for_ppc = first(adopa_region),
      .groups = "drop"
    )

  if (nrow(surveys_data) == 0) {
    cat("No survey data for", sex_label, "\n")
    return(data.frame())
  }

  ppc_results <- list()

  for (i in seq_len(nrow(surveys_data))) {
    survey_info <- surveys_data[i, ]

    observed_data <- ordinal_data %>%
      filter(iso3 == survey_info$iso3, survey == survey_info$survey,
             age >= 13, age <= 15) %>%
      group_by(days_activity_numeric) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      arrange(days_activity_numeric)

    if (nrow(observed_data) == 0) next

    observed_counts <- observed_data$count
    total_ess       <- sum(observed_counts)

    predicted_counts_matrix <- matrix(0, nrow = n_draws, ncol = 8)
    ages_1315 <- 13:15

    age_sample_sizes <- ordinal_data %>%
      filter(iso3 == survey_info$iso3, survey == survey_info$survey,
             age >= 13, age <= 15) %>%
      group_by(age) %>%
      summarise(age_ess = sum(count), .groups = "drop")

    age_lookup      <- setNames(age_sample_sizes$age_ess, age_sample_sizes$age)
    age_ess_vec     <- as.numeric(age_lookup[as.character(ages_1315)])
    age_ess_vec[is.na(age_ess_vec)] <- 0
    age_centered_vec <- ages_1315 - 14

    if (sum(age_ess_vec) == 0) next

    country_idx <- match(as.character(survey_info$iso3), stan_results$country_names)
    if (is.na(country_idx)) next

    region_idx <- NA_integer_
    if (!is.null(stan_results$stan_data$country_to_region) &&
        length(stan_results$stan_data$country_to_region) >= country_idx) {
      ri <- stan_results$stan_data$country_to_region[country_idx]
      if (is.finite(ri) && !is.na(ri)) region_idx <- as.integer(ri)
    }
    if (is.na(region_idx) || region_idx < 1) {
      region_idx <- match(as.character(survey_info$region_for_ppc), stan_results$region_names)
    }

    region_slope_age  <- if (!is.na(region_idx))
      get_draw_vector(paste0(rslope_param, "[", region_idx, ",1]")) else rep(0, n_draws)
    region_slope_year <- if (!is.na(region_idx))
      get_draw_vector(paste0(rslope_param, "[", region_idx, ",2]")) else rep(0, n_draws)

    country_int  <- get_draw_vector(paste0(int_param,   "[", country_idx, "]"))
    country_year <- if (has_country_year_slopes)
      get_draw_vector(paste0(slope_param, "[", country_idx, "]")) else rep(0, n_draws)

    slope_age  <- base_b_age  + region_slope_age
    slope_year <- base_b_year + region_slope_year + country_year

    mu_mat   <- outer(slope_age, age_centered_vec)
    mu_base  <- slope_year * survey_info$year_scaled + country_int
    mu_mat   <- sweep(mu_mat, 1, mu_base, "+")

    threshold_cols <- paste0(thresh_param, "[", country_idx, ",", 1:7, "]")
    thresholds_mat <- get_threshold_matrix(threshold_cols)
    if (is.null(thresholds_mat)) next

    thresholds_arr <- array(thresholds_mat, dim = c(n_draws, 7, length(ages_1315)))
    thresholds_arr <- aperm(thresholds_arr, c(1, 3, 2))
    mu_arr <- array(mu_mat, dim = c(n_draws, length(ages_1315), 7))

    prob_arr        <- array(0, dim = c(n_draws, length(ages_1315), 8))
    prob_arr[,, 1]  <- plogis(thresholds_arr[,, 1] - mu_arr[,, 1])
    for (k in 2:7) {
      prob_arr[,, k] <- plogis(thresholds_arr[,, k] - mu_arr[,, k]) -
                        plogis(thresholds_arr[,, k - 1] - mu_arr[,, k - 1])
    }
    prob_arr[,, 8] <- 1 - plogis(thresholds_arr[,, 7] - mu_arr[,, 7])

    prob_arr <- pmax(prob_arr, 1e-10)
    prob_sums <- apply(prob_arr, c(1, 2), sum)
    prob_sums[prob_sums <= 0] <- 1
    prob_arr  <- sweep(prob_arr, c(1, 2), prob_sums, "/")

    for (age_idx in seq_along(ages_1315)) {
      size_now <- age_ess_vec[age_idx]
      if (is.na(size_now) || size_now <= 0) next
      prob_slice <- prob_arr[, age_idx, , drop = FALSE]
      prob_mat   <- matrix(prob_slice, nrow = n_draws, ncol = 8)
      counts     <- sample_multinom_counts(rep(size_now, n_draws), prob_mat)
      predicted_counts_matrix <- predicted_counts_matrix + counts
    }

    pred_means <- apply(predicted_counts_matrix, 2, mean, na.rm = TRUE)
    pred_q025  <- apply(predicted_counts_matrix, 2, quantile, 0.025, na.rm = TRUE)
    pred_q975  <- apply(predicted_counts_matrix, 2, quantile, 0.975, na.rm = TRUE)
    pred_q25   <- apply(predicted_counts_matrix, 2, quantile, 0.25,  na.rm = TRUE)
    pred_q75   <- apply(predicted_counts_matrix, 2, quantile, 0.75,  na.rm = TRUE)

    actual_whoname <- country_names_meta$whoname[country_names_meta$iso3 == survey_info$iso3]
    if (length(actual_whoname) == 0) actual_whoname <- as.character(survey_info$iso3)

    ppc_entry <- data.frame(
      iso3     = rep(as.character(survey_info$iso3), 8),
      sex      = rep(sex_value, 8),
      sex_label = rep(sex_label, 8),
      whoname  = rep(actual_whoname, 8),
      day_num  = 0:7,
      observed_count = observed_counts,
      pred_mean  = pred_means,
      pred_q025  = pred_q025,
      pred_q975  = pred_q975,
      pred_q25   = pred_q25,
      pred_q75   = pred_q75,
      total_ess  = rep(total_ess, 8),
      survey     = rep(as.character(survey_info$survey), 8),
      stringsAsFactors = FALSE
    )
    ppc_results[[paste0(survey_info$iso3, "_", survey_info$survey)]] <- ppc_entry
  }

  if (length(ppc_results) > 0) {
    result <- do.call(rbind, ppc_results)
    cat("PPC generated:", nrow(result), "rows for", sex_label, "\n")
    return(result)
  } else {
    return(data.frame())
  }
}

ppc_data <- data.frame(
  iso3 = character(0), sex = numeric(0), sex_label = character(0),
  whoname = character(0), day_num = numeric(0), observed_count = numeric(0),
  pred_mean = numeric(0), pred_q025 = numeric(0), pred_q975 = numeric(0),
  pred_q25 = numeric(0), pred_q75 = numeric(0), total_ess = numeric(0),
  survey = character(0)
)

tryCatch({
  boys_ppc_results  <- generate_real_ppc(boys_results,  boys_data,  1, "Boys")
  girls_ppc_results <- generate_real_ppc(girls_results, girls_data, 2, "Girls")
  ppc_combined <- rbind(boys_ppc_results, girls_ppc_results)
  if (nrow(ppc_combined) > 0) {
    ppc_data <- ppc_combined
    cat("Total PPC rows:", nrow(ppc_data), "\n")
  }
}, error = function(e) {
  cat("Error generating PPC data:", as.character(e), "\n")
})

################################################################################
# SAVE OUTPUTS
################################################################################

cat("\n=== SAVING OUTPUTS ===\n")

# Prediction summary CSV
generate_csv_from_pred_data(pred_data, config$name, config$timestamp, "results")

# PPC CSV
if (nrow(ppc_data) > 0) {
  ppc_file <- paste0("results/ppc_", prefix, ".csv")
  write.csv(ppc_data, ppc_file, row.names = FALSE)
  cat("PPC data saved to:", ppc_file, "\n")
}

################################################################################
# MCMC DRAWS CSV
# Produces a wide CSV: one row per (iso3, year, age, sex), columns for each draw.
################################################################################

cat("\n=== GENERATING MCMC DRAWS CSV ===\n")

tryCatch({
  thin_rows <- 1

  # Extract in-sample trend probability arrays [n_draws, N_trends, 8]
  boys_trend_draws  <- rstan::extract(boys_results$fit,  "trend_category_probs")$trend_category_probs
  girls_trend_draws <- rstan::extract(girls_results$fit, "trend_category_probs")$trend_category_probs

  # Extract new-country probability arrays [n_draws, N_new, 8]
  boys_new_draws  <- rstan::extract(boys_results$fit,  "new_category_probs")$new_category_probs
  girls_new_draws <- rstan::extract(girls_results$fit, "new_category_probs")$new_category_probs

  n_draws_total <- nrow(boys_trend_draws)
  cat("Total draws available:", n_draws_total, "\n")

  thin_idx  <- seq(1, n_draws_total, by = thin_rows)
  cat("Using", length(thin_idx), "draws after thinning (thin_rows =", thin_rows, ")\n")

  boys_trend_thin  <- boys_trend_draws[thin_idx,, ]
  girls_trend_thin <- girls_trend_draws[thin_idx,, ]
  boys_new_thin    <- boys_new_draws[thin_idx,,  ]
  girls_new_thin   <- girls_new_draws[thin_idx,, ]

  # Reconstruct trend grid (same order as Stan)
  trend_countries <- boys_results$country_names
  trend_ages  <- seq(11, 17, by = 1)
  trend_years <- seq(2008, 2025, by = 1)
  trend_grid_draws <- expand.grid(
    iso3 = trend_countries, age = trend_ages, year = trend_years,
    stringsAsFactors = FALSE
  )

  # Reconstruct new-country grid
  covar_countries    <- sort(unique(as.character(covars$iso3)))
  new_countries_list <- setdiff(covar_countries, trend_countries)
  covar_region_map_draws <- covars %>%
    filter(iso3 %in% new_countries_list) %>%
    select(iso3, adopa_region) %>%
    distinct()

  new_grid_draws <- expand.grid(
    iso3 = new_countries_list, age = trend_ages, year = trend_years,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      adopa_region = covar_region_map_draws$adopa_region[match(iso3, covar_region_map_draws$iso3)]
    )

  cat("Trend grid:", nrow(trend_grid_draws), "rows |",
      "New-country grid:", nrow(new_grid_draws), "rows\n")

  boys_trend_df  <- make_prob_df(boys_trend_thin,  trend_grid_draws, 1)
  girls_trend_df <- make_prob_df(girls_trend_thin, trend_grid_draws, 2)
  boys_new_df    <- make_prob_df(boys_new_thin,    new_grid_draws,   1)
  girls_new_df   <- make_prob_df(girls_new_thin,   new_grid_draws,   2)

  combined_draws <- bind_rows(boys_trend_df, boys_new_df, girls_trend_df, girls_new_df) %>%
    arrange(iso3, year, age, sex)

  draws_file <- paste0("results/mcmc_draws_", prefix, "_all_countries.csv")
  write.csv(combined_draws, draws_file, row.names = FALSE)

  draws_size_mb <- file.info(draws_file)$size / 1024 / 1024
  cat("MCMC draws CSV saved to:", draws_file, "\n")
  cat("  Size:", round(draws_size_mb, 1), "MB |",
      "Rows:", nrow(combined_draws), "| Columns:", ncol(combined_draws), "\n")

}, error = function(e) {
  cat("Error generating MCMC draws CSV:", conditionMessage(e), "\n")
  cat("Prediction summaries are still available in results/predictions_*.csv\n")
})

################################################################################
# RUN SUMMARY
################################################################################

cat("\n=== RUN SUMMARY ===\n")
cat("Model:", config$name, "\n")
cat("Stan file:", config$stan_file, "\n")
cat("Chains:", config$chains, "| Iter:", config$iter, "| Warmup:", config$warmup,
    "| Thin:", config$thin, "\n")
cat("Adapt delta:", config$control$adapt_delta,
    "| Max treedepth:", config$control$max_treedepth, "\n")
cat("Timestamp:", config$timestamp, "\n")
cat("Elapsed time:", round(elapsed_time, 2), attr(elapsed_time, "units"), "\n")
cat("Output prefix:", prefix, "\n")
cat("Results saved to: results/\n")
cat("===================\n")
