
# Single Imputation with SuperLearner-Guided Truncated Normal Sampling
# Updated with prediction quality metrics output

library(SuperLearner)
library(dplyr)
library(truncnorm)

# Main function for single imputation with SuperLearner-guided truncated normal
single_imputation_sl_guided <- function(
  biomarker_df,        # Data with: subject_id, analyte_code, result, visit, analyte_lod, flg_lod, dataset
  cov_pred_df,         # Covariate data with subject_id and predictor variables
  lower_bounds = 0.0001,  # Lower bound for truncation
  sl_library = c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.mean"),
  seed = 123,
  log_transform = TRUE,
  weight_sl = 0.7,     # Weight for SuperLearner prediction (0-1)
  verbose = TRUE
) {

  set.seed(seed)

  # Load required packages
  required_packages <- c("SuperLearner", "dplyr", "truncnorm")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  # Input validation
  required_cols <- c("subject_id", "analyte_code", "result", "visit", "analyte_lod", "flg_lod")
  if (!all(required_cols %in% names(biomarker_df))) {
    stop(paste("biomarker_df must contain columns:", paste(required_cols, collapse = ", ")))
  }

  if (!"subject_id" %in% names(cov_pred_df)) {
    stop("cov_pred_df must contain subject_id column")
  }

  # Create analyte_visit combination for processing
  biomarker_df <- biomarker_df %>%
    mutate(analyte_visit = paste(analyte_code, visit, sep = "_"))

  # Get unique analyte codes for lower bounds
  unique_analytes <- unique(biomarker_df$analyte_code)

  # Handle lower_bounds
  if (length(lower_bounds) == 1) {
    lower_bounds <- setNames(rep(lower_bounds, length(unique_analytes)), unique_analytes)
  } else if (!all(unique_analytes %in% names(lower_bounds))) {
    stop("If lower_bounds is a vector, it must have names matching all analytes")
  }

  # Store original results and create working copy
  if (log_transform) {
    if (verbose) cat("Log-transforming biomarker values...\n")
    biomarker_df <- biomarker_df %>%
      mutate(
        result_original = result,
        result_log = ifelse(!is.na(result) & result > 0, log(result), result),
        analyte_lod_log = log(analyte_lod)
      )

    # Create log-transformed lower bounds
    log_lower_bounds <- log(lower_bounds)
  } else {
    biomarker_df <- biomarker_df %>%
      mutate(
        result_original = result,
        result_log = result,
        analyte_lod_log = analyte_lod
      )
    log_lower_bounds <- lower_bounds
  }

  # Get predictor variable names
  predictor_vars <- setdiff(names(cov_pred_df), "subject_id")

  if (verbose) {
    cat(paste("Number of predictor variables:", length(predictor_vars), "\n"))
    cat(paste("Number of analyte_code levels:", length(unique(biomarker_df$analyte_code)), "\n"))
    cat(paste("Number of visit levels:", length(unique(biomarker_df$visit)), "\n"))
    cat(paste("Number of analyte_code × visit combinations:", 
              length(unique(biomarker_df$analyte_visit)), "\n"))
    cat("\n=== Starting Single Imputation ===\n")
  }

  # Initialize output
  imputed_data <- biomarker_df
  imputation_details <- list()

  # Initialize quality metrics dataframe
  quality_metrics_list <- list()

  # Get unique analyte_visit combinations
  unique_analyte_visits <- unique(biomarker_df$analyte_visit)

  # Loop through each analyte_visit combination
  for (av_level in unique_analyte_visits) {

    if (verbose) cat(paste("\nProcessing:", av_level, "\n"))

    # Get analyte_code and visit
    parts <- strsplit(av_level, "_")[[1]]
    analyte_code <- paste(parts[-length(parts)], collapse = "_")  # Handle analyte codes with underscores
    current_visit <- parts[length(parts)]

    # Subset data for this analyte_visit
    av_indices <- which(imputed_data$analyte_visit == av_level)
    av_data <- imputed_data[av_indices, ]

    # Get LOD threshold and lower bound for this analyte
    lod_threshold <- unique(av_data$analyte_lod_log)[1]  # Should be same for all rows
    lower_bound <- log_lower_bounds[analyte_code]

    # Merge with covariates
    merged_data <- av_data %>%
      inner_join(cov_pred_df, by = "subject_id")

    # Skip if insufficient data
    if (nrow(merged_data) < 10) {
      if (verbose) cat("  Skipping - insufficient data (n =", nrow(merged_data), ")\n")

      # Record skipped in quality metrics
      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code,
        visit = current_visit,
        analyte_visit = av_level,
        n_total = nrow(merged_data),
        n_train = NA,
        n_imputed = NA,
        status = "insufficient_data",
        mse = NA,
        rmse = NA,
        mae = NA,
        r_squared = NA,
        adj_r_squared = NA,
        cv_risk = NA,
        best_algorithm = NA,
        n_predictors = NA,
        stringsAsFactors = FALSE
      )
      next
    }

    # Identify observations to impute (flg_lod == 1)
    below_lod_mask <- merged_data$flg_lod == 1

    if (sum(below_lod_mask) == 0) {
      if (verbose) cat("  No values below LOD to impute\n")

      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code,
        visit = current_visit,
        analyte_visit = av_level,
        n_total = nrow(merged_data),
        n_train = sum(!below_lod_mask & !is.na(merged_data$result_log)),
        n_imputed = 0,
        status = "no_imputation_needed",
        mse = NA,
        rmse = NA,
        mae = NA,
        r_squared = NA,
        adj_r_squared = NA,
        cv_risk = NA,
        best_algorithm = NA,
        n_predictors = NA,
        stringsAsFactors = FALSE
      )
      next
    }

    # Prepare training data (observed values only, flg_lod == 0)
    train_mask <- (merged_data$flg_lod == 0) & !is.na(merged_data$result_log)

    if (sum(train_mask) < 5) {
      if (verbose) cat("  Skipping - insufficient training data (n =", sum(train_mask), ")\n")

      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code,
        visit = current_visit,
        analyte_visit = av_level,
        n_total = nrow(merged_data),
        n_train = sum(train_mask),
        n_imputed = sum(below_lod_mask),
        status = "insufficient_training_data",
        mse = NA,
        rmse = NA,
        mae = NA,
        r_squared = NA,
        adj_r_squared = NA,
        cv_risk = NA,
        best_algorithm = NA,
        n_predictors = NA,
        stringsAsFactors = FALSE
      )
      next
    }

    # Prepare X (predictors) and Y (outcome)
    X_train <- merged_data[train_mask, predictor_vars, drop = FALSE]
    Y_train <- merged_data$result_log[train_mask]
    X_predict <- merged_data[below_lod_mask, predictor_vars, drop = FALSE]

    # Handle missing values in predictors
    X_train <- as.data.frame(X_train)
    X_predict <- as.data.frame(X_predict)

    # Simple imputation for predictors
    for (col in names(X_train)) {
      if (any(is.na(X_train[[col]]))) {
        if (is.numeric(X_train[[col]])) {
          med_val <- median(X_train[[col]], na.rm = TRUE)
          X_train[[col]][is.na(X_train[[col]])] <- ifelse(is.na(med_val), 0, med_val)
          X_predict[[col]][is.na(X_predict[[col]])] <- ifelse(is.na(med_val), 0, med_val)
        } else {
          X_train[[col]] <- as.numeric(as.factor(X_train[[col]]))
          X_predict[[col]] <- as.numeric(as.factor(X_predict[[col]]))
        }
      }
    }

    # Convert factors to numeric
    X_train <- X_train %>% 
      mutate(across(where(is.character), as.factor)) %>%
      mutate(across(where(is.factor), as.numeric))
    X_predict <- X_predict %>% 
      mutate(across(where(is.character), as.factor)) %>%
      mutate(across(where(is.factor), as.numeric))

    # Remove constant columns
    col_vars <- sapply(X_train, function(x) var(x, na.rm = TRUE))
    non_constant <- col_vars > 0 & !is.na(col_vars)

    if (sum(non_constant) == 0) {
      if (verbose) cat("  Skipping - no valid predictors\n")

      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code,
        visit = current_visit,
        analyte_visit = av_level,
        n_total = nrow(merged_data),
        n_train = sum(train_mask),
        n_imputed = sum(below_lod_mask),
        status = "no_valid_predictors",
        mse = NA,
        rmse = NA,
        mae = NA,
        r_squared = NA,
        adj_r_squared = NA,
        cv_risk = NA,
        best_algorithm = NA,
        n_predictors = 0,
        stringsAsFactors = FALSE
      )
      next
    }

    X_train <- X_train[, non_constant, drop = FALSE]
    X_predict <- X_predict[, non_constant, drop = FALSE]
    n_predictors_used <- ncol(X_train)

    # Step 1: Fit SuperLearner model and get predictions
    tryCatch({
      sl_fit <- SuperLearner(
        Y = Y_train,
        X = X_train,
        family = gaussian(),
        SL.library = sl_library,
        cvControl = list(V = min(5, sum(train_mask)))
      )

      # Get SuperLearner predictions for training data (for metrics)
      sl_pred_train <- predict(sl_fit, newdata = X_train, onlySL = TRUE)$pred

      # Calculate prediction quality metrics
      residuals <- Y_train - sl_pred_train
      mse <- mean(residuals^2)
      rmse <- sqrt(mse)
      mae <- mean(abs(residuals))

      # R-squared
      ss_tot <- sum((Y_train - mean(Y_train))^2)
      ss_res <- sum(residuals^2)
      r_squared <- 1 - (ss_res / ss_tot)

      # Adjusted R-squared
      n <- length(Y_train)
      p <- n_predictors_used
      adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

      # Cross-validated risk (from SuperLearner)
      cv_risk <- min(sl_fit$cvRisk)

      # Best performing algorithm
      best_algo <- names(sl_fit$coef)[which.max(sl_fit$coef)]

      # Get SuperLearner predictions for below-LOD observations
      sl_pred <- predict(sl_fit, newdata = X_predict, onlySL = TRUE)$pred

      # Cap predictions at LOD threshold and ensure above lower bound
      sl_pred_capped <- pmin(sl_pred, lod_threshold)
      sl_pred_capped <- pmax(sl_pred_capped, lower_bound)

      # Step 2: Calculate conventional truncated normal parameters
      conventional_mean <- lod_threshold / sqrt(2)
      conventional_sd <- sd(Y_train, na.rm = TRUE)

      # Step 3: Blend SuperLearner predictions with conventional approach
      guided_mean <- weight_sl * sl_pred_capped + (1 - weight_sl) * conventional_mean

      # Estimate residual SD from SuperLearner model
      residual_sd <- sd(residuals)

      # Use a weighted SD
      guided_sd <- weight_sl * residual_sd + (1 - weight_sl) * conventional_sd

      # Ensure SD is positive and reasonable
      if (is.na(guided_sd) || guided_sd <= 0) {
        guided_sd <- conventional_sd
        if (is.na(guided_sd) || guided_sd <= 0) {
          guided_sd <- (lod_threshold - lower_bound) / 4
        }
      }

      # Step 4: Generate truncated normal samples using guided parameters
      n_to_impute <- sum(below_lod_mask)
      imputed_values <- rtruncnorm(
        n = n_to_impute,
        a = lower_bound,
        b = lod_threshold,
        mean = guided_mean,
        sd = guided_sd
      )

      # Update imputed values in dataset (log scale)
      imputed_data$result_log[av_indices[below_lod_mask]] <- imputed_values

      # Store quality metrics
      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code,
        visit = current_visit,
        analyte_visit = av_level,
        n_total = nrow(merged_data),
        n_train = sum(train_mask),
        n_imputed = n_to_impute,
        pct_below_lod = round(100 * n_to_impute / nrow(merged_data), 2),
        status = "success",
        mse = round(mse, 6),
        rmse = round(rmse, 6),
        mae = round(mae, 6),
        r_squared = round(r_squared, 4),
        adj_r_squared = round(adj_r_squared, 4),
        cv_risk = round(cv_risk, 6),
        best_algorithm = best_algo,
        n_predictors = n_predictors_used,
        weight_sl_used = weight_sl,
        guided_mean_min = round(min(guided_mean), 4),
        guided_mean_max = round(max(guided_mean), 4),
        guided_sd = round(guided_sd, 4),
        conventional_mean = round(conventional_mean, 4),
        stringsAsFactors = FALSE
      )

      # Store detailed imputation info
      imputation_details[[av_level]] <- list(
        analyte_code = analyte_code,
        visit = current_visit,
        n_train = sum(train_mask),
        n_imputed = n_to_impute,
        sl_pred_range = range(sl_pred),
        sl_pred_capped_range = range(sl_pred_capped),
        guided_mean_range = range(guided_mean),
        guided_sd = guided_sd,
        conventional_mean = conventional_mean,
        conventional_sd = conventional_sd,
        lod_threshold = lod_threshold,
        lower_bound = lower_bound,
        sl_weights = sl_fit$coef,
        quality_metrics = list(
          mse = mse,
          rmse = rmse,
          mae = mae,
          r_squared = r_squared,
          adj_r_squared = adj_r_squared,
          cv_risk = cv_risk
        )
      )

      if (verbose) {
        cat(paste("  Analyte:", analyte_code, "| Visit:", current_visit, "\n"))
        cat(paste("  Training samples:", sum(train_mask), "\n"))
        cat(paste("  Imputed samples:", n_to_impute, "\n"))
        cat(paste("  R-squared:", round(r_squared, 3), "\n"))
        cat(paste("  RMSE:", round(rmse, 4), "\n"))
        cat(paste("  Best algorithm:", best_algo, "\n"))
      }

    }, error = function(e) {
      if (verbose) cat(paste("  Error:", e$message, "\n"))

      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code,
        visit = current_visit,
        analyte_visit = av_level,
        n_total = nrow(merged_data),
        n_train = sum(train_mask),
        n_imputed = sum(below_lod_mask),
        status = paste("error:", substr(e$message, 1, 50)),
        mse = NA,
        rmse = NA,
        mae = NA,
        r_squared = NA,
        adj_r_squared = NA,
        cv_risk = NA,
        best_algorithm = NA,
        n_predictors = n_predictors_used,
        stringsAsFactors = FALSE
      )
    })
  }

  # Combine quality metrics into a dataframe
  quality_metrics_df <- bind_rows(quality_metrics_list)

  # Back-transform if needed and prepare final output
  if (log_transform) {
    imputed_data <- imputed_data %>%
      mutate(
        result_imputed = ifelse(flg_lod == 1, exp(result_log), NA),
        result = ifelse(flg_lod == 1, exp(result_log), result_original)
      ) %>%
      select(-result_original, -result_log, -analyte_lod_log, -analyte_visit)
  } else {
    imputed_data <- imputed_data %>%
      mutate(
        result_imputed = ifelse(flg_lod == 1, result_log, NA),
        result = ifelse(flg_lod == 1, result_log, result_original)
      ) %>%
      select(-result_original, -result_log, -analyte_lod_log, -analyte_visit)
  }

  if (verbose) {
    cat("\n=== Single Imputation Completed ===\n")
    cat(paste("Total combinations processed:", nrow(quality_metrics_df), "\n"))
    cat(paste("Successful imputations:", 
              sum(quality_metrics_df$status == "success", na.rm = TRUE), "\n"))
    cat("\nQuality Metrics Summary:\n")
    successful <- quality_metrics_df %>% filter(status == "success")
    if (nrow(successful) > 0) {
      cat(paste("  Mean R-squared:", round(mean(successful$r_squared, na.rm = TRUE), 3), "\n"))
      cat(paste("  Mean RMSE:", round(mean(successful$rmse, na.rm = TRUE), 4), "\n"))
      cat(paste("  Mean MAE:", round(mean(successful$mae, na.rm = TRUE), 4), "\n"))
    }
  }

  # Return results
  result <- list(
    imputed_data = imputed_data,
    quality_metrics = quality_metrics_df,
    imputation_details = imputation_details,
    parameters = list(
      lower_bounds = lower_bounds,
      sl_library = sl_library,
      seed = seed,
      log_transform = log_transform,
      weight_sl = weight_sl
    )
  )

  class(result) <- "single_imputation_sl"
  return(result)
}


# Alternative: Adaptive weight based on SuperLearner prediction quality
single_imputation_sl_adaptive <- function(
  biomarker_df,
  cov_pred_df,
  lower_bounds = 0.0001,
  sl_library = c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.mean"),
  seed = 123,
  log_transform = TRUE,
  verbose = TRUE
) {

  set.seed(seed)

  # Load required packages
  required_packages <- c("SuperLearner", "dplyr", "truncnorm")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  # Input validation
  required_cols <- c("subject_id", "analyte_code", "result", "visit", "analyte_lod", "flg_lod")
  if (!all(required_cols %in% names(biomarker_df))) {
    stop(paste("biomarker_df must contain columns:", paste(required_cols, collapse = ", ")))
  }

  # Create analyte_visit combination
  biomarker_df <- biomarker_df %>%
    mutate(analyte_visit = paste(analyte_code, visit, sep = "_"))

  unique_analytes <- unique(biomarker_df$analyte_code)

  # Handle lower_bounds
  if (length(lower_bounds) == 1) {
    lower_bounds <- setNames(rep(lower_bounds, length(unique_analytes)), unique_analytes)
  }

  # Transform if needed
  if (log_transform) {
    biomarker_df <- biomarker_df %>%
      mutate(
        result_original = result,
        result_log = ifelse(!is.na(result) & result > 0, log(result), result),
        analyte_lod_log = log(analyte_lod)
      )
    log_lower_bounds <- log(lower_bounds)
  } else {
    biomarker_df <- biomarker_df %>%
      mutate(
        result_original = result,
        result_log = result,
        analyte_lod_log = analyte_lod
      )
    log_lower_bounds <- lower_bounds
  }

  predictor_vars <- setdiff(names(cov_pred_df), "subject_id")
  imputed_data <- biomarker_df
  imputation_details <- list()
  quality_metrics_list <- list()

  if (verbose) cat("\n=== Starting Adaptive Single Imputation ===\n")

  unique_analyte_visits <- unique(biomarker_df$analyte_visit)

  for (av_level in unique_analyte_visits) {

    if (verbose) cat(paste("\nProcessing:", av_level, "\n"))

    parts <- strsplit(av_level, "_")[[1]]
    analyte_code <- paste(parts[-length(parts)], collapse = "_")
    current_visit <- parts[length(parts)]

    av_indices <- which(imputed_data$analyte_visit == av_level)
    av_data <- imputed_data[av_indices, ]

    lod_threshold <- unique(av_data$analyte_lod_log)[1]
    lower_bound <- log_lower_bounds[analyte_code]

    merged_data <- av_data %>% inner_join(cov_pred_df, by = "subject_id")

    if (nrow(merged_data) < 10) {
      if (verbose) cat("  Skipping - insufficient data\n")
      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code, visit = current_visit, analyte_visit = av_level,
        n_total = nrow(merged_data), n_train = NA, n_imputed = NA,
        status = "insufficient_data", mse = NA, rmse = NA, mae = NA,
        r_squared = NA, adj_r_squared = NA, cv_risk = NA, best_algorithm = NA,
        n_predictors = NA, stringsAsFactors = FALSE
      )
      next
    }

    below_lod_mask <- merged_data$flg_lod == 1
    if (sum(below_lod_mask) == 0) {
      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code, visit = current_visit, analyte_visit = av_level,
        n_total = nrow(merged_data), n_train = sum(!below_lod_mask), n_imputed = 0,
        status = "no_imputation_needed", mse = NA, rmse = NA, mae = NA,
        r_squared = NA, adj_r_squared = NA, cv_risk = NA, best_algorithm = NA,
        n_predictors = NA, stringsAsFactors = FALSE
      )
      next
    }

    train_mask <- (merged_data$flg_lod == 0) & !is.na(merged_data$result_log)
    if (sum(train_mask) < 5) {
      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code, visit = current_visit, analyte_visit = av_level,
        n_total = nrow(merged_data), n_train = sum(train_mask), 
        n_imputed = sum(below_lod_mask), status = "insufficient_training_data",
        mse = NA, rmse = NA, mae = NA, r_squared = NA, adj_r_squared = NA,
        cv_risk = NA, best_algorithm = NA, n_predictors = NA, stringsAsFactors = FALSE
      )
      next
    }

    X_train <- as.data.frame(merged_data[train_mask, predictor_vars, drop = FALSE])
    Y_train <- merged_data$result_log[train_mask]
    X_predict <- as.data.frame(merged_data[below_lod_mask, predictor_vars, drop = FALSE])

    # Handle missing values
    for (col in names(X_train)) {
      if (any(is.na(X_train[[col]]))) {
        if (is.numeric(X_train[[col]])) {
          med_val <- median(X_train[[col]], na.rm = TRUE)
          X_train[[col]][is.na(X_train[[col]])] <- ifelse(is.na(med_val), 0, med_val)
          X_predict[[col]][is.na(X_predict[[col]])] <- ifelse(is.na(med_val), 0, med_val)
        } else {
          X_train[[col]] <- as.numeric(as.factor(X_train[[col]]))
          X_predict[[col]] <- as.numeric(as.factor(X_predict[[col]]))
        }
      }
    }

    X_train <- X_train %>% 
      mutate(across(where(is.character), as.factor)) %>%
      mutate(across(where(is.factor), as.numeric))
    X_predict <- X_predict %>% 
      mutate(across(where(is.character), as.factor)) %>%
      mutate(across(where(is.factor), as.numeric))

    col_vars <- sapply(X_train, function(x) var(x, na.rm = TRUE))
    non_constant <- col_vars > 0 & !is.na(col_vars)
    if (sum(non_constant) == 0) {
      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code, visit = current_visit, analyte_visit = av_level,
        n_total = nrow(merged_data), n_train = sum(train_mask), 
        n_imputed = sum(below_lod_mask), status = "no_valid_predictors",
        mse = NA, rmse = NA, mae = NA, r_squared = NA, adj_r_squared = NA,
        cv_risk = NA, best_algorithm = NA, n_predictors = 0, stringsAsFactors = FALSE
      )
      next
    }

    X_train <- X_train[, non_constant, drop = FALSE]
    X_predict <- X_predict[, non_constant, drop = FALSE]
    n_predictors_used <- ncol(X_train)

    tryCatch({
      # Fit SuperLearner
      sl_fit <- SuperLearner(
        Y = Y_train,
        X = X_train,
        family = gaussian(),
        SL.library = sl_library,
        cvControl = list(V = min(5, sum(train_mask)))
      )

      sl_pred <- predict(sl_fit, newdata = X_predict, onlySL = TRUE)$pred
      sl_pred_train <- predict(sl_fit, newdata = X_train, onlySL = TRUE)$pred

      # Calculate metrics
      residuals <- Y_train - sl_pred_train
      mse <- mean(residuals^2)
      rmse <- sqrt(mse)
      mae <- mean(abs(residuals))

      ss_tot <- sum((Y_train - mean(Y_train))^2)
      ss_res <- sum(residuals^2)
      r_squared <- 1 - (ss_res / ss_tot)

      n <- length(Y_train)
      p <- n_predictors_used
      adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

      cv_risk <- min(sl_fit$cvRisk)
      best_algo <- names(sl_fit$coef)[which.max(sl_fit$coef)]

      # Adaptive weight based on R-squared
      weight_sl <- max(0, min(1, r_squared))^0.5

      # Cap predictions
      sl_pred_capped <- pmin(sl_pred, lod_threshold)
      sl_pred_capped <- pmax(sl_pred_capped, lower_bound)

      # Conventional parameters
      conventional_mean <- lod_threshold / sqrt(2)
      conventional_sd <- sd(Y_train, na.rm = TRUE)

      # Adaptive blending
      guided_mean <- weight_sl * sl_pred_capped + (1 - weight_sl) * conventional_mean

      residual_sd <- sd(residuals)
      guided_sd <- weight_sl * residual_sd + (1 - weight_sl) * conventional_sd

      if (is.na(guided_sd) || guided_sd <= 0) {
        guided_sd <- conventional_sd
        if (is.na(guided_sd) || guided_sd <= 0) {
          guided_sd <- (lod_threshold - lower_bound) / 4
        }
      }

      # Impute
      n_to_impute <- sum(below_lod_mask)
      imputed_values <- rtruncnorm(
        n = n_to_impute,
        a = lower_bound,
        b = lod_threshold,
        mean = guided_mean,
        sd = guided_sd
      )

      imputed_data$result_log[av_indices[below_lod_mask]] <- imputed_values

      # Store quality metrics
      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code,
        visit = current_visit,
        analyte_visit = av_level,
        n_total = nrow(merged_data),
        n_train = sum(train_mask),
        n_imputed = n_to_impute,
        pct_below_lod = round(100 * n_to_impute / nrow(merged_data), 2),
        status = "success",
        mse = round(mse, 6),
        rmse = round(rmse, 6),
        mae = round(mae, 6),
        r_squared = round(r_squared, 4),
        adj_r_squared = round(adj_r_squared, 4),
        cv_risk = round(cv_risk, 6),
        best_algorithm = best_algo,
        n_predictors = n_predictors_used,
        adaptive_weight = round(weight_sl, 4),
        guided_mean_min = round(min(guided_mean), 4),
        guided_mean_max = round(max(guided_mean), 4),
        guided_sd = round(guided_sd, 4),
        conventional_mean = round(conventional_mean, 4),
        stringsAsFactors = FALSE
      )

      imputation_details[[av_level]] <- list(
        analyte_code = analyte_code,
        visit = current_visit,
        n_train = sum(train_mask),
        n_imputed = n_to_impute,
        r_squared = r_squared,
        adaptive_weight = weight_sl,
        guided_mean_range = range(guided_mean),
        guided_sd = guided_sd,
        sl_weights = sl_fit$coef,
        quality_metrics = list(
          mse = mse, rmse = rmse, mae = mae,
          r_squared = r_squared, adj_r_squared = adj_r_squared, cv_risk = cv_risk
        )
      )

      if (verbose) {
        cat(paste("  Analyte:", analyte_code, "| Visit:", current_visit, "\n"))
        cat(paste("  Training samples:", sum(train_mask), "\n"))
        cat(paste("  Imputed samples:", n_to_impute, "\n"))
        cat(paste("  R-squared:", round(r_squared, 3), "\n"))
        cat(paste("  Adaptive weight:", round(weight_sl, 3), "\n"))
        cat(paste("  Best algorithm:", best_algo, "\n"))
      }

    }, error = function(e) {
      if (verbose) cat(paste("  Error:", e$message, "\n"))
      quality_metrics_list[[av_level]] <- data.frame(
        analyte_code = analyte_code, visit = current_visit, analyte_visit = av_level,
        n_total = nrow(merged_data), n_train = sum(train_mask), 
        n_imputed = sum(below_lod_mask), 
        status = paste("error:", substr(e$message, 1, 50)),
        mse = NA, rmse = NA, mae = NA, r_squared = NA, adj_r_squared = NA,
        cv_risk = NA, best_algorithm = NA, n_predictors = n_predictors_used,
        stringsAsFactors = FALSE
      )
    })
  }

  # Combine quality metrics
  quality_metrics_df <- bind_rows(quality_metrics_list)

  # Back-transform
  if (log_transform) {
    imputed_data <- imputed_data %>%
      mutate(
        result_imputed = ifelse(flg_lod == 1, exp(result_log), NA),
        result = ifelse(flg_lod == 1, exp(result_log), result_original)
      ) %>%
      select(-result_original, -result_log, -analyte_lod_log, -analyte_visit)
  } else {
    imputed_data <- imputed_data %>%
      mutate(
        result_imputed = ifelse(flg_lod == 1, result_log, NA),
        result = ifelse(flg_lod == 1, result_log, result_original)
      ) %>%
      select(-result_original, -result_log, -analyte_lod_log, -analyte_visit)
  }

  if (verbose) {
    cat("\n=== Adaptive Single Imputation Completed ===\n")
    cat(paste("Total combinations processed:", nrow(quality_metrics_df), "\n"))
    successful <- quality_metrics_df %>% filter(status == "success")
    if (nrow(successful) > 0) {
      cat(paste("  Mean R-squared:", round(mean(successful$r_squared, na.rm = TRUE), 3), "\n"))
      cat(paste("  Mean Adaptive Weight:", round(mean(successful$adaptive_weight, na.rm = TRUE), 3), "\n"))
    }
  }

  result <- list(
    imputed_data = imputed_data,
    quality_metrics = quality_metrics_df,
    imputation_details = imputation_details,
    parameters = list(
      lower_bounds = lower_bounds,
      sl_library = sl_library,
      seed = seed,
      log_transform = log_transform,
      method = "adaptive"
    )
  )

  class(result) <- "single_imputation_sl_adaptive"
  return(result)
}


# Example usage:
# 
# result <- single_imputation_sl_guided(
#   biomarker_df = all_exp,
#   cov_pred_df = cov_pred,
#   weight_sl = 0.7,
#   verbose = TRUE
# )
# 
# # Extract results
# imputed_data <- result$imputed_data
# quality_metrics <- result$quality_metrics
# 
# # View quality metrics
# View(quality_metrics)
# 
# # Filter successful imputations
# successful_metrics <- quality_metrics %>% 
#   filter(status == "success") %>%
#   arrange(desc(r_squared))
# 
# # Summary statistics
# summary(successful_metrics[, c("r_squared", "rmse", "mae", "pct_below_lod")])
