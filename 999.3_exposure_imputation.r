# =============================================================================
# Truncated Single Imputation with SuperLearner — Robust + Full Diagnostics
# =============================================================================
# Every observation with flg_lod == 1 receives an imputed value.
# When SuperLearner cannot be fitted (insufficient data, insufficient training
# observations, or no valid predictors), a conventional fallback is used:
# values are drawn from TruncNormal(mean = LOD/sqrt(2), sd = fallback_sd,
# a = lower_bound, b = LOD).
#
# Required columns in biomarker_df:
#   subject_id, result, flg_lod, analyte_lod, analyte_visit
# =============================================================================

library(SuperLearner)
library(dplyr)
library(truncnorm)


# =============================================================================
# INTERNAL HELPER — conventional truncated normal fallback imputation
# Draws n values from TruncNormal bounded by [lower_bound, lod_threshold].
# mean  = lod_threshold / sqrt(2)   (standard LOD/sqrt(2) convention)
# sd    = (lod_threshold - lower_bound) / 4  if no training data available,
#         otherwise sd of observed training values
# =============================================================================
.conventional_impute <- function(n, lower_bound, lod_threshold,
                                 Y_train = NULL) {
  conv_mean <- lod_threshold / sqrt(2)
  conv_sd   <- if (!is.null(Y_train) && length(Y_train) >= 2) {
    sd(Y_train, na.rm = TRUE)
  } else {
    (lod_threshold - lower_bound) / 4
  }
  if (is.na(conv_sd) || conv_sd <= 0) conv_sd <- (lod_threshold - lower_bound) / 4
  if (is.na(conv_sd) || conv_sd <= 0) conv_sd <- 1e-6
  
  rtruncnorm(n    = n,
             a    = lower_bound,
             b    = lod_threshold,
             mean = conv_mean,
             sd   = conv_sd)
}


# =============================================================================
# INTERNAL HELPER — skeleton metric row for any analyte_visit outcome
# =============================================================================
.make_metric_row <- function(av_level, status,
                             n_total = NA, n_train = NA, n_imputed = NA,
                             n_predictors = NA,
                             lod_upper = NA, lower_bound_val = NA) {
  analyte_code  <- sub("_[^_]+$", "", av_level)
  current_visit <- sub(".*_",      "", av_level)
  data.frame(
    analyte_code      = analyte_code,
    visit             = current_visit,
    analyte_visit     = av_level,
    lod_upper         = lod_upper,
    lower_bound       = lower_bound_val,
    n_total           = n_total,
    n_train           = n_train,
    n_imputed         = n_imputed,
    pct_below_lod     = NA_real_,
    status            = status,
    imputation_method = NA_character_,
    mse               = NA_real_,
    rmse              = NA_real_,
    mae               = NA_real_,
    r_squared         = NA_real_,
    adj_r_squared     = NA_real_,
    cv_risk           = NA_real_,
    best_algorithm    = NA_character_,
    n_predictors      = n_predictors,
    adaptive_weight   = NA_real_,
    guided_mean_min   = NA_real_,
    guided_mean_max   = NA_real_,
    guided_sd         = NA_real_,
    conventional_mean = NA_real_,
    conventional_sd   = NA_real_,
    imp_mean_log      = NA_real_,
    imp_sd_log        = NA_real_,
    stringsAsFactors  = FALSE
  )
}


# =============================================================================
# MAIN FUNCTION
# =============================================================================

truncated_single_imputation_sl <- function(
    biomarker_df,              # Tibble/df: subject_id, result, flg_lod,
    #            analyte_lod, analyte_visit
    cov_pred_df,               # Covariates: subject_id + predictor columns
    lower_bounds    = 0.0001,  # Lower truncation bound — scalar or named vector
    #   by analyte_code
    maxit_truncated = 100,     # Number of Gibbs sampling iterations
    sl_library      = c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.mean"),
    seed            = 123,
    log_transform   = TRUE,    # Log-transform result and bounds before modelling
    adaptive_weight = TRUE,    # TRUE  -> w = sqrt(max(0, R^2)) per analyte_visit
    # FALSE -> use fixed weight_sl for all
    weight_sl       = 0.7,     # Fixed SL weight (only when adaptive_weight=FALSE)
    track_convergence = TRUE,  # Record imp mean/SD per Gibbs iter
    verbose         = TRUE
) {
  
  # ---------------------------------------------------------------------------
  # 0. Package checks
  # ---------------------------------------------------------------------------
  for (pkg in c("SuperLearner", "dplyr", "truncnorm")) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE))
      stop(paste("Package", pkg, "is required but not installed."))
  }
  
  # ---------------------------------------------------------------------------
  # 1. Input validation
  # ---------------------------------------------------------------------------
  required_bio_cols <- c("subject_id", "result", "flg_lod",
                         "analyte_lod", "analyte_visit")
  missing_cols <- setdiff(required_bio_cols, names(biomarker_df))
  if (length(missing_cols) > 0)
    stop(paste("biomarker_df is missing columns:",
               paste(missing_cols, collapse = ", ")))
  
  if (!"subject_id" %in% names(cov_pred_df))
    stop("cov_pred_df must contain a 'subject_id' column.")
  
  if (!is.numeric(weight_sl) || weight_sl < 0 || weight_sl > 1)
    stop("weight_sl must be numeric between 0 and 1.")
  
  if (maxit_truncated < 1)
    stop("maxit_truncated must be >= 1.")
  
  if (any(is.na(biomarker_df$analyte_lod)))
    stop("analyte_lod contains NA values. All rows must have a valid LOD.")
  
  if (any(biomarker_df$analyte_lod <= 0))
    stop("analyte_lod must be > 0 for all rows.")
  
  # ---------------------------------------------------------------------------
  # 2. Parse analyte_code and visit; expand lower_bounds
  # ---------------------------------------------------------------------------
  biomarker_df <- biomarker_df %>%
    mutate(
      analyte_code = sub("_[^_]+$", "", analyte_visit),
      visit        = sub(".*_",      "", analyte_visit)
    )
  
  unique_analytes       <- unique(biomarker_df$analyte_code)
  unique_analyte_visits <- unique(biomarker_df$analyte_visit)
  
  if (length(lower_bounds) == 1) {
    lower_bounds <- setNames(rep(lower_bounds, length(unique_analytes)),
                             unique_analytes)
  } else {
    missing_lb <- setdiff(unique_analytes, names(lower_bounds))
    if (length(missing_lb) > 0)
      stop(paste("lower_bounds missing entries for:",
                 paste(missing_lb, collapse = ", ")))
  }
  
  # Sanity check: lower_bound < analyte_lod for every row
  bound_violation <- biomarker_df %>%
    mutate(lb = lower_bounds[analyte_code]) %>%
    filter(lb >= analyte_lod) %>%
    distinct(analyte_code, analyte_lod) %>%
    mutate(lower_bound_val = lower_bounds[analyte_code])
  
  if (nrow(bound_violation) > 0)
    stop(paste0(
      "lower_bounds must be strictly < analyte_lod. Violations:\n",
      paste(sprintf("  %s: lower=%.6f, lod=%.6f",
                    bound_violation$analyte_code,
                    bound_violation$lower_bound_val,
                    bound_violation$analyte_lod),
            collapse = "\n")
    ))
  
  # ---------------------------------------------------------------------------
  # 3. Log-transform
  # ---------------------------------------------------------------------------
  if (log_transform) {
    if (verbose) cat("Log-transforming result, analyte_lod, and lower_bounds...\n")
    biomarker_df <- biomarker_df %>%
      mutate(
        result_original = result,
        result_log      = ifelse(!is.na(result) & result > 0,
                                 log(result), NA_real_),
        analyte_lod_log = log(analyte_lod)
      )
    log_lower_bounds <- log(lower_bounds)
  } else {
    biomarker_df <- biomarker_df %>%
      mutate(result_original = result,
             result_log      = result,
             analyte_lod_log = analyte_lod)
    log_lower_bounds <- lower_bounds
  }
  
  # ---------------------------------------------------------------------------
  # 4. Setup
  # ---------------------------------------------------------------------------
  predictor_vars <- setdiff(names(cov_pred_df), "subject_id")
  
  # Count total flg_lod == 1 observations to verify at the end
  total_to_impute <- sum(biomarker_df$flg_lod == 1, na.rm = TRUE)
  
  if (verbose) {
    cat(paste("Predictor variables        :", length(predictor_vars),             "\n"))
    cat(paste("Unique analyte codes       :", length(unique_analytes),            "\n"))
    cat(paste("Analyte x visit combos     :", length(unique_analyte_visits),      "\n"))
    cat(paste("Total flg_lod=1 to impute  :", total_to_impute,                   "\n"))
    cat(paste("Gibbs iterations           :", maxit_truncated,                    "\n"))
    cat(paste("Adaptive SL weighting      :", adaptive_weight,                    "\n"))
    cat("\n=== Starting Truncated Single Imputation ===\n")
  }
  
  set.seed(seed)
  
  current_data         <- biomarker_df
  quality_metrics_list <- list()
  imputation_details   <- list()
  convergence_tracker  <- list()
  
  for (av_level in unique_analyte_visits) {
    quality_metrics_list[[av_level]] <- .make_metric_row(
      av_level, status = "not_yet_processed"
    )
  }
  
  # ---------------------------------------------------------------------------
  # 5. Gibbs sampling loop
  # ---------------------------------------------------------------------------
  for (iter in seq_len(maxit_truncated)) {
    
    if (verbose) cat(paste0("\n--- Gibbs iteration ", iter,
                            " / ", maxit_truncated, " ---\n"))
    
    iter_convergence <- list()
    
    for (av_level in unique_analyte_visits) {
      
      if (verbose) cat(paste("  Processing:", av_level, "\n"))
      
      analyte_code  <- sub("_[^_]+$", "", av_level)
      current_visit <- sub(".*_",      "", av_level)
      
      # Get bounds for this analyte_visit
      lod_threshold <- current_data$analyte_lod_log[
        current_data$analyte_visit == av_level & 
          !is.na(current_data$analyte_lod_log)
      ][1]
      lower_bound <- log_lower_bounds[analyte_code]
      
      # Original-scale bounds for reporting
      lod_orig <- if (log_transform) exp(lod_threshold) else lod_threshold
      lb_orig  <- if (log_transform) exp(lower_bound)   else lower_bound
      
      av_data     <- current_data %>% filter(analyte_visit == av_level)
      merged_data <- av_data %>% inner_join(cov_pred_df, by = "subject_id")
      
      below_lod_mask <- merged_data$flg_lod == 1
      n_to_impute    <- sum(below_lod_mask)
      
      # If nothing to impute in this analyte_visit, record and move on
      if (n_to_impute == 0) {
        if (verbose) cat("    No below-LOD observations — nothing to impute\n")
        if (iter == maxit_truncated) {
          quality_metrics_list[[av_level]] <- .make_metric_row(
            av_level, "no_imputation_needed",
            n_total          = nrow(merged_data),
            n_train          = sum(!below_lod_mask & !is.na(merged_data$result_log)),
            n_imputed        = 0L,
            lod_upper        = lod_orig,
            lower_bound_val  = lb_orig
          )
        }
        next
      }
      
      train_mask <- (merged_data$flg_lod == 0) & !is.na(merged_data$result_log)
      
      # -----------------------------------------------------------------------
      # Determine imputation path:
      #   "superlearner" — enough data to fit SL
      #   "conventional" — fall back to LOD/sqrt(2) truncated normal
      #
      # Fallback conditions:
      #   (a) fewer than 10 total observations after covariate merge
      #   (b) fewer than 5 above-LOD training observations
      #   (c) no non-constant predictors (detected later inside tryCatch)
      # -----------------------------------------------------------------------
      
      use_fallback    <- FALSE
      fallback_reason <- NULL
      
      if (nrow(merged_data) < 10) {
        use_fallback    <- TRUE
        fallback_reason <- paste0("insufficient_total_data (n=", nrow(merged_data), ")")
        if (verbose) cat("    FALLBACK: insufficient total data (n =",
                         nrow(merged_data), ") — using conventional imputation\n")
      } else if (sum(train_mask) < 5) {
        use_fallback    <- TRUE
        fallback_reason <- paste0("insufficient_training_data (n_train=",
                                  sum(train_mask), ")")
        if (verbose) cat("    FALLBACK: insufficient training data (n =",
                         sum(train_mask), ") — using conventional imputation\n")
      }
      
      # -----------------------------------------------------------------------
      # FALLBACK PATH — conventional truncated normal
      # -----------------------------------------------------------------------
      if (use_fallback) {
        
        Y_train_avail  <- if (sum(train_mask) >= 2) merged_data$result_log[train_mask] else NULL
        imputed_values <- .conventional_impute(
          n             = n_to_impute,
          lower_bound   = lower_bound,
          lod_threshold = lod_threshold,
          Y_train       = Y_train_avail
        )
        
        below_lod_subjects <- merged_data$subject_id[below_lod_mask]
        update_indices     <- which(
          current_data$analyte_visit == av_level &
            current_data$subject_id   %in% below_lod_subjects
        )
        
        if (length(update_indices) == length(imputed_values)) {
          current_data$result_log[update_indices] <- imputed_values
        } else {
          warning(paste0("[", av_level, "] Fallback assignment mismatch (",
                         length(update_indices), " vs ",
                         length(imputed_values), ")."))
        }
        
        if (track_convergence) {
          iter_convergence[[av_level]] <- data.frame(
            iter              = iter,
            analyte_visit     = av_level,
            imp_mean          = mean(imputed_values),
            imp_sd            = sd(imputed_values),
            imp_min           = min(imputed_values),
            imp_max           = max(imputed_values),
            r_squared         = NA_real_,
            adaptive_weight   = NA_real_,
            imputation_method = "conventional",
            stringsAsFactors  = FALSE
          )
        }
        
        if (iter == maxit_truncated) {
          conv_mean_val <- lod_threshold / sqrt(2)
          quality_metrics_list[[av_level]] <- data.frame(
            analyte_code      = analyte_code,
            visit             = current_visit,
            analyte_visit     = av_level,
            lod_upper         = lod_orig,
            lower_bound       = lb_orig,
            n_total           = nrow(merged_data),
            n_train           = sum(train_mask),
            n_imputed         = n_to_impute,
            pct_below_lod     = round(100 * n_to_impute / nrow(merged_data), 2),
            status            = paste0("conventional_fallback: ", fallback_reason),
            imputation_method = "conventional",
            mse               = NA_real_,
            rmse              = NA_real_,
            mae               = NA_real_,
            r_squared         = NA_real_,
            adj_r_squared     = NA_real_,
            cv_risk           = NA_real_,
            best_algorithm    = NA_character_,
            n_predictors      = NA_real_,
            adaptive_weight   = NA_real_,
            guided_mean_min   = round(conv_mean_val, 4),
            guided_mean_max   = round(conv_mean_val, 4),
            guided_sd         = NA_real_,
            conventional_mean = round(conv_mean_val, 4),
            conventional_sd   = NA_real_,
            imp_mean_log      = round(mean(imputed_values), 4),
            imp_sd_log        = round(sd(imputed_values),   4),
            stringsAsFactors  = FALSE
          )
          
          imputation_details[[av_level]] <- list(
            analyte_code      = analyte_code,
            visit             = current_visit,
            imputation_method = "conventional",
            fallback_reason   = fallback_reason,
            n_train           = sum(train_mask),
            n_imputed         = n_to_impute,
            lod_threshold_log = lod_threshold,
            lower_bound_log   = lower_bound,
            imputed_range_log = range(imputed_values)
          )
          
          if (verbose) cat(paste0(
            "    [fallback] n_imp=", n_to_impute,
            " | mean_imp=", round(mean(imputed_values), 3),
            " | reason=", fallback_reason, "\n"
          ))
        }
        
        next   # move to next av_level — fallback imputation done
      }
      
      # -----------------------------------------------------------------------
      # SUPERLEARNER PATH
      # -----------------------------------------------------------------------
      X_train   <- as.data.frame(merged_data[train_mask,     predictor_vars, drop = FALSE])
      X_predict <- as.data.frame(merged_data[below_lod_mask, predictor_vars, drop = FALSE])
      Y_train   <- merged_data$result_log[train_mask]
      
      # Predictor cleaning
      for (col in names(X_train)) {
        if (is.numeric(X_train[[col]])) {
          med_val <- median(X_train[[col]], na.rm = TRUE)
          fill    <- ifelse(is.na(med_val), 0, med_val)
          X_train[[col]][is.na(X_train[[col]])]     <- fill
          X_predict[[col]][is.na(X_predict[[col]])] <- fill
        } else {
          lvls <- union(levels(as.factor(X_train[[col]])),
                        levels(as.factor(X_predict[[col]])))
          X_train[[col]]   <- as.integer(factor(X_train[[col]],   levels = lvls))
          X_predict[[col]] <- as.integer(factor(X_predict[[col]], levels = lvls))
          X_train[[col]][is.na(X_train[[col]])]     <- 0L
          X_predict[[col]][is.na(X_predict[[col]])] <- 0L
        }
      }
      
      X_train <- X_train %>%
        mutate(across(where(is.character), as.factor)) %>%
        mutate(across(where(is.factor),    as.numeric))
      X_predict <- X_predict %>%
        mutate(across(where(is.character), as.factor)) %>%
        mutate(across(where(is.factor),    as.numeric))
      
      col_vars     <- sapply(X_train, function(x) var(x, na.rm = TRUE))
      non_constant <- !is.na(col_vars) & col_vars > 0
      
      # If no non-constant predictors, fall back to conventional
      if (sum(non_constant) == 0) {
        if (verbose) cat("    FALLBACK: no non-constant predictors",
                         "— using conventional imputation\n")
        
        imputed_values <- .conventional_impute(
          n             = n_to_impute,
          lower_bound   = lower_bound,
          lod_threshold = lod_threshold,
          Y_train       = Y_train
        )
        
        below_lod_subjects <- merged_data$subject_id[below_lod_mask]
        update_indices     <- which(
          current_data$analyte_visit == av_level &
            current_data$subject_id   %in% below_lod_subjects
        )
        
        if (length(update_indices) == length(imputed_values))
          current_data$result_log[update_indices] <- imputed_values
        
        if (iter == maxit_truncated) {
          conv_mean_val <- lod_threshold / sqrt(2)
          quality_metrics_list[[av_level]] <- data.frame(
            analyte_code      = analyte_code,
            visit             = current_visit,
            analyte_visit     = av_level,
            lod_upper         = lod_orig,
            lower_bound       = lb_orig,
            n_total           = nrow(merged_data),
            n_train           = sum(train_mask),
            n_imputed         = n_to_impute,
            pct_below_lod     = round(100 * n_to_impute / nrow(merged_data), 2),
            status            = "conventional_fallback: no_valid_predictors",
            imputation_method = "conventional",
            mse               = NA_real_, rmse = NA_real_, mae = NA_real_,
            r_squared         = NA_real_, adj_r_squared = NA_real_,
            cv_risk           = NA_real_, best_algorithm = NA_character_,
            n_predictors      = 0L,       adaptive_weight = NA_real_,
            guided_mean_min   = round(conv_mean_val, 4),
            guided_mean_max   = round(conv_mean_val, 4),
            guided_sd         = NA_real_,
            conventional_mean = round(conv_mean_val, 4),
            conventional_sd   = round(sd(Y_train, na.rm = TRUE), 4),
            imp_mean_log      = round(mean(imputed_values), 4),
            imp_sd_log        = round(sd(imputed_values),   4),
            stringsAsFactors  = FALSE
          )
        }
        next
      }
      
      X_train      <- X_train[,   non_constant, drop = FALSE]
      X_predict    <- X_predict[, non_constant, drop = FALSE]
      n_preds_used <- ncol(X_train)
      
      # SuperLearner fit — if this errors, fall back to conventional
      tryCatch({
        
        sl_fit <- SuperLearner(
          Y          = Y_train,
          X          = X_train,
          family     = gaussian(),
          SL.library = sl_library,
          cvControl  = list(V = min(5L, sum(train_mask)))
        )
        
        sl_pred_train <- as.numeric(
          predict(sl_fit, newdata = X_train,   onlySL = TRUE)$pred)
        sl_pred       <- as.numeric(
          predict(sl_fit, newdata = X_predict, onlySL = TRUE)$pred)
        
        # Quality metrics
        residuals     <- Y_train - sl_pred_train
        mse           <- mean(residuals^2)
        rmse          <- sqrt(mse)
        mae           <- mean(abs(residuals))
        ss_tot        <- sum((Y_train - mean(Y_train))^2)
        ss_res        <- sum(residuals^2)
        r_squared     <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA_real_)
        n_tr          <- length(Y_train)
        adj_r_squared <- ifelse(
          !is.na(r_squared) && n_tr > n_preds_used + 1,
          1 - (1 - r_squared) * (n_tr - 1) / (n_tr - n_preds_used - 1),
          NA_real_
        )
        cv_risk   <- min(sl_fit$cvRisk, na.rm = TRUE)
        best_algo <- names(sl_fit$coef)[which.max(sl_fit$coef)]
        sl_coef   <- sl_fit$coef
        
        # Adaptive blend weight
        w_sl <- if (adaptive_weight) {
          sqrt(max(0, min(1, ifelse(is.na(r_squared), 0, r_squared))))
        } else {
          weight_sl
        }
        
        # Cap SL predictions within [lower_bound, lod_threshold]
        sl_pred_capped <- pmin(pmax(sl_pred, lower_bound), lod_threshold)
        
        # Conventional parameters
        conventional_mean <- lod_threshold / sqrt(2)
        conventional_sd   <- sd(Y_train, na.rm = TRUE)
        
        # Blend SL with conventional
        guided_mean <- w_sl * sl_pred_capped + (1 - w_sl) * conventional_mean
        residual_sd <- sd(residuals)
        guided_sd   <- w_sl * residual_sd + (1 - w_sl) * conventional_sd
        
        # Cascading SD fallback
        if (is.na(guided_sd) || guided_sd <= 0) guided_sd <- conventional_sd
        if (is.na(guided_sd) || guided_sd <= 0) guided_sd <- (lod_threshold - lower_bound) / 4
        if (is.na(guided_sd) || guided_sd <= 0) guided_sd <- 1e-6
        
        # Draw truncated normal samples
        imputed_values <- rtruncnorm(
          n    = n_to_impute,
          a    = lower_bound,
          b    = lod_threshold,
          mean = guided_mean,
          sd   = guided_sd
        )
        
        # Safe subject-ID-based assignment
        below_lod_subjects <- merged_data$subject_id[below_lod_mask]
        update_indices     <- which(
          current_data$analyte_visit == av_level &
            current_data$subject_id   %in% below_lod_subjects
        )
        
        if (length(update_indices) != length(imputed_values)) {
          warning(paste0("[", av_level, "] Assignment mismatch (",
                         length(update_indices), " vs ",
                         length(imputed_values), "). Skipping."))
        } else {
          current_data$result_log[update_indices] <- imputed_values
        }
        
        # Convergence tracking
        if (track_convergence) {
          iter_convergence[[av_level]] <- data.frame(
            iter              = iter,
            analyte_visit     = av_level,
            imp_mean          = mean(imputed_values),
            imp_sd            = sd(imputed_values),
            imp_min           = min(imputed_values),
            imp_max           = max(imputed_values),
            r_squared         = r_squared,
            adaptive_weight   = w_sl,
            imputation_method = "superlearner",
            stringsAsFactors  = FALSE
          )
        }
        
        # Full diagnostics on final iteration
        if (iter == maxit_truncated) {
          quality_metrics_list[[av_level]] <- data.frame(
            analyte_code      = analyte_code,
            visit             = current_visit,
            analyte_visit     = av_level,
            lod_upper         = lod_orig,
            lower_bound       = lb_orig,
            n_total           = nrow(merged_data),
            n_train           = sum(train_mask),
            n_imputed         = n_to_impute,
            pct_below_lod     = round(100 * n_to_impute / nrow(merged_data), 2),
            status            = "success",
            imputation_method = "superlearner",
            mse               = round(mse,             6),
            rmse              = round(rmse,             6),
            mae               = round(mae,              6),
            r_squared         = round(r_squared,        4),
            adj_r_squared     = round(adj_r_squared,    4),
            cv_risk           = round(cv_risk,          6),
            best_algorithm    = best_algo,
            n_predictors      = n_preds_used,
            adaptive_weight   = round(w_sl,             4),
            guided_mean_min   = round(min(guided_mean), 4),
            guided_mean_max   = round(max(guided_mean), 4),
            guided_sd         = round(guided_sd,        4),
            conventional_mean = round(conventional_mean, 4),
            conventional_sd   = round(conventional_sd,   4),
            imp_mean_log      = round(mean(imputed_values), 4),
            imp_sd_log        = round(sd(imputed_values),   4),
            stringsAsFactors  = FALSE
          )
          
          imputation_details[[av_level]] <- list(
            analyte_code         = analyte_code,
            visit                = current_visit,
            imputation_method    = "superlearner",
            lod_threshold_log    = lod_threshold,
            lower_bound_log      = lower_bound,
            lod_original         = lod_orig,
            lower_bound_original = lb_orig,
            n_train              = sum(train_mask),
            n_imputed            = n_to_impute,
            adaptive_weight      = w_sl,
            sl_coef              = sl_coef,
            sl_pred_range        = range(sl_pred),
            sl_pred_capped_range = range(sl_pred_capped),
            guided_mean_range    = range(guided_mean),
            guided_sd            = guided_sd,
            conventional_mean    = conventional_mean,
            conventional_sd      = conventional_sd,
            imputed_range_log    = range(imputed_values),
            quality_metrics      = list(
              mse = mse, rmse = rmse, mae = mae,
              r_squared = r_squared, adj_r_squared = adj_r_squared,
              cv_risk = cv_risk
            )
          )
        }
        
        if (verbose) cat(paste0(
          "    [SL] n_train=", sum(train_mask),
          " | n_imp=", n_to_impute,
          " | R2=",    round(r_squared, 3),
          " | w_SL=",  round(w_sl, 3),
          " | best=",  best_algo, "\n"
        ))
        
      }, error = function(e) {
        # SL errored — fall back to conventional so no observation is left unimputed
        if (verbose) cat(paste("    FALLBACK (SL error):", e$message,
                               "— using conventional imputation\n"))
        
        imputed_values <- .conventional_impute(
          n             = n_to_impute,
          lower_bound   = lower_bound,
          lod_threshold = lod_threshold,
          Y_train       = Y_train
        )
        
        below_lod_subjects <- merged_data$subject_id[below_lod_mask]
        update_indices     <- which(
          current_data$analyte_visit == av_level &
            current_data$subject_id   %in% below_lod_subjects
        )
        
        if (length(update_indices) == length(imputed_values))
          current_data$result_log[update_indices] <<- imputed_values
        
        if (iter == maxit_truncated) {
          conv_mean_val <- lod_threshold / sqrt(2)
          quality_metrics_list[[av_level]] <<- data.frame(
            analyte_code      = analyte_code,
            visit             = current_visit,
            analyte_visit     = av_level,
            lod_upper         = lod_orig,
            lower_bound       = lb_orig,
            n_total           = nrow(merged_data),
            n_train           = sum(train_mask),
            n_imputed         = n_to_impute,
            pct_below_lod     = round(100 * n_to_impute / nrow(merged_data), 2),
            status            = paste0("conventional_fallback: sl_error: ",
                                       substr(e$message, 1, 60)),
            imputation_method = "conventional",
            mse               = NA_real_, rmse = NA_real_, mae = NA_real_,
            r_squared         = NA_real_, adj_r_squared = NA_real_,
            cv_risk           = NA_real_, best_algorithm = NA_character_,
            n_predictors      = n_preds_used,
            adaptive_weight   = NA_real_,
            guided_mean_min   = round(conv_mean_val, 4),
            guided_mean_max   = round(conv_mean_val, 4),
            guided_sd         = NA_real_,
            conventional_mean = round(conv_mean_val, 4),
            conventional_sd   = round(sd(Y_train, na.rm = TRUE), 4),
            imp_mean_log      = round(mean(imputed_values), 4),
            imp_sd_log        = round(sd(imputed_values),   4),
            stringsAsFactors  = FALSE
          )
        }
      })
      
    }  # end av_level loop
    
    if (track_convergence && length(iter_convergence) > 0)
      convergence_tracker[[iter]] <- iter_convergence
    
  }  # end Gibbs loop
  
  # ---------------------------------------------------------------------------
  # 6. Back-transform
  # ---------------------------------------------------------------------------
  if (log_transform) {
    current_data <- current_data %>%
      mutate(
        result_imputed = ifelse(flg_lod == 1, exp(result_log), NA_real_),
        result         = ifelse(flg_lod == 1, exp(result_log), result_original)
      ) %>%
      select(-result_original, -result_log, -analyte_lod_log,
             -analyte_code, -visit)
  } else {
    current_data <- current_data %>%
      mutate(
        result_imputed = ifelse(flg_lod == 1, result_log, NA_real_),
        result         = ifelse(flg_lod == 1, result_log, result_original)
      ) %>%
      select(-result_original, -result_log, -analyte_lod_log,
             -analyte_code, -visit)
  }
  
  # ---------------------------------------------------------------------------
  # 7. Verify all flg_lod == 1 observations received an imputed value
  # ---------------------------------------------------------------------------
  n_still_missing <- sum(
    current_data$flg_lod == 1 & is.na(current_data$result_imputed),
    na.rm = TRUE
  )
  
  if (n_still_missing > 0) {
    warning(paste(n_still_missing,
                  "flg_lod=1 observations still have NA after imputation.",
                  "Check imputation_details for these analyte_visit levels."))
  }
  
  # ---------------------------------------------------------------------------
  # 8. Assemble outputs
  # ---------------------------------------------------------------------------
  quality_metrics_df <- bind_rows(quality_metrics_list)
  
  convergence_df <- NULL
  if (track_convergence && length(convergence_tracker) > 0)
    convergence_df <- bind_rows(lapply(convergence_tracker, bind_rows))
  
  if (verbose) {
    cat("\n=== Imputation Complete ===\n")
    n_combos   <- nrow(quality_metrics_df)
    n_sl       <- sum(quality_metrics_df$imputation_method == "superlearner", na.rm = TRUE)
    n_conv     <- sum(quality_metrics_df$imputation_method == "conventional", na.rm = TRUE)
    n_none     <- sum(quality_metrics_df$status == "no_imputation_needed",   na.rm = TRUE)
    cat(paste("Total analyte_visit combos       :", n_combos,          "\n"))
    cat(paste("  Imputed via SuperLearner        :", n_sl,              "\n"))
    cat(paste("  Imputed via conventional fallback:", n_conv,           "\n"))
    cat(paste("  No below-LOD obs (no imputation):", n_none,            "\n"))
    cat(paste("Total flg_lod=1 obs imputed      :",
              sum(quality_metrics_df$n_imputed, na.rm = TRUE), "/",
              total_to_impute,                                          "\n"))
    cat(paste("Obs still missing after imputation:", n_still_missing,   "\n"))
    
    sl_metrics <- quality_metrics_df %>%
      filter(imputation_method == "superlearner")
    if (nrow(sl_metrics) > 0) {
      cat("\nSuperLearner model quality (final iteration):\n")
      cat(paste("  Mean R-squared       :",
                round(mean(sl_metrics$r_squared,       na.rm = TRUE), 3), "\n"))
      cat(paste("  Mean RMSE            :",
                round(mean(sl_metrics$rmse,            na.rm = TRUE), 4), "\n"))
      cat(paste("  Mean Adaptive Weight :",
                round(mean(sl_metrics$adaptive_weight, na.rm = TRUE), 3), "\n"))
    }
  }
  
  result <- list(
    imputed_data       = current_data,
    quality_metrics    = quality_metrics_df,
    imputation_details = imputation_details,
    convergence        = convergence_df,
    parameters         = list(
      lower_bounds      = lower_bounds,
      maxit_truncated   = maxit_truncated,
      sl_library        = sl_library,
      seed              = seed,
      log_transform     = log_transform,
      adaptive_weight   = adaptive_weight,
      weight_sl         = weight_sl,
      track_convergence = track_convergence
    )
  )
  
  class(result) <- "truncated_si_sl"
  return(result)
}


# =============================================================================
# EXPORTED HELPERS
# =============================================================================

get_imputed_data <- function(si_result) {
  if (!inherits(si_result, "truncated_si_sl"))
    stop("Input must be class 'truncated_si_sl'.")
  si_result$imputed_data
}

# status_filter = NULL returns all rows
get_quality_metrics <- function(si_result, status_filter = "success") {
  if (!inherits(si_result, "truncated_si_sl"))
    stop("Input must be class 'truncated_si_sl'.")
  df <- si_result$quality_metrics
  if (!is.null(status_filter)) df <- df %>% filter(status == status_filter)
  df
}

get_convergence_trace <- function(si_result, av_level) {
  if (is.null(si_result$convergence))
    stop("Convergence tracking was not enabled (track_convergence = FALSE).")
  si_result$convergence %>% filter(analyte_visit == av_level)
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================
#
# result <- truncated_single_imputation_sl(
#   biomarker_df      = biomarker_df,
#   cov_pred_df       = cov_pred,
#   lower_bounds      = 0.0001,
#   maxit_truncated   = 50,
#   sl_library        = c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.mean"),
#   seed              = 42,
#   log_transform     = TRUE,
#   adaptive_weight   = TRUE,
#   track_convergence = TRUE,
#   verbose           = TRUE
# )
#
# imputed     <- get_imputed_data(result)
# metrics     <- get_quality_metrics(result, NULL)   # all combos
# sl_only     <- get_quality_metrics(result, "success")
# fallbacks   <- metrics %>% filter(imputation_method == "conventional")
#
# # Confirm all flg_lod=1 obs are imputed
# sum(is.na(imputed$result_imputed[imputed$flg_lod == 1]))  # should be 0
