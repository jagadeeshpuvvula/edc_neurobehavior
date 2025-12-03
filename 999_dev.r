# Function to run SuperLearner for each analyte_visit level
run_sl_by_analyte <- function(biomarker_df, cov_pred_df) {
  
  # Get unique analyte_visit levels
  analyte_levels <- unique(biomarker_df$analyte_visit)
  
  # Store results
  results_list <- list()
  
  for (analyte in analyte_levels) {
    
    cat("\n=== Analyzing:", analyte, "===\n")
    
    # Filter biomarker data for current analyte_visit
    bio_subset <- biomarker_df %>%
      filter(analyte_visit == analyte) %>%
      select(subject_id, result)
    
    # Merge with predictors
    merged_data <- bio_subset %>%
      inner_join(cov_pred_df, by = "subject_id")
    
    # Check if we have enough data
    if (nrow(merged_data) < 10) {
      cat("Skipping - insufficient data (n =", nrow(merged_data), ")\n")
      next
    }
    
    # Prepare data for SuperLearner
    # Remove subject_id and result to get predictors
    X <- merged_data %>%
      select(-subject_id, -result) %>%
      as.data.frame()
    
    Y <- merged_data$result
    
    # Remove rows with missing Y values first
    missing_y <- is.na(Y)
    if (any(missing_y)) {
      cat("Removing", sum(missing_y), "rows with missing outcome values\n")
      X <- X[!missing_y, , drop = FALSE]
      Y <- Y[!missing_y]
    }
    
    # Check for columns with all missing values and remove them
    all_na_cols <- sapply(X, function(col) all(is.na(col)))
    if (any(all_na_cols)) {
      cat("Removing", sum(all_na_cols), "columns with all missing values\n")
      X <- X[, !all_na_cols, drop = FALSE]
    }
    
    # Remove columns with near-zero variance or all same values
    constant_cols <- sapply(X, function(col) {
      if (is.numeric(col)) {
        length(unique(na.omit(col))) <= 1
      } else {
        FALSE
      }
    })
    if (any(constant_cols)) {
      cat("Removing", sum(constant_cols), "constant columns\n")
      X <- X[, !constant_cols, drop = FALSE]
    }
    
    # Convert factors/characters to numeric before imputation
    X <- X %>% mutate(across(where(is.character), as.factor))
    X <- X %>% mutate(across(where(is.factor), as.numeric))
    
    # Handle missing values - imputation with median for numeric
    for (col in names(X)) {
      if (is.numeric(X[[col]])) {
        na_count <- sum(is.na(X[[col]]))
        if (na_count > 0) {
          median_val <- median(X[[col]], na.rm = TRUE)
          # If median is still NA (all values were NA), use 0
          if (is.na(median_val)) {
            X[[col]][is.na(X[[col]])] <- 0
          } else {
            X[[col]][is.na(X[[col]])] <- median_val
          }
        }
      }
    }
    
    # Final check: remove any remaining rows with NA
    complete_rows <- complete.cases(X)
    if (!all(complete_rows)) {
      cat("Removing", sum(!complete_rows), "rows with remaining missing values\n")
      X <- X[complete_rows, , drop = FALSE]
      Y <- Y[complete_rows]
    }
    
    # Final validation
    if (any(is.na(X)) || any(is.na(Y))) {
      cat("Still have missing values - skipping this analyte\n")
      cat("NA in X:", sum(is.na(X)), "NA in Y:", sum(is.na(Y)), "\n")
      next
    }
    
    # Define SuperLearner library
    # Adjust based on your needs and computational resources
    SL.library <- c("SL.glm", "SL.glmnet", "SL.randomForest", 
                    "SL.xgboost", "SL.mean")
    
    # Run SuperLearner with cross-validation
    set.seed(123)
    
    tryCatch({
      sl_fit <- SuperLearner(
        Y = Y,
        X = X,
        family = gaussian(),
        SL.library = SL.library,
        cvControl = list(V = 10)  # 10-fold cross-validation
      )
      
      # Extract results
      cat("\nCross-validated Risk (MSE):\n")
      print(sl_fit$cvRisk)
      
      cat("\nAlgorithm Weights:\n")
      print(sl_fit$coef)
      
      cat("\nBest Algorithm:", 
          names(sl_fit$coef)[which.max(sl_fit$coef)], "\n")
      
      # Variable importance using permutation (for Random Forest)
      if ("SL.randomForest" %in% SL.library) {
        rf_fit <- randomForest::randomForest(X, Y)
        var_imp <- randomForest::importance(rf_fit)
        cat("\nTop 10 Important Variables:\n")
        print(head(var_imp[order(var_imp[,1], decreasing = TRUE), , drop = FALSE], 10))
      }
      
      # Store results
      results_list[[analyte]] <- list(
        sl_fit = sl_fit,
        n_samples = nrow(merged_data),
        n_predictors = ncol(X),
        cv_risk = sl_fit$cvRisk,
        weights = sl_fit$coef
      )
      
    }, error = function(e) {
      cat("Error fitting model:", e$message, "\n")
    })
  }
  
  return(results_list)
}

# Main function for truncated multiple imputation with SuperLearner and parallel processing
truncated_multiple_imputation_sl_parallel <- function(
  biomarker_df,           # Biomarker dataset with subject_id, analyte_visit, result, flg_lod
  cov_pred_df,            # Covariate dataset with subject_id and predictor variables
  analyte_lod_map,        # Named vector: names = analyte names, values = LOD thresholds
  lower_bounds = 0.0001,  # Lower bound for truncation (single value or named vector)
  m = 5,                  # Number of multiple imputation datasets
  maxit_conventional = 100, # Max iterations for conventional MI
  maxit_truncated = 100,   # Max iterations for truncated MI
  sl_library = c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.mean"),
  seed = 123,             # Random seed
  log_transform = TRUE,   # Whether to log-transform values
  n_cores = NULL,         # Number of cores (NULL = detect automatically)
  verbose = TRUE          # Print progress messages
) {
  
  # Load required packages
  required_packages <- c("mice", "truncnorm", "SuperLearner", "dplyr", "parallel", "doParallel", "foreach")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }
  
  # Setup parallel backend
  if (is.null(n_cores)) {
    n_cores <- max(1, detectCores() - 1)  # Leave one core free
  }
  
  if (verbose) cat(paste("Setting up parallel processing with", n_cores, "cores\n"))
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export necessary objects to cluster
  clusterExport(cl, c("biomarker_df", "cov_pred_df", "analyte_lod_map", 
                      "lower_bounds", "log_transform", "sl_library", 
                      "maxit_truncated", "verbose"),
                envir = environment())
  
  # Load packages on each worker
  clusterEvalQ(cl, {
    library(SuperLearner)
    library(dplyr)
    library(truncnorm)
  })
  
  # Input validation
  required_cols <- c("subject_id", "analyte_visit", "result", "flg_lod")
  if (!all(required_cols %in% names(biomarker_df))) {
    stopCluster(cl)
    stop(paste("biomarker_df must contain columns:", paste(required_cols, collapse = ", ")))
  }
  
  if (!"subject_id" %in% names(cov_pred_df)) {
    stopCluster(cl)
    stop("cov_pred_df must contain subject_id column")
  }
  
  # Get unique analyte names from analyte_visit
  biomarker_df <- biomarker_df %>%
    mutate(analyte = sub("_.*", "", analyte_visit))
  
  unique_analytes <- unique(biomarker_df$analyte)
  
  if (!all(unique_analytes %in% names(analyte_lod_map))) {
    stopCluster(cl)
    stop(paste("analyte_lod_map must contain LOD values for all analytes:",
               paste(unique_analytes, collapse = ", ")))
  }
  
  # Handle lower_bounds
  if (length(lower_bounds) == 1) {
    lower_bounds <- setNames(rep(lower_bounds, length(unique_analytes)), unique_analytes)
  } else if (!all(unique_analytes %in% names(lower_bounds))) {
    stopCluster(cl)
    stop("If lower_bounds is a vector, it must have names matching all analytes")
  }
  
  # Transform LOD values if requested
  if (log_transform) {
    if (verbose) cat("Log-transforming biomarker values and LOD thresholds...\n")
    biomarker_df <- biomarker_df %>%
      mutate(result_original = result,
             result = ifelse(!is.na(result) & result > 0, log(result), result))
    log_lod_map <- log(analyte_lod_map)
    log_lower_bounds <- log(lower_bounds)
  } else {
    log_lod_map <- analyte_lod_map
    log_lower_bounds <- lower_bounds
  }
  
  # Get predictor variable names
  predictor_vars <- setdiff(names(cov_pred_df), "subject_id")
  
  if (verbose) {
    cat(paste("Number of predictor variables:", length(predictor_vars), "\n"))
    cat(paste("Number of analyte_visit levels:", length(unique(biomarker_df$analyte_visit)), "\n"))
  }
  
  # Define function to process one imputation dataset
  process_imputation <- function(imp_num, biomarker_df, cov_pred_df, 
                                 log_lod_map, log_lower_bounds, 
                                 predictor_vars, sl_library, 
                                 maxit_truncated, seed, log_transform) {
    
    set.seed(seed + imp_num)
    
    # Start with original biomarker data
    current_imputed_data <- biomarker_df
    
    # Gibbs sampling iterations
    for (iter in 1:maxit_truncated) {
      
      # Loop through each analyte_visit level
      for (av_level in unique(biomarker_df$analyte_visit)) {
        
        # Get analyte name and corresponding LOD
        analyte_name <- sub("_.*", "", av_level)
        lod_threshold <- log_lod_map[analyte_name]
        lower_bound <- log_lower_bounds[analyte_name]
        
        # Subset data for this analyte_visit
        av_indices <- which(current_imputed_data$analyte_visit == av_level)
        av_data <- current_imputed_data[av_indices, ]
        
        # Merge with covariates
        merged_data <- av_data %>%
          inner_join(cov_pred_df, by = "subject_id")
        
        # Skip if insufficient data
        if (nrow(merged_data) < 10) {
          next
        }
        
        # Identify observations to impute (flg_lod == 1)
        below_lod_mask <- merged_data$flg_lod == 1
        
        if (sum(below_lod_mask) == 0) {
          next
        }
        
        # Prepare training data (observed values only)
        train_mask <- !below_lod_mask & !is.na(merged_data$result)
        
        if (sum(train_mask) < 5) {
          next
        }
        
        # Prepare X (predictors) and Y (outcome)
        X_train <- merged_data[train_mask, predictor_vars, drop = FALSE]
        Y_train <- merged_data$result[train_mask]
        
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
        X_train <- X_train %>% mutate(across(where(is.character), as.factor))
        X_train <- X_train %>% mutate(across(where(is.factor), as.numeric))
        X_predict <- X_predict %>% mutate(across(where(is.character), as.factor))
        X_predict <- X_predict %>% mutate(across(where(is.factor), as.numeric))
        
        # Remove constant columns
        col_vars <- sapply(X_train, function(x) var(x, na.rm = TRUE))
        non_constant <- col_vars > 0 & !is.na(col_vars)
        X_train <- X_train[, non_constant, drop = FALSE]
        X_predict <- X_predict[, non_constant, drop = FALSE]
        
        # Skip if no valid predictors
        if (ncol(X_train) == 0) {
          next
        }
        
        # Fit SuperLearner model
        tryCatch({
          sl_fit <- SuperLearner(
            Y = Y_train,
            X = X_train,
            family = gaussian(),
            SL.library = sl_library,
            cvControl = list(V = min(5, sum(train_mask)))
          )
          
          # Get predictions
          sl_pred <- predict(sl_fit, newdata = X_predict, onlySL = TRUE)$pred
          
          # Estimate residual standard deviation
          sl_pred_train <- predict(sl_fit, newdata = X_train, onlySL = TRUE)$pred
          residual_sd <- sd(Y_train - sl_pred_train)
          
          # Generate truncated normal samples
          n_to_impute <- sum(below_lod_mask)
          imputed_values <- rtruncnorm(
            n = n_to_impute,
            a = lower_bound,
            b = lod_threshold,
            mean = sl_pred,
            sd = residual_sd
          )
          
          # Update imputed values in current dataset
          current_imputed_data$result[av_indices[below_lod_mask]] <- imputed_values
          
        }, error = function(e) {
          # Silent error handling in parallel
        })
      }
    }
    
    # Back-transform if needed and create result_imputed only for flg_lod == 1
    if (log_transform) {
      current_imputed_data <- current_imputed_data %>%
        mutate(
          result_imputed = ifelse(flg_lod == 1, result, NA),  # Only for imputed values
          result = exp(result)  # Back-transform all results
        )
    } else {
      current_imputed_data <- current_imputed_data %>%
        mutate(result_imputed = ifelse(flg_lod == 1, result, NA))
    }
    
    return(current_imputed_data)
  }
  
  # Run imputations in parallel
  if (verbose) cat(paste("\n=== Running", m, "imputations in parallel ===\n"))
  
  imputed_datasets_list <- foreach(
    imp_num = 1:m,
    .packages = c("SuperLearner", "dplyr", "truncnorm"),
    .errorhandling = "pass"
  ) %dopar% {
    process_imputation(
      imp_num = imp_num,
      biomarker_df = biomarker_df,
      cov_pred_df = cov_pred_df,
      log_lod_map = log_lod_map,
      log_lower_bounds = log_lower_bounds,
      predictor_vars = predictor_vars,
      sl_library = sl_library,
      maxit_truncated = maxit_truncated,
      seed = seed,
      log_transform = log_transform
    )
  }
  
  # Stop cluster
  stopCluster(cl)
  
  if (verbose) cat("\n=== Parallel imputation completed ===\n")
  
  # Create summary (run on first imputation for model info)
  model_info <- list()
  
  if (verbose) {
    cat("\n=== Creating summary information ===\n")
    # Run a quick check on first dataset to get counts
    for (av_level in unique(biomarker_df$analyte_visit)) {
      av_data <- biomarker_df %>%
        filter(analyte_visit == av_level) %>%
        inner_join(cov_pred_df, by = "subject_id")
      
      if (nrow(av_data) >= 10) {
        n_imputed <- sum(av_data$flg_lod == 1, na.rm = TRUE)
        n_train <- sum(av_data$flg_lod == 0 & !is.na(av_data$result))
        
        if (n_imputed > 0 && n_train >= 5) {
          model_info[[av_level]] <- list(
            n_train = n_train,
            n_imputed = n_imputed
          )
        }
      }
    }
    
    cat("\n=== Imputation Summary ===\n")
    cat(paste("Total analyte_visit levels processed:", length(model_info), "\n"))
    for (av_level in names(model_info)) {
      cat(paste("\n", av_level, ":\n"))
      cat(paste("  Training samples:", model_info[[av_level]]$n_train, "\n"))
      cat(paste("  Imputed samples:", model_info[[av_level]]$n_imputed, "\n"))
    }
  }
  
  # Return results
  result <- list(
    imputed_datasets = imputed_datasets_list,
    model_info = model_info,
    parameters = list(
      analyte_lod_map = analyte_lod_map,
      lower_bounds = lower_bounds,
      m = m,
      maxit_truncated = maxit_truncated,
      sl_library = sl_library,
      seed = seed,
      log_transform = log_transform,
      n_cores = n_cores
    )
  )
  
  class(result) <- "truncated_mi_sl"
  return(result)
}

# Helper function to extract a specific imputed dataset
get_imputed_dataset <- function(mi_result, dataset_number = 1) {
  if (dataset_number > length(mi_result$imputed_datasets)) {
    stop(paste("Only", length(mi_result$imputed_datasets), "datasets available"))
  }
  return(mi_result$imputed_datasets[[dataset_number]])
}

# Helper function to pool results across imputed datasets
pool_imputed_results <- function(mi_result, analysis_function) {
  results <- lapply(mi_result$imputed_datasets, analysis_function)
  return(results)
}