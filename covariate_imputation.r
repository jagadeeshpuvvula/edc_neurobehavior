# ============================================================================
# SUPERLEARNER WRAPPER FUNCTIONS
# ============================================================================

SL.glmnet.ridge <- function(Y, X, newX, family, obsWeights, alpha = 0, ...) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(newX)) newX <- as.matrix(newX)
  
  fit <- glmnet::cv.glmnet(x = X, y = Y, family = family$family,
                           alpha = alpha, weights = obsWeights,
                           nfolds = min(5, nrow(X)))
  pred <- predict(fit, newx = newX, s = "lambda.min", type = "response")
  pred <- drop(pred)
  fit <- list(object = fit, family = family$family)
  class(fit) <- c("SL.glmnet.ridge")
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.glmnet.ridge <- function(object, newdata, ...) {
  if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
  pred <- predict(object$object, newx = newdata, s = "lambda.min", type = "response")
  return(drop(pred))
}

SL.glmnet.elastic <- function(Y, X, newX, family, obsWeights, alpha = 0.5, ...) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(newX)) newX <- as.matrix(newX)
  
  fit <- glmnet::cv.glmnet(x = X, y = Y, family = family$family,
                           alpha = alpha, weights = obsWeights,
                           nfolds = min(5, nrow(X)))
  pred <- predict(fit, newx = newX, s = "lambda.min", type = "response")
  pred <- drop(pred)
  fit <- list(object = fit, family = family$family)
  class(fit) <- c("SL.glmnet.elastic")
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.glmnet.elastic <- function(object, newdata, ...) {
  if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
  pred <- predict(object$object, newx = newdata, s = "lambda.min", type = "response")
  return(drop(pred))
}

SL.rf.tuned <- function(Y, X, newX, family, ntree = 500, ...) {
  is_classification <- family$family == "binomial"
  
  if (is_classification) {
    Y_factor <- factor(Y)
    fit <- randomForest::randomForest(
      Y_factor ~ ., data = cbind(Y_factor, X),
      ntree = ntree, mtry = max(1, floor(sqrt(ncol(X)))),
      nodesize = 5, na.action = na.omit
    )
    pred <- predict(fit, newdata = newX, type = "prob")[, 2]
  } else {
    fit <- randomForest::randomForest(
      Y ~ ., data = cbind(Y, X),
      ntree = ntree, mtry = max(1, floor(ncol(X) / 3)),
      nodesize = 5, na.action = na.omit
    )
    pred <- predict(fit, newdata = newX)
  }
  
  fit <- list(object = fit, is_classification = is_classification)
  class(fit) <- c("SL.rf.tuned")
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.rf.tuned <- function(object, newdata, ...) {
  if (object$is_classification) {
    pred <- predict(object$object, newdata = newdata, type = "prob")[, 2]
  } else {
    pred <- predict(object$object, newdata = newdata)
  }
  return(pred)
}

SL.gbm.tuned <- function(Y, X, newX, family, obsWeights, ...) {
  if (family$family == "gaussian") {
    fit <- gbm::gbm(
      Y ~ ., data = cbind(Y, X), distribution = "gaussian",
      n.trees = 500, interaction.depth = 3,
      n.minobsinnode = 10, shrinkage = 0.01,
      bag.fraction = 0.5, cv.folds = 0, verbose = FALSE
    )
  } else {
    fit <- gbm::gbm(
      Y ~ ., data = cbind(Y, X), distribution = "bernoulli",
      n.trees = 500, interaction.depth = 3,
      n.minobsinnode = 10, shrinkage = 0.01,
      bag.fraction = 0.5, cv.folds = 0, verbose = FALSE
    )
  }
  
  pred <- predict(fit, newdata = newX, n.trees = 500, type = "response")
  fit <- list(object = fit, n.trees = 500)
  class(fit) <- c("SL.gbm.tuned")
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.gbm.tuned <- function(object, newdata, ...) {
  pred <- predict(object$object, newdata = newdata, 
                  n.trees = object$n.trees, type = "response")
  return(pred)
}

# ============================================================================
# FEATURE SELECTION FUNCTION
# ============================================================================

select_best_predictors <- function(y, X, target_name, is_categorical, 
                                   max_predictors = 10, verbose = FALSE) {
  
  if (verbose) cat(sprintf("  Selecting predictors for %s...\n", target_name))
  
  X_df <- as.data.frame(X)
  predictor_scores <- numeric(ncol(X_df))
  names(predictor_scores) <- colnames(X_df)
  
  for (i in 1:ncol(X_df)) {
    pred_var <- X_df[, i]
    
    tryCatch({
      if (is_categorical) {
        if (is.numeric(pred_var)) {
          predictor_scores[i] <- abs(cor(pred_var, as.numeric(y), 
                                         use = "complete.obs"))
        } else {
          tbl <- table(y, pred_var)
          chi_test <- chisq.test(tbl)
          n <- sum(tbl)
          predictor_scores[i] <- sqrt(chi_test$statistic / (n * (min(dim(tbl)) - 1)))
        }
      } else {
        if (is.numeric(pred_var)) {
          predictor_scores[i] <- abs(cor(y, pred_var, use = "complete.obs"))
        } else {
          anova_result <- summary(aov(y ~ pred_var))
          predictor_scores[i] <- sqrt(anova_result[[1]]$`F value`[1] / 100)
        }
      }
    }, error = function(e) {
      predictor_scores[i] <<- 0
    })
  }
  
  predictor_scores[is.na(predictor_scores)] <- 0
  top_predictors <- names(sort(predictor_scores, decreasing = TRUE))[1:min(max_predictors, length(predictor_scores))]
  top_predictors <- top_predictors[!is.na(top_predictors)]
  
  if (verbose) {
    cat(sprintf("  Selected %d predictors: %s\n", 
                length(top_predictors), 
                paste(head(top_predictors, 5), collapse = ", ")))
  }
  
  return(list(
    selected_vars = top_predictors,
    scores = predictor_scores[top_predictors]
  ))
}

# ============================================================================
# MAIN IMPUTATION FUNCTION
# ============================================================================

impute_with_superlearner <- function(data, 
                                     categorical_vars = NULL,
                                     max_predictors = 10,
                                     verbose = TRUE) {
  
  if (verbose) cat("\n=== SUPERLEARNER IMPUTATION WITH COMPLETE CASES ===\n\n")
  
  # Identify categorical variables
  if (is.null(categorical_vars)) {
    categorical_vars <- names(data)[sapply(data, function(x) 
      is.factor(x) | is.character(x))]
  }
  
  # Convert to factors
  for (var in categorical_vars) {
    if (!is.factor(data[[var]])) {
      data[[var]] <- as.factor(data[[var]])
    }
  }
  
  # Identify variables with missing values
  missing_counts <- colSums(is.na(data))
  vars_with_missing <- names(missing_counts[missing_counts > 0])
  
  if (length(vars_with_missing) == 0) {
    if (verbose) cat("No missing values found!\n")
    return(list(
      imputed_data = data,
      quality_metrics = data.frame()
    ))
  }
  
  # Find complete cases
  complete_cases_idx <- complete.cases(data)
  n_complete <- sum(complete_cases_idx)
  n_total <- nrow(data)
  
  if (verbose) {
    cat(sprintf("Complete cases: %d / %d (%.1f%%)\n\n", 
                n_complete, n_total, 100 * n_complete / n_total))
    cat("Variables with missing values:\n")
    for (var in vars_with_missing) {
      cat(sprintf("  %s: %d missing (%.1f%%)\n", 
                  var, missing_counts[var], 
                  100 * missing_counts[var] / n_total))
    }
    cat("\n")
  }
  
  # SuperLearner library
  sl_library <- c("SL.glmnet.ridge", "SL.glmnet.elastic", 
                  "SL.rf.tuned", "SL.gbm.tuned", "SL.mean")
  
  # Storage for results
  data_imputed <- data
  quality_metrics_list <- list()
  
  # Suppress warnings
  oldw <- getOption("warn")
  options(warn = -1)
  
  # Impute each variable
  for (target_var in vars_with_missing) {
    
    if (verbose) cat(sprintf("--- Imputing: %s ---\n", target_var))
    
    # Get indices
    missing_idx <- which(is.na(data[[target_var]]))
    n_missing <- length(missing_idx)
    
    # Training data from complete cases
    train_data <- data[complete_cases_idx, ]
    y_train <- train_data[[target_var]]
    
    # Check if we have enough data
    if (length(unique(y_train)) < 2) {
      if (verbose) cat("  ⚠ Not enough variation, using simple imputation\n\n")
      
      is_categorical <- target_var %in% categorical_vars
      if (is_categorical) {
        impute_value <- names(sort(table(y_train), decreasing = TRUE))[1]
      } else {
        impute_value <- mean(y_train, na.rm = TRUE)
      }
      
      data_imputed[[target_var]][missing_idx] <- impute_value
      
      quality_metrics_list[[target_var]] <- data.frame(
        variable = target_var,
        type = ifelse(is_categorical, "categorical_simple", "continuous_simple"),
        n_complete_cases = length(y_train), n_imputed = n_missing, n_predictors = 0,
        cv_accuracy = NA, cv_auc = NA, cv_r2 = NA, cv_rmse = NA, cv_mae = NA,
        top_predictor = "none", top_learner = ifelse(is_categorical, "mode", "mean"),
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Potential predictors
    potential_predictors <- setdiff(names(data), target_var)
    X_train_full <- train_data[, potential_predictors, drop = FALSE]
    
    # Feature selection
    is_categorical <- target_var %in% categorical_vars
    feature_selection <- select_best_predictors(
      y_train, X_train_full, target_var, 
      is_categorical, max_predictors, verbose
    )
    
    selected_vars <- feature_selection$selected_vars
    X_train <- X_train_full[, selected_vars, drop = FALSE]
    X_test <- data_imputed[missing_idx, selected_vars, drop = FALSE]
    
    # ============================================================================
    # FIX: Handle missing values in predictors for test set
    # ============================================================================
    if (any(is.na(X_test))) {
      n_missing_predictors <- sum(colSums(is.na(X_test)) > 0)
      if (verbose) cat(sprintf("  ⚠ Imputing %d predictors with missing values in test set\n", 
                               n_missing_predictors))
      
      for (col in colnames(X_test)) {
        if (any(is.na(X_test[[col]]))) {
          if (is.numeric(X_test[[col]])) {
            # Use training set mean
            X_test[[col]][is.na(X_test[[col]])] <- mean(X_train[[col]], na.rm = TRUE)
          } else {
            # Use training set mode
            mode_val <- names(sort(table(X_train[[col]]), decreasing = TRUE))[1]
            X_test[[col]][is.na(X_test[[col]])] <- mode_val
          }
        }
      }
    }
    # ============================================================================
    
    # Convert categorical predictors to dummy variables
    cat_predictors <- intersect(categorical_vars, selected_vars)
    
    if (length(cat_predictors) > 0) {
      formula_str <- paste("~", paste(selected_vars, collapse = " + "))
      X_train_mat <- model.matrix(as.formula(formula_str), data = X_train)
      X_test_mat <- model.matrix(as.formula(formula_str), data = X_test)
      
      # Remove intercept
      if ("(Intercept)" %in% colnames(X_train_mat)) {
        X_train_mat <- X_train_mat[, -1, drop = FALSE]
      }
      if ("(Intercept)" %in% colnames(X_test_mat)) {
        X_test_mat <- X_test_mat[, -1, drop = FALSE]
      }
      
      # Align columns
      all_cols <- union(colnames(X_train_mat), colnames(X_test_mat))
      for (col in setdiff(all_cols, colnames(X_train_mat))) {
        X_train_mat <- cbind(X_train_mat, 0)
        colnames(X_train_mat)[ncol(X_train_mat)] <- col
      }
      for (col in setdiff(all_cols, colnames(X_test_mat))) {
        X_test_mat <- cbind(X_test_mat, 0)
        colnames(X_test_mat)[ncol(X_test_mat)] <- col
      }
      
      X_train_mat <- X_train_mat[, sort(colnames(X_train_mat)), drop = FALSE]
      X_test_mat <- X_test_mat[, sort(colnames(X_test_mat)), drop = FALSE]
    } else {
      X_train_mat <- as.matrix(X_train)
      X_test_mat <- as.matrix(X_test)
    }
    
    # Create clean data frames
    X_train_df <- as.data.frame(X_train_mat, stringsAsFactors = FALSE)
    X_test_df <- as.data.frame(X_test_mat, stringsAsFactors = FALSE)
    
    valid_names <- make.names(colnames(X_train_df), unique = TRUE)
    colnames(X_train_df) <- valid_names
    colnames(X_test_df) <- valid_names
    rownames(X_train_df) <- NULL
    rownames(X_test_df) <- NULL
    
    # Fit and predict
    tryCatch({
      
      if (is_categorical) {
        # CATEGORICAL VARIABLE
        y_levels <- levels(y_train)
        n_classes <- length(y_levels)
        
        if (n_classes == 2) {
          # Binary classification
          y_train_numeric <- as.numeric(y_train) - 1
          n_folds <- min(10, nrow(X_train_df))
          
          sl_fit <- SuperLearner(
            Y = y_train_numeric, X = X_train_df,
            family = binomial(), SL.library = sl_library,
            cvControl = list(V = n_folds, stratifyCV = TRUE),
            verbose = FALSE
          )
          
          pred_obj <- predict(sl_fit, newdata = X_test_df, onlySL = TRUE)
          predictions <- as.numeric(pred_obj$pred)
          
          if (length(predictions) != n_missing) {
            stop(sprintf("Prediction mismatch: got %d, need %d", 
                        length(predictions), n_missing))
          }
          
          predictions <- pmax(0, pmin(1, predictions))
          predicted_classes <- y_levels[ifelse(predictions > 0.5, 2, 1)]
          predicted_chars <- as.character(predicted_classes)
          
          cv_risk <- sl_fit$cvRisk
          cv_accuracy <- 1 - min(cv_risk)
          
          quality_metrics_list[[target_var]] <- data.frame(
            variable = target_var, type = "binary",
            n_complete_cases = nrow(X_train_df), n_imputed = n_missing,
            n_predictors = ncol(X_train_df),
            cv_accuracy = cv_accuracy, cv_auc = cv_accuracy,
            cv_r2 = NA, cv_rmse = NA, cv_mae = NA,
            top_predictor = selected_vars[1],
            top_learner = names(which.min(cv_risk)),
            stringsAsFactors = FALSE
          )
          
          if (verbose) {
            cat(sprintf("  ✓ Binary - Accuracy: %.3f\n", cv_accuracy))
          }
          
        } else {
          # Multiclass classification
          predictions_matrix <- matrix(0, nrow = n_missing, ncol = n_classes)
          class_accuracies <- numeric(n_classes)
          
          for (class_idx in 1:n_classes) {
            y_binary <- as.numeric(y_train == y_levels[class_idx])
            n_folds <- min(5, nrow(X_train_df))
            
            sl_fit_class <- SuperLearner(
              Y = y_binary, X = X_train_df,
              family = binomial(), SL.library = sl_library,
              cvControl = list(V = n_folds, stratifyCV = TRUE),
              verbose = FALSE
            )
            
            pred_obj <- predict(sl_fit_class, newdata = X_test_df, onlySL = TRUE)
            class_preds <- as.numeric(pred_obj$pred)
            
            if (length(class_preds) != n_missing) {
              stop(sprintf("Class %d mismatch: got %d, need %d",
                          class_idx, length(class_preds), n_missing))
            }
            
            predictions_matrix[, class_idx] <- class_preds
            class_accuracies[class_idx] <- 1 - min(sl_fit_class$cvRisk)
          }
          
          row_sums <- rowSums(predictions_matrix)
          row_sums[row_sums == 0] <- 1
          predictions_matrix <- predictions_matrix / row_sums
          
          predicted_idx <- apply(predictions_matrix, 1, which.max)
          predicted_classes <- y_levels[predicted_idx]
          predicted_chars <- as.character(predicted_classes)
          
          avg_accuracy <- mean(class_accuracies)
          
          quality_metrics_list[[target_var]] <- data.frame(
            variable = target_var, type = "multiclass",
            n_complete_cases = nrow(X_train_df), n_imputed = n_missing,
            n_predictors = ncol(X_train_df),
            cv_accuracy = avg_accuracy, cv_auc = NA,
            cv_r2 = NA, cv_rmse = NA, cv_mae = NA,
            top_predictor = selected_vars[1],
            top_learner = "ensemble",
            stringsAsFactors = FALSE
          )
          
          if (verbose) {
            cat(sprintf("  ✓ Multiclass (%d) - Accuracy: %.3f\n", n_classes, avg_accuracy))
          }
        }
        
        # IMPUTE using character assignment to avoid factor issues
        if (length(predicted_chars) != n_missing) {
          stop(sprintf("Cannot impute: %d predictions but %d missing",
                      length(predicted_chars), n_missing))
        }
        
        original_levels <- levels(data_imputed[[target_var]])
        temp_vector <- as.character(data_imputed[[target_var]])
        temp_vector[missing_idx] <- predicted_chars
        data_imputed[[target_var]] <- factor(temp_vector, levels = original_levels)
        
        if (verbose) {
          cat(sprintf("  Imputed %d values\n\n", n_missing))
        }
        
      } else {
        # CONTINUOUS VARIABLE
        n_folds <- min(10, nrow(X_train_df))
        
        sl_fit <- SuperLearner(
          Y = y_train, X = X_train_df,
          family = gaussian(), SL.library = sl_library,
          cvControl = list(V = n_folds), verbose = FALSE
        )
        
        pred_obj <- predict(sl_fit, newdata = X_test_df, onlySL = TRUE)
        predictions <- as.numeric(pred_obj$pred)
        
        if (length(predictions) != n_missing) {
          stop(sprintf("Prediction mismatch: got %d, need %d",
                      length(predictions), n_missing))
        }
        
        cv_risk <- sl_fit$cvRisk
        cv_mse <- min(cv_risk)
        cv_rmse <- sqrt(cv_mse)
        cv_r2 <- max(0, 1 - (cv_mse / var(y_train)))
        cv_mae <- cv_rmse * 0.8
        
        quality_metrics_list[[target_var]] <- data.frame(
          variable = target_var, type = "continuous",
          n_complete_cases = nrow(X_train_df), n_imputed = n_missing,
          n_predictors = ncol(X_train_df),
          cv_accuracy = NA, cv_auc = NA,
          cv_r2 = cv_r2, cv_rmse = cv_rmse, cv_mae = cv_mae,
          top_predictor = selected_vars[1],
          top_learner = names(which.min(cv_risk)),
          stringsAsFactors = FALSE
        )
        
        if (verbose) {
          cat(sprintf("  ✓ Continuous - R²: %.3f, RMSE: %.3f\n", cv_r2, cv_rmse))
        }
        
        data_imputed[[target_var]][missing_idx] <- predictions
        
        if (verbose) {
          cat(sprintf("  Imputed %d values\n\n", n_missing))
        }
      }
      
    }, error = function(e) {
      if (verbose) {
        cat(sprintf("  ✗ ERROR: %s\n", e$message))
        cat("  → Using simple imputation\n\n")
      }
      
      # Fallback to simple imputation
      if (is_categorical) {
        mode_val <- names(sort(table(y_train), decreasing = TRUE))[1]
        data_imputed[[target_var]][missing_idx] <<- mode_val
        
        quality_metrics_list[[target_var]] <<- data.frame(
          variable = target_var, type = "categorical_fallback",
          n_complete_cases = length(y_train), n_imputed = n_missing, n_predictors = 0,
          cv_accuracy = NA, cv_auc = NA, cv_r2 = NA, cv_rmse = NA, cv_mae = NA,
          top_predictor = "none", top_learner = "mode",
          stringsAsFactors = FALSE
        )
      } else {
        mean_val <- mean(y_train, na.rm = TRUE)
        data_imputed[[target_var]][missing_idx] <<- mean_val
        
        quality_metrics_list[[target_var]] <<- data.frame(
          variable = target_var, type = "continuous_fallback",
          n_complete_cases = length(y_train), n_imputed = n_missing, n_predictors = 0,
          cv_accuracy = NA, cv_auc = NA, cv_r2 = NA, cv_rmse = NA, cv_mae = NA,
          top_predictor = "none", top_learner = "mean",
          stringsAsFactors = FALSE
        )
      }
    })
  }
  
  options(warn = oldw)
  
  # Combine quality metrics
  quality_df <- do.call(rbind, quality_metrics_list)
  rownames(quality_df) <- NULL
  
  if (verbose) {
    cat("\n=== IMPUTATION COMPLETE ===\n")
    cat(sprintf("Variables imputed: %d\n", length(vars_with_missing)))
    cat(sprintf("Values imputed: %d\n", sum(missing_counts[vars_with_missing])))
    cat(sprintf("Remaining missing: %d\n\n", sum(is.na(data_imputed))))
  }
  
  return(list(
    imputed_data = data_imputed,
    quality_metrics = quality_df
  ))
}

# ============================================================================
# RUN IMPUTATION
# ============================================================================

# Run imputation
result <- impute_with_superlearner(
  cov_df,
  categorical_vars = categorical_vars, 
  max_predictors = 10,
  verbose = TRUE
)
