# ── Helpers ───────────────────────────────────────────────────────────────────

.count_sig <- function(res, fdr_thresholds, pval_thresholds) {
  counts <- list()
  for (thr in fdr_thresholds) {
    key <- gsub("\\.", "_", formatC(thr, format = "fg"))
    counts[[paste0("fdr_", key)]] <- sum(res$FDR < thr, na.rm = TRUE)
  }
  for (thr in pval_thresholds) {
    key <- formatC(thr, format = "e", digits = 2)
    counts[[paste0("raw_p_", key)]] <- sum(res$P.value < thr, na.rm = TRUE)
  }
  counts
}

.run_ewas_once <- function(bVals, pheno, variable, covariates,
                           fdr_thresholds, pval_thresholds) {
  
  ewas <- cpg.assoc(
    bVals,
    indep           = pheno[[variable]],
    covariates      = as.data.frame(pheno[, covariates, drop = FALSE]),
    logit.transform = TRUE
  )
  
  if (!"FDR" %in% colnames(ewas$results)) {
    warning("'FDR' not found in cpg.assoc results.")
    ewas$results$FDR <- NA_real_
  }
  
  if (!"P.value" %in% colnames(ewas$results)) {
    warning("'P.value' not found in cpg.assoc results.")
    ewas$results$P.value <- NA_real_
  }
  
  # ── Lambda ───────────────────────────────────────────────────────────────────
  lambda_val <- tryCatch({
    lam <- QCEWAS::P_lambda(ewas$results$P.value)
    if (!is.na(lam) && lam >= 2)
      warning("Lambda = ", round(lam, 4), " (>= 2): substantial oversignificance detected. ",
              "Consider genomic control correction on p-values.")
    if (!is.na(lam) && lam <= 0.8)
      warning("Lambda = ", round(lam, 4), " (<= 0.8): results are less significant than ",
              "expected from a random p-value distribution.")
    lam
  }, error = function(e) {
    warning("Lambda computation failed: ", e$message)
    NA
  })
  
  list(
    ewas   = ewas,
    lambda = lambda_val,
    counts = .count_sig(ewas$results, fdr_thresholds, pval_thresholds)
  )
}

# ── Main function ──────────────────────────────────────────────────────────────

ewas_loop_pll <- function(exposures,
                          covariates,
                          bVals_list,
                          bVals_names,
                          dat_pheno,
                          outputFolder,
                          sex_var         = "sex",
                          fdr_thresholds  = c(0.20, 0.05),
                          pval_thresholds = 7.058787e-08) {
  
  # ── Package check ─────────────────────────────────────────────────────────────
  if (!requireNamespace("CpGassoc", quietly = TRUE))
    stop("Package 'CpGassoc' is not installed. Install it with: install.packages('CpGassoc')")
  library(CpGassoc)
  
  if (!requireNamespace("QCEWAS", quietly = TRUE))
    stop("Package 'QCEWAS' is not installed. Install it with: install.packages('QCEWAS')")
  
  # ── Input validation ──────────────────────────────────────────────────────────
  if (!is.list(bVals_list))
    stop("'bVals_list' must be a list.")
  if (length(bVals_list) != length(bVals_names))
    stop("'bVals_list' and 'bVals_names' must be the same length (",
         length(bVals_list), " vs ", length(bVals_names), ").")
  if (!"basename" %in% colnames(dat_pheno))
    stop("'dat_pheno' must contain a 'basename' column.")
  
  bad_exp <- setdiff(exposures,  colnames(dat_pheno))
  bad_cov <- setdiff(covariates, colnames(dat_pheno))
  if (length(bad_exp) > 0) stop("Exposure(s) not in dat_pheno: ",  paste(bad_exp, collapse = ", "))
  if (length(bad_cov) > 0) stop("Covariate(s) not in dat_pheno: ", paste(bad_cov, collapse = ", "))
  
  if (!sex_var %in% colnames(dat_pheno))
    stop("'sex_var' column '", sex_var, "' not found in dat_pheno.")
  
  if (!is.numeric(fdr_thresholds)  || any(fdr_thresholds  <= 0 | fdr_thresholds  >= 1))
    stop("'fdr_thresholds' must be numeric values strictly between 0 and 1.")
  if (!is.numeric(pval_thresholds) || any(pval_thresholds <= 0 | pval_thresholds >= 1))
    stop("'pval_thresholds' must be numeric values strictly between 0 and 1.")
  
  fdr_thresholds  <- sort(fdr_thresholds)
  pval_thresholds <- sort(pval_thresholds)
  
  sex_levels <- sort(unique(na.omit(dat_pheno[[sex_var]])))
  if (length(sex_levels) == 0)
    stop("No non-NA levels found in sex variable '", sex_var, "'.")
  message("Sex strata: ", paste(sex_levels, collapse = ", "))
  
  # ── Build ordered column-name vectors ─────────────────────────────────────────
  fdr_keys       <- gsub("\\.", "_", formatC(fdr_thresholds,  format = "fg"))
  pval_keys      <- formatC(pval_thresholds, format = "e", digits = 2)
  fdr_cols       <- paste0("fdr_",   fdr_keys)
  raw_p_cols     <- paste0("raw_p_", pval_keys)
  all_count_cols <- c(fdr_cols, raw_p_cols)
  
  # ── Output folder ─────────────────────────────────────────────────────────────
  if (!dir.exists(outputFolder)) {
    message("Creating output folder: ", outputFolder)
    dir.create(outputFolder, recursive = TRUE)
  }
  
  # ── Accumulators ──────────────────────────────────────────────────────────────
  acc <- new.env(parent = emptyenv())
  acc$failed      <- list()
  acc$run_summary <- list()
  
  # ── Console count reporter ────────────────────────────────────────────────────
  .msg_counts <- function(cts, label = "") {
    pfx <- if (nchar(label)) paste0("[", label, "] ") else ""
    for (thr in fdr_thresholds) {
      key <- gsub("\\.", "_", formatC(thr, format = "fg"))
      message(pfx, "FDR <", thr, ": ", cts[[paste0("fdr_", key)]])
    }
    for (thr in pval_thresholds) {
      key <- formatC(thr, format = "e", digits = 2)
      message(pfx, "raw.p <", thr, ": ", cts[[paste0("raw_p_", key)]])
    }
  }
  
  # ── Main loop ─────────────────────────────────────────────────────────────────
  for (variable in exposures) {
    for (i in seq_along(bVals_list)) {
      
      bVals      <- bVals_list[[i]]
      bVals_name <- bVals_names[i]
      
      message("\n=== variable: '", variable, "' | dataset: '", bVals_name, "' ===")
      
      tryCatch({
        
        # ── Align basenames ──────────────────────────────────────────────────────
        bVals_basenames <- colnames(bVals)
        
        n_miss_pheno <- length(setdiff(bVals_basenames, dat_pheno$basename))
        n_miss_bvals <- length(setdiff(dat_pheno$basename, bVals_basenames))
        if (n_miss_pheno > 0)
          warning(n_miss_pheno, " basename(s) in '", bVals_name,
                  "' have no match in dat_pheno and will be dropped.")
        if (n_miss_bvals > 0)
          message(n_miss_bvals, " participant(s) in dat_pheno have no methylation data in '",
                  bVals_name, "' and will be excluded.")
        
        shared <- intersect(bVals_basenames, dat_pheno$basename)
        if (length(shared) == 0)
          stop("No shared basenames between '", bVals_name, "' and dat_pheno.")
        
        bVals_aligned <- bVals[, shared, drop = FALSE]
        pheno_aligned <- dat_pheno[match(shared, dat_pheno$basename), ]
        
        if (!identical(colnames(bVals_aligned), pheno_aligned$basename))
          stop("Sample alignment failed for '", bVals_name, "'.")
        
        message("Aligned ", length(shared), " samples.")
        
        # ── Drop samples with missing exposure ───────────────────────────────────
        non_miss <- !is.na(pheno_aligned[[variable]])
        if (sum(non_miss) == 0)
          stop("No non-missing observations for '", variable, "' in '", bVals_name, "'.")
        if (sum(!non_miss) > 0)
          message(sum(!non_miss), " sample(s) with missing '", variable, "' excluded.")
        
        bVals_use <- bVals_aligned[, non_miss, drop = FALSE]
        pheno_use <- pheno_aligned[non_miss, ]
        message(sum(non_miss), " samples in analysis.")
        
        # ── Overall ──────────────────────────────────────────────────────────────
        message("\n-- Overall --")
        res_overall <- .run_ewas_once(bVals_use, pheno_use, variable, covariates,
                                      fdr_thresholds, pval_thresholds)
        message("Lambda: ", round(res_overall$lambda, 4))
        .msg_counts(res_overall$counts)
        
        rds_file <- file.path(outputFolder,
                              paste0("ewas_", variable, "_", bVals_name, ".rds"))
        saveRDS(res_overall$ewas, file = rds_file)
        message("Saved: ", basename(rds_file))
        
        # Build summary row as a proper named list with scalar values only
        overall_row <- data.frame(
          exposure = variable,
          dataset  = bVals_name,
          stratum  = "Overall",
          n        = sum(non_miss),
          lambda   = round(res_overall$lambda, 4),
          stringsAsFactors = FALSE
        )
        for (col in all_count_cols)
          overall_row[[col]] <- res_overall$counts[[col]]
        
        acc$run_summary <- c(acc$run_summary, list(overall_row))
        
        # ── Sex-stratified ───────────────────────────────────────────────────────
        for (lvl in sex_levels) {
          
          message("\n-- Sex stratum: ", sex_var, " = ", lvl, " --")
          
          sex_idx   <- !is.na(pheno_use[[sex_var]]) & pheno_use[[sex_var]] == lvl
          n_stratum <- sum(sex_idx)
          
          if (n_stratum == 0) {
            warning("No samples for ", sex_var, " = '", lvl, "' — skipping.")
            next
          }
          message(n_stratum, " samples in stratum.")
          
          res_strat <- tryCatch(
            .run_ewas_once(bVals_use[, sex_idx, drop = FALSE],
                           pheno_use[sex_idx, ],
                           variable, covariates,
                           fdr_thresholds, pval_thresholds),
            error = function(e) {
              warning("Stratified EWAS failed for ", sex_var, "='", lvl,
                      "', variable='", variable, "': ", e$message)
              NULL
            }
          )
          
          if (is.null(res_strat)) next
          
          message("Lambda: ", round(res_strat$lambda, 4))
          .msg_counts(res_strat$counts, label = paste0(sex_var, "=", lvl))
          
          rds_strat <- file.path(outputFolder,
                                 paste0("ewas_", variable, "_", bVals_name,
                                        "_", sex_var, "_", lvl, ".rds"))
          saveRDS(res_strat$ewas, file = rds_strat)
          message("Saved: ", basename(rds_strat))
          
          strat_row <- data.frame(
            exposure = variable,
            dataset  = bVals_name,
            stratum  = paste0(sex_var, "=", lvl),
            n        = n_stratum,
            lambda   = round(res_strat$lambda, 4),
            stringsAsFactors = FALSE
          )
          for (col in all_count_cols)
            strat_row[[col]] <- res_strat$counts[[col]]
          
          acc$run_summary <- c(acc$run_summary, list(strat_row))
        }
        
        message("\nCompleted: variable='", variable, "' | dataset='", bVals_name, "'")
        
      }, error = function(e) {
        msg <- paste0("variable='", variable, "', dataset='", bVals_name, "': ", e$message)
        warning("Processing failed — ", msg)
        acc$failed <- c(acc$failed, list(msg))
      })
    }
  }
  
  # ── Summary CSV ───────────────────────────────────────────────────────────────
  if (length(acc$run_summary) > 0) {
    
    summary_df <- do.call(rbind, acc$run_summary)
    
    # Rename count columns to human-readable headers
    fdr_hdrs  <- paste0("FDR<",   fdr_thresholds)
    pval_hdrs <- paste0("raw.p<", formatC(pval_thresholds, format = "e", digits = 2))
    names(summary_df)[names(summary_df) %in% fdr_cols]   <- fdr_hdrs
    names(summary_df)[names(summary_df) %in% raw_p_cols] <- pval_hdrs
    
    csv_file <- file.path(outputFolder, "ewas_summary.csv")
    write.csv(summary_df, file = csv_file, row.names = FALSE)
    message("\nSummary saved: ", csv_file)
    
    # Print to console as well
    message("\n--- EWAS Run Summary ---")
    print(summary_df, row.names = FALSE)
  }
  
  # ── Failed iterations log ─────────────────────────────────────────────────────
  if (length(acc$failed) > 0) {
    failed_file <- file.path(outputFolder, "ewas_failed.txt")
    writeLines(c("Failed iterations:", unlist(acc$failed)), con = failed_file)
    message("\n", length(acc$failed), " iteration(s) failed — see: ", failed_file)
    for (f in acc$failed) message("  - ", f)
  } else {
    message("\nAll iterations completed successfully.")
  }
  
  invisible(list(failed = acc$failed, summary = acc$run_summary))
}