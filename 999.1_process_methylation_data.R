process_methylation_data <- function(
  basenames_df,
  pheno_data = NULL,
  output_folder,
  array_type = c("450k", "EPIC", "EPICv2"),
  detection_pval = 1e-6,
  beadcount_threshold = 3,
  sample_cutoff = 0.05,
  probe_detection_rate = 0.95,
  bisulfite_sd_threshold = 3,
  cpgs_to_keep_file = "~/Documents/k99_home/arzu_cpgs_to_drop/cpgs_to_keep.csv",
  compute_mvalues = TRUE,
  apply_bmiq = TRUE,
  apply_combat = TRUE,
  batch_variable = "batch",
  estimate_cell_proportions = TRUE,
  cell_reference = "FlowSorted.CordBlood.450k",
  cell_types = NULL,
  verbose = TRUE
) {
  suppressWarnings({
    options(matrixStats.useNames.NA = "deprecated")
  })
  
  # Load required packages
  require(minfi)
  require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  require(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  require(wateRmelon)
  require(readr)
  require(SummarizedExperiment)
  require(S4Vectors)
  require(sva)
  require(ENmix)
  require(FlowSorted.CordBlood.450k)
  
  array_type <- match.arg(array_type)
  probe_counts <- list()
  
  processing_log <- list(
    mixed_arrays_detected = FALSE,
    array_versions = NULL,
    samples_removed_qc = c(),
    probes_failed_detection = 0,
    probes_failed_beadcount = 0,
    samples_failed_bisulfite = c(),
    beadcount_available = TRUE,
    bmiq_applied = FALSE,
    bmiq_failed_samples = c(),
    bmiq_success_count = 0,
    combat_applied = FALSE,
    n_batches = 0,
    cell_proportions_estimated = FALSE,
    cell_reference_used = NULL
  )
  
  # ===== INPUT VALIDATION =====
  if (!is.data.frame(basenames_df)) stop("basenames_df must be a data frame")
  
  required_cols <- c("Sample_ID", "folder_location")
  if (!all(required_cols %in% colnames(basenames_df))) {
    stop("basenames_df must contain columns: Sample_ID, folder_location")
  }
  
  if (any(duplicated(basenames_df$Sample_ID))) {
    dup_ids <- basenames_df$Sample_ID[duplicated(basenames_df$Sample_ID)]
    stop("Duplicate Sample_IDs found: ", paste(unique(dup_ids), collapse = ", "))
  }
  
  basenames <- file.path(basenames_df$folder_location, basenames_df$Sample_ID)
  missing_files <- c()
  for (i in seq_along(basenames)) {
    red_file <- paste0(basenames[i], "_Red.idat")
    grn_file <- paste0(basenames[i], "_Grn.idat")
    if (!file.exists(red_file) || !file.exists(grn_file)) {
      missing_files <- c(missing_files, basenames_df$Sample_ID[i])
    }
  }
  if (length(missing_files) > 0) {
    stop("IDAT files not found for samples: ", paste(missing_files, collapse = ", "))
  }
  
  array_versions <- sapply(basenames, function(bn) {
    red_file <- paste0(bn, "_Red.idat")
    file_size <- file.info(red_file)$size
    if (file_size > 15000000) "EPICv2"
    else if (file_size > 7000000) "EPIC/EPICv1"
    else "450k"
  })
  names(array_versions) <- basenames_df$Sample_ID
  processing_log$array_versions <- array_versions
  
  if (length(unique(array_versions)) > 1) {
    processing_log$mixed_arrays_detected <- TRUE
    stop("Mixed array types detected. Process each array type separately.")
  }
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  if (is.null(pheno_data)) {
    pheno_data <- as.data.frame(basenames_df)
    rownames(pheno_data) <- basenames_df$Sample_ID
  } else {
    if (nrow(pheno_data) != nrow(basenames_df)) {
      stop("pheno_data must have the same number of rows as basenames_df")
    }
    pheno_data <- as.data.frame(pheno_data)
    rownames(pheno_data) <- basenames_df$Sample_ID
  }
  
  if (apply_combat) {
    if (!batch_variable %in% colnames(pheno_data)) {
      warning("Batch variable '", batch_variable, "' not found in pheno_data. ComBat will be skipped.")
      apply_combat <- FALSE
    } else if (verbose) {
      cat("ComBat batch correction will be applied using variable:", batch_variable, "\n\n")
    }
  }
  
  # ===== STEP 1: READ IDAT FILES =====
  if (verbose) cat("Step 1: Reading IDAT files (extended for bead count)...\n")
  
  rgSet <- tryCatch({
    read.metharray(basenames = basenames, extended = TRUE)
  }, error = function(e) {
    if (grepl("different array size", e$message)) {
      read.metharray(basenames = basenames, extended = TRUE, force = TRUE)
    } else stop(e)
  })
  colnames(rgSet) <- basenames_df$Sample_ID
  
  # CORRECTED: Convert to MethylSet to get actual CpG probe count (not M+U channels)
  mSet_temp <- preprocessRaw(rgSet)
  probe_counts$initial <- nrow(mSet_temp)
  rm(mSet_temp)
  gc(verbose = FALSE)
  
  if (verbose) cat("  Initial CpG probes:", probe_counts$initial, "\n")
  
  # ===== STEP 2: CALCULATE DETECTION P-VALUES =====
  if (verbose) cat("Step 2: Calculating detection p-values...\n")
  
  detP <- detectionP(rgSet)
  colnames(detP) <- colnames(rgSet)
  
  # ===== STEP 3: IDENTIFY AND REMOVE FAILED PROBES =====
  if (verbose) cat("Step 3: Identifying and removing failed probes...\n")
  if (verbose) cat("  Criteria: detection p >", detection_pval, "OR beads <", beadcount_threshold,
                  "in >=", sample_cutoff * 100, "% of samples\n", sep = "")
  
  bc <- tryCatch({
    beadcount(rgSet)
  }, error = function(e) {
    if (verbose) cat("  Warning: Bead count data not available\n")
    processing_log$beadcount_available <<- FALSE
    matrix(NA, nrow = nrow(rgSet), ncol = ncol(rgSet),
           dimnames = list(rownames(rgSet), colnames(rgSet)))
  })
  
  # Identify failed measurements
  failed_detection <- detP > detection_pval
  if (processing_log$beadcount_available) {
    failed_beadcount <- bc < beadcount_threshold | is.na(bc)
    failed_measurements <- failed_detection | failed_beadcount
  } else {
    failed_measurements <- failed_detection
  }
  
  # Calculate failure rate per probe
  n_samples <- ncol(rgSet)
  failure_count_per_probe <- rowSums(failed_measurements)
  failure_rate_per_probe <- failure_count_per_probe / n_samples
  
  # Identify probes that failed in >= sample_cutoff % of samples
  failed_probes <- failure_rate_per_probe >= sample_cutoff
  n_failed_probes <- sum(failed_probes)
  
  if (verbose) {
    cat("  Total samples:", n_samples, "\n")
    cat("  Total probes:", nrow(detP), "\n")
    cat("  Failed measurements:", sum(failed_measurements), "\n")
    cat("  Probes failing in >=", sample_cutoff * 100, "% of samples:", n_failed_probes, "\n", sep = "")
  }
  
  # Remove failed probes from rgSet, detP, and bc
  if (n_failed_probes > 0) {
    rgSet <- rgSet[!failed_probes, ]
    detP <- detP[!failed_probes, , drop = FALSE]
    bc <- bc[!failed_probes, , drop = FALSE]
    if (verbose) cat("  Probes remaining:", nrow(detP), "\n")
  }
  
  probe_counts$after_step3 <- nrow(detP)
  processing_log$probes_failed_qc <- n_failed_probes
  
  # Clear memory
  rm(failed_detection, failed_measurements, failure_count_per_probe, failure_rate_per_probe, failed_probes)
  if (processing_log$beadcount_available) rm(failed_beadcount)
  gc(verbose = FALSE)
  
  # ===== STEP 4: SAMPLE QC =====
  if (verbose) cat("Step 4: Sample QC (remove samples with poor overall quality)...\n")
  
  # Recalculate failed measurements on remaining probes
  failed_detection_step4 <- detP > detection_pval
  if (processing_log$beadcount_available) {
    failed_beadcount_step4 <- bc < beadcount_threshold | is.na(bc)
    failed_measurements_step4 <- failed_detection_step4 | failed_beadcount_step4
  } else {
    failed_measurements_step4 <- failed_detection_step4
  }
  
  # Calculate failure rate per sample
  failure_count_per_sample <- colSums(failed_measurements_step4)
  failure_rate_per_sample <- failure_count_per_sample / nrow(detP)
  
  # Identify samples with > probe_detection_rate failure
  failed_samples <- failure_rate_per_sample > (1 - probe_detection_rate)
  n_failed_samples <- sum(failed_samples)
  
  if (verbose) {
    cat("  Samples with >", (1 - probe_detection_rate) * 100, "% failed probes:", n_failed_samples, "\n", sep = "")
  }
  
  # Remove failed samples
  if (n_failed_samples > 0) {
    failed_sample_ids <- colnames(rgSet)[failed_samples]
    processing_log$samples_removed_qc <- failed_sample_ids
    rgSet <- rgSet[, !failed_samples]
    detP <- detP[, !failed_samples, drop = FALSE]
    bc <- bc[, !failed_samples, drop = FALSE]
    pheno_data <- pheno_data[!failed_samples, , drop = FALSE]
    if (verbose) {
      cat("  Removed samples:", paste(failed_sample_ids, collapse = ", "), "\n")
      cat("  Samples remaining:", ncol(rgSet), "\n")
    }
  }
  
  probe_counts$after_step4 <- nrow(detP)
  
  # Clear memory
  rm(failed_detection_step4, failed_measurements_step4, failure_count_per_sample, failure_rate_per_sample, failed_samples)
  if (processing_log$beadcount_available) rm(failed_beadcount_step4)
  gc(verbose = FALSE)
  
  # ===== STEP 5: BISULFITE CONVERSION QC =====
  if (verbose) cat("Step 5: Bisulfite conversion QC...\n")
  
  qc <- minfi::getQC(preprocessRaw(rgSet))
  
  # Calculate median and MAD (more robust than mean and SD)
  mMed <- median(qc$mMed)
  uMed <- median(qc$uMed)
  mMad <- mad(qc$mMed)
  uMad <- mad(qc$uMed)
  
  # Identify samples outside threshold
  outliers <- (qc$mMed < (mMed - bisulfite_sd_threshold * mMad)) |
    (qc$uMed < (uMed - bisulfite_sd_threshold * uMad))
  
  n_bisulfite_failed <- sum(outliers)
  
  if (verbose) {
    cat("  Median mMed:", round(mMed, 1), "± MAD:", round(mMad, 1), "\n")
    cat("  Median uMed:", round(uMed, 1), "± MAD:", round(uMad, 1), "\n")
    cat("  Samples failing bisulfite QC:", n_bisulfite_failed, "\n")
  }
  
  # Remove failed samples
  if (n_bisulfite_failed > 0) {
    failed_bisulfite_ids <- colnames(rgSet)[outliers]
    processing_log$samples_failed_bisulfite <- failed_bisulfite_ids
    rgSet <- rgSet[, !outliers]
    detP <- detP[, !outliers, drop = FALSE]
    bc <- bc[, !outliers, drop = FALSE]
    pheno_data <- pheno_data[!outliers, , drop = FALSE]
    if (verbose) {
      cat("  Removed samples:", paste(failed_bisulfite_ids, collapse = ", "), "\n")
      cat("  Samples remaining:", ncol(rgSet), "\n")
    }
  }
  
  probe_counts$after_step5 <- nrow(detP)
  rm(qc, outliers)
  gc(verbose = FALSE)
  
  # ===== STEP 6: NOOB NORMALIZATION =====
  if (verbose) cat("Step 6: Noob normalization (background correction & dye bias)...\n")
  
  mSet <- preprocessNoob(rgSet, dyeMethod = "single")
  beta <- getBeta(mSet)
  
  # Update detP to match beta probes
  common_probes <- intersect(rownames(beta), rownames(detP))
  beta <- beta[common_probes, , drop = FALSE]
  detP <- detP[common_probes, , drop = FALSE]
  
  probe_counts$after_step6 <- nrow(beta)
  if (verbose) cat("  Probes after Noob:", nrow(beta), "\n")
  
  rm(mSet)
  gc(verbose = FALSE)
  
  # ===== STEP 7: BMIQ NORMALIZATION =====
  if (apply_bmiq) {
    if (verbose) cat("Step 7: BMIQ probe-type bias correction...\n")
    
    # Aggressive garbage collection
    for(i in 1:2) gc(verbose = FALSE, full = TRUE)
    
    tryCatch({
      if (array_type == "450k") {
        ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      } else if (array_type == "EPIC") {
        ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      } else if (array_type == "EPICv2") {
        ann <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
      }
      
      common_probes_ann <- intersect(rownames(beta), rownames(ann))
      beta_for_bmiq <- beta[common_probes_ann, , drop = FALSE]
      probe_design <- ann[common_probes_ann, "Type"]
      design_vector <- ifelse(probe_design == "I", 1, 2)
      
      if (verbose) {
        cat("  Type I probes:", sum(design_vector == 1), "\n")
        cat("  Type II probes:", sum(design_vector == 2), "\n")
        cat("  Processing", ncol(beta_for_bmiq), "samples...\n")
      }
      
      beta_bmiq <- matrix(NA, nrow = nrow(beta_for_bmiq), ncol = ncol(beta_for_bmiq),
                          dimnames = dimnames(beta_for_bmiq))
      bmiq_success_count <- 0
      
      for (i in 1:ncol(beta_for_bmiq)) {
        if (verbose && i %% 10 == 0) {
          cat("  Sample", i, "of", ncol(beta_for_bmiq), "\n")
        }
        
        # Periodic garbage collection
        if (i %% 50 == 0) {
          gc(verbose = FALSE)
        }
        
        sample_beta <- beta_for_bmiq[, i]
        valid_idx <- !is.na(sample_beta)
        sample_beta_clean <- sample_beta[valid_idx]
        design_clean <- design_vector[valid_idx]
        
        if (sum(design_clean == 1) < 10 || sum(design_clean == 2) < 10) {
          if (verbose) cat("  Warning: Insufficient probes -",
                          colnames(beta_for_bmiq)[i], "\n")
          processing_log$bmiq_failed_samples <- c(processing_log$bmiq_failed_samples,
                                                   colnames(beta_for_bmiq)[i])
          beta_bmiq[, i] <- sample_beta
          next
        }
        
        sample_beta_clean[sample_beta_clean <= 0] <- 0.0001
        sample_beta_clean[sample_beta_clean >= 1] <- 0.9999
        
        bmiq_result <- tryCatch({
          BMIQ(sample_beta_clean, design_clean, plots = FALSE, nfit = 10000)
        }, error = function(e) {
          if (verbose) cat("  BMIQ failed:", colnames(beta_for_bmiq)[i], "\n")
          processing_log$bmiq_failed_samples <<- c(processing_log$bmiq_failed_samples,
                                                   colnames(beta_for_bmiq)[i])
          list(nbeta = sample_beta_clean)
        })
        
        corrected_values <- sample_beta
        corrected_values[valid_idx] <- bmiq_result$nbeta
        beta_bmiq[, i] <- corrected_values
        
        if (!identical(bmiq_result$nbeta, sample_beta_clean)) {
          bmiq_success_count <- bmiq_success_count + 1
        }
      }
      
      beta <- beta_bmiq
      processing_log$bmiq_applied <- TRUE
      processing_log$bmiq_success_count <- bmiq_success_count
      
      if (verbose) {
        cat("  BMIQ completed:", bmiq_success_count, "samples corrected\n")
        if (length(processing_log$bmiq_failed_samples) > 0) {
          cat("  Failed:", length(processing_log$bmiq_failed_samples), "samples\n")
        }
      }
      
    }, error = function(e) {
      if (verbose) cat("  BMIQ failed:", e$message, "\n")
      processing_log$bmiq_applied <- FALSE
    })
  } else {
    if (verbose) cat("Step 7: BMIQ skipped\n")
  }
  
  probe_counts$after_step7 <- nrow(beta)
  
  # ===== STEP 8: COMBAT BATCH CORRECTION =====
  cell_proportions <- NULL
  
  if (apply_combat) {
    if (verbose) cat("Step 8: ComBat batch correction...\n")
    
    batch <- as.factor(pheno_data[[batch_variable]])
    processing_log$n_batches <- length(unique(batch))
    
    if (verbose) cat("  Number of batches:", processing_log$n_batches, "\n")
    
    if (processing_log$n_batches > 1) {
      # Create model matrix (no covariates in this version)
      mod <- model.matrix(~ 1, data = pheno_data)
      
      # Apply ComBat
      beta_combat <- ComBat(dat = beta, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
      
      # Ensure beta values stay in [0,1]
      beta_combat[beta_combat < 0] <- 0
      beta_combat[beta_combat > 1] <- 1
      
      beta <- beta_combat
      processing_log$combat_applied <- TRUE
      
      if (verbose) cat("  ComBat correction complete\n")
      
      rm(beta_combat, mod, batch)
      gc(verbose = FALSE)
    } else {
      if (verbose) cat("  Only one batch detected, skipping ComBat\n")
      processing_log$combat_applied <- FALSE
    }
  } else {
    if (verbose) cat("Step 8: ComBat skipped\n")
  }
  
  probe_counts$after_step8 <- nrow(beta)
  
  # ===== STEP 9: ESTIMATE CELL TYPE PROPORTIONS =====
  # Cell type estimation using estimateCellProp from ENmix (BEFORE removing sex probes)
  cell_proportions <- NULL
  if (estimate_cell_proportions) {
    if (verbose) cat("Step 9: Estimating cell type proportions using estimateCellProp...\n")
    
    tryCatch({
      if (verbose) cat("  Using reference:", cell_reference, "\n")
      if (!is.null(cell_types) && length(cell_types) > 0) {
        if (verbose) cat("  Cell types to estimate:", paste(cell_types, collapse = ", "), "\n")
      } else {
        if (verbose) cat("  Cell types to estimate: all available cell types\n")
      }
      
      # Call estimateCellProp with reference name as string
      # Beta values are NOT manipulated - passed as-is
      # normalize=FALSE to avoid errors about rgDataSet/methDataSet requirement
      cell_counts <- estimateCellProp(
        userdata = beta,
        refdata = cell_reference,
        cellTypes = cell_types,
        nonnegative = TRUE,
        normalize = FALSE,
        nProbes = 50,
        refplot = FALSE
      )
      
      # Handle the result - extract only numeric columns (cell type proportions)
      # estimateCellProp may return mixed types, so we need to be careful
      if (is.matrix(cell_counts)) {
        # If it's a matrix, convert directly to data frame
        cell_proportions <- as.data.frame(cell_counts)
      } else if (is.data.frame(cell_counts)) {
        # If it's a data frame with mixed types, extract only numeric columns
        numeric_cols <- sapply(cell_counts, is.numeric)
        if (sum(numeric_cols) == 0) {
          stop("No numeric columns found in cell proportion result")
        }
        cell_proportions <- cell_counts[, numeric_cols, drop = FALSE]
      } else {
        # Try to coerce to data frame and handle non-numeric columns
        cell_proportions <- as.data.frame(cell_counts)
        numeric_cols <- sapply(cell_proportions, is.numeric)
        cell_proportions <- cell_proportions[, numeric_cols, drop = FALSE]
      }
      
      # Ensure rownames are sample identifiers from beta matrix
      rownames(cell_proportions) <- colnames(beta)
      
      # Add Sample_ID column
      cell_proportions$Sample_ID <- rownames(cell_proportions)
      
      # Reorder columns: cell types first, then Sample_ID
      cell_cols <- setdiff(colnames(cell_proportions), "Sample_ID")
      cell_proportions <- cell_proportions[, c(cell_cols, "Sample_ID")]
      
      processing_log$cell_proportions_estimated <- TRUE
      processing_log$cell_reference_used <- cell_reference
      
      if (verbose) cat("  Cell type estimation completed\n")
      if (verbose) cat("  Cell types estimated:", paste(cell_cols, collapse = ", "), "\n")
      if (verbose) cat("  Dimensions:", nrow(cell_proportions), "samples ×", length(cell_cols), "cell types\n")
      
      if (nrow(cell_proportions) > 0 && length(cell_cols) > 0) {
        if (verbose) cat("  Min proportion:", round(min(cell_proportions[, cell_cols], na.rm = TRUE), 4), "\n")
        if (verbose) cat("  Max proportion:", round(max(cell_proportions[, cell_cols], na.rm = TRUE), 4), "\n")
      }
      
    }, error = function(e) {
      if (verbose) cat("  Warning: Cell type estimation failed:", e$message, "\n")
      if (verbose) cat("  Attempting to continue without cell proportions...\n")
      processing_log$cell_proportions_estimated <<- FALSE
    })
  } else {
    if (verbose) cat("Step 9: Cell type estimation skipped\n")
  }
  
  # ===== STEP 10: REMOVE SEX CHROMOSOME PROBES =====
  if (verbose) cat("Step 10: Removing sex chromosome probes...\n")
  
  if (array_type == "450k") {
    anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } else if (array_type == "EPIC") {
    anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  } else if (array_type == "EPICv2") {
    anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  }
  
  common_probes <- intersect(rownames(beta), rownames(anno))
  anno_subset <- anno[common_probes, ]
  autosomes <- anno_subset$chr %in% paste0("chr", 1:22)
  
  beta <- beta[common_probes[autosomes], , drop = FALSE]
  
  # Update detP to match remaining probes
  detP <- detP[rownames(beta), , drop = FALSE]
  
  probe_counts$after_step10 <- nrow(beta)
  if (verbose) cat("  Removed", sum(!autosomes), "sex chromosome probes\n")
  
  rm(anno, anno_subset, common_probes, autosomes)
  gc(verbose = FALSE)
  
  # ===== STEP 11: FILTER TO KEEP-LIST PROBES ONLY =====
  # Only probes present in cpgs_to_keep_file are retained.
  if (verbose) cat("Step 11: Filtering to keep-list probes (SNPs, cross-reactive probes excluded)...\n")
  
  tryCatch({
    cpgs_to_keep <- read.csv(cpgs_to_keep_file, header = TRUE, stringsAsFactors = FALSE)$cpg
    
    if (verbose) cat("  Keep list contains", length(cpgs_to_keep), "probes\n")
    
    probes_before <- nrow(beta)
    keep_probes <- rownames(beta) %in% cpgs_to_keep   # <-- CHANGED: keep probes IN the list
    removed_count <- sum(!keep_probes)
    
    beta <- beta[keep_probes, , drop = FALSE]
    
    # Update detP to match remaining probes
    detP <- detP[rownames(beta), , drop = FALSE]
    
    probe_counts$after_step11 <- nrow(beta)
    if (verbose) cat("  Probes retained:", nrow(beta), "of", probes_before,
                     "(removed", removed_count, "not in keep list)\n")
    
  }, error = function(e) {
    stop("Error reading keep list: ", e$message)
  })
  
  probe_counts$final <- nrow(beta)
  
  # ===== STEP 12: COMPUTE M-VALUES =====
  if (compute_mvalues) {
    if (verbose) cat("Step 12: Computing M-values...\n")
    
    beta_bounded <- beta
    beta_bounded[beta_bounded <= 0] <- 1e-6
    beta_bounded[beta_bounded >= 1] <- 1 - 1e-6
    
    mvalues <- log2(beta_bounded / (1 - beta_bounded))
    
    rm(beta_bounded)
    gc(verbose = FALSE)
  } else {
    mvalues <- NULL
    if (verbose) cat("Step 12: M-values skipped\n")
  }
  
  # ===== CREATE GENOMIC RATIO SET =====
  if (verbose) cat("Creating GenomicRatioSet...\n")
  
  beta_final <- beta
  detP_final <- detP
  
  if (array_type == "450k") {
    grSet <- makeGenomicRatioSetFromMatrix(beta_final, array = "IlluminaHumanMethylation450k",
                                           annotation = "ilmn12.hg19", what = "Beta")
  } else if (array_type == "EPIC") {
    grSet <- makeGenomicRatioSetFromMatrix(beta_final, array = "IlluminaHumanMethylationEPIC",
                                           annotation = "ilm10b4.hg19", what = "Beta")
  } else if (array_type == "EPICv2") {
    grSet <- makeGenomicRatioSetFromMatrix(beta_final, array = "IlluminaHumanMethylationEPICv2",
                                           annotation = "20a1.hg38", what = "Beta")
  }
  
  pheno_data <- pheno_data[colnames(grSet), , drop = FALSE]
  colData(grSet) <- S4Vectors::DataFrame(pheno_data)
  
  # ===== SAVE RESULTS =====
  if (verbose) cat("Saving results...\n")
  
  save(grSet, file = file.path(output_folder, "grSet.rda"))
  save(beta_final, file = file.path(output_folder, "beta.rda"))
  if (!is.null(mvalues)) save(mvalues, file = file.path(output_folder, "mvalues.rda"))
  save(detP_final, file = file.path(output_folder, "detP.rda"))
  save(pheno_data, file = file.path(output_folder, "pheno_data.rda"))
  
  if (!is.null(cell_proportions)) {
    save(cell_proportions, file = file.path(output_folder, "cell_proportions.rda"))
    write.csv(cell_proportions, file = file.path(output_folder, "cell_proportions.csv"), row.names = TRUE)
  }
  
  save(probe_counts, file = file.path(output_folder, "probe_counts.rda"))
  save(processing_log, file = file.path(output_folder, "processing_log.rda"))
  
  # ===== GENERATE SUMMARY REPORT =====
  summary_file <- file.path(output_folder, "processing_summary.txt")
  sink(summary_file)
  
  cat("DNA METHYLATION PROCESSING SUMMARY\n")
  cat("===================================\n\n")
  cat("Processing Date:", as.character(Sys.time()), "\n\n")
  
  cat("PARAMETERS\n")
  cat("----------\n")
  cat("Array Type:", array_type, "\n")
  cat("Detection p-value threshold:", detection_pval, "\n")
  cat("Bead count threshold:", beadcount_threshold, "\n")
  cat("Bead count available:", processing_log$beadcount_available, "\n")
  cat("Sample QC cutoff:", sample_cutoff * 100, "%\n")
  cat("Probe detection rate:", probe_detection_rate * 100, "%\n")
  cat("Bisulfite SD threshold:", bisulfite_sd_threshold, "\n")
  cat("BMIQ applied:", processing_log$bmiq_applied, "\n")
  if (processing_log$bmiq_applied) {
    cat("  Success:", processing_log$bmiq_success_count, "samples\n")
    cat("  Failed:", length(processing_log$bmiq_failed_samples), "samples\n")
  }
  cat("ComBat applied:", processing_log$combat_applied, "\n")
  if (processing_log$combat_applied) {
    cat("  Number of batches:", processing_log$n_batches, "\n")
  }
  cat("Cell proportions estimated:", processing_log$cell_proportions_estimated, "\n")
  if (processing_log$cell_proportions_estimated) {
    cat("  Reference:", processing_log$cell_reference_used, "\n")
  }
  cat("M-values computed:", compute_mvalues, "\n\n")
  
  cat("PROCESSING STEPS\n")
  cat("----------------\n")
  cat("1. Read IDAT files (extended for bead count)\n")
  cat("2. Calculate detection p-values\n")
  cat("3. Remove failed probes (detection p >", detection_pval, "OR beads <", beadcount_threshold, "in >=", sample_cutoff * 100, "% of samples)\n", sep = "")
  cat("4. Sample QC - Remove samples with poor overall quality\n")
  cat("5. Bisulfite conversion QC (Mean -", bisulfite_sd_threshold, "× SD)\n")
  cat("6. Noob normalization (background correction & dye bias)\n")
  cat("7. BMIQ probe-type correction (Type I vs II)\n")
  cat("8. ComBat batch correction\n")
  cat("9. Cell type proportion estimation (AFTER ComBat - using batch-corrected data)\n")
  cat("10. Remove sex chromosome probes\n")
  cat("11. Filter to keep-list probes only (retain probes in cpgs_to_keep_file)\n")
  cat("12. Compute M-values: M = log2(beta / (1-beta))\n\n")
  
  cat("PROBE COUNTS\n")
  cat("------------\n")
  cat("Initial probes (CpG sites):", probe_counts$initial, "\n")
  cat("After step 3 (failed probes):", probe_counts$after_step3,
      " (-", probe_counts$initial - probe_counts$after_step3, ")\n")
  cat("After step 4 (sample QC):", probe_counts$after_step4, "\n")
  cat("After step 5 (bisulfite QC):", probe_counts$after_step5, "\n")
  cat("After step 6 (Noob):", probe_counts$after_step6, "\n")
  cat("After step 7 (BMIQ):", probe_counts$after_step7, "\n")
  cat("After step 8 (ComBat):", probe_counts$after_step8, "\n")
  cat("After step 10 (sex chr):", probe_counts$after_step10,
      " (-", probe_counts$after_step8 - probe_counts$after_step10, ")\n")
  cat("After step 11 (keep list):", probe_counts$after_step11,
      " (-", probe_counts$after_step10 - probe_counts$after_step11, "not in keep list)\n")
  cat("Final probes:", probe_counts$final, "\n")
  cat("Total removed:", probe_counts$initial - probe_counts$final,
      " (", round((probe_counts$initial - probe_counts$final) / probe_counts$initial * 100, 1), "%)\n\n")
  
  cat("SAMPLE COUNTS\n")
  cat("-------------\n")
  cat("Initial samples:", nrow(basenames_df), "\n")
  cat("Removed (QC):", length(processing_log$samples_removed_qc), "\n")
  cat("Removed (bisulfite):", length(processing_log$samples_failed_bisulfite), "\n")
  cat("Final samples:", ncol(grSet), "\n")
  cat("Total removed:", nrow(basenames_df) - ncol(grSet), "\n\n")
  
  if (length(processing_log$samples_removed_qc) > 0) {
    cat("SAMPLES REMOVED (QC)\n")
    cat("--------------------\n")
    cat(paste(processing_log$samples_removed_qc, collapse = ", "), "\n\n")
  }
  
  if (length(processing_log$samples_failed_bisulfite) > 0) {
    cat("SAMPLES REMOVED (BISULFITE)\n")
    cat("---------------------------\n")
    cat(paste(processing_log$samples_failed_bisulfite, collapse = ", "), "\n\n")
  }
  
  cat("FINAL DATASET\n")
  cat("-------------\n")
  cat("Dimensions:", nrow(grSet), "probes ×", ncol(grSet), "samples\n")
  cat("Probes retained:", round(nrow(grSet) / probe_counts$initial * 100, 1), "%\n")
  cat("Samples retained:", round(ncol(grSet) / nrow(basenames_df) * 100, 1), "%\n\n")
  
  cat("BETA VALUE STATISTICS\n")
  cat("---------------------\n")
  cat("Min:", round(min(beta_final, na.rm = TRUE), 4), "\n")
  cat("Max:", round(max(beta_final, na.rm = TRUE), 4), "\n")
  cat("Mean:", round(mean(beta_final, na.rm = TRUE), 4), "\n")
  cat("Median:", round(median(beta_final, na.rm = TRUE), 4), "\n\n")
  
  if (!is.null(mvalues)) {
    cat("M-VALUE STATISTICS\n")
    cat("------------------\n")
    cat("Min:", round(min(mvalues, na.rm = TRUE), 2), "\n")
    cat("Max:", round(max(mvalues, na.rm = TRUE), 2), "\n")
    cat("Mean:", round(mean(mvalues, na.rm = TRUE), 2), "\n")
    cat("Median:", round(median(mvalues, na.rm = TRUE), 2), "\n")
  }
  
  sink()
  
  if (verbose) cat("\nProcessing complete!\n")
  if (verbose) cat("Summary saved to:", summary_file, "\n")
  if (!is.null(cell_proportions)) {
    if (verbose) cat("Cell proportions saved to: cell_proportions.rda and cell_proportions.csv\n")
  }
  
  invisible(list(
    grSet = grSet,
    beta = beta_final,
    mvalues = mvalues,
    detP = detP_final,
    pheno_data = pheno_data,
    cell_proportions = cell_proportions,
    output_folder = output_folder,
    probe_counts = probe_counts,
    processing_log = processing_log
  ))
}
