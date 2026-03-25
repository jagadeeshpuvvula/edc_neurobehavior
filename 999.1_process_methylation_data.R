
#' DNA Methylation Preprocessing Pipeline - Publication Ready (No Sample Removal)
#' Removes failed probes only. Logs QC issues (including bisulfite outliers).
#' Compatible with 450k/EPIC/EPICv2. Follows minfi/ENmix best practices (2026).
#'
#' NOTE: This version keeps all samples for maximum statistical power and
#'       documents outliers in `processing_log$bisulfite_outliers`.
#'
#' @param basenames_df Data frame with Sample_ID, folder_location
#' @param pheno_data Optional phenotype data (rows aligned to basenames_df)
#' @param output_folder Output directory
#' @param array_type "450k", "EPIC", or "EPICv2"
#' @param detection_pval Detection P-value threshold (default 1e-6)
#' @param beadcount_threshold Minimum bead count per probe (default 3)
#' @param sample_cutoff Fraction of samples failing for probe removal (default 0.05)
#' @param bisulfite_mad_threshold Bisulfite outlier threshold in MAD units (default 15)
#'   [FIX #6: renamed from bisulfite_sd_threshold - this threshold uses MAD, not SD]
#' @param cpgs_to_keep_file CSV with column `cpg` listing CpG IDs to keep
#' @param compute_mvalues Logical, compute M-values (default FALSE)
#' @param apply_bmiq Logical, apply BMIQ (default TRUE)
#' @param apply_combat Logical, apply ComBat (default TRUE)
#' @param batch_variable Batch column name in pheno_data (default "batch")
#' @param estimate_cell_proportions Logical, estimate cell proportions (default TRUE)
#' @param cell_reference ENmix reference dataset name (default "FlowSorted.CordBlood.450k")
#' @param cell_types Vector of cell types to estimate (NULL = all)
#' @param verbose Logical, print progress (default TRUE)
#'
#' @return List with grSet, beta, mvalues, detP, pheno_data, cell_proportions,
#'         probe_counts, processing_log.
#'
#' FIXES APPLIED:
#'   FIX #1 (Bug)   - Step 4: colnames(qc)[outliers] -> rownames(qc)[outliers]
#'                    getQC() returns a DataFrame with rows=samples, cols=mMed/uMed.
#'   FIX #2 (Bug)   - Step 2: beadcount fallback matrix now uses nrow/ncol/dimnames
#'                    from detP instead of rgSet. rgSet has raw R/G channel probe
#'                    dimensions which do not match detP (CpG-level), causing a
#'                    dimension mismatch in the failed_bc computation.
#'   FIX #3 (Bug)   - Step 5: colnames(cell_proportions)[numeric_cols] ->
#'                    colnames(cell_proportions). cell_proportions is already subset
#'                    to numeric columns; re-indexing with the original logical vector
#'                    (wrong length) would throw an error or return NAs.
#'   FIX #4 (Risk)  - Package loading: require() -> library() with explicit
#'                    requireNamespace() check. require() silently returns FALSE on
#'                    missing packages, causing cryptic downstream errors.
#'   FIX #5 (Minor) - Step 2 verbose message: "> X% samples" -> ">= X% samples"
#'                    to match the actual >= operator used in failed_probes.
#'   FIX #6 (Minor) - Parameter rename: bisulfite_sd_threshold -> bisulfite_mad_threshold
#'                    to accurately reflect that MAD (not SD) is used internally.

process_methylation_data <- function(
    basenames_df,
    pheno_data = NULL,
    output_folder,
    array_type = c("450k", "EPIC", "EPICv2"),
    detection_pval = 1e-6,
    beadcount_threshold = 3,
    sample_cutoff = 0.05,
    bisulfite_mad_threshold = 15,           # FIX #6: renamed from bisulfite_sd_threshold
    cpgs_to_keep_file,
    compute_mvalues = FALSE,
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
  
  # ==== FIX #4: Load packages with library() + explicit missing-package check =
  # require() silently returns FALSE on missing packages; library() stops immediately
  # with an informative error, which is safer in a production pipeline.
  pkgs <- c(
    "minfi",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    "wateRmelon",
    "readr",
    "SummarizedExperiment",
    "S4Vectors",
    "sva",
    "ENmix",
    "FlowSorted.CordBlood.450k"
  )
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Required package '", pkg, "' is not installed. ",
        "Install it with: BiocManager::install('", pkg, "')",
        call. = FALSE
      )
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  array_type <- match.arg(array_type)
  
  probe_counts <- list()
  processing_log <- list(
    array_type = array_type,
    n_samples_initial = nrow(basenames_df),
    probes_failed_qc = 0,
    beadcount_available = TRUE,
    bisulfite_outliers = character(0),
    bmiq_applied = FALSE,
    bmiq_failed_samples = character(0),
    bmiq_success_count = 0,
    combat_applied = FALSE,
    n_batches = 0,
    cell_proportions_estimated = FALSE,
    cell_reference_used = NA_character_
  )
  
  # ==== INPUT VALIDATION ======================================================
  required_cols <- c("Sample_ID", "folder_location")
  if (!is.data.frame(basenames_df)) {
    stop("basenames_df must be a data frame.")
  }
  if (!all(required_cols %in% colnames(basenames_df))) {
    stop("basenames_df must contain columns: Sample_ID, folder_location")
  }
  if (anyDuplicated(basenames_df$Sample_ID)) {
    dup_ids <- unique(basenames_df$Sample_ID[duplicated(basenames_df$Sample_ID)])
    stop("Duplicate Sample_IDs found: ", paste(dup_ids, collapse = ", "))
  }
  
  basenames <- file.path(basenames_df$folder_location, basenames_df$Sample_ID)
  
  # Check IDAT files
  missing_files <- character(0)
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
  
  # Quick array-type check by file size (defensive only)
  array_versions <- sapply(basenames, function(bn) {
    red_file <- paste0(bn, "_Red.idat")
    fs <- file.info(red_file)$size
    if (fs > 15000000) "EPICv2" else if (fs > 7000000) "EPIC" else "450k"
  })
  if (length(unique(array_versions)) > 1) {
    stop("Mixed array types detected. Process each type separately.")
  }
  
  # Output directory
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  # Phenotype data
  if (is.null(pheno_data)) {
    pheno_data <- as.data.frame(basenames_df)
    rownames(pheno_data) <- basenames_df$Sample_ID
  } else {
    pheno_data <- as.data.frame(pheno_data)
    rownames(pheno_data) <- basenames_df$Sample_ID
  }
  
  if (apply_combat && !batch_variable %in% colnames(pheno_data)) {
    warning("Batch variable '", batch_variable, "' not found in pheno_data. Skipping ComBat.")
    apply_combat <- FALSE
  }
  
  # ==== STEP 1: READ IDAT FILES ==============================================
  if (verbose) cat("Step 1: Reading IDAT files...\n")
  
  rgSet <- tryCatch({
    read.metharray(basenames = basenames, extended = TRUE)
  }, error = function(e) {
    if (grepl("different array size", e$message)) {
      read.metharray(basenames = basenames, extended = TRUE, force = TRUE)
    } else {
      stop(e)
    }
  })
  colnames(rgSet) <- basenames_df$Sample_ID
  
  # Initial probe count from annotation
  anno_pkg <- switch(
    array_type,
    "450k"   = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "EPIC"   = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "EPICv2" = "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
  )
  anno_obj <- get(anno_pkg)
  probe_counts$initial <- nrow(getAnnotation(anno_obj))
  
  if (verbose) {
    cat("  Initial probes:", probe_counts$initial, " | Samples:", ncol(rgSet), "\n")
  }
  
  # ==== STEP 2: PROBE QC (DETECTION P + BEADCOUNT), NO SAMPLE DROPS ==========
  if (verbose) cat("Step 2: Probe QC (detection/beadcount)...\n")
  
  detP <- detectionP(rgSet)
  colnames(detP) <- colnames(rgSet)
  
  # FIX #2: Beadcount fallback matrix must match detP dimensions (CpG-level),
  # NOT rgSet dimensions (raw R/G channel probes). The original code used
  # nrow/ncol/dimnames(rgSet), which has a different (larger) row count and
  # caused a dimension mismatch when computing failed_bc | failed_det below.
  bc <- tryCatch({
    beadcount(rgSet)
  }, error = function(e) {
    if (verbose) cat("  Warning: beadcount data not available; using detectionP only.\n")
    processing_log$beadcount_available <<- FALSE
    matrix(
      NA,
      nrow     = nrow(detP),
      ncol     = ncol(detP),
      dimnames = dimnames(detP)           # FIX #2: was dimnames(rgSet)
    )
  })
  
  failed_det <- detP > detection_pval
  if (processing_log$beadcount_available) {
    failed_bc  <- bc < beadcount_threshold | is.na(bc)
    failed_any <- failed_det | failed_bc
  } else {
    failed_any <- failed_det
  }
  
  fail_rate_per_probe <- rowSums(failed_any) / ncol(rgSet)
  failed_probes <- fail_rate_per_probe >= sample_cutoff
  
  n_failed <- sum(failed_probes)
  processing_log$probes_failed_qc <- n_failed
  
  if (n_failed > 0) {
    rgSet <- rgSet[!failed_probes, ]
    detP  <- detP[!failed_probes, , drop = FALSE]
    bc    <- bc[!failed_probes,   , drop = FALSE]
  }
  
  probe_counts$after_probe_qc <- nrow(detP)
  
  if (verbose) {
    # FIX #5: message now says ">=" to match the >= operator used in failed_probes
    cat("  Removed", n_failed, "failed probes (>=", sample_cutoff * 100, "% samples)\n")
    cat("  Probes remaining:", nrow(detP), "\n")
  }
  
  rm(failed_det, failed_any, fail_rate_per_probe, failed_probes)
  if (processing_log$beadcount_available) rm(failed_bc)
  gc(verbose = FALSE)
  
  # ==== STEP 3: NOOB NORMALIZATION ===========================================
  if (verbose) cat("Step 3: Noob normalization...\n")
  
  mSet <- preprocessNoob(rgSet, dyeMethod = "single")
  beta <- getBeta(mSet)
  
  # Sync detP with beta
  common_probes <- intersect(rownames(beta), rownames(detP))
  beta <- beta[common_probes, , drop = FALSE]
  detP <- detP[common_probes, , drop = FALSE]
  
  probe_counts$after_noob <- nrow(beta)
  
  # NOTE: We keep mSet for bisulfite QC; do NOT rm(mSet) yet.
  rm(rgSet)
  gc(verbose = FALSE)
  
  # ==== STEP 4: BISULFITE QC (LOGGING ONLY) ==================================
  if (verbose) cat("Step 4: Bisulfite QC (logging outliers)...\n")
  
  qc   <- getQC(mSet)
  mMed <- median(qc$mMed)
  uMed <- median(qc$uMed)
  mMad <- mad(qc$mMed)
  uMad <- mad(qc$uMed)
  
  # FIX #6: use bisulfite_mad_threshold (renamed parameter)
  outliers <- (qc$mMed < (mMed - bisulfite_mad_threshold * mMad)) |
    (qc$uMed < (uMed - bisulfite_mad_threshold * uMad))
  
  # FIX #1: getQC() returns a DataFrame with rows = samples and cols = mMed/uMed.
  # The original code used colnames(qc), which returns c("mMed","uMed"), not sample IDs.
  # rownames(qc) correctly returns the sample identifiers.
  processing_log$bisulfite_outliers <- rownames(qc)[outliers]  # FIX #1: was colnames(qc)
  
  if (verbose) {
    cat("  Median mMed:", round(mMed, 1), "+/-", round(mMad, 1), "\n")
    cat("  Median uMed:", round(uMed, 1), "+/-", round(uMad, 1), "\n")
    cat("  Bisulfite outliers (kept):", length(processing_log$bisulfite_outliers), "\n")
  }
  
  rm(qc, mSet)
  gc(verbose = FALSE)
  
  # ==== STEP 5: CELL TYPE PROPORTION ESTIMATION (POST-NOOB) ==================
  cell_proportions <- NULL
  
  if (estimate_cell_proportions) {
    if (verbose) cat("Step 5: Estimating cell type proportions...\n")
    
    tryCatch({
      cell_counts <- estimateCellProp(
        userdata    = beta,
        refdata     = cell_reference,
        cellTypes   = cell_types,
        nonnegative.g = TRUE,
        normalize   = TRUE,
        nProbes     = 50,
        refplot     = FALSE
      )
      
      numeric_cols     <- sapply(cell_counts, is.numeric)
      cell_proportions <- as.data.frame(cell_counts[, numeric_cols, drop = FALSE])
      rownames(cell_proportions) <- colnames(beta)
      cell_proportions$Sample_ID <- rownames(cell_proportions)
      
      processing_log$cell_proportions_estimated <- TRUE
      processing_log$cell_reference_used        <- cell_reference
      
      if (verbose) {
        # FIX #3: cell_proportions is already subset to numeric columns above.
        # Indexing with numeric_cols (length = ncol(cell_counts)) against a shorter
        # vector (ncol(cell_proportions) = sum(numeric_cols)) would fail or return NAs.
        # Use colnames(cell_proportions) directly, excluding the appended Sample_ID col.
        numeric_col_names <- setdiff(colnames(cell_proportions), "Sample_ID")  # FIX #3
        cat("  Cell types estimated:", paste(numeric_col_names, collapse = ", "), "\n")
      }
    }, error = function(e) {
      warning("Cell type estimation failed: ", e$message)
      processing_log$cell_proportions_estimated <<- FALSE
    })
  } else if (verbose) {
    cat("Step 5: Cell type estimation skipped.\n")
  }
  
  # ==== STEP 6: BMIQ NORMALIZATION ===========================================
  if (apply_bmiq) {
    if (verbose) cat("Step 6: BMIQ normalization...\n")
    
    ann          <- getAnnotation(anno_obj)
    common_ann   <- intersect(rownames(beta), rownames(ann))
    beta_for_bmiq <- beta[common_ann, , drop = FALSE]
    probe_design <- ann[common_ann, "Type"]
    design_v     <- ifelse(probe_design == "I", 1, 2)
    
    beta_bmiq    <- beta
    bmiq_success <- 0L
    
    for (i in seq_len(ncol(beta_for_bmiq))) {
      if (verbose && i %% 20 == 0) {
        cat("  BMIQ sample", i, "of", ncol(beta_for_bmiq), "\n")
      }
      
      sample_beta <- beta_for_bmiq[, i]
      valid_idx   <- !is.na(sample_beta)
      if (sum(design_v[valid_idx] == 1) < 10 || sum(design_v[valid_idx] == 2) < 10) {
        next
      }
      
      sample_clean <- sample_beta[valid_idx]
      sample_clean[sample_clean <= 0] <- 0.0001
      sample_clean[sample_clean >= 1] <- 0.9999
      design_clean <- design_v[valid_idx]
      
      bmiq_out <- tryCatch({
        BMIQ(sample_clean, design_clean, plots = FALSE, nfit = 10000)$nbeta
      }, error = function(e) {
        processing_log$bmiq_failed_samples <<- c(
          processing_log$bmiq_failed_samples,
          colnames(beta_for_bmiq)[i]
        )
        sample_clean
      })
      
      beta_bmiq[common_ann, i][valid_idx] <- bmiq_out
      if (!identical(unname(bmiq_out), unname(sample_clean))) {
        bmiq_success <- bmiq_success + 1L
      }
      if (i %% 50 == 0) gc(verbose = FALSE)
    }
    
    beta <- beta_bmiq
    processing_log$bmiq_applied      <- TRUE
    processing_log$bmiq_success_count <- bmiq_success
    
    if (verbose) {
      cat("  BMIQ completed;", bmiq_success, "samples successfully corrected\n")
      if (length(processing_log$bmiq_failed_samples) > 0) {
        cat("  BMIQ failures:", length(processing_log$bmiq_failed_samples), "samples\n")
      }
    }
    
    rm(beta_bmiq, ann, common_ann, probe_design, design_v)
    gc(verbose = FALSE)
  } else if (verbose) {
    cat("Step 6: BMIQ skipped.\n")
  }
  
  probe_counts$after_bmiq <- nrow(beta)
  
  # ==== STEP 7: COMBAT BATCH CORRECTION ======================================
  if (apply_combat) {
    if (verbose) cat("Step 7: ComBat batch correction...\n")
    
    batch <- as.factor(pheno_data[[batch_variable]])
    processing_log$n_batches <- length(unique(batch))
    
    if (processing_log$n_batches > 1) {
      mod         <- model.matrix(~ 1, data = pheno_data)
      beta_combat <- ComBat(dat = beta, batch = batch, mod = mod,
                            par.prior = TRUE, prior.plots = FALSE)
      beta_combat[beta_combat < 0] <- 0
      beta_combat[beta_combat > 1] <- 1
      beta <- beta_combat
      processing_log$combat_applied <- TRUE
      if (verbose) cat("  ComBat correction applied on", processing_log$n_batches, "batches\n")
      rm(beta_combat, mod, batch)
      gc(verbose = FALSE)
    } else if (verbose) {
      cat("  Only one batch detected; ComBat skipped.\n")
    }
  } else if (verbose) {
    cat("Step 7: ComBat skipped.\n")
  }
  
  probe_counts$after_combat <- nrow(beta)
  
  # ==== STEP 8: REMOVE SEX CHROMOSOME PROBES =================================
  if (verbose) cat("Step 8: Removing sex chromosome probes...\n")
  
  ann2          <- getAnnotation(anno_obj)
  common_probes2 <- intersect(rownames(beta), rownames(ann2))
  autosomes     <- ann2[common_probes2, "chr"] %in% paste0("chr", 1:22)
  
  beta <- beta[common_probes2[autosomes], , drop = FALSE]
  detP <- detP[rownames(beta), , drop = FALSE]
  
  probe_counts$after_sex_remove <- nrow(beta)
  if (verbose) {
    cat("  Removed", sum(!autosomes), "sex chromosome probes\n")
  }
  
  rm(ann2, common_probes2, autosomes)
  gc(verbose = FALSE)
  
  # ==== STEP 9: KEEP-LIST FILTERING ==========================================
  if (verbose) cat("Step 9: Filtering to keep-list probes...\n")
  
  cpgs_to_keep_df <- readr::read_csv(cpgs_to_keep_file, show_col_types = FALSE)
  if (!"cpg" %in% colnames(cpgs_to_keep_df)) {
    stop("cpgs_to_keep_file must contain a column named 'cpg'")
  }
  cpgs_to_keep <- cpgs_to_keep_df$cpg
  
  keep_mask    <- rownames(beta) %in% cpgs_to_keep
  removed_keep <- sum(!keep_mask)
  
  beta <- beta[keep_mask, , drop = FALSE]
  detP <- detP[keep_mask, , drop = FALSE]
  
  probe_counts$final <- nrow(beta)
  
  if (verbose) {
    cat("  Probes retained:", nrow(beta), "of", probe_counts$after_sex_remove,
        "(removed", removed_keep, "not in keep-list)\n")
  }
  
  # ==== STEP 10: M-VALUES (OPTIONAL) =========================================
  mvalues <- NULL
  if (compute_mvalues) {
    if (verbose) cat("Step 10: Computing M-values...\n")
    beta_bounded <- beta
    beta_bounded[beta_bounded <= 0] <- 1e-6
    beta_bounded[beta_bounded >= 1] <- 1 - 1e-6
    mvalues <- log2(beta_bounded / (1 - beta_bounded))
    rm(beta_bounded)
    gc(verbose = FALSE)
  } else if (verbose) {
    cat("Step 10: M-values skipped.\n")
  }
  
  # ==== STEP 11: CREATE GENOMIC RATIO SET & SAVE =============================
  if (verbose) cat("Step 11: Creating GenomicRatioSet and saving outputs...\n")
  
  ann_tag <- switch(array_type,
                    "450k"   = "ilmn12.hg19",
                    "EPIC"   = "ilm10b4.hg19",
                    "EPICv2" = "20a1.hg38"
  )
  
  grSet <- makeGenomicRatioSetFromMatrix(
    beta,
    array      = paste0("IlluminaHumanMethylation", array_type),
    annotation = ann_tag,
    what       = "Beta"
  )
  
  pheno_data <- pheno_data[colnames(grSet), , drop = FALSE]
  colData(grSet) <- S4Vectors::DataFrame(pheno_data)
  
  # Save core objects
  save(grSet,       file = file.path(output_folder, "grSet.rda"))
  save(beta,        file = file.path(output_folder, "beta.rda"))
  if (!is.null(mvalues)) save(mvalues, file = file.path(output_folder, "mvalues.rda"))
  save(detP,        file = file.path(output_folder, "detP.rda"))
  save(pheno_data,  file = file.path(output_folder, "pheno_data.rda"))
  
  if (!is.null(cell_proportions)) {
    save(cell_proportions, file = file.path(output_folder, "cell_proportions.rda"))
    utils::write.csv(cell_proportions,
                     file = file.path(output_folder, "cell_proportions.csv"),
                     row.names = TRUE)
  }
  
  save(probe_counts,   file = file.path(output_folder, "probe_counts.rda"))
  save(processing_log, file = file.path(output_folder, "processing_log.rda"))
  
  # Summary report
  summary_file <- file.path(output_folder, "processing_summary.txt")
  sink(summary_file)
  cat("DNA METHYLATION PROCESSING SUMMARY (NO SAMPLE REMOVAL)\n")
  cat("====================================================\n\n")
  cat("Date:", as.character(Sys.time()), "\n\n")
  
  cat("PARAMETERS\n")
  cat("----------\n")
  cat("Array type:", array_type, "\n")
  cat("Detection p-value threshold:", detection_pval, "\n")
  cat("Beadcount threshold:", beadcount_threshold, "\n")
  cat("Probe failure cutoff (fraction of samples):", sample_cutoff, "\n")
  cat("Bisulfite MAD threshold:", bisulfite_mad_threshold, "\n")  # FIX #6: renamed
  cat("BMIQ applied:", processing_log$bmiq_applied, "\n")
  cat("ComBat applied:", processing_log$combat_applied, "\n")
  cat("Cell proportions estimated:", processing_log$cell_proportions_estimated, "\n")
  cat("M-values computed:", compute_mvalues, "\n\n")
  
  cat("PROBE COUNTS\n")
  cat("------------\n")
  cat("Initial probes:", probe_counts$initial, "\n")
  cat("After probe QC:", probe_counts$after_probe_qc, "\n")
  cat("After Noob:", probe_counts$after_noob, "\n")
  cat("After BMIQ:", probe_counts$after_bmiq, "\n")
  cat("After ComBat:", probe_counts$after_combat, "\n")
  cat("After sex removal:", probe_counts$after_sex_remove, "\n")
  cat("Final (keep-list):", probe_counts$final, "\n")
  cat("Total removed:", probe_counts$initial - probe_counts$final, " (",
      round((probe_counts$initial - probe_counts$final) / probe_counts$initial * 100, 1), "%)\n\n")
  
  cat("SAMPLE COUNTS\n")
  cat("-------------\n")
  cat("Initial samples:", processing_log$n_samples_initial, "\n")
  cat("Final samples (no removals):", ncol(grSet), "\n")
  cat("Bisulfite outliers flagged (retained):",
      length(processing_log$bisulfite_outliers), "\n\n")
  
  cat("FINAL DATASET\n")
  cat("-------------\n")
  cat("Dimensions:", nrow(grSet), "probes x", ncol(grSet), "samples\n")
  cat("Probes retained:", round(nrow(grSet) / probe_counts$initial * 100, 1), "%\n")
  cat("Samples retained:",
      round(ncol(grSet) / processing_log$n_samples_initial * 100, 1), "%\n\n")
  
  sink()
  
  if (verbose) {
    cat("\nProcessing complete. Summary saved to:", summary_file, "\n")
  }
  
  invisible(list(
    grSet            = grSet,
    beta             = beta,
    mvalues          = mvalues,
    detP             = detP,
    pheno_data       = pheno_data,
    cell_proportions = cell_proportions,
    probe_counts     = probe_counts,
    processing_log   = processing_log
  ))
}


# ==== USAGE EXAMPLE ===========================================================
#
# process_methylation_data(
#   basenames_df              = basenames_df_t,
#   pheno_data                = pheno_data_matched_t,
#   output_folder             = "~/Documents/k99_home/new",
#   array_type                = "EPIC",
#   detection_pval            = 1e-6,
#   beadcount_threshold       = 3,
#   sample_cutoff             = 0.05,
#   bisulfite_mad_threshold   = 15,   # FIX #6: renamed from bisulfite_sd_threshold
#   cpgs_to_keep_file         = "~/Documents/k99_home/arzu_cpgs_to_drop/cpgs_to_keep.csv",
#   compute_mvalues           = FALSE,
#   apply_bmiq                = TRUE,
#   apply_combat              = TRUE,
#   batch_variable            = "batch",
#   estimate_cell_proportions = TRUE,
#   verbose                   = TRUE
# )
