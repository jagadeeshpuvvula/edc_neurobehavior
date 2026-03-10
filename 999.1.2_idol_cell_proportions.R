estimate_idol_cell_proportions <- function(
  basenames_df,
  pheno_data        = NULL,
  output_folder,
  composite_cell_type = "CordBloodCombined",
  reference_platform  = "IlluminaHumanMethylationEPIC",
  cell_types          = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),
  verbose             = TRUE
) {

  # ===== PACKAGES =====
  require(minfi)
  require(FlowSorted.Blood.EPIC)
  require(ExperimentHub)

  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

  # ===== READ IDAT FILES =====
  if (verbose) cat("Reading IDAT files...\n")

  basenames <- file.path(basenames_df$folder_location, basenames_df$Sample_ID)

  rgSet <- tryCatch({
    read.metharray(basenames = basenames, extended = FALSE, verbose = verbose)
  }, error = function(e) {
    if (grepl("different array size", e$message)) {
      read.metharray(basenames = basenames, extended = FALSE, force = TRUE, verbose = verbose)
    } else stop(e)
  })
  colnames(rgSet) <- basenames_df$Sample_ID

  if (verbose) cat("  Samples loaded:", ncol(rgSet), "\n")
  if (verbose) cat("  Array annotation:", annotation(rgSet)[["array"]], "\n")

  # ===== AGGRESSIVE GC BEFORE LOADING REFERENCE =====
  for (i in 1:3) gc(verbose = FALSE, full = TRUE)

  # ===== RUN estimateCellCounts2 =====
  if (verbose) cat("Running estimateCellCounts2 with IDOL probes...\n")
  if (verbose) cat("  Composite cell type:", composite_cell_type, "\n")
  if (verbose) cat("  Reference platform:", reference_platform, "\n")
  if (verbose) cat("  Cell types:", paste(cell_types, collapse = ", "), "\n\n")

  idol_raw <- estimateCellCounts2(
    rgSet,
    compositeCellType = composite_cell_type,
    processMethod     = "preprocessNoob",
    probeSelect       = "IDOL",
    cellTypes         = cell_types,
    referencePlatform = reference_platform,
    referenceset      = NULL,
    CustomCpGs        = NULL,
    returnAll         = TRUE,   
    meanPlot          = FALSE,
    verbose           = verbose,
    lessThanOne       = TRUE
  )

  # ===== EXTRACT PROPORTIONS =====
  # estimateCellCounts2 with returnAll=TRUE returns a list with:
  #   $prop       — the proportions matrix (samples x cell types)  [primary target]
  #   $normalizedData — the combined normalised beta matrix
  # With returnAll=FALSE it returns the matrix directly.
  # Handle both cases defensively.

  if (verbose) cat("\nExtracting proportions from result object...\n")
  if (verbose) cat("  Result class:", class(idol_raw), "\n")

  if (is.list(idol_raw) && !is.data.frame(idol_raw)) {
    if (verbose) cat("  Result is a list; names:", paste(names(idol_raw), collapse = ", "), "\n")

    if ("prop" %in% names(idol_raw)) {
      prop_mat <- idol_raw$prop
    } else if ("counts" %in% names(idol_raw)) {
      prop_mat <- idol_raw$counts
    } else {
      # Fall back to first numeric matrix element
      mat_idx <- which(sapply(idol_raw, function(x) is.matrix(x) && is.numeric(x)))
      if (length(mat_idx) == 0) stop("Could not find a numeric proportions matrix in estimateCellCounts2 output.")
      prop_mat <- idol_raw[[mat_idx[1]]]
      if (verbose) cat("  Using list element:", names(idol_raw)[mat_idx[1]], "\n")
    }
  } else {
    # Returned directly as a matrix
    prop_mat <- idol_raw
  }

  if (verbose) cat("  Proportions matrix dimensions:", nrow(prop_mat), "x", ncol(prop_mat), "\n")

  # ===== BUILD OUTPUT DATA FRAME =====
  cell_proportions_idol <- as.data.frame(prop_mat)

  # Ensure rownames match sample IDs
  if (nrow(cell_proportions_idol) == ncol(rgSet)) {
    rownames(cell_proportions_idol) <- colnames(rgSet)
  } else {
    warning("Row count mismatch between proportions (", nrow(cell_proportions_idol),
            ") and samples (", ncol(rgSet), "). Rownames not reassigned.")
  }

  # Append pheno_data columns if provided
  if (!is.null(pheno_data)) {
    pheno_df <- as.data.frame(pheno_data)
    rownames(pheno_df) <- basenames_df$Sample_ID
    # Only merge rows that exist in both
    common_ids <- intersect(rownames(cell_proportions_idol), rownames(pheno_df))
    cell_proportions_idol <- cbind(
      cell_proportions_idol[common_ids, , drop = FALSE],
      pheno_df[common_ids, , drop = FALSE]
    )
  }

  # ===== SUMMARY =====
  cell_cols <- colnames(prop_mat)
  if (verbose) {
    cat("\nIDOL cell proportion estimation complete\n")
    cat("  Samples:", nrow(cell_proportions_idol), "\n")
    cat("  Cell types:", paste(cell_cols, collapse = ", "), "\n")
    cat("\n  Mean proportions per cell type:\n")
    for (ct in cell_cols) {
      cat("   ", ct, ":", round(mean(cell_proportions_idol[[ct]], na.rm = TRUE), 4), "\n")
    }
  }

  # ===== SAVE =====
  save(cell_proportions_idol, file = file.path(output_folder, "cell_proportions_idol.rda"))
  write.csv(cell_proportions_idol, file = file.path(output_folder, "cell_proportions_idol.csv"),
            row.names = TRUE)

  if (verbose) {
    cat("\nSaved:\n")
    cat("  ", file.path(output_folder, "cell_proportions_idol.rda"), "\n")
    cat("  ", file.path(output_folder, "cell_proportions_idol.csv"), "\n")
  }

  invisible(cell_proportions_idol)
}


# ===== RUN =====
cell_proportions_idol <- estimate_idol_cell_proportions(
  basenames_df        = basenames_df_t,
  pheno_data          = pheno_data_matched_t,
  output_folder       = "~/Documents/k99_home/dnam_prep_data",
  composite_cell_type = "CordBloodCombined",
  reference_platform  = "IlluminaHumanMethylationEPIC",
  cell_types          = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"),
  verbose             = TRUE
)
