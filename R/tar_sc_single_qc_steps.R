#' @title Single-Cell Quality Control Steps
#' @description A collection of functions to perform quality control (QC) steps 
#' for single-cell RNA sequencing (scRNA-seq) data. These steps include reading 
#' 10x Genomics data, detecting and removing empty droplets, calculating per-cell 
#' QC metrics, applying dataset-sensitive and custom filters, and preparing the 
#' filtered single-cell experiment (SCE) object.
#'
#' @section Functions:
#' \describe{
#'   \item{\code{tar_sc_single_qc_step_read_10x_counts(path)}}{
#'     Reads 10x Genomics data from the specified path and assigns barcodes as column names.
#'   }
#'   \item{\code{tar_sc_single_qc_step_detect_empty_droplets(sce_raw, empty_lower_droplets, BPPARAM)}}{
#'     Detects empty droplets using the \code{DropletUtils::emptyDrops} function.
#'   }
#'   \item{\code{tar_sc_single_qc_step_remove_empty_drop(sce_raw, empty_droplets, empty_droplets_fdr_threshold)}}{
#'     Removes empty droplets based on the FDR threshold and updates the SCE object.
#'   }
#'   \item{\code{tar_sc_single_qc_step_cal_per_cell_qc_metrics(sce_no_empty_drop, BPPARAM)}}{
#'     Calculates per-cell QC metrics, including mitochondrial and ribosomal gene percentages.
#'   }
#'   \item{\code{tar_sc_single_qc_step_create_unfiltered_sce(sce_no_empty_drop, per_cell_qc_metrics, replace)}}{
#'     Creates an unfiltered SCE object by applying dataset-sensitive and custom filters.
#'   }
#'   \item{\code{tar_sc_single_qc_tep_make_sensitive_filter(sce_unfiltered)}}{
#'     Applies dataset-sensitive filters to the unfiltered SCE object.
#'   }
#'   \item{\code{tar_sc_single_step_make_custom_filter(sce_unfiltered)}}{
#'     Applies custom filters to the unfiltered SCE object.
#'   }
#'   \item{\code{tar_sc_single_qc_step_after_filter_rowsums(sce_filter)}}{
#'     Computes row sums of the counts matrix for the filtered SCE object.
#'   }
#'   \item{\code{tar_sc_single_qc_step_apply_filter(sce_sensitive_filter, sce_custom_filter, save_dataset_sensitive_filtering)}}{
#'     Applies the selected filter (dataset-sensitive or custom) and removes genes 
#'     with low expression across cells.
#'   }
#' }
#'
#' @param path Character. Path to the 10x Genomics data directory.
#' @param sce_raw SingleCellExperiment. The raw SCE object.
#' @param empty_lower_droplets Numeric. Lower threshold for detecting empty droplets.
#' @param BPPARAM BiocParallelParam. Parallelization parameters.
#' @param empty_droplets DataFrame. Results from \code{DropletUtils::emptyDrops}.
#' @param empty_droplets_fdr_threshold Numeric. FDR threshold for identifying cells.
#' @param sce_no_empty_drop SingleCellExperiment. SCE object after removing empty droplets.
#' @param per_cell_qc_metrics DataFrame. Per-cell QC metrics.
#' @param replace Logical. Whether to replace existing columns in \code{colData}.
#' @param sce_unfiltered SingleCellExperiment. Unfiltered SCE object.
#' @param sce_filter SingleCellExperiment. Filtered SCE object.
#' @param sce_sensitive_filter SingleCellExperiment. SCE object after dataset-sensitive filtering.
#' @param sce_custom_filter SingleCellExperiment. SCE object after custom filtering.
#' @param save_dataset_sensitive_filtering Logical. Whether to save dataset-sensitive filtering.
#' @param min_umi Numeric. Minimum UMI count for gene filtering.
#' @param min_ratio_cells Numeric. Minimum ratio of cells expressing a gene for it to be retained.
#'
#' @return Various outputs depending on the function, including SCE objects, QC metrics, 
#' and filtered data.
#' 
#' @examples
#' # Example usage of the functions:
#' sce_raw <- tar_sc_single_qc_step_read_10x_counts("path/to/data")
#' empty_droplets <- tar_sc_single_qc_step_detect_empty_droplets(sce_raw, 100, BPPARAM = BiocParallel::SerialParam())
#' sce_no_empty_drop <- tar_sc_single_qc_step_remove_empty_drop(sce_raw, empty_droplets, 0.01)
#' per_cell_qc_metrics <- tar_sc_single_qc_step_cal_per_cell_qc_metrics(sce_no_empty_drop, BPPARAM = BiocParallel::SerialParam())
#' sce_unfiltered <- tar_sc_single_qc_step_create_unfiltered_sce(sce_no_empty_drop, per_cell_qc_metrics, replace = TRUE)
#' sce_sensitive_filter <- tar_sc_single_qc_tep_make_sensitive_filter(sce_unfiltered)
#' sce_custom_filter <- tar_sc_single_step_make_custom_filter(sce_unfiltered)
#' filtered_counts <- tar_sc_single_qc_step_after_filter_rowsums(sce_sensitive_filter)
#' sce_final <- tar_sc_single_qc_step_apply_filter(sce_sensitive_filter, sce_custom_filter, save_dataset_sensitive_filtering = TRUE)

# 1
tar_sc_single_qc_step_read_10x_counts <- function(path) {
  sce_raw <- DropletUtils::read10xCounts(path)
  colnames(sce_raw) <- colData(sce_raw)$Barcode
  sce_raw
}

# 2
tar_sc_single_qc_step_detect_empty_droplets <- function(sce_raw, empty_lower_droplets, BPPARAM) {
  DropletUtils::emptyDrops(
    m = counts(sce_raw), lower = empty_lower_droplets, BPPARAM = BPPARAM
  )
}

# 3
tar_sc_single_qc_step_remove_empty_drop <- function(sce_raw, empty_droplets, empty_droplets_fdr_threshold) {

  stopifnot(
    ncol(sce_raw) > 0
  )
  if (is.null(empty_droplets)) {
    sce_raw$is_empty_fdr <- NA
  } else {
    is_cell_index <- which(empty_droplets$FDR <= empty_droplets_fdr_threshold)
    colnames_orig <- colnames(sce_raw)
    sce_raw <- sce_raw[, is_cell_index]
    sce_raw$is_empty_fdr <- empty_droplets$FDR[is_cell_index]
  }

  sce_raw
}

# 4
tar_sc_single_qc_step_cal_per_cell_qc_metrics <- function(sce_no_empty_drop, BPPARAM) {
  mito_genes <- stringr::str_which(rowData(sce_no_empty_drop)[["Symbol"]], stringr::regex("MT-", ignore_case = TRUE))
  ribo_genes <- stringr::str_which(rowData(sce_no_empty_drop)[["Symbol"]], stringr::regex("^RP[SL]", ignore_case = TRUE))

  scater::perCellQCMetrics(
    sce_no_empty_drop,
    subsets = list(mito = mito_genes, ribo = ribo_genes), 
    BPPARAM = ignore(BPPARAM)
  )
}

# 5
tar_sc_single_qc_step_create_unfiltered_sce <- function(sce_no_empty_drop, per_cell_qc_metrics, replace) {
  # ### Dataset-sensitive cell filtering ##########################################
  # MAD_THRESHOLD: 3
  # DATASET_SENSITIVE_FILTERS_OPERATOR: "&"
  # ###############################################################################

  # ### Custom cell filtering #####################################################
  # MIN_UMI_CF: 1000
  # MAX_UMI_CF: 50000
  # MIN_FEATURES: 1000
  # MAX_MITO_RATIO: 0.2
  # CUSTOM_FILTERS_OPERATOR: "&"
  # ###############################################################################

  qc_filter <- list(
      qc_lib = scater::isOutlier(per_cell_qc_metrics$total, log = TRUE, nmads = 3, type = "lower"),
      qc_nexprs = scater::isOutlier(per_cell_qc_metrics$detected, nmads = 3, log = TRUE, type = "lower"),
      qc_mito = scater::isOutlier(per_cell_qc_metrics$subsets_mito_percent, nmads = 3, type = "higher")
      # qc_ribo = isOutlier(per_cell_qc_metrics$subsets_ribo_percent, nmads = !!cfg$MAD_THRESHOLD, type = "higher")
    ) %>%
    purrr::map(~ as.logical(.) %>% tidyr::replace_na(replace = FALSE)) %>%
    purrr::reduce(`&`)


  ## -- Custom filters.
  custom_filter <- list(
      low_count = per_cell_qc_metrics$total <= 1000,
      high_count = per_cell_qc_metrics$total >= 50000,
      low_expression = per_cell_qc_metrics$detected <= 1000,
      high_mito = per_cell_qc_metrics$subsets_mito_percent >= 0.2 * 100
      # low_ribo = per_cell_qc_metrics$subsets_ribo_percent <= !!cfg$MIN_RIBO_RATIO * 100
    ) %>%
  purrr::map( ~ as.logical(.) %>% tidyr::replace_na(replace = FALSE)) %>%
  purrr::reduce(`&`)

  ## -- Add filters to sce and create Seurat object.

  df = data.frame(per_cell_qc_metrics, discard_qc = qc_filter, discard_custom = custom_filter)

  stopifnot(
    nrow(df) == ncol(sce_no_empty_drop)
  )

  if (replace) {
    existing_columns <- intersect(colnames(colData(sce_no_empty_drop)), colnames(df))
    colData(sce_no_empty_drop) <- colData(sce_no_empty_drop)[, !colnames(colData(sce_no_empty_drop)) %in% existing_columns]
  }

  colData(sce_no_empty_drop) <- cbind(colData(sce_no_empty_drop), df)
  
  sce_no_empty_drop
}

# 6
tar_sc_single_qc_tep_make_sensitive_filter <- function(sce_unfiltered) {
  sce_unfiltered[, !sce_unfiltered$discard_qc]
}

# 7
tar_sc_single_step_make_custom_filter <- function(sce_unfiltered) {
  sce_unfiltered[, !sce_unfiltered$discard_custom]
}

# 8
tar_sc_single_qc_step_after_filter_rowsums <- function(sce_filter) {
  sce_filtered %>% 
    counts() %>% 
    rowSums()
}

# 9
tar_sc_single_qc_step_apply_filter <- function(sce_sensitive_filter, sce_custom_filter, save_dataset_sensitive_filtering) {

  min_umi <- 1
  min_ratio_cells <- 0.01

  if (save_dataset_sensitive_filtering) {
    sce_filter <- sce_sensitive_filter
  } else {
    sce_filter <- sce_custom_filter
  }

  num_cells <- min_ratio_cells * ncol(sce_filter)
  sce_gene_filter <- !rowSums(counts(sce_filter) >= min_umi) >= num_cells

  sce_filter[!sce_gene_filter, ]
}



