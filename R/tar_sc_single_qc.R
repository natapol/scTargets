
# 2
#' Detect Empty Droplets in Single-Cell Experiment Data
#'
#' This function identifies empty droplets in single-cell experiment data using the `DropletUtils::emptyDrops` method.
#'
#' @param sce_raw A `SingleCellExperiment` object containing raw count data.
#' @param empty_lower_droplets An integer specifying the lower bound for the total UMI count to consider a droplet as potentially empty.
#' @param BPPARAM A `BiocParallelParam` object specifying the parallelization parameters for computation.
#'
#' @return A `DataFrame` object with the results of the empty droplet detection, including p-values and other statistics.
#'
#' @details This function applies the `emptyDrops` method from the `DropletUtils` package to identify droplets that are likely to be empty based on their UMI counts. The `lower` parameter sets the threshold for the minimum UMI count to consider a droplet for testing.
#'
#' @importFrom DropletUtils emptyDrops
#' @importFrom SummarizedExperiment counts
#'
#' @examples
#' # Example usage:
#' library(SingleCellExperiment)
#' library(DropletUtils)
#' library(BiocParallel)
#'
#' # Create a mock SingleCellExperiment object
#' sce_raw <- SingleCellExperiment(assays = list(counts = matrix(rpois(100, lambda = 1), ncol = 10)))
#'
#' # Detect empty droplets
#' result <- tar_sc_single_qc_step_detect_empty_droplets(
#'   sce_raw = sce_raw,
#'   empty_lower_droplets = 100,
#'   BPPARAM = SerialParam()
#' )
#' 
#' @export
empty_droplets <- function(
  sce, 
  BPPARAM = BiocParallel::SerialParam()
) {
  DropletUtils::emptyDrops(
    m = counts(sce), lower = TRUE, BPPARAM = BPPARAM
  )
}


#' Remove Empty Droplets Based on FDR Threshold
#'
#' This function filters out empty droplets from a SingleCellExperiment (SCE) object
#' based on a specified false discovery rate (FDR) threshold. It also adds an 
#' `is_empty_fdr` column to the SCE object to indicate the FDR values for the retained cells.
#'
#' @param sce_raw A `SingleCellExperiment` object containing raw single-cell data. 
#'   Must have non-zero columns.
#' @param empty_droplets A data frame containing information about droplets, 
#'   including their FDR values. Can be `NULL` if no empty droplet information is available.
#' @param empty_droplets_fdr_threshold A numeric value specifying the FDR threshold 
#'   for identifying cells. Droplets with FDR values less than or equal to this 
#'   threshold are retained.
#'
#' @return A `SingleCellExperiment` object with empty droplets removed and an 
#'   additional `is_empty_fdr` column indicating the FDR values for the retained cells.
#'
#' @details
#' - If `empty_droplets` is `NULL`, the function adds an `is_empty_fdr` column 
#'   with `NA` values to the SCE object.
#' - If `empty_droplets` is provided, the function filters the SCE object to 
#'   retain only the cells with FDR values less than or equal to the specified threshold.
#' - The function ensures that the input `sce_raw` has at least one column.
#'
#' @examples
#' # Example usage:
#' sce_filtered <- tar_sc_single_qc_step_remove_empty_drop(
#'   sce_raw = sce_object,
#'   empty_droplets = droplet_data,
#'   empty_droplets_fdr_threshold = 0.01
#' )
#'
#' @export
remove_empty_drop <- function(
  sce,
  empty_droplets,
  empty_droplets_fdr_threshold = 0.01
) {

  if (is.null(empty_droplets)) {
    sce$is_empty_fdr <- NA
  } else {
    is_cell_index <- which(empty_droplets$FDR <= empty_droplets_fdr_threshold)
    sce <- sce[, is_cell_index]
    sce$is_empty_fdr <- empty_droplets$FDR[is_cell_index]
  }
  sce
}

#' Create an Unfiltered SingleCellExperiment Object with QC Metrics
#'
#' This function applies dataset-sensitive and custom cell filtering criteria 
#' to a SingleCellExperiment (SCE) object and appends the filtering results 
#' as metadata columns. The function is designed to handle quality control (QC) 
#' metrics and allows for the replacement of existing metadata columns if specified.
#'
#' @param sce_no_empty_drop A `SingleCellExperiment` object containing the input data 
#'   without empty droplets.
#' @param per_cell_qc_metrics A `data.frame` containing per-cell quality control metrics. 
#'   Expected columns include `total`, `detected`, and `subsets_mito_percent`.
#' @param replace A logical value indicating whether to replace existing metadata 
#'   columns in the SCE object if they overlap with the new QC metrics. Default is `FALSE`.
#'
#' @details
#' The function performs the following steps:
#' 1. **Dataset-sensitive cell filtering**: Identifies outliers based on library size, 
#'    number of detected features, and mitochondrial content using the Median Absolute 
#'    Deviation (MAD) method.
#' 2. **Custom cell filtering**: Applies user-defined thresholds for library size, 
#'    number of detected features, and mitochondrial content.
#' 3. Combines the results of the dataset-sensitive and custom filters using logical 
#'    AND operations.
#' 4. Appends the filtering results as new metadata columns (`discard_qc` and 
#'    `discard_custom`) to the SCE object.
#'
#' @return A `SingleCellExperiment` object with updated metadata columns containing 
#'   the QC filtering results.
#'
#' @examples
#' # Example usage:
#' sce <- tar_sc_single_qc_step_create_unfiltered_sce(
#'   sce_no_empty_drop = sce_object,
#'   per_cell_qc_metrics = qc_metrics,
#'   replace = TRUE
#' )
#'
#' @importFrom scater isOutlier
#' @importFrom purrr map reduce
#' @importFrom tidyr replace_na
#' @importFrom SingleCellExperiment colData
#' @export
create_unfiltered_sce <- function(
  sce,
  mad_threshold = 3,
  min_umi_cf = 1000,
  max_umi_cf = 50000,
  min_feature = 1000,
  max_mito_ratio = 0.2,
  replace_unfiltered = TRUE,
  BPPARAM = BiocParallel::SerialParam()
) {
  
  mito_genes <- stringr::str_which(rowData(sce)[["Symbol"]], stringr::regex("MT-", ignore_case = TRUE))
  ribo_genes <- stringr::str_which(rowData(sce)[["Symbol"]], stringr::regex("^RP[SL]", ignore_case = TRUE))
  
  per_cell_qc_metrics <- scater::perCellQCMetrics(
    sce,
    subsets = list(mito = mito_genes, ribo = ribo_genes),
    BPPARAM = BPPARAM
  )
  
  qc_filter <- list(
      qc_lib = scater::isOutlier(per_cell_qc_metrics$total, log = TRUE, nmads = mad_threshold, type = "lower"),
      qc_nexprs = scater::isOutlier(per_cell_qc_metrics$detected, nmads = mad_threshold, log = TRUE, type = "lower"),
      qc_mito = scater::isOutlier(per_cell_qc_metrics$subsets_mito_percent, nmads = mad_threshold, type = "higher")
      # qc_ribo = isOutlier(per_cell_qc_metrics$subsets_ribo_percent, nmads = !!cfg$MAD_THRESHOLD, type = "higher")
    ) %>%
    purrr::map(~ as.logical(.) %>% tidyr::replace_na(replace = FALSE)) %>%
    purrr::reduce(`&`)


  ## -- Custom filters.
  custom_filter <- list(
      low_count = per_cell_qc_metrics$total <= min_umi_cf,
      high_count = per_cell_qc_metrics$total >= max_umi_cf,
      low_expression = per_cell_qc_metrics$detected <= min_feature,
      high_mito = per_cell_qc_metrics$subsets_mito_percent >= max_mito_ratio * 100
      # low_ribo = per_cell_qc_metrics$subsets_ribo_percent <= !!cfg$MIN_RIBO_RATIO * 100
    ) %>%
  purrr::map( ~ as.logical(.) %>% tidyr::replace_na(replace = FALSE)) %>%
  purrr::reduce(`&`)

  ## -- Add filters to sce and create Seurat object.

  df = data.frame(per_cell_qc_metrics, discard_qc = qc_filter, discard_custom = custom_filter)

  stopifnot(
    nrow(df) == ncol(sce)
  )

  if (replace_unfiltered) {
    existing_columns <- intersect(colnames(colData(sce)), colnames(df))
    colData(sce) <- colData(sce)[, !colnames(colData(sce)) %in% existing_columns]
  }

  colData(sce) <- cbind(colData(sce), df)
  
  sce
}

#' Apply Filtering to Single-Cell Experiment Object
#'
#' This function applies gene filtering to a SingleCellExperiment (SCE) object 
#' based on either a sensitive filtering approach or a custom filtering approach. 
#' The filtering removes genes that do not meet the minimum expression threshold 
#' across a specified proportion of cells.
#'
#' @param sce_sensitive_filter A SingleCellExperiment object to be used for 
#'   sensitive filtering.
#' @param sce_custom_filter A SingleCellExperiment object to be used for 
#'   custom filtering.
#' @param save_dataset_sensitive_filtering A logical value indicating whether 
#'   to use the sensitive filtering dataset (`TRUE`) or the custom filtering 
#'   dataset (`FALSE`).
#'
#' @return A filtered SingleCellExperiment object with genes that meet the 
#'   minimum expression threshold retained.
#'
#' @details The function filters genes based on two criteria:
#'   - A gene must have at least `min_umi` counts in at least `min_ratio_cells` 
#'     proportion of cells.
#'   - Genes that do not meet this criterion are removed from the dataset.
#'
#' @examples
#' # Example usage:
#' # filtered_sce <- tar_sc_single_qc_step_apply_filter(
#' #   sce_sensitive_filter = sce1,
#' #   sce_custom_filter = sce2,
#' #   save_dataset_sensitive_filtering = TRUE
#' # )
#'
#' @importFrom SingleCellExperiment counts
#' @export
apply_filter <- function(
  sce_unfiltered, 
  save_dataset_sensitive_filtering = TRUE,
  min_umi = 1,
  min_ratio_cells = 0.01
) {

  if (save_dataset_sensitive_filtering) {
    sce_filtered <- sce_unfiltered[, !sce_unfiltered$discard_qc]
  } else {
    sce_filtered <- sce_unfiltered[, !sce_unfiltered$discard_custom]
  }

  # num_cells <- min_ratio_cells * ncol(sce_filtered)
  # sce_gene_filter <- !rowSums(counts(sce_filtered) >= min_umi) >= num_cells
  # 
  # sce_filtered[!sce_gene_filter, ]
  # sce_filtered
}

#' Perform Quality Control Step After Filtering Rowsums
#'
#' This function calculates the row sums of the counts matrix from a 
#' SingleCellExperiment object after filtering.
#'
#' @param sce_filter A `SingleCellExperiment` object that has been filtered.
#' 
#' @return A numeric vector containing the row sums of the counts matrix.
#'
#' @examples
#' # Assuming `sce_filtered` is a filtered SingleCellExperiment object:
#' row_sums <- tar_sc_single_qc_step_after_filter_rowsums(sce_filtered)
#'
#' @importFrom SingleCellExperiment counts
#' @export
tar_sc_single_qc_step_after_filter_rowsums <- function(sce_filtered) {
  sce_filtered %>% 
    counts() %>% 
    rowSums()
}


