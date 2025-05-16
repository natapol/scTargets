#' Read 10x Genomics Count Matrix
#'
#' This function reads a 10x Genomics count matrix from the specified path and
#' returns a SingleCellExperiment object. The column names of the resulting
#' object are set to the barcodes from the cell metadata.
#'
#' @param path A string specifying the path to the 10x Genomics data directory.
#'             This directory should contain the matrix, features, and barcodes files.
#'
#' @return A \code{SingleCellExperiment} object containing the count matrix
#'         and associated metadata.
#'
#' @importFrom DropletUtils read10xCounts
#' @importFrom SingleCellExperiment colData
#'
#' @examples
#' \dontrun{
#'   # Example usage:
#'   sce <- tar_sc_single_qc_step_read_10x_counts("path/to/10x/data")
#' }
#'
#' @export
scTarget_single_qc <- function(
  name,
  x10_file_path,
  empty_droplets_fdr_threshold = 0.01,
  replace_unfiltered = TRUE,
  save_dataset_sensitive_filtering = TRUE,
  BPPARAM = BiocParallel::SerialParam()
) {

  targets::tar_assert_path(
    x10_file_path,
    "Input path for 10x data was not found."
  )
  
  read_10x_counts <- function(path) {
    sce_raw <- DropletUtils::read10xCounts(path)
    SummarizedExperiment::colnames(sce_raw) <- SummarizedExperiment::colData(sce_raw)$Barcode
    targets::tar_assert_ge(ncol(sce_raw), 0)
    sce_raw
  }
  
  empty_droplets <- function(sce_raw, BPPARAM) {
    DropletUtils::emptyDrops(
      m = SummarizedExperiment::counts(sce_raw), lower = TRUE, BPPARAM = BPPARAM
    )
  }
  
  remove_empty_droplets <- function(sce_raw, empty_droplets, empty_droplets_fdr_threshold) {
    
    if (is.null(empty_droplets)) {
      sce_raw$is_empty_fdr <- NA
    } else {
      is_cell_index <- which(empty_droplets$FDR <= empty_droplets_fdr_threshold)
      sce_raw <- sce_raw[, is_cell_index]
      sce_raw$is_empty_fdr <- empty_droplets$FDR[is_cell_index]
    }
    
    sce_raw
  }

  list(
    targets::tar_target_raw("path", x10_file_path, format = "file"),
    targets::tar_target_raw(
      name = "sce_raw",
      command = quote(read_10x_counts(path))
    ),
    targets::tar_target_raw(
      name = "empty_droplets", 
      command = substitute(
        empty_droplets(sce_raw, BPPARAM),
        env = list(BPPARAM = BPPARAM)
      )
    ),
    targets::tar_target_raw(
      name = "sce_no_empty_drop", 
      command = substitute(
        remove_empty_droplets(sce_raw, empty_droplets, empty_droplets_fdr_threshold),
        env = list(empty_droplets_fdr_threshold = empty_droplets_fdr_threshold)
      )
    )
  )
}