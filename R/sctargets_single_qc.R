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
  empty_lower_droplets = TRUE,
  replace_unfiltered = TRUE,
  save_dataset_sensitive_filtering = TRUE,
  BPPARAM = BiocParallel::SerialParam()
) {

  read_10x_counts <- function(path) {
    sce_raw <- DropletUtils::read10xCounts(path)
    colnames(sce_raw) <- SummarizedExperiment::colData(sce_raw)$Barcode
    sce_raw
  }

  targets::tar_assert_path(
    x10_file_path,
    "Input path for 10x data was not found."
  )

  list(
    targets::tar_target_raw("path", x10_file_path, format = "file"),
    targets::tar_target_raw(
      name = "sce_raw",
      command = quote(read_10x_counts(path))
    )
  )
}