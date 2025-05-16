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
read_10x_counts <- function(path) {
  sce_raw <- DropletUtils::read10xCounts(path)
  targets::tar_assert_expr(ncol(sce_raw) > 0)
  colnames(sce_raw) <- colData(sce_raw)$Barcode
  sce_raw
}

#' Read 10X Genomics Data
#'
#' This function reads 10X Genomics single-cell RNA-seq data from the specified path.
#'
#' @param path Character. The file path to the 10X Genomics data directory.
#'
#' @return None. (Function implementation is currently incomplete.)
#'
#' @examples
#' \dontrun{
#' read10x("/path/to/10x/data")
#' }
#'
#' @export
read10x <- function(path) {

  targets::tar_assert_path(path)

  # sce <- SingleCellExperiment(
  #   list(counts=Matrix::readMM("/Users/natapolpornputtapong/Projects/singlecell/raw_feature_bc_matrix/matrix.mtx.gz")),
  #   colData=read.table(file = '/Users/natapolpornputtapong/Projects/singlecell/raw_feature_bc_matrix/barcodes.tsv.gz', sep = '\t', col.names = c("Barcode")),
  #   rowData=read.table(file = '/Users/natapolpornputtapong/Projects/singlecell/raw_feature_bc_matrix/features.tsv.gz', sep = '\t', col.names = c("ID", "Symbol", "Type")),
  #   metadata=list(study="GSE111111")
  # )
  # colnames(sce) <- colData(sce_raw)$Barcode
  # rownames(sce_test) <- rowData(sce_test)$ID
}