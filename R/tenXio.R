
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