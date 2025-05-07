read_10x_counts <- function(path) {
  sce_raw <- read10xCounts(path)
  colnames(sce_raw) <- colData(sce_raw)$Barcode
  sce_raw
}

detect_empty_droplets <- function(sce_raw, empty_lower_droplets = TRUE, BPPARAM = BiocParallel::SerialParam()) {
  DropletUtils::emptyDrops(
    m = counts(sce_raw), lower = empty_lower_droplets, BPPARAM = BPPARAM
  )
}

remove_empty_drop <- function(sce_raw, empty_droplets, empty_droplets_fdr_threshold = 0.01) {
  # empty_droplets_fdr_threshold == !!cfg$EMPTY_DROPLETS_FDR_THRESHOLD

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

cal_per_cell_qc_metrics <- function(sce_no_empty_drop) {
  mito_genes <- stringr::str_which(rowData(sce_no_empty_drop)[["Symbol"]], stringr::regex("MT-", ignore_case = TRUE))
  ribo_genes <- stringr::str_which(rowData(sce_no_empty_drop)[["Symbol"]], stringr::regex("^RP[SL]", ignore_case = TRUE))

  scater::perCellQCMetrics(
    sce_no_empty_drop,
    subsets = list(mito = mito_genes, ribo = ribo_genes), BPPARAM = ignore(BiocParallel::bpparam())
  )
}

create_unfiltered_sce <- function(sce_no_empty_drop, per_cell_qc_metrics, replace = TRUE) {
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

make_sensitive_filter <- function(sce_unfiltered) {
  sce_unfiltered[, !sce_unfiltered$discard_qc]
}

make_custom_filter <- function(sce_unfiltered) {
  sce_unfiltered[, !sce_unfiltered$discard_custom]
}

after_filter_rowsums <- function(sce_filter) {
  sce_filtered %>% 
    counts() %>% 
    rowSums()
}

apply_filter <- function(sce_sensitive_filter, sce_custom_filter, save_dataset_sensitive_filtering=TRUE) {

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
