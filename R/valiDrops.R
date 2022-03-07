#' @title valiDrops
#'
#' @description Wrapper function that runs the entire valiDrops flow with default parameters.
#'
#' @param counts A matrix containing counts for barcodes passing the barcode-rank threshold.
#' @param rank_barcodes A boolean (TRUE or FALSE) indicating whether or not to rank barcodes [default = TRUE].
#' @param label_apoptotic A boolean (TRUE or FALSE) indicating whether or not to label putative apoptotic cells [default = FALSE].
#' @param ... Pass parameters to functions within valiDrops. See \link{rank_barcodes} \link{quality_metrics} \link{quality_filter} \link{expression_metrics} \link{expression_filter} \link{label_apoptotic}
#'
#' @return A data frame containing quality metrics, as well as quality control labels (and if requested, apoptotic labels) for all barcodes passing the rank threshold.
#' @export
#' @import R.utils
#' @import Matrix
#' @importFrom Seurat GetAssayData
#' @importFrom SingleCellExperiment counts

valiDrops = function(counts, rank_barcodes = TRUE, label_apoptotic = FALSE, ...) {
  ## Check the rank_barcodes parameter
  if (!isTRUE(rank_barcodes) & !isFALSE(rank_barcodes)) { stop("rank_barcodes must be either TRUE or FALSE") }

  ## Check the label_apoptotic parameter
  if (!isTRUE(label_apoptotic) & !isFALSE(label_apoptotic)) { stop("label_apoptotic must be either TRUE or FALSE") }

  ## Extract counts from Seurat or SCE objects
  if (any(class(counts) %in% c("SingleCellExperiment", "Seurat"))) {
    if (class(counts) == "SingleCellExperiment") {
      counts <- counts(counts)
    } else {
      counts <- Seurat::GetAssayData(counts, slot = "counts")
    }
  }

  ## Validate the counts object
  if(missing(counts)) {
    stop('No count matrix was provided', call. = FALSE)
  } else {
    if (!any(class(counts) == c("dgTMatrix", "Matrix","matrix", "dgCMatrix"))) { stop('Count matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix', call. = FALSE) }
  }

  ## convert the counts into dgCMatrix if its class() is not dgCMatrix
  if(class(counts) != "dgCMatrix") { counts = as(counts, "dgCMatrix") }

  ## Run rank_barcodes
  if (rank_barcodes) {
    threshold <- R.utils::doCall(valiDrops::rank_barcodes, args = ..., alwaysArgs = list(counts = counts))
	rank.pass <- rownames(threshold$ranks[ threshold$ranks$counts >= threshold$lower.threshold,])
  } else {
	rank.pass <- colnames(counts)
  }

  ## Subset the counts
  counts.subset <- counts[, colnames(counts) %in% rank.pass]

  ## Run quality_metrics
  metrics <- R.utils::doCall(valiDrops::quality_metrics, args = ..., alwaysArgs = list(counts = counts.subset))

  ## Run quality_filter
  qc.pass <- R.utils::doCall(valiDrops::quality_filter, args = ..., alwaysArgs = list(metrics = metrics$metrics))

  ## Run expression_metrics
  counts.subset.filtered <- counts.subset[ rownames(counts.subset) %in% metrics$protein_coding, colnames(counts.subset) %in% qc.pass$final]
  expr.metrics <- R.utils::doCall(valiDrops::expression_metrics, args = ..., alwaysArgs = list(counts = counts.subset.filtered))

  ## Run expression_filter
  valid <- R.utils::doCall(valiDrops::expression_filter, args = ..., alwaysArgs = list(stats = expr.metrics$stats, clusters = expr.metrics$clusters))

  ## Setup the results
  met <- metrics$metrics
  met$qc.pass <- "fail"
  met[ met$barcode %in% valid, "qc.pass"] <- "pass"

  ## Run label_apoptotic
  if (label_apoptotic) {
	rownames(metrics$metrics) <- metrics$metrics$barcode
    apoptotic <- R.utils::doCall(valiDrops::label_apoptotic, args = ..., alwaysArgs = list(counts = counts.subset, metrics = metrics$metrics))
	met$apoptotic <- "healthy"
	met[ met$barcode %in% apoptotic, "apoptotic"] <- "apoptotic"
  }

  ## Return
  return(met)
}