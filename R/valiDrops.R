#' @title valiDrops
#'
#' @description Wrapper function that runs the entire valiDrops flow with default parameters.
#'
#' @param counts A matrix containing counts for barcodes passing the barcode-rank threshold.
#' @param rank_barcodes A boolean (TRUE or FALSE) indicating whether or not to rank barcodes [default = TRUE].
#' @param mitochondrial_clusters NULL or an interger indicating how many deviations a cluster must be to be removed [default = 3]. 
#' @param ribosomal_clusters NULL or an interger indicating how many deviations a cluster must be to be removed [default = 3]. 
#' @param label_dead A boolean (TRUE or FALSE) indicating whether or not to label putative dead cells [default = FALSE].
#' @param status A boolean (TRUE or FALSE) indicating whether or not to print the progress of valiDrops [default = TRUE].
#' @param stageThree A boolean (TRUE or FALSE) indicating whether or not to run stage 3 of valiDrops [default = TRUE].
#' @param ... Pass parameters to functions within valiDrops. See \link{rank_barcodes} \link{quality_metrics} \link{quality_filter} \link{expression_metrics} \link{expression_filter} \link{label_dead}
#'
#' @return A data frame containing quality metrics, as well as quality control labels (and if requested, dead labels) for all barcodes passing the rank threshold.
#' @export
#' @import R.utils
#' @import Matrix
#' @import Seurat
#' @import SingleCellExperiment

valiDrops = function(counts, rank_barcodes = TRUE, mitochondrial_clusters = 3, ribosomal_clusters = 3, label_dead = FALSE, status = TRUE, stageThree = TRUE, ...) {
  ## Check the rank_barcodes parameter
  if (!isTRUE(rank_barcodes) & !isFALSE(rank_barcodes)) { stop("rank_barcodes must be either TRUE or FALSE") }
  
  ## Check the label_dead parameter
  if (!isTRUE(label_dead) & !isFALSE(label_dead)) { stop("label_dead must be either TRUE or FALSE") }
  
  ## Check the stageThree parameter
  if (!isTRUE(stageThree) & !isFALSE(stageThree)) { stop("stageThree must be either TRUE or FALSE") }
  
  ## Extract counts from Seurat or SCE objects
  if (any(class(counts) %in% c("SingleCellExperiment", "Seurat"))) {
    if (class(counts) == "SingleCellExperiment") {
      counts <- SingleCellExperiment::counts(counts)
    } else {
      counts <- Seurat::GetAssayData(counts, slot = "counts")
    }
  }
  
  ## Validate the mitochondrial_clusters argument
  if(!is.null(mitochondrial_clusters) & !is.numeric(mitochondrial_clusters) & !is.integer(mitochondrial_clusters)) {
    stop("The mitochondrial_clusters argument must be either NULL or a numeric", call. = FALSE) 	  
  }
  
  ## Validate the ribosomal_clusters argument
  if(!is.null(ribosomal_clusters) & !is.numeric(ribosomal_clusters) & !is.integer(ribosomal_clusters)) {
    stop("The ribosomal_clusters argument must be either NULL or a numeric", call. = FALSE) 	  
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
    if (status) { message("Step 1: Filtering on the barcode-rank plot.")}
    threshold <- R.utils::doCall(valiDrops::rank_barcodes, args = ..., alwaysArgs = list(counts = counts))
    rank.pass <- rownames(threshold$ranks[ threshold$ranks$counts >= threshold$lower.threshold,])
  } else {
    if (status) { message("Step 1: Removing barcodes with zero counts.")}
    rank.pass <- colnames(counts)[which(colSums(counts) > 0)]
  }
  
  ## Check the number of barcodes that pass initial filtering
  if (length(rank.pass) > 20000) {
    message("More than 20.000 barcodes passed initial filtering. It is like that breakpoint estimation did not work satisfactorily. If it didn't, you can try to increase alpha, alpha.max and/or psi.max in the rank_barcodes() function.")
  }
  
  ## Subset the counts
  counts.subset <- counts[, colnames(counts) %in% rank.pass]
  
  ## Run quality_metrics
  if (status) { message("Step 2: Collecting quality metrics.")}
  metrics <- R.utils::doCall(valiDrops::quality_metrics, args = ..., alwaysArgs = list(counts = counts.subset))
  
  ## Run quality_filter
  if (status) { message("Step 3: Filtering on quality metrics.")}
  qc.pass <- R.utils::doCall(valiDrops::quality_filter, args = ..., alwaysArgs = list(metrics = metrics$metrics))
  
  if (stageThree) {
    ## Run expression_metrics
    if (status) { message("Step 4: Collecting expression-based metrics.")}
    counts.subset.filtered <- counts.subset[ rownames(counts.subset) %in% metrics$protein_coding, colnames(counts.subset) %in% qc.pass$final]
    expr.metrics <- R.utils::doCall(valiDrops::expression_metrics, args = ..., alwaysArgs = list(counts = counts.subset.filtered, mito = metrics$mitochondrial, ribo = metrics$ribosomal))
    
    ## Run expression_filter
    if (status) { message("Step 5: Filtering on expression-based metrics.")}
    valid <- R.utils::doCall(valiDrops::expression_filter, args = ..., alwaysArgs = list(stats = expr.metrics$stats, clusters = expr.metrics$clusters, mito = mitochondrial_clusters, ribo = ribosomal_clusters))
  } else {
    valid <- qc.pass$final
  }
  
  ## Setup the results
  met <- metrics$metrics
  met$qc.pass <- "fail"
  met[ met$barcode %in% valid, "qc.pass"] <- "pass"
  
  ## Run label_dead
  if (label_dead) {
    if (status) { 
      if (stageThree) { message("Step 6: Predicting dead cells.") } else { message("Step 4: Predicting dead cells without running stage 3. CAUTION this has not been tested.") }
    }
    dead <- R.utils::doCall(valiDrops::label_dead, args = ..., alwaysArgs = list(counts = counts, metrics = met, qc.labels = setNames(as.character(met$qc), met$barcode)))	  
    dead <- dead$metrics
    met$label <- "live"
    met[ met$barcode %in% dead[ dead$label == "dead","barcode"], "label"] <- "dead"
    met[ met$barcode %in% dead[ dead$label == "uncertain","barcode"], "label"] <- "uncertain"
  }
  
  ## Output to user
  if (status) { message(paste("\t", nrow(met[ met$qc.pass == "pass",]), " barcodes passed quality control.", sep=""))}
  if (label_dead) { if (status) { message(paste("\t", nrow(met[ met$qc.pass == "pass" & met$label == "dead",]), " barcodes that passed quality control are predicted to be dead.", sep=""))} }
  
  ## Return
  return(met)
}
