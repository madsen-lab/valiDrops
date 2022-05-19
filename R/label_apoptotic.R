#' @title label_apoptotic
#'
#' @description Function that uses Isolation Forest to label apoptotic cells.
#'
#' @param counts A matrix containing counts for barcodes passing the barcode-rank threshold. See \link{rank_barcodes}
#' @param metrics A data frame containing quality metrics. See \link{quality_metrics} and \link{quality_filter}
#' @param mitochondrial.threshold A mitochondrial threshold identified using quality_filter method. See \link{quality_filter}
#' @param threshold A number indicating the score to use for soft labeling of apoptotic cells [default = 0.2].
#' @param nfeats A number indicating the number of variable features to use for clustering [default = 5000].
#' @param npcs A number indicating the number of singular values to use for clustering [default = 100].
#' @param seed Set a seed for reproducible results. Set to NULL to run in non-deterministic mode [default = 42].
#' @param verbose Set to TRUE for verbose mode. Set to FALSE for silent mode [default = TRUE].
#'
#' @return A data frame object that contains the quality metrics of the data and results indicating the apoptotic cells.
#' @export
#' @import Matrix
#' @import mixtools 
#' @import solitude
#' @import scry
#' @import classInt
#' @importFrom Seurat GetAssayData
#' @importFrom SingleCellExperiment counts
#' @importFrom sparseMatrixStats colSds
#' 
label_apoptotic = function(counts, metrics, mitochondrial.threshold, threshold = 0.2, nfeats = 5000, npcs = 20, seed = 42, verbose = FALSE){
  
  ## Check the verbose parameter
  if (!isTRUE(verbose) & !isFALSE(verbose)) { stop("verbose must be either TRUE or FALSE") }
  
  ## Check the threshold parameter
  if (mitochondrial.threshold <= 0) { stop("zscore must be larger than 0") }
  
  ## Check the threshold parameter
  if (threshold <= 0) { stop("zscore must be larger than 0") }
  
  # nfeats argument
  if(class(nfeats) != "numeric" | nfeats <= 0 | nfeats < npcs) stop('nfeats needs to be a numeric greater than 0 and than npcs', call. = FALSE)
  
  # npcs argument
  if(class(npcs) != "numeric" | npcs <= 0) stop('npcs needs to be a numeric greater than 0', call. = FALSE)
  
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
  
  ## Check metric dataframe if it contains barcodes, the fraction of mitochondrial genes,
  ##the fraction of ribosomal genes and coding gene fraction
  if(is.null(metrics$barcode) && is.null(metrics$mitochondrial_fraction) && is.null(metrics$ribosomal_fraction) && is.null(metrics$coding_fraction)) {
    stop("Provided metrics dataframe does not contain mitochondrial fraction, ribosomal fraction, coding fraction or barcodes.")
  }
  
  # Seed
  if (!is.null(seed)) { set.seed(seed) }
  
  ## Start message
  if (verbose) { message("Data formats and parameters succesfully checked. Starting the apoptotic cell run.\n") }
  start = Sys.time()
  
  metrics <- metrics[metrics$mitochondrial_fraction <= mitochondrial.threshold,]
  
  # Transformation the fractions
  metrics$ribosomal_fraction <- asin(sqrt((metrics$ribosomal_fraction)))/(pi/2)
  metrics$coding_fraction <- asin(sqrt((metrics$coding_fraction)))/(pi/2)
  metrics$mitochondrial_fraction <- asin(sqrt(metrics$mitochondrial_fraction))/(pi/2)
  
  # Calculate a score
  metrics$score <- metrics$mitochondrial_fraction * 6.1 - metrics$mitochondrial_fraction * metrics$ribosomal_fraction * 37 - metrics$mitochondrial_fraction * metrics$coding_fraction * 5 + metrics$mitochondrial_fraction * metrics$coding_fraction * metrics$ribosomal_fraction * 37
  
  # Threshold to get initial labels
  metrics$s <- 0
  metrics[metrics$score > threshold, "s"] <- 1
  
  counts = counts[ , colnames(counts) %in% metrics$barcode]
  nonzero = counts[ Matrix::rowSums(counts) > 0,]
  
  # Calculate size factors
  sf <- 10000 / Matrix::colSums(nonzero)
  
  # Normalization and log1p transformation
  norm_transform <- Matrix::t(Matrix::t(nonzero) * sf)
  norm_transform@x <- log1p(norm_transform@x)
  
  # Variable features for the full dataset
  dev <- scry::devianceFeatureSelection(as.matrix(nonzero))
  var.feats <- names(which(rank(-dev) <= nfeats))
  
  # Feature selection
  data <- Matrix::t(norm_transform[ rownames(norm_transform) %in% var.feats,])
  
  # Feature scaling
  means <- Matrix::colMeans(data)
  sds <- sparseMatrixStats::colSds(data)
  sds[ sds == 0 ] <- 1
  data.scaled <- Matrix::t((Matrix::t(data) - means)/sds)
  
  # SVD
  svd <- irlba::irlba(data.scaled, nv=npcs, nu=npcs)
  svd_data = cbind(as.numeric(as.character(metrics$s)), svd$u)
  svd_data = as.data.frame(svd_data)
  colnames(svd_data) = c("s", paste("SVD_", seq(1,npcs,1), sep=""))
  rownames(svd_data) = rownames(data.scaled)
  
  #subsetting for Isolation Forest
  iF_dead = subset(svd_data, s == 1)
  iF_dead = iF_dead[,-1]
  
  iF_healthy = subset(svd_data, s == 0)
  iF_healthy = iF_healthy[,-1]
  
  index = sample(ceiling(nrow(iF_dead) * 0.8))
  
  #Positive labeling
  iforest = isolationForest$new(sample_size = length(index),
                                num_trees = 500)
  iforest$fit(iF_dead)
  
  # Calculate anomaly scores for 'healthy' labels
  scores_unlabeled <- as.data.frame(iforest$predict(data = iF_healthy))
  scores_unlabeled <- scores_unlabeled[ order(scores_unlabeled$anomaly_score, decreasing = TRUE),]
  
  #Identify the breaks on the anomaly scores
  intervals = classIntervals(scores_unlabeled$anomaly_score, style = "quantile", thr = 0)
  intervals = intervals$brks[3]
  
  #Subset the data based on the biggest jump
  midpoint = subset(scores_unlabeled, anomaly_score < intervals)
  
  #break into intervals and define 
  decision_boundary = classIntervals(midpoint$anomaly_score, style = "quantile", thr = 0)
  decision_boundary = unique(decision_boundary$brks)
  
  #Identifying plateau point
  plateau_point <- c(decision_boundary[5:length(decision_boundary)], rep(NA, 4)) - decision_boundary
  plateua_index <- decision_boundary %in% decision_boundary[which.min(plateau_point):(which.min(plateau_point)+4)]

  #Final cut-off
  decision_boundary = decision_boundary[plateua_index == TRUE][1]
  
  # Re-label based on the anomaly scores
  svd_data[ rownames(svd_data) %in% rownames(iF_healthy)[ scores_unlabeled$anomaly_score <= decision_boundary],"s"] <- 1
  svd_data$barcode <- rownames(svd_data)
  
  #Select barcodes that were marked as apoptotic
  apoptotic <- svd_data[ svd_data$s == 1, "barcode"]
  
  #return apoptotic cells
  return(apoptotic)
}






