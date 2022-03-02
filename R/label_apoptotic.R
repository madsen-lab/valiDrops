#' @title label_apoptotic
#'
#' @description Function that uses Isolation Forest to label apoptotic cells.
#'
#' @param counts A matrix containing counts for barcodes passing the barcode-rank threshold. See \link{rank_barcodes}
#' @param metrics A data frame containing quality metrics. See \link{quality_metrics} and \link{quality_filter}
#' @param zscore A number indicating the Z score to use for soft labeling of apoptotic cells [default = 3].
#' @param nfeats A number indicating the number of variable features to use for clustering [default = 5000].
#' @param npcs A number indicating the number of singular values to use for clustering [default = 100].
#' @param anomaly.quantile A number indicating the quantile to use for detecting anomalous barcodes [default = 0.1].
#' @param plot A boolean (TRUE or FALSE) indicating whether or not to produce plots [default = TRUE].
#' @param seed Set a seed for reproducible results. Set to NULL to run in non-deterministic mode [default = 42].
#' @param verbose Set to TRUE for verbose mode. Set to FALSE for silent mode [default = TRUE].
#'
#' @return A data frame object that contains the quality metrics of the data and results indicating the apoptotic cells.
#' @export
#' @import Matrix
#' @import mixtools 
#' @import solitude
#' @import scry
#' @importFrom Seurat GetAssayData
#' @importFrom SingleCellExperiment counts
#' @importFrom sparseMatrixStats colSds

label_apoptotic = function(counts, metrics, zscore = 3, nfeats = 5000, npcs = 20, anomaly.quantile = 0.15, plot = TRUE, seed = 42, verbose = FALSE) {
  
  ## Check the plot parameter
  if (!isTRUE(plot) & !isFALSE(plot)) { stop("plot must be either TRUE or FALSE") }
  
  ## Check the verbose parameter
  if (!isTRUE(verbose) & !isFALSE(verbose)) { stop("verbose must be either TRUE or FALSE") }
  
  ## Check the anomaly.quantile parameter
  if (anomaly.quantile < 0 | anomaly.quantile > 1) { stop("anomaly.quantile must be larger than 0 and less than 1") }
  
  ## Check the zscore parameter
  if (zscore <= 0) { stop("zscore must be larger than 0") }
  
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
  
  #uses quality_counts.metrics to identify anomalies based on the z-score for soft labeling 
  metrics$fraction = (1-metrics$ribosomal_fraction) + (1-metrics$coding_fraction) + metrics$mitochondrial_fraction
  metrics$fraction_Z = (metrics$fraction - mean(metrics$fraction)) / sd(metrics$fraction)
  metrics$s <- 0
  metrics[ metrics$fraction_Z > zscore, "s"] <- 1
 
  ##Return a plot showing the fraction of z scores and selected threshold 
  if(plot){
    plot(metrics$fraction_Z, col = ifelse(metrics$fraction_Z > zscore,"red", "black"), pch = 16, ylab = "Z score", las = 1)
    abline(h = zscore, col = "black", lwd = 2, lty=2)
  }
  
  ####### TODO: Check that *enough* barcodes pass the threshold, else stop

  ## Extract feature sets
  counts = counts[ , colnames(counts) %in% rownames(metrics)]
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
  
  # Define groups for isolationForest
  iF_dead = subset(svd_data, s == 1)
  iF_dead = iF_dead[,-1]
  iF_healthy = subset(svd_data, s == 0)
  iF_healthy = iF_healthy[,-1]
  index = sample(ceiling(nrow(iF_dead) * 0.9))
  
  # Positive labeling
  iforest = isolationForest$new(sample_size = length(index), num_trees = 200)
  iforest$fit(iF_dead)
  
  # Calculate anomaly scores for 'healthy' labels
  scores_unlabeled <- as.data.frame(iforest$predict(data = iF_healthy))
  scores_unlabeled <- scores_unlabeled[ order(scores_unlabeled$anomaly_score, decreasing = TRUE),]
  quantile_cut <- quantile(scores_unlabeled$anomaly_score, anomaly.quantile)

  ####### POTENTIAL IMPROVEMENT: Create better function to detect the breakpoint in anomaly scores, instead of using a quantile-based cutoff

  # Re-label based on the anomaly scores
  svd_data[ rownames(svd_data) %in% rownames(iF_healthy)[ scores_unlabeled$anomaly_score < quantile_cut],"s"] <- 1
  
  # Return
  return(rownames(svd_data[ svd_data$s == 1,]))
 }
  
 