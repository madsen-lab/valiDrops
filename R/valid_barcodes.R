#' @title valid_barcodes
#'
#' @description Function to detect valid barcodes 
#'
#' @param counts a (dense or sparse) matrix containing gene counts
#' @param barcodes a vector of barcodes passing the lower threshold on the barcode-rank plot
#' @param ncpus a number of CPUs to be used (default: will be detected identifying user's computer specs)
#'
#' @return A list of valid barcodes
#' @export
#' @importFrom scuttle librarySizeFactors
#' @importFrom irlba irlba
#' @importFrom scry devianceFeatureSelection
#' @importFrom Seurat FindClusters FindNeighbors
#' @importFrom coop pcor
#' @importFrom sparseMatrixStats colSds
#' @import Matrix
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import BiocParallel

#finding valid barcodes
valid_barcodes = function(counts, barcodes, ncpus = "detect"){

  #library dependencies
  if (!require(Matrix)){
    stop("Matrix library not installed")
  } else if (!require(scuttle)) {
    stop("scuttle library not installed")
  } else if (!require(Seurat)) {
    stop("Seurat library not installed")
  } else if (!require(foreach)) {
    stop("foreach library not installed")
  } else if (!require(doParallel)) {
    stop("doParallel library not installed")
  } else if (!require(parallel)) {
    stop("parallel library not installed")
  } else if (!require(BiocParallel)) {
    stop("BiocParallel library not installed")
  } else if (!require(scry)) {
    stop("scry library not installed")
  } else if (!require(irlba)) {
    stop("irlba library not installed")
  } else if (!require(coop)) {
    stop("coop library not installed")
  } else if (!require(sparseMatrixStats)) {
    stop("sparseMatrixStats library not installed")
  }

  ## evaluate arguments
  try(if(!any(class(counts) == c("dgTMatrix", "Matrix", "dgCMatrix"))) stop('Incorrect input format. Accepted formats: matrix, dgTMatrix, dgCMatrix', call. = FALSE))
  stopifnot(sum(colnames(counts) %in% barcodes) == length(barcodes))
  
  if(ncpus != "detect"){
    try(if(class(ncpus) != "numeric") stop("Please specify the number of CPUs to be used or leave as a default to be detected", call. = FALSE))
  }
   
  #convert the given matrix into dgCMatrix if its class() Matrix or dgTMatrix
  if(class(counts) == "Matrix"){
    counts = as(counts, "dgCMatrix")
  } else if (class(counts) == "dgTMatrix") {
    counts = as(counts, "dgCMatrix")
  }

  #check the quality of the matrix
  try(if(nrow(counts) > ncol(counts))
    stop('You have more genes than cells. This is unexpected. Is your count matrix was used unfiltered?', call. = TRUE))
	
  # Start a cluster
  if(ncpus == "detect") {
    ncpus <- detectCores() - 1
    cl <- makePSOCKcluster(ncpus)
    registerDoParallel(cl)
	bp <- BiocParallel::DoparParam()
  } else if ( ncpus > 1 )  {
    cl <- makePSOCKcluster(ncpus)
    registerDoParallel(cl)
	bp <- BiocParallel::DoparParam()
  } else {
	registerDoSEQ()
	bp <- BiocParallel::SerialParam()
  }
  
  # Subset the count matrix
  nonzero <- counts[, colnames(counts) %in% barcodes,]
  nonzero <- nonzero[ Matrix::rowSums(nonzero) > 0,]

  # Calculate and center size factors
  sf <- scuttle::librarySizeFactors(nonzero, BPPARAM=bp)
  sf <- sf/exp(mean(log(sf)))

  # Normalization and acosh transformation (# Source: transformGamPoi)
  norm_transform <- Matrix::t(Matrix::t(nonzero) / sf)
  norm_transform@x <- 1/sqrt(0.01) * acosh((2 * 0.01 * norm_transform@x) + 1)

  # Variable features for the full dataset
  dev <- scry::devianceFeatureSelection(as.matrix(nonzero))
  var.feats <- names(which(rank(-dev) <= 5000))

  ## Clustering
  # Feature selection
  data <- Matrix::t(norm_transform[ rownames(norm_transform) %in% var.feats,])

  # Feature scaling
  means <- Matrix::colMeans(data)
  sds <- sparseMatrixStats::colSds(data)
  sds[ sds == 0 ] <- 1
  data.scaled <- Matrix::t((Matrix::t(data) - means)/sds)

  # SVD
  svd <- irlba::irlba(data.scaled, nv=10, nu=10)
  sv <- svd$u %*% diag(svd$d)
  rownames(sv) <- rownames(data.scaled)

  # Louvaine
  snn <- Seurat::FindNeighbors(sv, verbose=F)$snn
  clusters <- Seurat::FindClusters(snn, verbose=F, res=20)

  ## Process clusters
  # Cluster-to-cluster correlation
  mean.cl <- as.data.frame(matrix(nrow=ncol(sv), ncol=length(unique(clusters[,1]))))
  for (cl.idx in seq(0,length(unique(clusters[,1]))-1,1)) {
      mean.cl[,cl.idx+1] <- colMeans(sv[ which(clusters[,1] == cl.idx),])
  }
  cor.mat <- coop::pcor(mean.cl)
  colnames(cor.mat) <- seq(0,length(unique(clusters[,1]))-1,1)

  # Cluster specificity
  stats <- foreach(cl.idx=seq(0,length(unique(clusters[,1]))-1,1), .combine="rbind", .packages="Matrix") %dopar% {
      exprs <- (Matrix::rowSums(norm_transform[ , which(clusters[,1] == cl.idx)]) / sum(norm_transform[ , which(clusters[,1] == cl.idx)])) - (Matrix::rowSums(norm_transform[ , which(clusters[,1] %in% names(which(cor.mat[cl.idx+1,] <= 0.8)))]) / sum(norm_transform[ , which(clusters[,1] %in% names(which(cor.mat[cl.idx+1,] <= 0.8)))]))
      freq <- (Matrix::rowSums(norm_transform[ , which(clusters[,1] == cl.idx)] > 0) / sum(norm_transform[ , which(clusters[,1] == cl.idx)] > 0)) - (Matrix::rowSums(norm_transform[ , which(clusters[,1] %in% names(which(cor.mat[cl.idx+1,] <= 0.8)))] > 0) / sum(norm_transform[ , which(clusters[,1] %in% names(which(cor.mat[cl.idx+1,] <= 0.8)))] > 0))
      c(cl.idx, mean(exprs* freq))
  }
  
  # Filtering
  valid.barcodes <- rownames(clusters)[ clusters[,1] %in% stats[ stats[,2] < median(stats[,2])+mad(stats[,2]),1]]
  
  #stop the cluster
  if (ncpus > 1) { stopCluster(cl) }
  
  #output
  return(valid.barcodes)
}
