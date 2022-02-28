#' @title expression_metrics
#'
#' @description Function to calculate expression metrics 
#'
#' @param counts A matrix containing counts for barcodes passing quality filtering. See \link{quality_filter}
#' @param nfeats A number indicating the number of variable features to use for clustering [default = 5000].
#' @param npcs A number indicating the number of singular values to use for clustering [default = 10].
#' @param k.min A number indicating the lowest number barcodes in a cluster [default = 5].
#' @param res.shallow A number indicating the resolution to use for shallow clustering [default = 0.1].
#' @param top.n A number indicating the top genes to use for analysis [default = 10].
#'
#' @return A data frame containing expression metrics
#' @export
#' @importFrom irlba irlba
#' @importFrom scry devianceFeatureSelection
#' @importFrom Seurat FindClusters FindNeighbors
#' @importFrom sparseMatrixStats colSds
#' @import Matrix
#' @import matrixTests

#calculate expression metrics
expression_metrics = function(counts, nfeats = 5000, npcs = 10, k.min = 5, res.shallow = 0.1, top.n = 10) {
  ## evaluate arguments
  # count matrix
  if(missing(counts)) {
    stop('No count matrix was provided', call. = FALSE)
  } else {
    if (!any(class(counts) == c("dgTMatrix", "Matrix","matrix", "dgCMatrix"))) { stop('Count matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix', call. = FALSE) }
  }
  
  #check the quality of the matrix
  if (nrow(counts) < ncol(counts)) { stop('You have more cells than genes. This is unexpected. Did you supply an unfiltered count matrix?', call. = TRUE) }
  
  ## convert the counts into dgCMatrix if its class() is not dgCMatrix
  if(class(counts) != "dgCMatrix") { counts = as(counts, "dgCMatrix") }
  
  # nfeats argument
  if(class(nfeats) != "numeric" | nfeats <= 0 | nfeats < npcs) stop('nfeats needs to be a numeric greater than 0 and than npcs', call. = FALSE)
  
  # npcs argument
  if(class(npcs) != "numeric" | npcs <= 0) stop('npcs needs to be a numeric greater than 0', call. = FALSE)
  
  # res.shallow argument
  if(class(res.shallow) != "numeric" | res.shallow <= 0) stop('res.shallow needs to be a numeric greater than 0', call. = FALSE)
  
  # k.min argument
  if(class(k.min) != "numeric" | k.min <= 0) stop('k.min needs to be a numeric greater than 0', call. = FALSE)
  
  # top.n argument
  if(class(top.n) != "numeric" | top.n <= 0) stop('top.n needs to be a numeric greater than 0', call. = FALSE)

    # Subset the count matrix
  nonzero <- counts[ Matrix::rowSums(counts) > 0,]
  
  # Calculate size factors
  sf <- 10000 / Matrix::colSums(nonzero)
  
  # Normalization and log1p transformation
  norm_transform <- Matrix::t(Matrix::t(nonzero) * sf)
  norm_transform@x <- log1p(norm_transform@x)
  
  # Variable features for the full dataset
  dev <- scry::devianceFeatureSelection(as.matrix(nonzero))
  var.feats <- names(which(rank(-dev) <= nfeats))
  
  ## Clustering
  # Feature selection
  data <- Matrix::t(norm_transform[ rownames(norm_transform) %in% var.feats,])
  
  # Feature scaling
  means <- Matrix::colMeans(data)
  sds <- sparseMatrixStats::colSds(data)
  sds[ sds == 0 ] <- 1
  data.scaled <- Matrix::t((Matrix::t(data) - means)/sds)
  
  # SVD
  svd <- irlba::irlba(data.scaled, nv=npcs, nu=npcs)
  sv <- svd$u %*% diag(svd$d)
  rownames(sv) <- rownames(data.scaled)
  colnames(sv) <- paste("PC_", seq(1,npcs,1), sep="")
  
  ## Neighborhood graph
  snn <- Seurat::FindNeighbors(sv, verbose=F)$snn
  
  ### Clustering
  ## Shallow
  clusters.shallow <- Seurat::FindClusters(snn, verbose=F, res=res.shallow)
  
  ## Deep
  # Coarse clustering and processing
  clusters.deep <- Seurat::FindClusters(snn, verbose=F, res=seq(1,20, by = 1))
  mins <- apply(clusters.deep,2,FUN=function(x) { min(table(x))})
  closest <- min(abs(mins - k.min))
  close.res <- as.numeric(substr(names(which(abs(mins - k.min) == closest)),5,100))
  fine.res <- c()
  for (res in 1:length(close.res)) { fine.res <- c(fine.res, seq(close.res[res] - 0.9,close.res[res] + 0.9, by=0.1)) }
  fine.res <- unique(fine.res)
  
  # Fine clustering and processing
  clusters.deep <- Seurat::FindClusters(snn, verbose=F, res=fine.res)
  res.deep <- as.numeric(substr(names(which.max(which(apply(clusters.deep,2,FUN=function(x) { min(table(x))}) == k.min))),5,100))
  if (identical(res.deep, numeric(0))) {
    mins <- apply(clusters.deep,2,FUN=function(x) { min(table(x))})
    k.min <- min(mins[ mins > k.min])
    res.deep <- as.numeric(substr(names(which.max(which(apply(clusters.deep,2,FUN=function(x) { min(table(x))}) == k.min))),5,100))
  }
  
  # Final clustering
  clusters.deep <- Seurat::FindClusters(snn, verbose=F, res=res.deep)

  ## Loop across clusters and calculate stats
  stats <- as.data.frame(matrix(ncol = 8, nrow=length(unique(clusters.deep[,1]))))
  counter <- 1
  
  for (cl.idx in unique(clusters.deep[,1])) {
    # Define groups
    target <- rownames(clusters.deep)[which(clusters.deep[,1] == cl.idx)]
    rest <- rownames(clusters.shallow)[which(clusters.shallow[,1] != names(which.max(table(clusters.shallow[ which(rownames(clusters.shallow) %in% target),1]))))]
    
    # Calculate percentages and differences in percentages for all genes
    pct.1 <- round(x = rowSums(x = norm_transform[, target, drop = FALSE] > 0) / length(x = target), digits = 3)
    pct.2 <- round(x = rowSums(x = norm_transform[, rest, drop = FALSE] > 0) / length(x = rest),  digits = 3 )
    pct.diff <- (pct.1 - pct.2)/pct.1
    
    # Calculate log fold change
    fc <- log(x = rowMeans(x = expm1(x = norm_transform[, target, drop = FALSE])) + 1, base = 2) - log(x = rowMeans(x = expm1(x = norm_transform[, rest, drop = FALSE])) + 1, base = 2)
    
    # Define genes to include for testing
    alpha.min <- pmax(pct.1, pct.2)
    features <- names(x = which(x = alpha.min >= 0.1))
    features.diff <- names(x = which(x = fc >= 0.25))
    features <- intersect(x = features, y = features.diff)
    
    # Run test
    wilcox.res <- suppressWarnings(matrixTests::row_wilcoxon_twosample(as.matrix(norm_transform[features, target, drop = FALSE]),as.matrix(norm_transform[features, rest, drop = FALSE])))
    wilcox.res <- wilcox.res[ order(wilcox.res$pvalue),]
    wilcox.res$FDR <- p.adjust(wilcox.res$pvalue, method="bonferroni", n = nrow(nonzero))
    
    stats[counter,1] <- cl.idx
    stats[counter,2] <- mean(pct.diff[ rownames(wilcox.res[1:min(c(sum(wilcox.res$FDR <= 0.05), top.n)),])])
    stats[counter,3] <- mean(pct.1[ rownames(wilcox.res[1:min(c(sum(wilcox.res$FDR <= 0.05), top.n)),])])
    stats[counter,4] <- mean(pct.2[ rownames(wilcox.res[1:min(c(sum(wilcox.res$FDR <= 0.05), top.n)),])])
    stats[counter,5] <- sum(wilcox.res$FDR <= 0.05)
    stats[counter,6] <- nrow(wilcox.res)
    stats[counter,7] <- sum(pct.diff[ rownames(wilcox.res[1:min(c(sum(wilcox.res$FDR <= 0.05), top.n)),])] < -0.01)
    stats[counter,8] <- min(wilcox.res$FDR)
    counter <- counter + 1
  }
  
  # Remove barcodes with a negative percentual difference
  stats[,9] <- stats[,5] / stats[,6]
  colnames(stats) <- c("cluster","pct.diff","pct.1","pct.2","n_de","n_total","n_negative","min_fdr","de_fraction")
  
  #output
  output <- list(stats = stats, clusters = clusters.deep)
  return(output)
}