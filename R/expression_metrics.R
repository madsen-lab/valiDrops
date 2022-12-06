#' @title expression_metrics
#'
#' @description Function to calculate expression metrics 
#'
#' @param counts A matrix containing counts for barcodes passing quality filtering. See \link{quality_filter}
#' @param mito A vector containing names of mitochondrial genes. See \link{quality_metrics}
#' @param ribo A vector containing names of ribosomal genes. See \link{quality_metrics}
#' @param nfeats A number indicating the number of variable features to use for clustering [default = 5000].
#' @param npcs A number indicating the number of singular values to use for clustering [default = 10].
#' @param k.min A number indicating the lowest number barcodes in a cluster [default = 5].
#' @param res.shallow A number indicating the resolution to use for shallow clustering [default = 0.1].
#' @param top.n A number indicating the top genes to use for analysis [default = 10].
#'
#' @return A data frame containing expression metrics
#' @export
#' @import irlba
#' @import scry
#' @import Seurat
#' @import Matrix

expression_metrics = function(counts, mito, ribo, nfeats = 5000, npcs = 10, k.min = 5, res.shallow = 0.1, top.n = 10) {
  ## evaluate arguments
  # count matrix
  if(missing(counts)) {
    stop('No count matrix was provided', call. = FALSE)
  } else {
    if (!any(class(counts) == c("dgTMatrix", "Matrix","matrix", "dgCMatrix"))) { stop('Count matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix', call. = FALSE) }
  }
		  
  ## convert the counts into dgCMatrix if its class() is not dgCMatrix
  if(class(counts) != "dgCMatrix") { counts = as(counts, "dgCMatrix") }
  
  # mito argument
  if (is.null(mito)) {
    stop("Mitochondrial genes needs to be supplied.", call. = FALSE)
  }

  # ribo argument
  if (is.null(ribo)) {
    stop("Ribosomal genes needs to be supplied.", call. = FALSE)
  }
	
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
  nr <- nrow(data)
  sds <- sqrt((Matrix::colMeans(data*data) - Matrix::colMeans(data)^2) * (nr / (nr-1)))
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
    k.min <- mins[ which.min(abs(mins - k.min)) ]
    res.deep <- as.numeric(substr(names(which.max(which(apply(clusters.deep,2,FUN=function(x) { min(table(x))}) == k.min))),5,100))
  }
  
  # Final clustering
  clusters.deep <- Seurat::FindClusters(snn, verbose=F, res=res.deep)

  ## Setup to catch stats across loop
  stats <- as.data.frame(matrix(ncol = 10, nrow=length(unique(clusters.deep[,1]))))
  counter <- 1
	
  ## Check for presto
  presto.flag <- require("presto", quietly = TRUE)
  if (!presto.flag) { message("The 'presto' package is not availible. Consider installing 'presto' from immunogenomics/presto on GitHub for a 80 - 100x speed-up." ) }

  ## Execute loop
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
    
    # Only proceed if at least two genes passed filtering
    if (length(features) >= 2) {    
      # Setup test
      y <- rep("excluded", ncol(norm_transform))
      y[which(colnames(norm_transform) %in% target)] <- "target"
      y[which(colnames(norm_transform) %in% rest)] <- "rest"

      # Run test depending on method
      if (presto.flag) {
	  wilcox.res <- presto::wilcoxauc(X = norm_transform[features,], y = y, groups_use = c("target","rest"))
	  wilcox.res <- wilcox.res[ wilcox.res$group == "target",]
	  wilcox.res <- wilcox.res[ order(wilcox.res$pval),]
	  wilcox.res$FDR <- p.adjust(wilcox.res$pval, method="bonferroni", n = nrow(nonzero))
	  rownames(wilcox.res) <- wilcox.res$feature
      } else {
	  wilcox.res <- as.data.frame(matrix(ncol=2, nrow=length(features)))
	  rownames(wilcox.res) <- features
	  colnames(wilcox.res) <- c("Pvalue","FDR")
	  for (feat in features) {
	  	wilcox.res[ feat,1] <- wilcox.test(norm_transform[ feat, which(y=="target")], norm_transform[ feat, which(y=="rest")])$p.value
	  }
	  wilcox.res[,2] <- p.adjust(wilcox.res[,1], method="bonferroni", n = nrow(nonzero))
      }

      # Collect stats
      stats[counter,1] <- cl.idx
      stats[counter,2] <- mean(pct.diff[ rownames(wilcox.res[1:min(c(sum(wilcox.res$FDR <= 0.05), top.n)),])])
      stats[counter,3] <- mean(pct.1[ rownames(wilcox.res[1:min(c(sum(wilcox.res$FDR <= 0.05), top.n)),])])
      stats[counter,4] <- mean(pct.2[ rownames(wilcox.res[1:min(c(sum(wilcox.res$FDR <= 0.05), top.n)),])])
      stats[counter,5] <- sum(wilcox.res$FDR <= 0.05)
      stats[counter,6] <- nrow(wilcox.res)
      stats[counter,7] <- sum(pct.diff[ rownames(wilcox.res[1:min(c(sum(wilcox.res$FDR <= 0.05), top.n)),])] < -0.01)
      stats[counter,8] <- min(wilcox.res$FDR)
      stats[counter,9] <- 0
      if (length(mito[mito %in% rownames(nonzero[,target])]) > 0) {
    	  stats[counter, 10] <- median((colSums(nonzero[mito[mito %in% rownames(nonzero[,target])],target]) / colSums(nonzero[,target])))
      } else {
	  stats[counter, 10] <- 0
      }
      if (length(ribo[ribo %in% rownames(nonzero[,target])]) > 0) {
    	  stats[counter, 11] <- median((colSums(nonzero[ribo[ribo %in% rownames(nonzero[,target])],target]) / colSums(nonzero[,target])))
      } else {
  	  stats[counter, 11] <- 0
      }
      counter <- counter + 1
    }
  }

  # Calculate fractions and set colnames
  stats <- stats[!is.na(stats[,2]),]	
  stats[,9] <- stats[,5] / stats[,6]
  colnames(stats) <- c("cluster","pct.diff","pct.1","pct.2","n_de","n_total","n_negative","min_fdr","de_fraction","mito_fraction", "ribo_fraction")
  
  # Output
  output <- list(stats = stats, clusters = clusters.deep)
  return(output)
}
