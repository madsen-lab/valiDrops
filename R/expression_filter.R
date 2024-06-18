#' @title expression_filter
#'
#' @description Function to filter on expression metrics 
#'
#' @param stats A data frame containing expression metrics for clusters. See \link{expression_metrics}
#' @param clusters A data frame mapping barcodes to clusters. See \link{expression_metrics}
#' @param mito NULL or an interger indicating how many deviations a cluster must be to be removed [default = 3].
#' @param ribo NULL or an interger indicating how many deviations a cluster must be to be removed [default = 3].
#' @param min.significant A number indicating the minimum number of significant genes a cluster must have to pass filtering [default = 1].
#' @param min.target.pct A number indicating the minimum mean fraction of barcodes in a cluster that expresses top.n marker genes [default = 0.3].
#' @param max.background.pct A number indicating the maximum mean fraction of barcodes outside a cluster that expresses top.n marker genes [default = 0.8].
#' @param min.diff.pct A number indicating the minimum difference in percent of barcodes within versus outside a cluster that expresses top.n marker genes [default = 0.2].
#' @param min.de.frac A number indicating the minimum fraction of tested genes that must be significant for a cluster to pass filtering [default = 0.01]
#' @param min.significance.level A number indicating the minimum significance level of the most signficiant marker gene for a cluster to pass filtering. Set to NULL to identify the threshold automatically [default = NULL].
#' @param plot A boolean (TRUE or FALSE), which indicates if a plot should be returned [default = FALSE]
#' @param tol A number indicating the tolerance parameter for segmentation [default = 1e-05]
#' @param maxit.glm A number indicating the maximum number of iterations for inner IWLS iterations in segmentation [default = NULL]
#' @param h A number indicating the positive factor by which to increment the breakpoint updates in segmentation [default = 1.25]
#' @param quant A boolean (TRUE or FALSE), which indicates if quantiles (TRUE) or equally spaced values (FALSE) should be used for starting values in segmentation [default = FALSE]
#'
#' @return A vector of valid barcodes
#' @import segmented
#' @export

#finding valid barcodes
expression_filter = function(stats, clusters, mito = 3, ribo = 3, min.significant = 1, min.target.pct = 0.3, max.background.pct = 0.7, min.diff.pct = 0.2, min.de.frac = 0.01, min.significance.level = NULL, plot = FALSE, tol = 1e-50, maxit.glm = 2500, h = 0.01, quant = TRUE) {
  ## evaluate arguments
  # min.significant argument
  if(class(min.significant) != "numeric" | min.significant < 0) stop('min.significant needs to be a numeric greater than or equal to 0', call. = FALSE)
  
  # Validate the mito argument
  if(!is.null(mito) & !is.numeric(mito) & !is.integer(mito)) {
    stop("The mito argument must be either NULL or a numeric", call. = FALSE) 	  
  }
  
  # Validate the ribo argument
  if(!is.null(ribo) & !is.numeric(ribo) & !is.integer(ribo)) {
    stop("The ribo argument must be either NULL or a numeric", call. = FALSE) 	  
  }
  
  # min.target.pct argument
  if(class(min.target.pct) != "numeric" | min.target.pct < 0 | min.target.pct > 1) stop('min.target.pct needs to be a numeric greater than or equal to 0 and less than or equal to 1', call. = FALSE)
  
  # max.background.pct argument
  if(class(max.background.pct) != "numeric" | max.background.pct < 0 | max.background.pct > 1) stop('max.background.pct needs to be a numeric greater than or equal to 0 and less than or equal to 1', call. = FALSE)
  
  # min.diff.pct argument
  if(class(min.diff.pct) != "numeric" | min.diff.pct < 0 | min.diff.pct > 1) stop('min.diff.pct needs to be a numeric greater than or equal to 0 and less than or equal to 1', call. = FALSE)
  
  # min.significance.level argument
  if(class(min.significance.level) != "NULL" & class(min.significance.level) != "numeric") { 
    stop('min.significance.level needs to be a numeric greater than 0', call. = FALSE)
  } else if (class(min.significance.level) == "numeric") {
    if (min.significance.level < 0) { stop('min.significance.level needs to be a numeric greater than 0', call. = FALSE) }
  }
  
  # plot argument
  if(class(plot) != "logical") stop('plot needs to be a boolean (TRUE or FALSE)', call. = FALSE)

  # h argument
  if(class(h) != "numeric" | h < 0) stop('h needs to be a numeric greater than 0', call. = FALSE)

  # maxit.glm argument
  if((class(maxit.glm) != "NULL" & class(maxit.glm) != "numeric") | (class(maxit.glm) == "numeric" & maxit.glm < 0)) { stop('maxit.glm needs to be NULL or a numeric greater than 0', call. = FALSE)  }
  
  # tol argument
  if(class(tol) != "numeric" | tol < 0) stop('tol needs to be a numeric greater than 0', call. = FALSE)
  
  # quant argument
  if(class(quant) != "logical") stop('quant needs to be a boolean (TRUE or FALSE)', call. = FALSE)

  # model function to return error if necessary
  model.significance.level.function <- function(model){result <- tryCatch({segmented::segmented(model, npsi = 1, control = segmented::seg.control(quant = quant, tol = tol, maxit.glm = maxit.glm, h = h))$psi[2]}, error = function(e){NA})}
  
  # Find threshold on significance level
  if (is.null(min.significance.level)) {
    subset <- stats[stats[, 8] > 0, ]
    y <- subset[, 2]
    x <- -log10(subset[, 8])
    threshold.significance.level <- median(x[y <= 0.4]) + (robustbase::Sn(x[y <= 0.4]) * 3)
    model <- lm(y ~ x)
    model.significance.level <- model.significance.level.function(model)
    min.significance.level <- min(threshold.significance.level, model.significance.level, na.rm=TRUE)
  }

  # Stop function if min.significance.level is 0
  if(min.significance.level == 0){
    message("Error fitting model. Skipping stage 3 filtering")
    valid.barcodes <- rownames(clusters)
  } else {
    # Filter
    stats.filtered <- stats[ stats[,7] == 0,]
    stats.filtered <- stats.filtered[ stats.filtered[,2] >= min.diff.pct,]
    stats.filtered <- stats.filtered[ stats.filtered[,3] >= min.target.pct,]
    stats.filtered <- stats.filtered[ stats.filtered[,4] <= max.background.pct,]
    stats.filtered <- stats.filtered[ stats.filtered[,5] >= min.significant,]
    stats.filtered <- stats.filtered[ -log10(stats.filtered[,8]) >= min.significance.level,]
    stats.filtered <- stats.filtered[ stats.filtered[,9] > min.de.frac,]
    if (!is.null(mito)) {
      stats.filtered <- stats.filtered[ stats.filtered[,1] %in% stats[ stats[,10] <= median(stats[,10]) + (mito * robustbase::Sn(stats[,10])),1],]
    }
    if (!is.null(ribo)) {
      stats.filtered <- stats.filtered[ stats.filtered[,1] %in% stats[ stats[,11] <= median(stats[,11]) + (ribo * robustbase::Sn(stats[,11])),1],]
    }
    
    # Plot
    if (plot) {
      plot(stats[,7], col = ifelse(stats[,7] == 0, "green","red"), pch=ifelse(stats[,1] %in% stats.filtered[,1], 16, 3), ylab="Number of genes with negative percentual enrichment", xlab="Cluster number", las = 1)
      plot(stats[,2], col = ifelse(stats[,2] >= min.diff.pct, "green","red"), pch=ifelse(stats[,1] %in% stats.filtered[,1], 16, 3), ylab="Average percentual difference between cluster and background", las = 1)
      abline(h = min.diff.pct)
      plot(stats[,3], col = ifelse(stats[,3] >= min.target.pct, "green","red"), pch=ifelse(stats[,1] %in% stats.filtered[,1], 16, 3), ylab="Average percentual expression frequency in cluster", las = 1, ylim=c(0,1))
      abline(h = min.target.pct)
      plot(stats[,4], col = ifelse(stats[,4] <= max.background.pct, "green","red"), pch=ifelse(stats[,1] %in% stats.filtered[,1], 16, 3), ylab="Average percentual expression frequency in background", las = 1, ylim = c(0,1))
      abline(h = max.background.pct)
      plot(stats[,5], col = ifelse(stats[,5] >= min.significant, "green","red"), pch=ifelse(stats[,1] %in% stats.filtered[,1], 16, 3), ylab="Number of significant genes", las = 1)
      abline(h = min.significant)
      plot(-log10(stats[,8]), col = ifelse(-log10(stats[,8]) >= min.significance.level, "green","red"), pch=ifelse(stats[,1] %in% stats.filtered[,1], 16, 3), ylab="Maximum significance level in cluster", las = 1)
      abline(h = min.significance.level)
      if (!is.null(mito)) {
      plot(stats[,10], col = ifelse(stats[,10] <= median(stats[,10]) + (mito * robustbase::Sn(stats[,10])), "green","red"), pch=ifelse(stats[,1] %in% stats.filtered[,1], 16, 3), ylab="Mitochondrial content", las = 1)
      abline(h = median(stats[,10]) + (mito * robustbase::Sn(stats[,10])))
      }
      if (!is.null(ribo)) {
      plot(stats[,11], col = ifelse(stats[,11] <= median(stats[,11]) + (mito * robustbase::Sn(stats[,11])), "green","red"), pch=ifelse(stats[,1] %in% stats.filtered[,1], 16, 3), ylab="Ribosomal content", las = 1)
      abline(h = median(stats[,11]) + (ribo * robustbase::Sn(stats[,11])))
      }
      plot(stats[,9], col = ifelse(stats[,9] >= min.de.frac, "green","red"), pch=ifelse(stats[,1] %in% stats.filtered[,1], 16, 3), ylab="Fraction of genes that are significant", las = 1, ylim=c(0,1))
      abline(h = min.de.frac)
      legend(x = "topright", legend = c("Kept by current filter", "Removed by current filter", "Removed by any filter"), pch = c(16,16,3), col=c("green","red","black"))
    }
    
    # Return good barcodes
    valid.barcodes <- rownames(clusters)[ which(clusters[,1] %in% stats.filtered[,1])]
  }
  #output
  return(valid.barcodes)
}
