#' @title rank_barcodes
#'
#' @description Function to rank barcodes according to the number of UMIs (or genes) and detect a cut-off point for putative cell- (or nuclei-) containing barcodes.
#'
#' @param counts A matrix containing counts for all barcodes prior to any filtering.
#' @param type A string ("UMI" or "Genes"), which indicates which feature to use for ranking barcodes [default = "UMI"]. See details for more information.
#' @param psi.min A number indicating the lowest number of breakpoints to test when approximating the curve [default = 1]
#' @param psi.max A number indicating the highest number of breakpoints to test when approximating the curve [default = 20]
#' @param threshold A boolean (TRUE or FALSE), which indicates if the threshold to be included in the output [default = TRUE]
#' @param plot A boolean (TRUE or FALSE), which indicates if a plot should be returned [default = TRUE]
#'
#' @details
#' \strong{Choosing the type of feature to use for ranking of barcodes}\cr
#' In our experience, using the number of UMIs for barcode ranking is the best approach for most single-cell RNA-seq datasets. It generally performs well when the ambient RNA contamination is low to medium. For more contaminated datasets, it may be more robust to use the number of genes for barcode ranking. Generally, we suggest to test using the number of UMIs for barcode ranking first and visually inspect the plot If the threshold is not satisfactory, you may get better result using the number of genes for barcode ranking.
#'
#' \strong{Parallel processing}\cr
#' This function is uses parallel processing through the foreach and doRNG packages. To use this functionality, the backend needs to be registered beforehand. This can be done by the user (using e.g. doFuture) or by the valiDrops wrapper. For more details, see \link{rank_barcodes}
#' @return A list or a data frame object that contains ranked barcodes and if chosen a threshold.
#' @export
#' @importFrom segmented segmented slope
#' @importFrom zoo rollmean
#' @import Matrix

rank_barcodes = function(counts, type = "UMI", psi.min = 1, psi.max = 20, threshold = TRUE, plot = TRUE) {
  ## evaluate arguments
  # count matrix
  if(missing(counts)) {
    stop('No count matrix was provided', call. = FALSE)
  } else {
    if (!any(class(counts) == c("dgTMatrix", "Matrix","matrix", "dgCMatrix"))) { stop('Count matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix', call. = FALSE) }
  }
  
  # type argument
  if(!any(type %in% c("Genes","gene","genes", "UMI", "umi","UMIS","umis","UMIs"))) stop('Incorrect input. Did you choose UMI or Genes?', call. = FALSE)
  
  # psi.min argument
  if(class(psi.min) != "numeric" | psi.min <= 0 | psi.min > psi.max) stop('psi.min needs to be a numeric greater than 0 and less than or equal to psi.max ', call. = FALSE)
  
  # psi.max argument
  if(class(psi.max) != "numeric" | psi.max <= 0 | psi.max < psi.min) stop('psi.max needs to be a numeric greater than 0 and greater than or equal to psi.min ', call. = FALSE)
  
  # threshold argument
  if(class(threshold) != "logical") stop('threshold needs to be a boolean (TRUE or FALSE)', call. = FALSE)
  
  # plot argument
  if(class(plot) != "logical") stop('plot needs to be a boolean (TRUE or FALSE)', call. = FALSE)
  
  ## convert the counts into dgCMatrix if its class() is not dgCMatrix
  if(class(counts) != "dgCMatrix") { counts = as(counts, "dgCMatrix") }
  
  ## get the feature type (allowing for spelling variants)
  feature_type <- "UMI"
  if (type %in% c("Genes","gene","genes")) { feature_type <- "Genes" }
  
  ## get barcode ranks using the appropriate type of feature
  switch(feature_type,
         UMI = {bcranks <- data.frame(counts = Matrix::colSums(counts))
         rownames(bcranks) <- colnames(counts)
         },
         Genes = {bcranks <- data.frame(counts = Matrix::colSums(counts > 0) )
         rownames(bcranks) <- colnames(counts)
         }
  )
  
  ## rank, filter and sort
  bcranks$rank <- rank(-bcranks$counts)
  bcranks <- bcranks[ bcranks$counts > 0,]
  bcranks <- bcranks[ order(-bcranks$counts, -bcranks$rank),]
  
  ## get unique counts, and log transform
  unique.counts <- bcranks[duplicated(bcranks$counts)==F,]
  unique.counts$counts <- log(unique.counts$counts)
  unique.counts$rank <- log(unique.counts$rank)
  
  ## Smooth the data using a moving average
  n <- ceiling(2*(nrow(unique.counts)^(1/3)))
  y <- zoo::rollmean(unique.counts$counts, k = n, align = "center")
  x <- zoo::rollmean(unique.counts$rank, k = n, align = "center")
  
  ## breakpoint analysis
  rmse <- as.data.frame(matrix(ncol=2, nrow = length(psi.min:psi.max)))
  counter <- 1
  for (psi in psi.min:psi.max) {
    model <- lm(y ~ x)
    out <- suppressWarnings(segmented::segmented(model, npsi=psi))
    if (class(out)[1] == "segmented") {
      rmse[counter,1] <- psi
      rmse[counter,2] <- sqrt(mean(out$residuals^2))
      counter <- counter + 1
    }
  }

  ## fit the best model (within a factor 1.5 of the smallest RMSE)
  rmse <- rmse[!is.na(rmse[,1]),]
  model <- lm(y ~ x)
  out <- segmented::segmented(model, npsi=min(rmse[ rmse[,2] <= min(rmse[,2])*1.5,1]))
  
  ## select lower threshold
  slope <- segmented::slope(out)$x[,1]
  diffs <- c()
  for (iter in 1:(length(slope)-1)) { diffs <- c(diffs, slope[iter] / slope[iter+1]) }
  best_bpt <- which.max(diffs[ -1 ])+1
  lower_rank <- unique.counts[ which.min(abs(unique.counts$rank - out$psi[best_bpt,2])),2]
  lower <- unique.counts[ which.min(abs(unique.counts$rank  - out$psi[best_bpt,2])),1]
  
  ## output a plot of ranks if requested
  if (plot){
    plot(y=unique.counts$counts, x=unique.counts$rank, xlab="log Rank", ylab=paste("log ", feature_type, sep=""), pch=16, las=1, col="#CDCDCD20")
    lines(y=c(0,lower),x=c(lower_rank, lower_rank), col = "red")
    lines(y=c(lower,lower),x=c(0, lower_rank), col = "red")
    legend("topright", box.lty=0, legend = as.expression(bquote("n"^"Lower" ~ " = " ~ .(nrow(bcranks[ bcranks$counts >= exp(lower),])))))
  }
  
  ## finalize results depending on arguments
  if(threshold){
    output = list(ranks = bcranks, lower.threshold = exp(lower))
  } else {
    output = bcranks
  }
  
  ## return results
  return(output)
}
