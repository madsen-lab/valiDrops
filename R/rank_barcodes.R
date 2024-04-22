#' @title rank_barcodes
#'
#' @description Function to rank barcodes according to the number of UMIs (or genes) and detect a cut-off point for putative cell- (or nuclei-) containing barcodes.
#'
#' @param counts A matrix containing counts for all barcodes prior to any filtering.
#' @param type A string ("UMI" or "Genes"), which indicates which feature to use for ranking barcodes [default = "UMI"]. See details for more information.
#' @param psi.min A number indicating the lowest number of breakpoints to test when approximating the curve [default = 2]
#' @param psi.max A number indicating the highest number of breakpoints to test when approximating the curve [default = 5]
#' @param alpha A number indicating the region to find breakpoints within [default = 0.001]. See details for more information.
#' @param alpha.max The maximum allowable number for indicating the region to find breakpoints within [default = 0.05]. See details for more information.
#' @param boot A number indicating the number of bootstrap replicates used to infer breakpoints [default = 10]. Maybe necessary to increase if setting psi.max to a large number.
#' @param factor A number indicating the number of folds above the error of the best model is allowed [default = 1.5]. See details for more information.
#' @param threshold A boolean (TRUE or FALSE), which indicates if the threshold to be included in the output [default = TRUE]
#' @param plot A boolean (TRUE or FALSE), which indicates if a plot should be returned [default = TRUE]
#'
#' @details
#' \strong{Choosing the type of feature to use for ranking of barcodes}\cr
#' In our experience, using the number of UMIs for barcode ranking is the best approach for most single-cell RNA-seq datasets. It generally performs well when the ambient RNA contamination is low to medium. For more contaminated datasets, it may be more robust to use the number of genes for barcode ranking. Generally, we suggest to test using the number of UMIs for barcode ranking first and visually inspect the plot If the threshold is not satisfactory, you may get better result using the number of genes for barcode ranking.
#'
#' \strong{Setting alpha}\cr
#' Breakpoints can only be found within the two vertical blacklines on the diagnostic plot. If the desired breakpoint is located outside of these lines, it is necessary to decrease alpha. If the desired breakpoint is far within the region, it is possible to increase alpha to improve convergence. If the selected alpha is too low, the function will automatically increment alpha untill it reaches alpha.max.
#' 
#' \strong{Setting boot}\cr
#' The function uses the root-mean-squared error to select the best segmentation model. The RMSE decreases with more breakpoints, therefore to choose a simple model that approximates the best model, the selected
#' @return A list or a data frame object that contains ranked barcodes and if chosen a threshold.
#' @export
#' @import segmented
#' @import zoo
#' @import Matrix

rank_barcodes = function(counts, type = "UMI", psi.min = 2, psi.max = 5, alpha = 0.001, alpha.max = 0.05, boot = 10, factor = 1.5, threshold = TRUE, plot = TRUE) {
  ## evaluate arguments
  # count matrix
  if(missing(counts)) {
    stop('No count matrix was provided', call. = FALSE)
  } else {
    if (!any(class(counts) == c("dgTMatrix", "Matrix","matrix", "dgCMatrix","DelayedMatrix"))) { stop('Count matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix, DelayedMatrix', call. = FALSE) }
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
  rmse <- as.data.frame(matrix(ncol=3, nrow = length(psi.min:psi.max)))
  models <- list()
  counter <- 1
  for (psi in psi.min:psi.max) {
    curr.alpha <- alpha
	  model <- lm(y ~ x)
    out <- tryCatch(suppressWarnings(segmented::segmented(model, psi = seq(quantile(x, prob = curr.alpha),quantile(x, prob = (1-curr.alpha)), length.out = psi), control = segmented::seg.control(alpha = (curr.alpha-(curr.alpha/1000)), n.boot = boot))), error = function(e) e)
    if (class(out)[1] == "segmented") {
      rmse[counter,1] <- psi
      rmse[counter,2] <- sqrt(mean(out$residuals^2))
      rmse[counter,3] <- counter
      models[[counter]] <- out
      counter <- counter + 1
    } else if (any(grepl("psi starting values too close", out[[1]]))) {
		  stop = 0
		  while (stop == 0) {
			  curr.alpha <- curr.alpha + alpha
			  if (curr.alpha > alpha.max) { stop = 1 }
			    out <- tryCatch(suppressWarnings(segmented::segmented(model, psi = seq(quantile(x, prob = curr.alpha),quantile(x, prob = (1-curr.alpha)), length.out = psi), control = segmented::seg.control(alpha = (curr.alpha-(curr.alpha/1000)), n.boot = boot))), error = function(e) e)
			    if (class(out)[1] == "segmented") {
				    rmse[counter,1] <- psi
				    rmse[counter,2] <- sqrt(mean(out$residuals^2))
				    rmse[counter,3] <- counter
				    models[[counter]] <- out
				    counter <- counter + 1
				    stop = 1
			    }
		    }
	    }
    }

  ## select the best model (within a factor of the smallest RMSE)
  rmse <- rmse[!is.na(rmse[,1]),]
  out <- models[[min(rmse[ rmse[,2] <= min(rmse[,2])*factor,3])]]
  
  ## select lower threshold
  slope <- segmented::slope(out)$x[,1]
  angles <- c()
  for (iter in 1:(length(slope) - 1)) { angles <- c(angles, atan((slope[iter]-slope[iter + 1]) / (1 + (slope[iter] * slope[iter + 1])))*(180/pi)) }
  best_bpt <- which.min(angles[ -1 ])+1
  lower_rank <- unique.counts[ which.min(abs(unique.counts$rank - out$psi[best_bpt,2])),2]
  lower <- unique.counts[ which.min(abs(unique.counts$rank  - out$psi[best_bpt,2])),1]
  
  ## output a plot of ranks if requested
  if (plot){
    plot(y=unique.counts$counts, x=unique.counts$rank, xlab="log Rank", ylab=paste("log ", feature_type, sep=""), pch=16, las=1, col="#CDCDCD20")
    lines(y=c(0,lower),x=c(lower_rank, lower_rank), col = "red")
    lines(y=c(lower,lower),x=c(0, lower_rank), col = "red")
    abline(v = quantile(x, prob = alpha))
    abline(v = quantile(x, prob = (1-alpha)))
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
