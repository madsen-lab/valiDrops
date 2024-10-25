#' @title quality_filter
#'
#' @description Function to filter barcodes based on quality metrics
#'
#' @param metrics A data frame containing quality metrics. See \link{quality_metrics}
#' @param mito A boolean (TRUE or FALSE) indicating whether or not to filter using a threshold on the fraction of UMIs derived from mitochondrial genes [default = TRUE].
#' @param distance A boolean (TRUE or FALSE) indicating whether or not to filter using a threshold on the residuals from robust regression between number of features and number of UMIs [default = TRUE].
#' @param coding A boolean (TRUE or FALSE) indicating whether or not to filter using a threshold on the fraction of UMIs derived from protein-coding genes [default = TRUE].
#' @param contrast A boolean (TRUE or FALSE) indicating whether or not to filter using a threshold on the contrast fraction [default = FALSE]. See \link{quality_metrics}
#' @param mito.nreps A numeric indicating the number of times to repeat threshold identification for mitochondrial filtering [default = 10].
#' @param mito.max The maximum allowable mitochondrial threshold to detect before falling back to an alternative strategy [default = 0.3].
#' @param npsi An integer indicating the number of breakpoints for feature to UMI fitting [default = 3].
#' @param dist.thres The maximum number of standard deviations below the mean that passes the QC [default = 5].
#' @param coding.threshold The maximum number of standard deviations around the mean that passes the QC [default = 5].
#' @param contrast.threshold The maximum number of standard deviations around the mean that passes the QC [default = 5].
#' @param plot A boolean (TRUE or FALSE) indicating whether or not to produce plots [default = TRUE].
#' @param tol A number indicating the tolerance parameter for segmentation [default = 1e-05]
#' @param h A number indicating the positive factor by which to increment the breakpoint updates in segmentation [default = 0.01]
#' @param quant A boolean (TRUE or FALSE), which indicates if quantiles (TRUE) or equally spaced values (FALSE) should be used for starting values in segmentation [default = FALSE]
#'
#' @return A list of vectors containing the mitochondrial threshold, number of barcodes filtered at each step and the final barcodes that pass QC filtering.
#' @export
#' @import Matrix
#' @import mixtools 
#' @import inflection
#' @import robustbase
#' @import segmented

quality_filter = function(metrics, mito = TRUE, distance = TRUE, coding = TRUE, contrast = FALSE, mito.nreps = 10, mito.max = 0.3, npsi = 3, dist.threshold = 5, coding.threshold = 3, contrast.threshold = 3, plot = TRUE, tol = 1e-05, h = 0.01, quant = FALSE) {
  ## evaluate arguments
  # metrics matrix
  if (missing(metrics)) { stop('No metrics data frame was provided', call. = FALSE) }
  
  # npsi argument
  if (floor(npsi) <= 0) stop('npsi needs to be a numeric greater than 0', call. = FALSE)
  
  # dist.threshold argument
  if (class(dist.threshold) != "numeric" | dist.threshold <= 0) stop('dist.threshold needs to be a numeric greater than 0', call. = FALSE)
  
  # coding.threshold argument
  if (class(coding.threshold) != "numeric" | coding.threshold <= 0) stop('coding.threshold needs to be a numeric greater than 0', call. = FALSE)

  # contrast.threshold argument
  if (class(contrast.threshold) != "numeric" | contrast.threshold <= 0) stop('contrast.threshold needs to be a numeric greater than 0', call. = FALSE)
  
  # mito.nreps argument
  if (class(mito.nreps) != "numeric" | mito.nreps <= 0) stop('mito.nreps needs to be a numeric greater than 0', call. = FALSE)
  
  # Create a list for the output
  output <- list()
  
  # Set rownames
  rownames(metrics) <- metrics$barcode
  
  ## Mitochondrial filtering	
  if (isTRUE(mito) | is.numeric(mito)) {
    if (sum(colnames(metrics) %in% c("mitochondrial_fraction","logFeatures")) == 2) {  
      if (isTRUE(mito)) {
      # Loop X times to protect against bad fits
      mito.thresholds <- c()
      for (rep in 1:mito.nreps) {
        # Find the group of highly sequenced barcodes
        log <- capture.output({ model <- mixtools::normalmixEM(metrics$logFeatures) })
        grp <- metrics[ which(apply(model$posterior, 1, FUN="which.max") == which.max(model$mu)),]
    
        # Find a threshold
        if (nrow(grp) > 0) {
          sequence <- seq(median(grp$mitochondrial_fraction),1,by=0.001)
        } else {
          sequence <- seq(median(metrics$mitochondrial_fraction),1,by=0.001)
        }
        cnts <- c()
        for (i in sequence) {
          cnts <- c(cnts, sum(grp$mitochondrial_fraction <= i))
        }
        mito.thresholds <- c(mito.thresholds, inflection::uik(sequence, cnts))
      }
     
	    
    # Check if mito threshold is above maximum and fall back to alternative technique using segmentation
    mito.threshold <- median(mito.thresholds)
    if (mito.threshold > mito.max) {
      sample.size <- min(5000, floor(nrow(metrics) * 0.8))
      mito.thresholds <- c()
      for (rep in 1:mito.nreps) {
        psi <- 1
        stop <- 0
        metrics.subsample <- metrics[ sample(1:nrow(metrics), sample.size),]
        model <- lm(logFeatures ~ mitochondrial_fraction, data = metrics.subsample)
        seg <- segmented::segmented(model, npsi = psi, control = segmented::seg.control(quant = TRUE, tol = tol, h = h, quant = quant))
        if (min(seg$psi[,2]) > mito.max) { stop <- 1 }
          while (stop == 1) {
            psi <- psi + 1
            seg <- segmented::segmented(model, npsi = psi, control = segmented::seg.control(quant = TRUE, tol = tol, h = h, quant = quant))
            if (psi >= 5 | min(seg$psi[,2]) <= mito.max) { stop <- 0 }
          }		
          mito.thresholds <- c(mito.thresholds, min(seg$psi[,2]))
        }
        mito.threshold <- median(mito.thresholds)
	  }
      } else {
	mito.threshold <- mito      
      }
  
      # Filter
      output$mitochondrial.threshold <- mito.threshold
      qc.pass <- metrics[ metrics$mitochondrial_fraction <= mito.threshold,"barcode"]
      output$pass.mitochondrial_filter <- qc.pass
  
      # Plot
      if (plot) {
        plot(y=metrics$mitochondrial_fraction, x=metrics$logFeatures, pch=16, xlab="log Total features", ylab="Mitochondrial fraction", col = ifelse(metrics$mitochondrial_fraction > mito.threshold, "red","black"), las=1)
        abline(h = mito.threshold)
        mtext(paste("Threshold = ", mito.threshold))
        mtext(paste("Kept ", length(qc.pass), " barcodes",sep=""), line = 1)
      }
    
    # Filter the metrics internally
    metrics <- metrics[ metrics$barcode %in% qc.pass,]
    } else {
      message("Columns named mitochondrial_fraction and logFeatures does not both exist. Skipping filtering using the mitochondrial fraction.")
    }
  }
  
  if (distance) {
    if (sum(colnames(metrics) %in% c("logUMIs","logFeatures")) == 2) {  
      # Segmented model
      model <- lm(logFeatures ~ logUMIs, data = metrics)
      while (floor(npsi) >= 1) {
	out <- suppressWarnings(segmented::segmented(model, npsi = floor(npsi), control = segmented::seg.control(quant = TRUE, tol = tol, h = h, quant = quant)))
	if (class(out)[1] == "segmented") {
	  break
	} else {
	  npsi <- floor(npsi) - 1
	}
      }
  
      # Define upper and lower thresholds
      upper.threshold <- median(resid(out)) + (robustbase::Sn(resid(out)) * dist.threshold)
      lower.threshold <- median(resid(out)) - (robustbase::Sn(resid(out)) * dist.threshold)
  
      # Select barcodes to keep
      qc.pass <- metrics[ which(resid(out) <= upper.threshold & resid(out) >= lower.threshold),1]
  
      # Plot it
      if (plot) {
        plot(y = metrics$logFeatures, x = metrics$logUMIs, pch=16, las=1, xlab="Total UMIs", ylab="Total features", col = ifelse(metrics$barcode %in% qc.pass, "grey","red"))
        if (class(out)[1] == "segmented") {
	      plot(out, add = TRUE, lwd = 4, col = "blue", rug = FALSE)
	} else {
		abline(out, lwd=4, col="blue")
	}
        mtext(paste("Kept ", length(qc.pass), " barcodes",sep=""))
      }
    
      # Save and filter
      output$pass.distance_filter <- qc.pass
      metrics <- metrics[ metrics$barcode %in% qc.pass,]
    } else {
      message("Columns named logUMIs and logFeatures does not both exist. Skipping filtering using the distance.")
    }
  }
  
  if (coding) {
    if (any(colnames(metrics) == "coding_fraction")) {
      # Define upper and lower thresholds
      lower.threshold <- median(metrics$coding_fraction) - (robustbase::Sn(metrics$coding_fraction) * coding.threshold)
      upper.threshold <- median(metrics$coding_fraction) + (robustbase::Sn(metrics$coding_fraction) * coding.threshold)

      # Select barcodes to keep
      qc.pass <- metrics[ metrics$coding_fraction >= lower.threshold & metrics$coding_fraction <= upper.threshold,"barcode"]
    
      # Plot it
      if (plot) {
        hist(metrics$coding_fraction, breaks = "fd", las = 1, xlab="Fraction of UMIs from protein-coding genes", main = "")
        if (lower.threshold > 0 & lower.threshold < 1) { abline(v = lower.threshold, col ="red", lty=2) } 
        if (upper.threshold > 0 & upper.threshold < 1) { abline(v = upper.threshold, col ="red", lty=2) } 
        mtext(paste("Kept ", length(qc.pass), " barcodes",sep=""))
      }
    
      # Save and filter
      output$pass.coding_filter <- qc.pass
      metrics <- metrics[ metrics$barcode %in% qc.pass,]
    } else {
      message("Column named coding_fraction does not exist. Skipping filtering using the coding fraction.")
    }
  }
  
  if (contrast) {
    if (any(colnames(metrics) == "contrast_fraction")) {
      # Define upper and lower thresholds
      lower.threshold <- median(metrics$contrast_fraction) - (robustbase::Sn(metrics$contrast_fraction) * contrast.threshold)
      upper.threshold <- median(metrics$contrast_fraction) + (robustbase::Sn(metrics$contrast_fraction) * contrast.threshold)
    
      # Select barcodes to keep
      qc.pass <- metrics[ metrics$contrast_fraction >= lower.threshold & metrics$contrast_fraction <= upper.threshold,"barcode"]
    
      # Plot it
      if (plot) {
        hist(metrics$contrast_fraction, breaks = "fd", las = 1, xlab="Contrast fraction", main = "")
        if (lower.threshold > 0 & lower.threshold < 1) { abline(v = lower.threshold, col ="red", lty=2) } 
        if (upper.threshold > 0 & upper.threshold < 1) { abline(v = upper.threshold, col ="red", lty=2) } 
        mtext(paste("Kept ", length(qc.pass), " barcodes",sep=""))
      }
    
      # Save and filter
      output$pass.contrast_filter <- qc.pass
      metrics <- metrics[ metrics$barcode %in% qc.pass,]
    } else {
      message("Column named contrast_fraction does not exist. Skipping filtering using the contrast fraction.")
    }
  }
  
  ## Final output
  output$final <- metrics$barcode
  
  ## return results
  return(output)
}
