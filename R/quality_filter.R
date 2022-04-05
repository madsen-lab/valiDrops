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
#' @param dist.degree A numeric indicating the degree for the polynomial for robust regression [default = 3].
#' @param dist.thres The maximum number of standard deviations below the mean that passes the QC. Set to ""auto" to automatically determine the threshold [default = "auto"]. See Details for more information
#' @param coding.threshold The maximum number of standard deviations around the mean that passes the QC [default = 5].
#' @param plot A boolean (TRUE or FALSE) indicating whether or not to produce plots [default = TRUE].
#'
#' @details
#' \strong{Thresholding regression residuals on distance between number of features and number of UMIs}\cr
#' In our experience, the residuals from polynomial regression analysis of the number of features and number of UMIs almost always produces a left-tailed normal distribution. The default behaviour is to find the smallest number of standard deviations above the mean that includes all data points and use to threshold below the mean.
#'
#' @return A list of vectors containing the mitochondrial threshold, number of barcodes filtered at each step and the final barcodes that pass QC filtering.
#' @export
#' @import Matrix
#' @import mixtools 
#' @import inflection
#' @import robustbase
#' @import MASS

quality_filter = function(metrics, mito = TRUE, distance = TRUE, coding = TRUE, contrast = FALSE, mito.nreps = 10, dist.degree = 3, dist.thres = "auto", coding.threshold = 5, plot = TRUE) {
  ## evaluate arguments
  # metrics matrix
  if(missing(metrics)) { stop('No metrics data frame was provided', call. = FALSE) }
  
  # dist.degree argument
  if(class(dist.degree) != "numeric" | dist.degree <= 0) stop('dist.degree needs to be a numeric greater than 0', call. = FALSE)
  
  # dist.thres argument
  if (dist.thres != "auto") {
    if(!(class(dist.thres) %in% "numeric")) stop('dist.thres needs to be "auto" or a numeric', call. = FALSE)
  }
  
  # mito.nreps argument
  if(class(mito.nreps) != "numeric" | mito.nreps <= 0) stop('mito.nreps needs to be a numeric greater than 0', call. = FALSE)
  
  # Create a list for the output
  output <- list()
  
  # Set rownames
  rownames(metrics) <- metrics$barcode
  
  ## Mitochondrial filtering
  if (mito) {
    if (sum(colnames(metrics) %in% c("mitochondrial_fraction","logFeatures")) == 2) {  
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
  
      # Select threshold and filter
      mito.threshold <- median(mito.thresholds)
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
      # Robust polynomial regression using a variety of methods
      reg <- robustbase::lmrob(logFeatures ~ poly(logUMIs,dist.degree), data = metrics, setting = "KS2014", method = "SMDM")
  
      # Fit a normal distribution to the residuals
      fit <- MASS::fitdistr(reg$residuals[which(metrics$logUMIs >= (max(metrics$logUMIs) + min(metrics$logUMIs)) / 2)], densfun = "normal")
  
      # Set (or get) the number of standard deviations from the mean to use for filtering
      if (dist.thres == "auto") {
        dist.thres <- (max(reg$residuals[which(metrics$logUMIs >= (max(metrics$logUMIs) + min(metrics$logUMIs)) / 2)]) - fit$estimate[1]) / fit$estimate[2]
      }
  
      # Define upper and lower thresholds
      upper.threshold <- fit$estimate[1] + (fit$estimate[2] * dist.thres)
      lower.threshold <- fit$estimate[1] - (fit$estimate[2] * dist.thres)
  
      # Select barcodes to keep
      qc.pass <- names(which(reg$residuals <= upper.threshold & reg$residuals >= lower.threshold))
  
      # Plot it
      if (plot) {
        plot(y = metrics$logFeatures, x = metrics$logUMIs, pch=16, las=1, xlab="Total UMIs", ylab="Total features", col = ifelse(metrics$barcode %in% qc.pass, "grey","red"))
        lines(x=sort(metrics$logUMIs), y=fitted(reg)[order(metrics$logUMIs)], col='blue', lwd=4) 
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
      # Fit a normal distribution to the fraction
      fit <- MASS::fitdistr(metrics$coding_fraction, densfun = "normal")
    
      # Define upper and lower thresholds
      lower.threshold <- fit$estimate[1] - (fit$estimate[2] * coding.threshold)
      upper.threshold <- fit$estimate[1] + (fit$estimate[2] * coding.threshold)

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
      # Fit a normal distribution to the fraction
      fit <- MASS::fitdistr(metrics$contrast_fraction, densfun = "normal")
    
      # Define upper and lower thresholds
      lower.threshold <- fit$estimate[1] - (fit$estimate[2] * coding.threshold)
      upper.threshold <- fit$estimate[1] + (fit$estimate[2] * coding.threshold)
    
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
      output$pass.coding_filter <- qc.pass
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
