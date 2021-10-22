#' @title tf_idf
#'
#' @description Function uses 'Term Frequancy - Inverse Document Frequency' statistical measure to select informative genes
#'
#' @param bcranks a dataframe of ranked genes or UMI
#' @param lower_threshold a threshold to define barcodes that are within acceptable range
#' @param max.frac a maximum fraction
#' @param diff.frac the difference for the fractions
#' @param nsteps.frac steps taken between the fraction range
#' @param plot a boolean
#' @param plot_format a format for the output plot (e.g., available choises: pdf, tiff, jpg, png or none)
#' @param directory a boolean
#' @param ncpus a number of CPUs to be used (default: will be detected identifying user's computer specs)
#'
#' @return A dataframe with statistical metrics from TF-IDF
#'
#' @export
#' @importFrom foreach
#' @importFrom doParallel
#' @importFrom parallel
#' @importFrom segmented
#' @importFrom Matrix
#' @importFrom Seurat
#' @importFrom sparseMatrixStats
#' @importFrom scry
#' @importFrom robustbase

tf_idf = function(bcranks, lower_threshold = 0.0, max.frac = 1.0, diff.frac = 0.5, nsteps.frac = 20, plot = TRUE, plot_format = "none", directory = FALSE, ncpus = "detect") {

  #library dependencies
  if (!require(foreach)){
    stop("foreach library not installed")
  } else if (!require(doParallel)) {
    stop("doParallel library not installed")
  } else if (!require(parallel)) {
    stop("parallel library not installed")
  } else if (!require(segmented)) {
    stop("segmented library not installed")
  } else if (!require(Matrix)) {
    stop("Matrix library not installed")
  } else if (!require(Seurat)) {
    stop("Seurat library not installed")
  } else if (!require(sparseMatrixStats)) {
    stop("sparseMatrixStats library not installed")
  } else if (!require(scry)) {
    stop("scry library not installed")
  } else if (!require(robustbase)) {
    stop("robustbase library not installed")
  }

  #checking if the user chose to use the plot format
  format = TRUE
  if(plot_format == "none"){
    format = FALSE
  }

  ## evaluate arguments
  stopifnot(any(plot %in% c(TRUE, FALSE)))
  stopifnot(any(directory %in% c(TRUE, FALSE)))
  stopifnot(class(max.frac) == "numeric" | class(diff.frac) == "numeric" | class(nsteps.frac) == "numeric" | class(lower_threshold) == "numeric")
  try(if(!any(plot_format %in% c("pdf", "jpg", "tiff", "png", "none"))) stop('Incorrect input. Supported formats: pdf, jpg, tiff, png or none', call. = FALSE))

  if(ncpus != "detect"){
    try(if(class(ncpus) != "numeric") stop("Please specify the number of CPUs to be used or leave as a default to be detected", call. = FALSE))
  }

  ### Identify and remove unspecific barcodes using TF-IDF
  ## Find threshold by injecting a variable number of noisy barcodes
  # Define which barcodes are within acceptable range from the barcode-rank plot
  if (!lower_threshold == 0.0){
    lower <- lower_threshold
    lower <- exp(lower)
    barcodes.above <- rownames(bcranks[ bcranks$counts >= lower,])
  }

  # Calculate frac_range
  # NOTE: Warn the user if min.frac becomes less than 1!!
  max.frac <- min(max.frac,nrow(bcranks)/nrow(bcranks[ bcranks$counts >= lower,]))
  min.frac <- max.frac - diff.frac
  frac_range <- seq(min.frac,max.frac, length.out=nsteps.frac)

  try(if(min.frac < 1) stop('Fraction is less than 1', call. = FALSE))

  # Precalculate the TF matrix
  # Inject max barcodes and select only genes with expression
  nonzero <- counts[, colnames(counts) %in% rownames(bcranks[ 1:ceiling(nrow(bcranks[ bcranks$counts >= lower,])*(1+max.frac)),])]
  nonzero <- nonzero[ Matrix::rowSums(nonzero) > 0,]
  tf = t(t(nonzero) / Matrix::colSums(nonzero)) #computing the TF matrix
  tf@x = log1p(tf@x * 100000) #tf(t,d)=log(1+freq(t,d)), normalization to total count * 100k

  # Loop across fractions (in parallel)
  # NOTE: Is the overhead is worth it for parallel processing?
  if(ncpus == "detect") {
    ncpus <- detectCores() - 1
    cl <- makeCluster(ncpus)
  } else {
    cl <- makeCluster(ncpus)
  }

  registerDoParallel(cl)
  thresholds <- foreach(i=1:length(frac_range), .packages=c("sparseMatrixStats","segmented"), .combine="c") %dopar% {
    # Inject barcodes and select only genes with expression
    frac <- frac_range[i]
    nonzero.subset <- nonzero[, colnames(nonzero) %in% rownames(bcranks[ 1:ceiling(nrow(bcranks[ bcranks$counts >= lower,])*(1+frac)),])]
    nonzero.subset <- nonzero.subset[ Matrix::rowSums(nonzero.subset) > 0,]

    # Calculate TF-IDF (see http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/)
    tf.subset <- tf[, colnames(tf) %in% rownames(bcranks[ 1:ceiling(nrow(bcranks[ bcranks$counts >= lower,])*(1+frac)),])]
    tf.subset <- tf.subset[ rownames(tf.subset) %in% rownames(nonzero.subset),]
    idf = log(1 + ncol(nonzero.subset) / Matrix::rowSums(nonzero.subset)) #idf(t,D)=log(N/count(d ∈ D:t ∈ d))
    tfidf = tf.subset * idf

    # Select informative genes
    binning <- hist(Matrix::rowSums(tfidf), breaks=1000,plot=FALSE)
    x = binning$counts
    y = binning$mids
    model <- lm(x ~ y)
    out <- segmented(model, npsi=1)
    tfidf <- tfidf[ Matrix::rowSums(tfidf) >= out$psi[2],]

    # Calculate the per barcode standard deviation of TF-IDF scores
    tfidf.sd <- data.frame(Mean = Matrix::colMeans(tfidf), SD = sparseMatrixStats::colSds(tfidf))
    tfidf.sd$barcode <- colnames(nonzero.subset)
    tfidf.sd$CV <- tfidf.sd$SD / tfidf.sd$Mean

    # Fit a model
    y <- sort(tfidf.sd$CV, decreasing=TRUE)
    x <- 1:nrow(tfidf.sd)
    model <- lm(y ~ x)
    out <- segmented(model, npsi=2)

    # Save the two threshold
    return(min(y[floor(out$psi[,2])]))
  }
  stopCluster(cl)
  #Time difference of 7.68257 mins
  #3 samples, Mac specs: 16 GB of ram, 1TB of ssd, 15 cpu

  # Select a threshold
  # NOTE: Warn the user if the solutions are (very) unstable (less than X% of solutions within 2 SD?)
  sum(thresholds >= median(thresholds) - (2*robustbase::Sn(thresholds)) & thresholds <= median(thresholds) + (2*robustbase::Sn(thresholds))) / length(thresholds)
  threshold <- median(thresholds)

  ## Define barcodes to keep (within max_frac)
  # Calculate TF-IDF (see http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/)
  idf = log(1 + ncol(nonzero) / Matrix::rowSums(nonzero))
  tfidf = tf * idf

  # Select informative genes
  binning <- hist(Matrix::rowSums(tfidf), breaks=1000,plot=FALSE)
  x = binning$counts
  y = binning$mids
  model <- lm(x ~ y)
  out <- segmented(model, npsi=1)
  tfidf <- tfidf[ Matrix::rowSums(tfidf) >= out$psi[2],]

  # Calculate the per barcode standard deviation of TF-IDF scores
  tfidf.sd <- data.frame(Mean = Matrix::colMeans(tfidf), SD = sparseMatrixStats::colSds(tfidf))
  tfidf.sd$barcode <- colnames(nonzero)
  tfidf.sd$CV <- tfidf.sd$SD / tfidf.sd$Mean
  barcodes.keep <- tfidf.sd[ tfidf.sd$CV <= threshold,3]

  #for plot
  tfidf.sd$col <- "#9999CC"
  tfidf.sd[ tfidf.sd$barcode %in% barcodes.keep,"col"] <- "#99CC66"
  tfidf.sd[ !(tfidf.sd$barcode %in% barcodes.above),"col"] <- "#CC9999"
  tfidf.sd$pch <- 16
  tfidf.sd[ !(tfidf.sd$barcode %in% barcodes.above),"pch"] <- 8

  # Generate data for plotting
  if(plot && directory && format){
    mainDir = getwd()
    subDir = "Quality_control"
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
    do.call(
      plot_format,
      args = list(paste(paste(input_type, "_TFIDF", sep=""), plot_format, sep = "."))
    )
    # Plot the mixes
    plot(tfidf.sd$CV, col = tfidf.sd$col, pch=tfidf.sd$pch, las=1, ylab="TFIDF [CV]")
    dev.off()
  } else if (plot && format) {
    do.call(
      plot_format,
      args = list(paste(paste(input_type, "_TFIDF", sep=""), plot_format, sep = "."))
    )

    # Plot the mixes
    plot(tfidf.sd$CV, col = tfidf.sd$col, pch=tfidf.sd$pch, las=1, ylab="TFIDF [CV]")
    dev.off()
  } else if (plot){
    # Plot the mixes
    plot(tfidf.sd$CV, col = tfidf.sd$col, pch=tfidf.sd$pch, las=1, ylab="TFIDF [CV]")
  }

  #output format
  drops <- c("col","pch")
  output <- tfidf.sd[ ,!(names(tfidf.sd) %in% drops)]
  return(output)

}
