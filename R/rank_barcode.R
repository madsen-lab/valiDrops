#' @title rank_barcode
#'
#' @description Function ranks UMI/genes according to their expression levels
#'
#' @param object a (dense or sparse) matrix containing gene counts
#' @param input_type an object to be ranked (e.g., UMI or genes)
#' @param psi.min a number indicating the lowest number of breakpoints to try when approximating the curve
#' @param psi.max a number indicating the highest number of breakpoints to try when approximating the curve
#' @param threshold a boolean for threshold to be included in the output
#' @param ncpus a number of CPUs to be used (default: will be detected identifying user's computer specs)
#' @param plot a boolean
#' @param plot_format a format for the output plot (e.g., available choises: pdf, tiff, jpg, png or none)
#' @param directory a boolean
#'
#' @return A list or a data frame object that contains ranked barcodes and if chosen a threshold.
#' @export
#' @importFrom segmented segmented
#' @importFrom zoo rollmean
#' @import Matrix
#' @import foreach
#' @import doParallel
#' @import parallel

#barcode ranking and quality filtering
rank_barcode = function(object, input_type = "UMI", psi.min = 1, psi.max = 20, threshold = TRUE, ncpus = "detect", plot = TRUE, plot_format = "none", directory = FALSE){

  #library dependencies
  if (!require(Matrix)){
    stop("Matrix library not installed")
  } else if (!require(segmented)) {
    stop("segmented library not installed")
  } else if (!require(zoo)) {
    stop("zoo library not installed")
  } else if (!require(foreach)) {
    stop("foreach library not installed")
  } else if (!require(doParallel)) {
    stop("doParallel library not installed")
  } else if (!require(parallel)) {
    stop("parallel library not installed")
  } else if (!require(BiocParallel)) {
    stop("BiocParallel library not installed")
  }

  #checking if the user chose to use the plot format
  format = TRUE
  if(plot_format == "none"){
    format = FALSE
  }

  ## evaluate arguments
  try(if(!any(class(counts) == c("dgTMatrix", "Matrix", "dgCMatrix"))) stop('Incorrect input format. Accepted formats: matrix, dgtMatrix, dgCMatrix', call. = FALSE))
  stopifnot(any(plot %in% c(TRUE, FALSE)))
  stopifnot(any(directory %in% c(TRUE, FALSE)))
  stopifnot(any(threshold %in% c(TRUE, FALSE)))
  try(if(!any(input_type %in% c("Genes", "UMI"))) stop('Incorrect input. Did you choose UMI or Genes?', call. = FALSE))
  try(if(!any(plot_format %in% c("pdf", "jpg", "tiff", "png", "none"))) stop('Incorrect input. Supported formats: pdf, jpg, tiff, png or none', call. = FALSE))
  
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
  } else if ( ncpus > 1 )  {
    cl <- makePSOCKcluster(ncpus)
    registerDoParallel(cl)
  } else {
	registerDoSEQ()
  }
  
  ## Get barcode ranks
  switch(input_type,
         UMI = {bcranks <- data.frame(counts = Matrix::colSums(counts))
         rownames(bcranks) <- colnames(counts)
         },
         Genes = {bcranks <- data.frame(counts = Matrix::colSums(counts > 0) )
         rownames(bcranks) <- colnames(counts)
         }
  )
  bcranks$rank <- rank(-bcranks$counts)
  bcranks <- bcranks[ bcranks$counts > 0,]
  bcranks <- bcranks[ order(-bcranks$counts, -bcranks$rank),]

  ## Get only unique counts, and log transform
  unique.counts <- bcranks[duplicated(bcranks$counts)==F,]
  unique.counts$counts <- log(unique.counts$counts)
  unique.counts$rank <- log(unique.counts$rank)

  ## Smooth the data using a moving average
  n <- ceiling(2*(nrow(unique.counts)^(1/3)))
  y <- zoo::rollmean(unique.counts$counts, k = n, align = "center")
  x <- zoo::rollmean(unique.counts$rank, k = n, align = "center")
				
  ## Breakpoint analysis
  rmse <- foreach(psi=psi.min:psi.max, .combine="rbind", .packages = "segmented", .errorhandling = 'remove') %dopar% {
    model <- lm(y ~ x)
	out <- segmented(model, npsi=psi)
    return(c(psi, sqrt(mean(out$residuals^2))))
    }
                
  ## Fit the best model (within a factor 1.5 of the smallest RMSE)
  model <- lm(y ~ x)
  out <- segmented(model, npsi=min(rmse[ rmse[,2] <= min(rmse[,2])*1.5,1]))
				
  ## Select lower threshold
  slope <- slope(out)$x[,1]
  diffs <- c()
  for (iter in 1:(length(slope)-1)) { diffs <- c(diffs, slope[iter] / slope[iter+1]) }
  best_bpt <- which.max(diffs[ -1 ])+1
  lower_rank <- unique.counts[ which.min(abs(unique.counts$rank - out$psi[best_bpt,2])),2]
  lower <- unique.counts[ which.min(abs(unique.counts$rank  - out$psi[best_bpt,2])),1]
				
  #check if counts are bellow the chosen threshold
  #later find a better way to evaluate this
  if(length(unique.counts[unique.counts >= lower]) < 3000) warning('counts are bellow 3000')

  #create plot of ranks with chosen arguments
  if(plot && directory && format){
    mainDir = getwd()
    subDir = "Quality_control"
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
    do.call(
      plot_format,
      args = list(paste(paste(input_type, "_Rank", sep=""), plot_format, sep = "."))
    )
	plot(y=unique.counts$counts, x=unique.counts$rank, xlab="log Rank", ylab="log count", pch=16, las=1, col="#CDCDCD20")
	lines(y=c(0,lower),x=c(lower_rank, lower_rank), col = "red")
	lines(y=c(lower,lower),x=c(0, lower_rank), col = "red")
	legend("topright", box.lty=0, legend = as.expression(bquote("n"^"Lower" ~ " = " ~ .(nrow(bcranks[ bcranks$counts >= exp(lower),])))))
    dev.off()
  } else if (plot && format) {
    do.call(
      plot_format,
      args = list(paste(paste(input_type, "_Rank", sep=""), plot_format, sep = "."))
    )
	plot(y=unique.counts$counts, x=unique.counts$rank, xlab="log Rank", ylab="log count", pch=16, las=1, col="#CDCDCD20")
	lines(y=c(0,lower),x=c(lower_rank, lower_rank), col = "red")
	lines(y=c(lower,lower),x=c(0, lower_rank), col = "red")
	legend("topright", box.lty=0, legend = as.expression(bquote("n"^"Lower" ~ " = " ~ .(nrow(bcranks[ bcranks$counts >= exp(lower),])))))
    dev.off()
  } else if (plot){
	plot(y=unique.counts$counts, x=unique.counts$rank, xlab="log Rank", ylab="log count", pch=16, las=1, col="#CDCDCD20")
	lines(y=c(0,lower),x=c(lower_rank, lower_rank), col = "red")
	lines(y=c(lower,lower),x=c(0, lower_rank), col = "red")
	legend("topright", box.lty=0, legend = as.expression(bquote("n"^"Lower" ~ " = " ~ .(nrow(bcranks[ bcranks$counts >= exp(lower),])))))
  }

  #output
  if(threshold){
    output = list(ranks = bcranks, lower.threshold = exp(lower))
  } else {
    output = bcranks
  }
  return(output)
}
