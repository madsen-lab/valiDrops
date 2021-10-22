#' @title rank_barcode
#'
#' @description Function ranks UMI/genes according to their expression levels
#'
#' @param counts a matrix containing counts of the genes or UMI
#' @param input_type an object to be ranked (e.g., UMI or genes)
#' @param threshold a boolean for threshold to be included in the output
#' @param plot a boolean
#' @param plot_format a format for the output plot (e.g., available choises: pdf, tiff, jpg, png or none)
#' @param directory a boolean
#'
#' @return A list or a data frame object that contains ranked barcodes and if chosen a threshold.
#' @export
#' @importFrom Matrix, Seurat, inflection, sparseMatrixStats, segmented

#barcode ranking and quality filtering
rank_barcode = function(counts, input_type = "UMI", threshold = TRUE, plot = TRUE, plot_format = "none", directory = FALSE){

  #library dependencies
  if (!require(Matrix)){
    stop("Matrix library not installed")
  } else if (!require(Seurat)) {
    stop("Seurat library not installed")
  } else if (!require(inflection)) {
    stop("inflection library not installed")
  } else if (!require(sparseMatrixStats)) {
    stop("sparseMatrixStats library not installed")
  } else if (!require(segmented)) {
    stop("segmented library not installed")
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

  #convert the given matrix into dgTMatrix if its class() Matrix or dgCMatrix
  if(class(counts) == "Matrix"){
    counts = as(counts, "dgtMatrix")
  } else if (class(counts) == "dgCMatrix") {
    counts = as(counts, "dgtMatrix")
  }

  #check the quality of the matrix
  try(if(nrow(counts) > ncol(counts))
    stop('you have more genes than cells. This might have happened because count matrix was used unfiltered', call. = FALSE))

  ## Get barcode ranks
  switch(input_type,
         UMI = {bcranks <- data.frame(counts = Matrix::colSums(counts))
         rownames(bcranks) <- colnames(counts)
         },
         Genes = {bcranks <- data.frame(counts = Matrix::colSums(counts) )
         rownames(bcranks) <- colnames(counts)
         bcranks <- subset(bcranks, bcranks$counts > 0)
         }
  )
  bcranks$rank <- rank(-bcranks$counts)
  bcranks <- bcranks[ bcranks$counts > 0,]
  bcranks <- bcranks[ order(-bcranks$counts, -bcranks$rank),]

  ## Get only unique counts, and log transform
  unique.counts <- bcranks[duplicated(bcranks$counts)==F,]
  unique.counts$counts <- log(unique.counts$counts)
  unique.counts$rank <- log(unique.counts$rank)

  ## Loop across n and find break points
  lower.threshold <- c()
  upper.threshold <- c()
  for (n in seq(100,500,by=5)) {
    # Setup to approximating the curve  using n bins
    tmp.counts <- unique.counts
    tmp.counts$cuts <- cut(tmp.counts$counts, breaks = n, labels = F)

    # Approximate the curve
    ap <- approx(y=tmp.counts$counts, x=tmp.counts$rank, xout = aggregate(tmp.counts$rank, by=list(tmp.counts$cuts), FUN="mean")$x)
    x <- ap$x
    y <- ap$y

    # Make a linear model
    model <- lm(y ~ x)

    # Predict 2 breakpoints
    out <- segmented(model, npsi=2)

    # Save them
    lower.threshold  <- c(lower.threshold , max(out$psi[,2]))
    upper.threshold <- c(upper.threshold, min(out$psi[,2]))
  }

  # Select the median solution
  # NOTE: Warn the user if the solutions are (very) unstable (less than X% of solutions within 2 SD?)
  # NOTE: If the solution is unstable (or otherwise poor (hard to define!!)) it might be worth retrying with the number of detected genes as a metric instead
  # NOTE: Trying with number of genes could be a user definable parameter
  upper_rank <- median(upper.threshold)
  upper <- min(unique.counts[ unique.counts$rank <= upper_rank,1])
  lower_rank <- median(lower.threshold)
  lower <- min(unique.counts[ unique.counts$rank <= lower_rank,1])

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
    plot(y=unique.counts$counts, x=unique.counts$rank, xlab="Rank",
         ylab= paste(input_type, " count", sep=""), main = paste(input_type, " rank plot", sep=""), pch=16, las=1) +
      abline(v=lower_rank, col="red") +
      abline(v=upper_rank, col="green") +
      abline(h=lower, col="red") +
      abline(h=upper, col="green")
    dev.off()
  } else if (plot && format) {
    do.call(
      plot_format,
      args = list(paste(paste(input_type, "_Rank", sep=""), plot_format, sep = "."))
    )
    plot(y=unique.counts$counts, x=unique.counts$rank, xlab="Rank",
         ylab= paste(input_type, " count", sep=""), main = paste(input_type, " rank plot", sep=""), pch=16, las=1) +
      abline(v=lower_rank, col="red") +
      abline(v=upper_rank, col="green") +
      abline(h=lower, col="red") +
      abline(h=upper, col="green")
    dev.off()
  } else if (plot){
    plot(y=unique.counts$counts, x=unique.counts$rank, xlab="Rank",
         ylab= paste(input_type, " count", sep=""), main = paste(input_type, " rank plot", sep=""), pch=16, las=1) +
      abline(v=lower_rank, col="red") +
      abline(v=upper_rank, col="green") +
      abline(h=lower, col="red") +
      abline(h=upper, col="green")
  }

  #output
  if(threshold){
    output = list(ranks = bcranks, lower.threshold = lower, upper.threshold = upper)
  } else {
    output = bcranks
  }
  return(output)
}
