#' @title cell_classifier
#'
#' @description Uses two steps classification method to separate healthy cells from damaged/dying cells.
#' The mixture of regressions is used to obtain the primary labels in identifying two populations.
#' Later identified decision boundary is used to label dead cells positively and train one-class Support Vector Machines on them to identify ambiguous cells from unlabeled (= healthy) data
#'
#' @param counts a matrix containing counts of the genes or UMI
#' @param barcodes.above ranked barcodes that passed the threshold
#' @param gamma a hypermeter used in Support Vector Machines to decide the decision boundaries (the higher the gamma the more curvative is the hyperplane) (default = 0.1)
#' @param nu a hypermeter used in SVM as an upper bound on the fraction of margin errors and a lower bound of the fraction of support vectors (0,1] (default = 0.01)
#' @param plot a boolean
#' @param plot_format a format for the output plot (e.g., available choises: pdf, tiff, jpg, png or none)
#' @param directory a boolean
#' @param ncpus a number of CPUs to be used (default: will be detected identifying user's computer specs)
#'
#' @return A dataframe with claasified cells into healthy and dying cells.
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
#' @importFrom e1071
#' @importFrom dplyr

cell_classifier = function(counts, barcodes.above, gamma = 0.1, nu = 0.01, plot = TRUE, plot_format = "none", directory = FALSE, ncpus = "detect"){

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
  } else if (!require(e1071)) {
    stop("e1017 library not installed")
  } else if (!require(magdplyrrittr)) {
    stop("dplyr library not installed")
  }

  #checking if the user chose to use the plot format
  format = TRUE
  if(plot_format == "none"){
    format = FALSE
  }

  ##Evaluate arguments
  stopifnot(any(plot %in% c(TRUE, FALSE)))
  stopifnot(any(directory %in% c(TRUE, FALSE)))
  stopifnot(class(gamma) == "numeric" | class(nu) == "numeric")
  try(if(!any(plot_format %in% c("pdf", "jpg", "tiff", "png", "none"))) stop('Incorrect input. Supported formats: pdf, jpg, tiff, png or none', call. = FALSE))

  if(ncpus != "detect"){
    try(if(class(ncpus) != "numeric") stop("Please specify the number of CPUs to be used or leave as a default to be detected", call. = FALSE))
  }


  ### Identify and remove dead barcodes (using mitochondrial fractions)
  # Extract stats about the number of counts and mitochondrial fraction
  nonzero <- counts[, colnames(counts) %in% barcodes.above]
  nonzero <- nonzero[ Matrix::rowSums(nonzero) > 0,]
  stats <- data.frame(total_features = log(Matrix::colSums(nonzero > 0)), mitochondrial_fraction = (Matrix::colSums(nonzero[ rownames(nonzero) %in% features[ grep("^MT-", features$V2),1], ]) / Matrix::colSums(nonzero)))

  # Loop across fractions (in parallel)
  # NOTE: Is the overhead is worth it for parallel processing?
  if(ncpus == "detect") {
    ncpus <- detectCores() - 1
  }

  # Setup for parallel processing
  if (nreps_mitothreshold > ncpus) {
    cl <- makeCluster(ncpus)
    registerDoParallel(cl)
    pernode.reps <- ceiling(nreps_mitothreshold/ncpus)
  } else {
    cl <- makeCluster(nreps_mitothreshold)
    registerDoParallel(cl)
    pernode.reps <- 1
  }

  ######inform the user if it fails to check the stats data
  #check the variables to be removed
  #potentially more variables than observation

  # Loop across nreps_mitothreshold (makes it stable against a few bad fits)
  mito.thresholds  <- foreach(i=1:nreps_mitothreshold, .packages="mixtools", .combine="c") %dopar% {
    # Setup on each node
    thresholds <- c()

    # Loop across replicates per node
    for (m in 1:pernode.reps) {
      # Use mixtools to fit two components to identify a decision boundary
      model <- mixtools::regmixEM(y=stats$mitochondrial_fraction,  x=stats$total_features)

      # Identify the mitochondrial fraction threshold
      high.mito <- rownames(stats)[which(apply(model$posterior,1,FUN="which.max") == which.max(model$beta[1,]))]
      low.mito <- rownames(stats)[which(apply(model$posterior,1,FUN="which.max") == which.min(model$beta[1,]))]
      high.mito <- high.mito[ !(stats[ rownames(stats) %in% high.mito,2] <= median(stats[ rownames(stats) %in% low.mito,2])) ]

      # Return the threshold
      thresholds <- c(thresholds, min(stats[high.mito ,2]))
    }
    # Return the result
    return(thresholds)
  }
  stopCluster(cl)

  # Select a threshold
  # NOTE: Warn the user if the solutions are (very) unstable (less than X% of solutions within 2 SD?)
  sum(mito.thresholds >= median(mito.thresholds) - (2*robustbase::Sn(mito.thresholds)) & mito.thresholds <= median(mito.thresholds) + (2*robustbase::Sn(mito.thresholds))) / length(mito.thresholds)
  mito.threshold <- median(mito.thresholds)

  ## Labels using mix regression
  PU_set = stats %>% mutate(y = if_else(mitochondrial_fraction > mito.threshold, 1, 0))
  PU_set = PU_set %>% mutate(s = if_else(y == 1, 1, 0))

  ### The combined way
  ## Extract feature sets
  features = counts[ , colnames(counts) %in% rownames(PU_set)]
  features = features[ Matrix::rowSums(features) > 0,]

  ## Normalization + Transformation
  features.normalized = Matrix::t(Matrix::t(features) / (Matrix::colSums(features) / 10000))
  features.normalized = log1p(features.normalized)

  ## Feature selection
  fvf = Seurat::FindVariableFeatures(features.normalized)
  fvf = fvf[ order(fvf$vst.variance.standardized),]

  ## Normalization using nullResiduals
  features <- features[ rownames(features) %in% rownames(fvf[ 1:1000,]),] #adjust 1000 genes
  features <- as.matrix(features)
  nr <- nullResiduals(features, fam = "poisson", type = "deviance")

  ## SVD (100 dims)
  svd = irlba::irlba(t(nr), nu = 100, nv = 100)
  unscaled.vectors = svd$u
  scaled.vectors = svd$u %*% diag(svd$d)

  #training_data = cbind(as.numeric(as.character(PU_set$s)), PU_set$mitochondrial_fraction, PU_set$total_features, unscaled.vectors)
  training_data = cbind(as.numeric(as.character(PU_set$s)),  unscaled.vectors)
  training_data = as.data.frame(training_data)
  colnames(training_data) = c("s",paste("SVD_", seq(1,100,1), sep=""))

  #Subset data
  rownames(training_data) = rownames(PU_set)
  trainPositive = subset(training_data, s == 1)
  positive_ = trainPositive[,-1]
  trainUnlabeled = subset(training_data, s == 0)
  unlabeled_ = trainUnlabeled[,-1]

  model_OCC = svm(x = positive_, type = "one-classification", kernel = "radial", gamma = gamma, nu = nu)

  #use the above model to predict outliers in the unlabeled_
  pred_oneclasssvm <- predict(model_OCC, unlabeled_)
  summary(pred_oneclasssvm) #identified 383 labels as bad
  pred_oneclasssvm = as.data.frame(pred_oneclasssvm)
  colnames(pred_oneclasssvm) = "prediction"

  #projection on the previous data
  PU_set_unlabeled = subset(PU_set, y == 0)
  PU_set_labeled = subset(PU_set, y == 1)
  PU_set_labeled$prediction = "TRUE"
  predicted_unlabeled = merge(PU_set_unlabeled, pred_oneclasssvm, by = 0, all = TRUE)
  predicted_unlabeled2 = predicted_unlabeled[,-1]
  rownames(predicted_unlabeled2) = predicted_unlabeled[,1]
  PU_set_new = rbind(predicted_unlabeled2, PU_set_labeled)
  PU_set_new = PU_set_new %>% mutate(s = if_else(prediction == FALSE, 0, 1))

  PU_set_new$col = "darkgrey"
  PU_set_new[ PU_set_new$s == 1,"col"] = "darkgrey"
  PU_set_new[ PU_set_new$s == 0,"col"] = "darkolivegreen3"

  total_count = sum(PU_set_new$s == 0)

  #conditional plotting
  if(plot && directory && format){
    mainDir = getwd()
    subDir = "Cell_classifier"
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
    do.call(
      plot_format,
      args = list(paste(paste("Cell_classifier", sep=""), plot_format, sep = "."))
    )
    plot(PU_set_new$total_features, PU_set_new$mitochondrial_fraction, col = PU_set_new$col, pch = 19,
         main = paste("One-class SVM classifier (gamma = ", gamma, " nu = ", nu, ")","\nHealthy cells: ", total_count), xlab = "Gene expression", ylab = "Mito%")
    dev.off()
  } else if (plot && format){
    do.call(
      plot_format,
      args = list(paste(paste("Cell_classifier", sep=""), plot_format, sep = "."))
    )
    plot(PU_set_new$total_features, PU_set_new$mitochondrial_fraction, col = PU_set_new$col, pch = 19,
         main = paste("One-class SVM classifier (gamma = ", gamma, " nu = ", nu, ")","\nHealthy cells: ", total_count), xlab = "Gene expression", ylab = "Mito%")
    dev.off()
  } else if (plot){
    plot(PU_set_new$total_features, PU_set_new$mitochondrial_fraction, col = PU_set_new$col, pch = 19,
         main = paste("One-class SVM classifier (gamma = ", gamma, " nu = ", nu, ")","\nHealthy cells: ", total_count), xlab = "Gene expression", ylab = "Mito%")
  }

  #output
  classified_cells = PU_set_new[,c(1,2,4)]
  classified_cells = classified_cells %>% mutate(type = if_else(s == 0, "healthy", "dead"))
  classified_cells = classified_cells[,c(1,2,4)]
  return(classified_cells)

}
