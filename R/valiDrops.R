#' @title valiDrops
#'
#' @description Wrapper function that runs the entire valiDrops flow with default parameters.
#'
#' @param counts A matrix containing counts for barcodes passing the barcode-rank threshold.
#' @param rank_barcodes A boolean (TRUE or FALSE) indicating whether or not to rank barcodes [default = TRUE].
#' @param mitochondrial_clusters NULL or an interger indicating how many deviations a cluster must be to be removed [default = 3]. 
#' @param ribosomal_clusters NULL or an interger indicating how many deviations a cluster must be to be removed [default = 3]. 
#' @param label_dead A boolean (TRUE or FALSE) indicating whether or not to label putative dead cells [default = FALSE].
#' @param status A boolean (TRUE or FALSE) indicating whether or not to print the progress of valiDrops [default = TRUE].
#' @param stageThree A boolean (TRUE or FALSE) indicating whether or not to run stage 3 of valiDrops [default = TRUE].
#' @param ... Pass parameters to functions within valiDrops. See \link{rank_barcodes} \link{quality_metrics} \link{quality_filter} \link{expression_metrics} \link{expression_filter} \link{label_dead}
#'
#' @return A data frame containing quality metrics, as well as quality control labels (and if requested, dead labels) for all barcodes passing the rank threshold.
#' @export
#' @import R.utils
#' @import Matrix
#' @import Seurat
#' @import SingleCellExperiment

valiDrops = function(counts, filtered_counts = NULL, rank_barcodes = TRUE, status = TRUE, mito_fixed = NULL, stageThree = TRUE, label_dead = FALSE, plot = TRUE, verbose = TRUE, tol = 1e-50, maxit.glm = 2500, h = 0.01, timeout = Inf,
                    type = "UMI", psi.min = 2, psi.max = 5, alpha = 0.001, alpha.max = 0.05, boot = 10, factor = 1.5, threshold = TRUE,
                    contrast = NULL, contrast_type = "denominator", species = "auto", annotation = "auto", mito = "auto", ribo = "auto", coding = "auto",
                    mitol = TRUE, distancel = TRUE, codingl = TRUE, contrastl = FALSE, mito.nreps = 10, mito.max = 0.3, npsi = 3, dist.threshold = 5, coding.threshold = 3, contrast.threshold = 3,
                    nfeats = 5000, npcs = 30, k.min = 5, res.shallow = 0.1, top.n = 10,
                    mitochondrial_clusters = 3, ribosomal_clusters = 3, min.significant = 1, min.target.pct = 0.3, max.background.pct = 0.7, min.diff.pct = 0.2, min.de.frac = 0.01, min.significance.level = NULL,
                    cor.threshold = NULL, train = TRUE, rep = 10, n.min = 8, n.relabel = 1, feature.try = 3, label.thrs = NULL, label.frac = 0.1, nfeats2 = 2000, alpha2 = 0, npcs2 = 100, weight = TRUE, epochs = 20, nfolds = 5, nrep = 10, fail.weight = 0.2, cor.min = 0.0001, cor.max = 0.005, cor.steps = 50, nrep.cor = 10, min.dead = 100, max.live = 500, bpparam = SerialParam()) {
  ## Check the rank_barcodes parameter
  if (!isTRUE(rank_barcodes) & !isFALSE(rank_barcodes)) { stop("rank_barcodes must be either TRUE or FALSE") }
  
  ## Check the label_dead parameter
  if (!isTRUE(label_dead) & !isFALSE(label_dead)) { stop("label_dead must be either TRUE or FALSE") }
  
  ## Check the stageThree parameter
  if (!isTRUE(stageThree) & !isFALSE(stageThree)) { stop("stageThree must be either TRUE or FALSE") }
  
  ## Extract counts from Seurat or SCE objects
  if (any(class(counts) %in% c("SingleCellExperiment", "Seurat"))) {
    if (class(counts) == "SingleCellExperiment") {
      counts <- SingleCellExperiment::counts(counts)
    } else {
      counts <- Seurat::GetAssayData(counts, slot = "counts")
    }
  }
  
  ## Validate the mitochondrial_clusters argument
  if (!is.null(mitochondrial_clusters) & !is.numeric(mitochondrial_clusters) & !is.integer(mitochondrial_clusters)) {
    stop("The mitochondrial_clusters argument must be either NULL or a numeric", call. = FALSE) 	  
  }
  
  ## Validate the ribosomal_clusters argument
  if (!is.null(ribosomal_clusters) & !is.numeric(ribosomal_clusters) & !is.integer(ribosomal_clusters)) {
    stop("The ribosomal_clusters argument must be either NULL or a numeric", call. = FALSE) 	  
  }
  
  ## Validate the counts object
  if (missing(counts)) {
    stop('No count matrix was provided', call. = FALSE)
  } else {
    if (!any(class(counts) == c("dgTMatrix", "Matrix","matrix", "dgCMatrix", "DelayedMatrix"))) { stop('Count matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix, and DelayedMatrix', call. = FALSE) }
  }
  
  ## convert the counts into dgCMatrix if its a dense matrix
  if (class(counts) == "matrix") { 
    counts = as(counts, "dgCMatrix")
  }

  ## Check for column and rownames
  if (is.null(colnames(counts))) {
    if (status) { message("Input check: No column names found. Setting to Setting cell names to cell_columnidx (e.g 'cell_1')") }
    colnames(counts) <- paste("cell", 1:ncol(counts), sep="_")
  }
  if (is.null(rownames(counts))) {
    if (status) { message("Input check: No row names found. Setting to Setting cell names to feature_rowidx (e.g 'feature_1')") }
    rownames(counts) <- paste("feature", 1:nrow(counts), sep="_")
  }
  
  ## Run rank_barcodes
  if (rank_barcodes) {
    if (status) { message("Step 1: Filtering on the barcode-rank plot.") }
    threshold <- R.utils::withTimeout( {valiDrops::rank_barcodes(counts = counts, type = type, psi.min = psi.min, psi.max = psi.max, alpha = alpha, alpha.max = alpha.max, boot = boot, factor = factor, threshold = threshold, plot = plot) },
                                      timeout = timeout,
                                      onTimeout = "error")
    rank.pass <- rownames(threshold$ranks[ threshold$ranks$counts >= threshold$lower.threshold,])
  } else {
    if (status) { message("Step 1: Removing barcodes with zero counts.") }
    rank.pass <- colnames(counts)[which(colSums(counts) > 0)]
  }
  
  ## Check the number of barcodes that pass initial filtering
  if (length(rank.pass) > 20000) {
    message("More than 20.000 barcodes passed initial filtering. It is like that breakpoint estimation did not work satisfactorily. If it didn't, you can try to increase alpha, alpha.max and/or psi.max in the rank_barcodes() function.")
  }
  if (length(rank.pass) < ncol(filtered_counts)) {
    to_continue <- readline("Less barcodes passed ranking than in filtered count matrix. Type [Y] if you would like to continue and [N] if not: ")
    if (to_continue %in% c("N", "No", "n", "no")) {
        stop("Pipeline terminated by user. Please try again on filtered counts and skip barcode ranking.")
      }
  }
  
  ## Subset the counts
  counts.subset <- counts[, colnames(counts) %in% rank.pass]

  ## Run quality_metrics
  if (status) { message("Step 2: Collecting quality metrics.") }
  metrics <- valiDrops::quality_metrics(counts = counts.subset, contrast = contrast, contrast_type = contrast_type, species = species, annotation = annotation, mito = mito, ribo = ribo, coding = coding, verbose = verbose)
  
  ## Run quality_filter
  if (status) { message("Step 3: Filtering on quality metrics.") }
  qc.pass <- R.utils::withTimeout( {valiDrops::quality_filter(metrics = metrics$metrics, mito = mitol, distance = distancel, coding = codingl, contrast = contrastl, mito.nreps = mito.nreps, mito.max = mito.max, npsi = npsi, dist.threshold = dist.threshold, coding.threshold = coding.threshold, contrast.threshold = contrast.threshold, plot = plot, tol = tol, maxit.glm = maxit.glm, h = h) },
                                   timeout = timeout,
                                   onTimeout = "error")
  if (!is.null(filtered_counts)) {
    if (length(qc.pass$final) >= 2*ncol(filtered_counts)) {
      to_continue <- readline(prompt = "More than twice as many cells detected than in the filtered matrix. Type [Y] if you would like to continue and [N] if not: ")
      if (to_continue %in% c("N", "No", "n", "no")) {
        stop("Pipeline terminated by user. Please try again on filtered counts and skip barcode ranking.")
      }
    }
    if (length(qc.pass$final) <= 0.5*ncol(filtered_counts)) {
      to_continue <- readline(prompt = "Less than half as many cells detected than in the filtered matrix.\nType [Y] if you would like to refilter with a fixed proportion of mitochondrial reads and [N] if not: ")
      if (to_continue %in% c("N", "No", "n", "no")) {
        message("Very well. Continuing pipeline.")
      } else{
        if (is.null(mito_fixed)) {
          mito_fixed <- as.numeric(readline(prompt = "Fixed proportion of mitochondrial reads not provided.\nWhat percentage of mitochondrial reads would you like to use: "))
        }
        if (mito_fixed < 1) {mito_fixed <- mito_fixed*100}
        message(paste0("Repeating Step 3: Filtering on quality metrics with ", mito_fixed, "% mitochondrial reads."))
        qc.pass <- R.utils::withTimeout( {valiDrops::quality_filter(metrics = metrics$metrics, mito = mito_fixed/100, distance = distancel, coding = codingl, contrast = contrastl, mito.nreps = mito.nreps, mito.max = mito.max, npsi = npsi, dist.threshold = dist.threshold, coding.threshold = coding.threshold, contrast.threshold = contrast.threshold, plot = plot, tol = 1e-100, maxit.glm = 10000, h = 1e-5) },
                                   timeout = timeout,
                                   onTimeout = "error")
      }
    }
  }
      
  if (stageThree) {
    ## Convert counts to Seurat capatible format
    if (status) { message("Step 4: Collecting expression-based metrics.") }
    counts.subset.filtered <- counts.subset[ rownames(counts.subset) %in% metrics$protein_coding, colnames(counts.subset) %in% qc.pass$final]
    if (class(counts.subset.filtered) != "dgCMatrix") { 
      counts.subset.filtered = as(counts.subset.filtered, "dgCMatrix")
    }

    ## Run expression_matrix
    expr.metrics <- R.utils::withTimeout( {valiDrops::expression_metrics(counts = counts.subset.filtered, mito = metrics$mitochondrial, ribo = metrics$ribosomal, nfeats = nfeats, npcs = min(npcs, ncol(counts.subset.filtered)), k.min = k.min, res.shallow = res.shallow, top.n = top.n) },
                                         timeout = timeout*3,
                                         onTimeout = "error")
    
    ## Run expression_filter
    if (status) { message("Step 5: Filtering on expression-based metrics.") }
    valid <- R.utils::withTimeout( {valiDrops::expression_filter(stats = expr.metrics$stats, clusters = expr.metrics$clusters, mito = mitochondrial_clusters, ribo = ribosomal_clusters, min.significant = min.significant, min.target.pct = min.target.pct, max.background.pct = max.background.pct, min.diff.pct = min.diff.pct, min.de.frac = min.de.frac, min.significance.level = min.significance.level, plot = plot, tol = tol, maxit.glm = maxit.glm, h = h) },
                                  timeout = timeout,
                                  onTimeout = "error")
  } else {
    valid <- qc.pass$final
  }
  
  ## Setup the results
  met <- metrics$metrics
  met$qc.pass <- "fail"
  met[ met$barcode %in% valid, "qc.pass"] <- "pass"
  
  ## Run label_dead
  if (label_dead) {
    if (status) { 
      if (stageThree) { message("Step 6: Predicting dead cells.") } else { message("Step 4: Predicting dead cells without running stage 3. CAUTION this has not been tested.") }
    }
    dead <- valiDrops::label_dead(counts = counts, metrics = met, qc.labels = setNames(as.character(met$qc), met$barcode), cor.threshold = cor.threshold, train = train, rep = rep, n.min = n.min, n.relabel = n.relabel, feature.try = feature.try, verbose = FALSE, label.thrs = label.thrs, label.frac = label.frac, nfeats = nfeats2, alpha = alpha2, npcs = npcs2, weight = weight, epochs = epochs, nfolds = nfolds, nrep = nrep, fail.weight = fail.weight, cor.min = cor.min, cor.max = cor.max, cor.steps = cor.steps, nrep.cor = nrep.cor, min.dead = min.dead, max.live = max.live, plot = plot, bpparam = bpparam)	  
    dead <- dead$metrics
    met$label <- "live"
    met[ met$barcode %in% dead[ dead$label == "dead","barcode"], "label"] <- "dead"
    met[ met$barcode %in% dead[ dead$label == "uncertain","barcode"], "label"] <- "uncertain"
  }
  
  ## Output to user
  if (status) { message(paste("\t", nrow(met[ met$qc.pass == "pass",]), " barcodes passed quality control.", sep="")) }
  if (label_dead) { if (status) { message(paste("\t", nrow(met[ met$qc.pass == "pass" & met$label == "dead",]), " barcodes that passed quality control are predicted to be dead.", sep="")) } }
  
  ## Return
  return(met)
}
