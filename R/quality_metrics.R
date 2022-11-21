#' @title quality_metrics
#'
#' @description Function to calculate various quality metrics for barcodes passing the barcode-rank threshold
#'
#' @param counts A matrix containing counts for barcodes passing the barcode-rank threshold. See \link{rank_barcodes}
#' @param contrast An optional matrix containing contrasting counts for barcodes passing the barcode-rank threshold [default = NULL]. See details for more information.
#' @param contrast_type A string ("denominator" or "numerator") indicating the type of contrast to calculate [default = "denominator"]. See details for more information.
#' @param mito A vector of the identifiers for mitochondrial genes or "auto". If set to "auto", valiDrops will automatically find the identifiers of mitochondrial genes [default = "auto"].
#' @param ribo A vector of the identifiers for ribosomal genes or "auto". If set to "auto", valiDrops will automatically find the identifiers of ribosomal genes [default = "auto"].
#' @param coding A vector of the identifiers for protein-coding genes or "auto". If set to "auto", valiDrops will automatically find the identifiers of protein-coding genes [default = "auto"].
#' @param species A character identifying the species of origin or "auto". See details for more information. If set to "auto", valiDrops will automatically determine the species if possible [default = "auto"].
#' @param annotation A character identifying the gene annotation or "auto". See details for more information. If set to "auto", valiDrops will automatically determine the gene annotation if possible [default = "auto"].
#' @param verbose A boolean (TRUE or FALSE) indicating whether or not to be verbose [default = TRUE].
#'
#' @details
#' \strong{Using a contrast matrix}\cr
#' The fraction of reads derived from exons is a useful quality metric for both single-cell RNA-seq and single-nucleus RNA-seq. To calculate this metric, it is necessary to provide an additional matrix (the contrast matrix) and specify the constrast type. The exact specification of these two parameters depends on the which reads are included in the main count matrix. This often depends on the type of experiment. As a guideline, for:
#' \itemize{
#'   \item \strong{Single-cell RNA-seq:} The main count matrix often contain only exon-derived reads. The contrast matrix should contains all reads (exons + intron) and the contrast_type should be "numerator".
#'   \item \strong{Single-nucleus RNA-seq:} The main count matrix often contains all reads (exons + intron). The contrast matrix should only contains exon-derived reads and the contrast_type should be "denominator".
#' }
#'
#' \strong{Available species}\cr
#' Current valiDrops supports the following species: "human", "mouse", "rat", "C.elegans", "drosophila" and "zebrafish"
#'
#' \strong{Available gene annotations}\cr
#' Current valiDrops supports the following gene annotations: "symbol", "ensembl", "entrez", "HGNC" and "MGI"
#' @return A data frame object that contains the name of each barcode along with various quality metrics
#' @export
#' @import Matrix

quality_metrics = function(counts, contrast = NULL, contrast_type = "denominator", species = "auto", annotation = "auto", mito = "auto", ribo = "auto", coding = "auto", verbose = TRUE) {
  ## evaluate arguments
  # count matrix
  if(missing(counts)) {
    stop('No count matrix was provided', call. = FALSE)
  } else {
    if (!any(class(counts) == c("dgTMatrix", "Matrix","matrix", "dgCMatrix"))) { stop('Count matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix', call. = FALSE) }
  }

  # species arugment
  if (!any(species %in% c("Human","human","sapiens","Sapiens","h.sapiens", "H.Sapiens","H.sapiens","Mouse","mouse","Musculus","musculus","m.musculus","M.Musculus","M.musculus", "rat","Rat","Norvegicus","norvegicus","r.norvegicus","R.Norvegicus","R.norvegicus","worm","Worm","elegans","Elegans","C.Elegans","c.elegens","C.elegans","fly","Fly","drosophila","Drosophila","D.melanogaster","d.melanogaster","D.Melanogaster","zebrafish","Zebrafish","D.rerio","d.rerio","D.Rerio", "auto"))) {
    stop('Species must be either "auto", "human", "mouse","rat", "C.elegans", "drosophila" or "zebrafish"', call. = FALSE)
  }

  # annotation arugment
  if (!any(annotation %in% c("Symbol", "symbol","entrez","Entrez","NCBI","Ensembl","ensembl","HGNC","MGI", "auto"))) {
    stop('annotation must be either "auto", "symbol", "ensembl" or "entrez"', call. = FALSE)
  }

  # contrast matrix
  if(!is.null(contrast)) {
    if (!any(class(contrast) == c("dgTMatrix", "Matrix","matrix", "dgCMatrix"))) { stop('Contrast matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix', call. = FALSE) }
  }

  # contrast_type argument
  if (!any(contrast_type == c("denominator", "numerator"))) { stop('contrast_type should be set to either "denomiator" or "numerator"', call. = FALSE) }

  # mito argument
  if (mito != "auto") {
    if (sum(rownames(counts) %in% mito) != length(mito)) { stop('The count matrix does not contain all of the entered mitochondrial genes. Continueing', call. = TRUE) }
  }

  # ribo argument
  if (ribo != "auto") {
    if (sum(rownames(counts) %in% ribo) != length(ribo)) { stop('The count matrix does not contain all of the entered ribosomal genes. Continueing', call. = TRUE) }
  }

  # coding argument
  if (coding != "auto") {
    if (sum(rownames(counts) %in% coding) != length(coding)) { stop('The count matrix does not contain all of the entered protein-coding genes. Continueing', call. = TRUE) }
  }

  ## Check the verbose parameter
  if (!isTRUE(verbose) & !isFALSE(verbose)) { stop("verbose must be either TRUE or FALSE") }

  ## convert the counts into dgCMatrix if its class() is not dgCMatrix
  if(class(counts) != "dgCMatrix") { counts = as(counts, "dgCMatrix") }

  ## create a list for holding the output
  output <- list()

  ## Create a list
  ds <- list()
  ds[[1]] <- as.data.frame(annotation[[1]])
  ds[[2]] <- as.data.frame(annotation[[2]])
  ds[[3]] <- as.data.frame(annotation[[3]])
  ds[[4]] <- as.data.frame(annotation[[6]])
  ds[[5]] <- as.data.frame(annotation[[5]])
  ds[[6]] <- as.data.frame(annotation[[4]])

  ## Determine dataset
  if (species == "auto") {
    ds.range <- 1:6
  } else if (any(species %in% c("Human","human","sapiens","Sapiens","h.sapiens", "H.Sapiens","H.sapiens"))) {
    ds.range <- 1
  } else if (any(species %in% c("Mouse","mouse","Musculus","musculus","m.musculus","M.Musculus","M.musculus"))) {
    ds.range <- 2
  } else if (any(species %in% c("rat","Rat","Norvegicus","norvegicus","r.norvegicus","R.Norvegicus","R.norvegicus"))) {
    ds.range <- 3
  } else if (any(species %in% c("zebrafish","Zebrafish","D.rerio","d.rerio","D.Rerio"))) {
    ds.range <- 4
  } else if (any(species %in% c("worm","Worm","elegans","Elegans","C.Elegans","c.elegens","C.elegans"))) {
    ds.range <- 5
  } else if (any(species %in% c("fly","Fly","drosophila","Drosophila","D.melanogaster","d.melanogaster","D.Melanogaster"))) {
    ds.range <- 6
  }

  ## extract the names of all expressed genes
  genes <- data.frame(name = rownames(counts))
  genes$clean <- genes$name

  ## strip dots if the input annotation is ensembl (dots to use versions of Ensembl genes is used by GENCODE)
  ens.idx <- c(grep("^ENSG00", genes$name), grep("^ENSMUSG00", genes$name))
  if (length(ens.idx) > 0) {
    ens.id <- genes[ ens.idx, "name"]
    lengths <- data.frame(dot = as.numeric(regexpr("\\.", ens.id))-1, total = nchar(as.character(ens.id)))
    lengths$use <- apply(lengths, 1, FUN = function(x) min(x[x > 0]))
    ens.clean <- substr(genes[ens.idx, "name"], 0, lengths$use)
    genes[ ens.idx, "clean"] <- ens.clean
  }

  ## Loop across selected datasets
  ann.stats <- as.data.frame(matrix(ncol=3, nrow=length(ds.range)))
  counter <- 1
  for (ds.sel in ds.range) {
    tmp <- ds[[ds.sel]]
    if (annotation == "auto") {
      for (col in 1:ncol(tmp)) {
        ann.stats[counter,1] <- ds.sel
        ann.stats[counter,2] <- col
        ann.stats[counter,3] <- sum(genes$clean %in% tmp[,col])
        counter <- counter + 1
      }
    } else {
      if (annotation %in% c("Symbol", "symbol")) { col <- which(colnames(tmp) == "Symbol") }
      if (annotation %in% c("entrez", "Entrez", "NCBI")) { col <- which(colnames(tmp) == "NCBI") }
      if (annotation %in% c("Ensembl", "ensembl")) { col <- which(colnames(tmp) == "Ensembl") }
      if (annotation %in% c("HGNC")) { col <- which(colnames(tmp) == "HGNC") }
      if (annotation %in% c("MGI")) { col <- which(colnames(tmp) == "MGI") }
      ann.stats[counter,1] <- ds.sel
      ann.stats[counter,2] <- col
      ann.stats[counter,3] <- sum(genes$clean %in% tmp[,col])
    }
  }

  ## Get the best hits
  best.ds <- ann.stats[ which.max(ann.stats[,3]),1]
  best.col <- ann.stats[ which.max(ann.stats[,3]),2]
  best.sum <- ann.stats[ which.max(ann.stats[,3]),3]

  ## Protein-coding genes
  if (coding == "auto") {
    pcgenes <- ds[[best.ds]][ ds[[best.ds]]$Type %in% "protein_coding",best.col]
    pcgenes <- genes[ genes$clean %in% pcgenes,1]
  } else {
    pcgenes <- coding
  }

  ## Mitochondrial genes
  if (mito == "auto") {
    mitogenes <- ds[[best.ds]][ ds[[best.ds]]$Chr %in% c("MtDNA","MT","mitochondrion_genome"),best.col]
    mitogenes <- genes[ genes$clean %in% mitogenes,1]
  } else {
    mitogenes <- mito
  }

  ## Ribosomal genes
  if (ribo == "auto") {
    ribogenes <- ds[[best.ds]][ c(grep("^rpl", tolower(ds[[best.ds]]$Symbol)),grep("^rps", tolower(ds[[best.ds]]$Symbol))),best.col]
    ribogenes <- genes[ genes$clean %in% ribogenes,1]
  } else {
    ribogenes <- ribo
  }

  ## Tell the user about the match
  if (verbose) {
    species <-  c("Human", "Mouse","Rat","Zebrafish","Worm","Fly")[best.ds]
    annotation <- colnames(ds[[best.ds]])[best.col]
    if (annotation == "NCBI") { annotation <- "Entrez" }
    message(paste("Detected sample origin: ", species, ". Detected gene annotation: ", annotation, ". Mapped ", best.sum, "/", nrow(genes), " (", signif((best.sum/nrow(genes))*100,3), "%) of input IDs.", sep=""))
    message(paste("Found ", length(mitogenes), " mitochondrial genes, ", length(ribogenes), " ribosomal genes, and ", length(pcgenes), " protein-coding genes.", sep=""))
  }

  ## calculate quality metrics
  metrics <- as.data.frame(matrix(ncol=6, nrow=ncol(counts)))
  colnames(metrics) <- c("barcode","logUMIs","logFeatures","mitochondrial_fraction","ribosomal_fraction","coding_fraction")
  metrics[,1] <- colnames(counts)
  metrics[,2] <- Matrix::colSums(counts)
  metrics[,3] <- log(Matrix::colSums(counts > 0))
  metrics[,4] <- Matrix::colSums(counts[ rownames(counts) %in% mitogenes,]) / metrics[,2]
  metrics[,5] <- Matrix::colSums(counts[ rownames(counts) %in% ribogenes,]) / metrics[,2]
  metrics[,6] <- Matrix::colSums(counts[ rownames(counts) %in% pcgenes,]) / metrics[,2]
  metrics[,2] <- log(metrics[,2])

  ## calculate metric using the contrast matrix
  if (!is.null(contrast)) {
    # Match barcodes between the input counts and the contrast matrix
    contrast <- contrast[, colnames(contrast) %in% colnames(counts)]
    contrast <- contrast[, match(colnames(counts), colnames(contrast))]
    if (contrast_type == "denominator") {
      metrics[,7] <- colSums(counts) / colSums(contrast)
    } else if (contrast_type == "numerator") {
      metrics[,7] <- colSums(contrast) / colSums(counts)
    }
    colnames(metrics)[7] <- "contrast_fraction"
  }

  ## setup the output
  output$metrics <- metrics
  output$mitochondrial <- mitogenes
  output$ribosomal <- ribogenes
  output$protein_coding <- pcgenes

  ## return results
  return(output)
}
