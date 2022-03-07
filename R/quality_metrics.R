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
#' @param species A character identifying the species of interest [default = "human"].
#'
#' @details
#' \strong{Using a contrast matrix}\cr
#' The fraction of reads derived from exons is a useful quality metric for both single-cell RNA-seq and single-nucleus RNA-seq. To calculate this metric, it is necessary to provide an additional matrix (the contrast matrix) and specify the constrast type. The exact specification of these two parameters depends on the which reads are included in the main count matrix. This often depends on the type of experiment. As a guideline, for:
#' \itemize{
#'   \item \strong{Single-cell RNA-seq:} The main count matrix often contain only exon-derived reads. The contrast matrix should contains all reads (exons + intron) and the contrast_type should be "numerator".
#'   \item \strong{Single-nucleus RNA-seq:} The main count matrix often contains all reads (exons + intron). The contrast matrix should only contains exon-derived reads and the contrast_type should be "denominator".
#' }
#' @return A data frame object that contains the name of each barcode along with various quality metrics
#' @export
#' @import Matrix
#' @import mygene

quality_metrics = function(counts, contrast = NULL, contrast_type = "denominator", mito = "auto", ribo = "auto", coding = "auto", species = "human") {
  ## evaluate arguments
  # count matrix
  if(missing(counts)) {
    stop('No count matrix was provided', call. = FALSE)
  } else {
    if (!any(class(counts) == c("dgTMatrix", "Matrix","matrix", "dgCMatrix"))) { stop('Count matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix', call. = FALSE) }
  }

  # contrast matrix
  if(!is.null(contrast)) {
    if (!any(class(contrast) == c("dgTMatrix", "Matrix","matrix", "dgCMatrix"))) { stop('Contrast matrix has an unacceptable format. Accepted formats: matrix, Matrix, dgTMatrix, dgCMatrix', call. = FALSE) }
  }

  # contrast_type argument
  if (!any(contrast_type == c("denominator", "numerator"))) { stop('contrast_type should be set to either "denomiator" or "numerator"', call. = FALSE) }

  # mito argument
  if (mito != "auto") {
    if (sum(rownames(counts) %in% mito) != length(mito)) { stop('The count matrix does not contain all of the entered mitochondrial genes', call. = FALSE) }
  }

  # ribo argument
  if (ribo != "auto") {
    if (sum(rownames(counts) %in% ribo) != length(ribo)) { stop('The count matrix does not contain all of the entered ribosomal genes', call. = FALSE) }
  }

  # coding argument
  if (coding != "auto") {
    if (sum(rownames(counts) %in% coding) != length(coding)) { stop('The count matrix does not contain all of the entered protein-coding genes', call. = FALSE) }
  }

  ## convert the counts into dgCMatrix if its class() is not dgCMatrix
  if(class(counts) != "dgCMatrix") { counts = as(counts, "dgCMatrix") }

  ## create a list for holding the output
  output <- list()

  ## extract the names of all expressed genes
  genes <- data.frame(name = rownames(counts))

  ## find species and scope if either ribo or mito is set to auto
  if (ribo == "auto" | mito == "auto") {
    ## strip version names of gene ids (dot separated)
    lengths <- data.frame(dot = as.numeric(regexpr("\\.", genes$name))-1, total = nchar(as.character(genes$name)))
    lengths$use <- apply(lengths, 1, FUN = function(x) min(x[x > 0]))
    genes$clean <- substr(genes$name, 0, lengths$use)

    ## determine the species using the first 100 genes to query mygene
    out <- capture.output(query <- mygene::queryMany(qterms = genes[1:100,"clean"], scopes = c("entrezgene", "ensembl.gene","symbol","name","alias","refseq","unigene","ensembl.transcript","accession","xenbase","zfin","wormbase","flybase","rgd","mgi","hgnc")))
    query <- query[ order(query$query, -query$X_score),]
    query <- query[ duplicated(query$query)==F,]
    species_id <- names(which.max(table(query$taxid)))

    ## determine the scope (i.e. annotator) using the first 100 genes to query mygene
    scopes <- data.frame()
    counter <- 1
    for (scope in c("entrezgene", "ensembl.gene","symbol","name","alias","refseq","unigene","ensembl.transcript","accession","xenbase","zfin","wormbase","flybase","rgd","mgi","hgnc")) {
      out <- capture.output(test <- mygene::queryMany(qterms = genes[1:100,"clean"], scopes = scope, species = species_id, return.as = "records"))
      scopes[counter,1] <- scope
      scopes[counter,2] <- sum(unlist(lapply(test,FUN="length")) > 2)
      counter <- counter + 1
    }
    scope <- scopes[ which.max(scopes[,2]),1]
  }

  ## process mitochondrial genes
  if (mito == "auto") {
    ## lookup mitochondrial genes
    mitoquery <- as.data.frame(mygene::query("'chrMT:1-100,000,000'", species=species_id, fields = scope, size = 1000)$hits)
    mitogenes <- genes[ genes$clean %in% unlist(mitoquery[,3]),1]
  }  else {
    mitogenes <- mito
  }

  if (length(mitogenes) == 0){
    species_data = get(data(list = species))
    mitogenes_species = subset(species_data, Chromosome.scaffold.name %in% c("mt", "MT", "MTDNA", "mitochondrion_genome"))
    for(i in colnames(mitogenes_species)){
      mitogenes_search <- genes[genes$clean %in% mitogenes_species[,i],]
      if(!(nrow(mitogenes_search) == 0)){
        mitogenes <- mitogenes_search[,1]
      }
    }
  }

  ## process ribosomal genes
  if (ribo == "auto") {
    ## lookup ribosomal genes
    ribogenes <- as.data.frame(mygene::query(q = 'name:ribosomal AND name:protein AND symbol:rp* NOT name:mitochondrial', species=species_id, size=1000, fields = scope)$hits)
    ribogenes <- genes[ genes$clean %in% unlist(ribogenes[,3]),1]
  } else {
    ribogenes <- ribo
  }

  ## process non-coding genes
  if (coding == "auto") {
    ## lookup all nonzero genes and retrieve the type of gene from Ensembl
    #nz <- names(which(rowMeans(counts) > 0)) #gives error: 'x' must be an array of at least two dimensions
    nz <- names(Matrix::rowMeans(counts > 0))
    out <- capture.output(pcgenes <- suppressMessages(as.data.frame(mygene::queryMany(genes[ genes$name %in% nz, "clean"], species=species_id, fields = c(scope, "ensembl.type_of_gene")))))
    pcgenes <- do.call(rbind.data.frame, pcgenes[,4])
    pcgenes <- genes[genes$clean %in% pcgenes[ pcgenes$type_of_gene == "protein_coding",1],1]
  } else {
    pcgenes <- coding
  }

  ## if mygene fails invoke provided data
  if(length(pcgenes) == 0){
    species_data = get(data(list = species))
    pc_species = subset(species_data, Gene.type %in% c("protein_coding"))
    for(i in colnames(pc_species)){
      pc_search <- genes[genes$clean %in% pc_species[,i],]
      if(!(nrow(pc_search) == 0)){
        pcgenes <- pc_search[,1]
      }
    }
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
