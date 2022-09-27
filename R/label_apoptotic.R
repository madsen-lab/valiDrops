#' @title label_apoptotic
#'
#' @description Function that uses Isolation Forest to label apoptotic cells.
#'
#' @param counts A matrix containing counts for barcodes passing the barcode-rank threshold. See \link{rank_barcodes}
#' @param metrics A data frame containing quality metrics. See \link{quality_metrics} and \link{quality_filter}
#' @param mitochondrial.threshold A mitochondrial threshold identified using quality_filter method [default = 0.5].
#' @param threshold A number indicating the score to use for soft labeling of apoptotic cells [default = 0.2].
#' @param nfeats A number indicating the number of variable features to use for clustering [default = 5000].
#' @param npcs A number indicating the number of singular values to use for clustering [default = 100].
#' @param seed Set a seed for reproducible results. Set to NULL to run in non-deterministic mode [default = 42].
#' @param verbose Set to TRUE for verbose mode. Set to FALSE for silent mode [default = TRUE].
#'
#' @return A data frame object that contains the quality metrics of the data and results indicating the apoptotic cells.
#' @export
#' @import Matrix
#' @import pROC
#' @import glmnet
#' @import pcaPP
#' 
label_apoptotic <- function(counts, metrics, qc.labels, cor.threshold = NULL, verbose = TRUE, label.thrs = 13.5, nfeats = 2000, alpha = 0, npcs = 100, weight = TRUE, epochs = 10, nfolds = 5, nrep = 50, fail.weight = 0.2, cor.min = 0.0001, cor.max = 0.005, cor.steps = 50, nrep.cor = 20, min.apoptotic = 100, max.healthy = 500) {
  # Soft label the dataset
  metrics$logUMIs <- scale(metrics$logUMIs, scale = FALSE)
  metrics$logFeatures <- scale(metrics$logFeatures, scale = FALSE)
  metrics$ribosomal_fraction <- asin(sqrt((metrics$ribosomal_fraction)))/(pi/2)
  metrics$coding_fraction <- asin(sqrt((metrics$coding_fraction)))/(pi/2)
  metrics$mitochondrial_fraction <- asin(sqrt(metrics$mitochondrial_fraction))/(pi/2)
  metrics$score <- metrics$logUMIs * -11.82 + metrics$logFeatures * 2.08 + metrics$ribosomal_fraction * 158.98 + metrics$logFeatures * metrics$coding_fraction * 18.87 + metrics$ribosomal_fraction * metrics$coding_fraction * -125.9
  metrics$label <- "healthy"
  metrics[ metrics$score <= label.thrs,"label"] <- "apoptotic"
  metrics$label <- factor(metrics$label, levels = c("healthy","apoptotic"))

  # Define tracking column
  metrics$track = metrics$label

  # Break if there are not enough observations
  if (nrow(metrics[ metrics$label == "apoptotic",]) < 50) {
    if (verbose) { message("Soft-labeling identified less than 50 apoptotic barcodes. Aborting", sep="") }
    metrics$label <- "healthy"
    return(metrics)
  } else {
    if (verbose) { message("Soft-labeling identified ", nrow(metrics[ metrics$label == "apoptotic",])," apoptotic barcodes.", sep="") }
  }

  # Normalize and log transform counts
  counts <- counts[ , colnames(counts) %in% metrics$barcode]
  nonzero <-  counts[ Matrix::rowSums(counts) > 0,]
  sf <- 10000 / Matrix::colSums(nonzero)
  norm_transform <- Matrix::t(Matrix::t(nonzero) * sf)
  norm_transform@x <- log1p(norm_transform@x)

  # Find variable features
  dev <- scry::devianceFeatureSelection(as.matrix(nonzero))
  var.feats <- names(which(rank(-dev) <= nfeats))

  # Scaling and SVD
  means <- Matrix::colMeans(nonzero)
  nr <- nrow(nonzero)
  sds <- sqrt((Matrix::colMeans(nonzero*nonzero) - Matrix::colMeans(nonzero)^2) * (nr / (nr-1)))
  sds[ sds == 0 ] <- 1
  data.scaled <- Matrix::t((Matrix::t(nonzero) - means)/sds)

  # Setup for training
  Y <- factor(metrics$label)
  fraction_0 <- rep(1 - sum(Y == "apoptotic") / length(Y), sum(Y == "apoptotic"))
  fraction_1 <- rep(1 - sum(Y == "healthy") / length(Y), sum(Y == "healthy"))
  weights <- numeric(length(Y))
  if (weight) {
    weights[Y == "apoptotic"] <- fraction_0
    weights[Y == "healthy"] <- fraction_1
    score_weights <- log1p(abs(metrics$score - label.thrs))
    score_weights <- (score_weights - min(score_weights)) / (max(score_weights) - min(score_weights))
    qc_weight <- rep(1, length(Y))
    qc_weight[metrics$qc == "fail"] <- fail.weight
    weights <- weights * score_weights * qc_weight
  } else {
    weights[1:length(weights)] <- 1
  }

  # Calculate SVDs
  svd <- irlba::irlba(Matrix::t(data.scaled), nv=npcs, nu=npcs)

  # Determine correlation threshold
  if (is.null(cor.threshold)) {
    if (verbose) { message("Determining feature selection threshold.") }
    stats <- as.data.frame(matrix(ncol=6, nrow = cor.steps))
    counter <- 1
    sequence <- 2^seq(log2(cor.min), log2(cor.max), length.out = cor.steps)
    for (cor.threshold in sequence) {
      # Set probabilities
      metrics$prob = 1
      metrics[ metrics$label == "healthy", "prob"] <- 0

      # Extract SVDs
      X = svd$u[,1:npcs]
      rownames(X) <- colnames(data.scaled)

      # Subset SVDs
      cor.coef <- c()
      for (dim in 1:npcs) { cor.coef <- c(cor.coef, pcaPP::cor.fk(X[,dim], metrics$label)) }
      X = X[,which(cor.coef^2 >= cor.threshold)]

      # Loop across replicates
      prob.mat <- as.data.frame(matrix(ncol = nrep.cor, nrow = nrow(metrics)))
      for (rep in 1:nrep.cor) {
        # Extract and combine the samples
        idx.apoptotic.fail <- as.numeric(sample(x=as.character(which(metrics$qc == "fail" & metrics$label == "apoptotic")), size=max(min.apoptotic, length(which(metrics$qc == "fail" & metrics$label == "apoptotic"))), replace=TRUE, prob=metrics[ which(metrics$qc == "fail" & metrics$label == "apoptotic"), "prob"]))
        idx.apoptotic.pass <- sample(x=which(metrics$qc == "pass" & metrics$label == "apoptotic"), size=max(min.apoptotic, length(which(metrics$qc == "pass" & metrics$label == "apoptotic"))), replace=TRUE, prob=metrics[ which(metrics$qc == "pass" & metrics$label == "apoptotic"), "prob"])
        idx.healthy.fail <- sample(x=which(metrics$qc == "fail" & metrics$label == "healthy"), size=min(max.healthy, length(which(metrics$qc == "fail" & metrics$label == "healthy"))), replace=TRUE, prob=abs(metrics[ which(metrics$qc == "fail" & metrics$label == "healthy"), "prob"]-1))
        idx.healthy.pass <- sample(x=which(metrics$qc == "pass" & metrics$label == "healthy"), size=min(max.healthy, length(which(metrics$qc == "pass" & metrics$label == "healthy"))), replace=TRUE, prob=abs(metrics[ which(metrics$qc == "pass" & metrics$label == "healthy"), "prob"]-1))
        new.X.fail <- X[c(idx.apoptotic.fail, idx.healthy.fail),]
        new.X.pass <- X[c(idx.apoptotic.pass, idx.healthy.pass),]
        new.X <- rbind(new.X.fail, new.X.pass)
        new.Y.fail <- metrics[c(idx.apoptotic.fail, idx.healthy.fail),"label"]
        new.Y.pass <- metrics[c(idx.apoptotic.pass, idx.healthy.pass),"label"]
        new.Y <- c(new.Y.fail, new.Y.pass)
        new.Y <- factor(new.Y, levels = c("healthy","apoptotic"))
        weights <- rep(1, length(new.Y))
        weights[ 1:length(new.Y.fail) ] <- fail.weight

        # Jitter and randomize the data (reduces overfitting, and makes sampling with replacement a smaller problem)
        for (dim in 1:ncol(X)) { new.X[,dim] <- jitter(new.X[,dim], amount = sd(new.X[,dim])/5) }
        random.order <- sample(1:length(new.Y))

        # Train the model
        net <- glmnet::cv.glmnet(x = new.X[random.order ,], y = new.Y[random.order], family = "binomial", alpha = alpha, nfolds = nfolds, weights = weights[random.order])
        probs <- predict(net, newx = X, type = "response", s = "lambda.1se")
        prob.mat[,rep] <- probs[,1]
      }
      metrics$prob <- apply(prob.mat, 1, FUN="median")
      ROC <-  suppressMessages(pROC::roc(response = factor(metrics[ metrics$qc == "pass", "label"], levels = c("healthy","apoptotic")), predictor = as.numeric(metrics[ metrics$qc == "pass", "prob"])))
      cords <- pROC::coords(ROC,"best")
      metrics$prediction <- "healthy"
      metrics[ metrics$prob > cords$threshold,"prediction"] <- "apoptotic"
      metrics$prediction <- factor(metrics$prediction, levels = c("apoptotic","healthy"))
      stats[counter,1] <- cor.threshold
      stats[counter,2] <- cords$specificity
      stats[counter,3] <- table(metrics$prediction)[1]
      stats[counter,4] <- table(metrics$prediction, metrics$label)[3]
      stats[counter,5] <- table(metrics$prediction, metrics$label)[4]
      stats[counter,6] <- table(metrics$label)[2]
      counter <- counter + 1

    }
    # Select best threshold
    if (nrow(stats[ stats[,2] >= 0.99,]) > 0) {
      stats.subset <- stats[ stats[,2] >= 0.99,]
      if (nrow(stats.subset[ stats.subset[,5] / stats.subset[,6] <= 0.5,]) > 0) {
        stats.subset <- stats.subset[ stats.subset[,5] / stats.subset[,6] <= 0.5,]
        if (nrow(stats.subset[ stats.subset[,3] / stats.subset[,6] >= 2,]) > 0) {
          stats.subset <- stats.subset[ stats.subset[,3] / stats.subset[,6] >= 2,]
          cor.threshold <- stats.subset[ which.min(stats.subset[,4]),1]
        } else {
          cor.threshold <- stats.subset[ which.max(stats.subset[,3] / stats.subset[,6]),1]
        }
      } else {
        cor.threshold <- stats.subset[ which.max(stats.subset[,5] / stats.subset[,6]),1]
      }
    } else {
      cor.threshold <- stats[ which.max(stats[,2]),1]
    }
  }

  # Set probabilities
  metrics$prob = 1
  metrics[ metrics$label == "healthy", "prob"] <- 0

  #initialize attributes
  track_flag = list()
  counter = 0
  limited_epochs = 0
  balance = c()

  # Loop across epochs
  if (verbose) { message("Training models.") }
  for (eph in 1:epochs) {
    # Get all features
    X = svd$u[,1:npcs]
    rownames(X) <- colnames(data.scaled)

    #conditional looping
    if(counter < 5) {
      #Training first 5 epochs with determined threshold
      cor.coef <- c()
      for (dim in 1:npcs) { cor.coef <- c(cor.coef, pcaPP::cor.fk(X[,dim], metrics$label)) }
      X = X[,which(cor.coef^2 >= cor.threshold)]
    } else if ((flag < mean(c(track_flag[[eph-1]], track_flag[[eph-2]])))) {
      #compares how many got re-labeled in this epoch comparison to the mean of the previous 2 epochs
      if((balance[counter-1]/balance[counter-2]) < 1.5) {
        #compares the class balance of relabels
        cor.coef <- c()
        cor.threshold = cor.threshold + 0.0001
        for (dim in 1:npcs) { cor.coef <- c(cor.coef, pcaPP::cor.fk(X[,dim], metrics$label)) }
        X = X[,which(cor.coef^2 >= cor.threshold)]
      }
    } else if (limited_epochs < 1){
      #if conditions are not met executes once with the base threshold
      cor.coef <- c()
      cor.threshold = cor.threshold
      for (dim in 1:npcs) { cor.coef <- c(cor.coef, pcaPP::cor.fk(X[,dim], metrics$label)) }
      X = X[,which(cor.coef^2 >= cor.threshold)]
      limited_epochs = limited_epochs + 1
    } else {
      if (verbose) { message("Epoch limit reached.") }
      break
    }

    # Loop across replicates
    prob.mat <- as.data.frame(matrix(ncol = nrep, nrow = nrow(metrics)))
    for (rep in 1:nrep) {

      # Extract and combine the samples
      idx.apoptotic.fail <- sample(x=which(metrics$qc == "fail" & metrics$label == "apoptotic"), size=max(min.apoptotic, length(which(metrics$qc == "fail" & metrics$label == "apoptotic"))), replace=TRUE, prob=metrics[ which(metrics$qc == "fail" & metrics$label == "apoptotic"), "prob"])
      idx.apoptotic.pass <- sample(x=which(metrics$qc == "pass" & metrics$label == "apoptotic"), size=max(min.apoptotic, length(which(metrics$qc == "pass" & metrics$label == "apoptotic"))), replace=TRUE, prob=metrics[ which(metrics$qc == "pass" & metrics$label == "apoptotic"), "prob"])
      idx.healthy.fail <- sample(x=which(metrics$qc == "fail" & metrics$label == "healthy"), size=min(max.healthy, length(which(metrics$qc == "fail" & metrics$label == "healthy"))), replace=TRUE, prob=abs(metrics[ which(metrics$qc == "fail" & metrics$label == "healthy"), "prob"]-1))
      idx.healthy.pass <- sample(x=which(metrics$qc == "pass" & metrics$label == "healthy"), size=min(max.healthy, length(which(metrics$qc == "pass" & metrics$label == "healthy"))), replace=TRUE, prob=abs(metrics[ which(metrics$qc == "pass" & metrics$label == "healthy"), "prob"]-1))
      new.X.fail <- X[c(idx.apoptotic.fail, idx.healthy.fail),]
      new.X.pass <- X[c(idx.apoptotic.pass, idx.healthy.pass),]
      new.X <- rbind(new.X.fail, new.X.pass)
      new.Y.fail <- metrics[c(idx.apoptotic.fail, idx.healthy.fail),"label"]
      new.Y.pass <- metrics[c(idx.apoptotic.pass, idx.healthy.pass),"label"]
      new.Y <- c(new.Y.fail, new.Y.pass)
      new.Y <- factor(new.Y, levels = c("healthy","apoptotic"))
      weights <- rep(1, length(new.Y))
      weights[ 1:length(new.Y.fail) ] <- fail.weight

      # Jitter and randomize the data (reduces overfitting, and makes sampling with replacement a smaller problem)
      for (dim in 1:ncol(X)) { new.X[,dim] <- jitter(new.X[,dim], amount = sd(new.X[,dim])/5) }
      random.order <- sample(1:length(new.Y))

      # Train the model
      net <- glmnet::cv.glmnet(x = new.X[random.order ,], y = new.Y[random.order], family = "binomial", alpha = alpha, nfolds = nfolds, weights = weights[random.order])
      probs <- predict(net, newx = X, type = "response", s = "lambda.1se")
      prob.mat[,rep] <- probs[,1]
    }
    metrics$prob <- apply(prob.mat, 1, FUN="median")
    ROC <-  suppressMessages(pROC::roc(response = factor(metrics[ metrics$qc == "pass", "label"], levels = c("healthy","apoptotic")), predictor = as.numeric(metrics[ metrics$qc == "pass", "prob"])))
    cords <- pROC::coords(ROC,"best")
    metrics$prediction <- "healthy"
    metrics[ metrics$prob > cords$threshold,"prediction"] <- "apoptotic"
    metrics$prediction <- factor(metrics$prediction, levels = c("apoptotic","healthy"))
    metrics$label <- factor(metrics$prediction, levels = c("healthy","apoptotic"))

    #flag to keep a tract of how many got relabeled
    flag = ifelse(metrics$track==metrics$label,"Yes","No")
    flag = sum(flag == "No")
    track_flag[[eph]] = flag

    #for the other epoch update track labels
    metrics$track = metrics$label

    #check the class relabel balance
    res = do.call('rbind', lapply(label_list, FUN = function(x) { matrix(table(x[,2]), ncol=2) }))
    balance = res[,2] / res[,1]

    counter = counter + 1
  }

  # Return
  return(metrics)
}
