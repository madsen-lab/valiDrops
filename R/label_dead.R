#' @title label_dead
#'
#' @description Identification of dead barcode using a simple heurestic to soft-label followed by label optimization using glmnet and adaptive resampling.
#'
#' @param counts A matrix containing counts for barcodes passing the barcode-rank threshold. See \link{rank_barcodes}
#' @param metrics A data frame containing quality metrics. See \link{quality_metrics} and \link{quality_filter}
#' @param qc.labels A named vector of barcodes quality filter labels ("pass" or "fail").
#' @param cor.threshold Threshold on correlation to use for feature selection. Set to NULL to automatically identify the threshold [default = NULL].
#' @param train A boolean (TRUE or FALSE) indicating whether or not to optimize labels [default = TRUE].
#' @param rep An integer indicating the number of times to optimize labels [default = 10].
#' @param n.min An integer indicating the number of times a predicted label must be consistent to be retained [default = 8].
#' @param n.relabel An integer indicating the number of dead cells to label if soft-labelling identifies no dead cells [default = 1].
#' @param feature.try An integer indicating the number of times to reduce the correlation threshold boundaries if no good solutions are found [default = 3].
#' @param verbose A boolean (TRUE or FALSE) indicating whether or not to be verbose [default = FALSE].
#' @param label.thrs Threshold to use for soft-labeling barcodes as either dead or live. If set to NULL, automatic threshold identification. [default = NULL].
#' @param label.frac Maximum fraction of barcodes to soft-label as dead [default = 0.1].
#' @param nfeats A number indicating the number of variable features to use for SVD [default = 2000].
#' @param alpha Alpha to use in glmnet [default = 0].
#' @param npcs A number indicating the number of singular values to calculate [default = 100].
#' @param weight A boolean (TRUE or FALSE) indicating whether or not to be weight barcodes during training [default = TRUE].
#' @param epochs The maximum number of epochs to train for [default = 20].
#' @param nfolds The number of folds to use in cross-validation [default = 5].
#' @param nrep The number of replicates to use during train time [default = 10].
#' @param fail.weight The weight to assign to barcodes that failed quality control [default = 0.2].
#' @param cor.min The minimum threshold to test for correlation-based feature selection [default = 0.0001].
#' @param cor.max The maximum threshold to test for correlation-based feature selection [default = 0.005].
#' @param cor.steps The number of steps between cor.min and cor.max to test for correlation-based feature selection [default = 50].
#' @param nrep.cor The number of replicates to perform when determining the correlation threshold [default = 10].
#' @param min.dead The minimum number of dead barcodes to use for training [default = 100].
#' @param max.live The maximum number of live barcodes to use for training [default = 500].
#' @param plot A boolean (TRUE or FALSE) indicating whether or not to plot the scoring threshold [default = TRUE].
#' @param bpparam The backend to use for processing. Set to MulticoreParam() or SnowParam() for parallel processing [default = SeriamParam()].
#'
#' @return A data frame object that contains the quality metrics of the data and results indicating the apoptotic cells.
#' @export
#' @import Matrix
#' @import pROC
#' @import glmnet
#' @import pcaPP
#' @import BiocParallel
#' @import inflection
#' 
label_dead <- function(counts, metrics, qc.labels, cor.threshold = NULL, train = TRUE, rep = 10, n.min = 8, n.relabel = 1, feature.try = 3, verbose = FALSE, label.thrs = NULL, label.frac = 0.1, nfeats = 2000, alpha = 0, npcs = 100, weight = TRUE, epochs = 20, nfolds = 5, nrep = 10, fail.weight = 0.2, cor.min = 0.0001, cor.max = 0.005, cor.steps = 50, nrep.cor = 10, min.dead = 100, max.live = 500, plot = TRUE, bpparam = SerialParam()) {
  # Soft label the dataset
  metrics$logUMIs <- scale(metrics$logUMIs, scale = FALSE)
  metrics$logFeatures <- scale(metrics$logFeatures, scale = FALSE)
  metrics$ribosomal_fraction <- asin(sqrt((metrics$ribosomal_fraction)))/(pi/2)
  metrics$coding_fraction <- asin(sqrt((metrics$coding_fraction)))/(pi/2)
  metrics$mitochondrial_fraction <- asin(sqrt(metrics$mitochondrial_fraction))/(pi/2)
  metrics$score <- metrics$logUMIs * -11.82 + metrics$logFeatures * 2.08 + metrics$ribosomal_fraction * 158.98 + metrics$logFeatures * metrics$coding_fraction * 18.87 + metrics$ribosomal_fraction * metrics$coding_fraction * -125.9
  metrics$qc <- qc.labels
  
  # Define a flag for a succesfull run
  flag <- "Succes" 
  
  # Determine the optimal threshold
  if (is.null(label.thrs)) {
	if (verbose) { message("Optimizing labeling threshold") }
	# Try default quantile cutoff
	max.quantile <- 0.1
	break.data <- as.data.frame(matrix(ncol=2, nrow = 200))
	cnt <- 1
	for (brk in seq(0.0001, max.quantile, 0.0001)) {
		break.data[cnt,1] <- brk
		break.data[cnt,2] <- quantile(metrics$score, brk)
		cnt <- cnt + 1
	}
	last.thrs <- quantile(metrics$score, inflection::uik(break.data$V1, break.data$V2))
	metrics$label <- "live"
	metrics[ metrics$score <= last.thrs,"label"] <- "dead"
	last.distribution <- table(metrics$label, metrics$qc)
	
	# Increase quantile cutoff and evaluate
	stop <- 0
	while (stop == 0) {
		# Label using a new cutoff
		max.quantile <- max.quantile + 0.1
		break.data <- as.data.frame(matrix(ncol=2, nrow = 200))
		cnt <- 1
		for (brk in seq(0.0001, max.quantile, 0.0001)) {
			break.data[cnt,1] <- brk
			break.data[cnt,2] <- quantile(metrics$score, brk)
			cnt <- cnt + 1
		}
		new.thrs <- quantile(metrics$score, inflection::uik(break.data$V1, break.data$V2))
		metrics$label <- "live"
		metrics[ metrics$score <= new.thrs,"label"] <- "dead"
		new.distribution <- table(metrics$label, metrics$qc)
		
		# Evaluate whether or not to stop
		if (min(new.distribution) > 0) {
			if (min(last.distribution) > 0) {
				if (last.distribution[1,2] == new.distribution[1,2]) {
					last.distribution <- new.distribution
					last.thrs <- new.thrs
				} else if (last.distribution[1,2] < new.distribution[1,2]) {
					label.thrs <- last.thrs
					stop <- 1
				}
			} else {
				label.thrs <- new.thrs
				stop <- 1
			}
		} else {
			if (max.quantile >= 0.85) {
				label.thrs <- new.thrs
				flag <- "Failed"
				stop <- 1
			} else {
				last.distribution <- new.distribution
				last.thrs <- new.thrs
			}
		}
	}
	
	# Label using final cutoff
	metrics$label <- "live"
	metrics[ metrics$score <= label.thrs,"label"] <- "dead"
	metrics$original.qc <- metrics$qc
  } else {
    # Label the barcodes
	metrics$label <- "live"
	metrics[ metrics$score <= label.thrs,"label"] <- "dead"
	metrics$original.qc <- metrics$qc
  }
  
  # Break if there are not enough observations or if there are too many
  if (nrow(metrics[ metrics$label == "dead",]) < 3) {
    if (verbose) { message("Soft-labeling identified less than 3 barcodes associated with dead cells. Returning results from soft-thresholding") }
	train <- FALSE
	flag <- "Failed"
  } else if (sum(metrics$qc == "pass" & metrics$label == "dead") == 0) {
    if (verbose) { message(paste("Soft-labeling labelled 0 barcode that pass QC as dead cells. Relabeling the top ", n.relabel, " least dead-like barcodes.", sep="")) }
	metrics[ metrics$barcode %in% metrics[ metrics$label == "dead","barcode"][order(metrics[ metrics$label == "dead","score"], decreasing = TRUE)][1:n.relabel],"qc"] <- "pass"
	flag <- "Caution"
  } else if (nrow(metrics[ metrics$label == "dead",]) / nrow(metrics) >= label.frac) {
    if (verbose) { message(paste("Soft-labeling identified more than ", label.frac*100, "% of barcodes as dead. Aborting.", sep="")) }
	train <- FALSE
	metrics$label <- "live"
	flag <- "Failed"
  } else {
	if (verbose) { message("Soft-labeling identified ", nrow(metrics[ metrics$label == "dead",])," dead barcodes.", sep="") }
  }
  
  # Plot if requested
  if (plot) {
	plot(sort(metrics$score), las = 1, ylab="Score", xlab="Rank", pch = 16)
	abline(h = label.thrs, col="red", lty=2)
  }
  
  # Define objects
  labels.df <- as.data.frame(matrix(ncol=rep, nrow=nrow(metrics)))
  
  # Start training loop
  if (train) {
	  if (verbose) { message("Running feature selection and dimensional reduction") }
	  # Normalize and log transform counts
	  counts <- counts[ , colnames(counts) %in% metrics$barcode]
	  counts <- counts[, match(metrics$barcode, colnames(counts))]
	  nonzero <-  counts[ Matrix::rowSums(counts) > 0,]
	  sf <- 10000 / Matrix::colSums(nonzero)
	  norm_transform <- Matrix::t(Matrix::t(nonzero) * sf)
	  norm_transform@x <- log1p(norm_transform@x)

	  # Find variable features
	  dev <- scry::devianceFeatureSelection(nonzero)
	  var.feats <- names(which(rank(-dev) <= nfeats))

	  # Scaling and SVD
	  means <- Matrix::colMeans(nonzero)
	  nr <- nrow(nonzero)
	  sds <- sqrt((Matrix::colMeans(nonzero*nonzero) - means^2) * (nr / (nr-1)))
	  sds[ sds == 0 ] <- 1
	  data.scaled <- Matrix::t((Matrix::t(nonzero) - means)/sds)

	  # Calculate SVDs
	  svd <- irlba::irlba(Matrix::t(data.scaled), nv=npcs, nu=npcs)
	  
	  # Optimize labels
	  X = svd$u[,1:npcs]
	  rownames(X) <- colnames(data.scaled)
	  Y <- factor(metrics$label, levels = c("live","dead"))
	  net <- glmnet::cv.glmnet(x = X, y = Y, alpha = 0, family = "binomial", nfolds = nfolds)
	  tmp <- metrics
	  tmp$prob <- predict(net, newx = X, type = "response", s = "lambda.1se")
	  ROC <-  suppressMessages(pROC::roc(response = factor(tmp[ tmp$qc == "pass" , "label"], levels = c("live","dead")), predictor = as.numeric(tmp[ tmp$qc == "pass", "prob"])))
	  cords <- pROC::coords(ROC)
	  threshold <- min(cords[ cords$specificity >= 0.99,1])
	  if (cords[ cords$threshold == threshold, 3] > 0) {
	  	metrics[ metrics$barcode %in% tmp[ tmp$qc == "pass" & tmp$label == "dead" & tmp$prob <= threshold,"barcode"],"label"] <- "live"
	  }
	  metrics$label <- factor(metrics$label, levels = c("live","dead"))
	  
	  # Pseudolist
pseudolist <- list()
for (outer in 1:rep) { pseudolist[[outer]] <- list(flag = flag, nrep = nrep, epochs = epochs, alpha = alpha, nfolds = nfolds, max.live = max.live, min.dead = min.dead, nrep.cor = nrep.cor, npcs = npcs, feature.try = feature.try, names = colnames(data.scaled), metrics = metrics, weight = weight, label.thrs = label.thrs, fail.weight = fail.weight, cor.max = cor.max, cor.min = cor.min, cor.steps = cor.steps, cor.threshold = cor.threshold) }

# Loop across the pseudolist
if (verbose) { message("Optimizing labels.") }
results <- bplapply(pseudolist, BPPARAM = bpparam, FUN = function(x) {
  # Extract information
  metrics <- x$metrics
  weight <- x$weight
  fail.weight <- x$fail.weight
  label.thrs <- x$label.thrs
  cor.max <- x$cor.max
  cor.min <- x$cor.min
  cor.steps <- x$cor.steps
  cor.threshold <- x$cor.threshold
  feature.try <- x$feature.try
  npcs <- x$npcs
  nrep.cor <- x$nrep.cor
  min.dead <- x$min.dead
  max.live <- x$max.live
  nfolds <- x$nfolds
  alpha <- x$alpha
  epochs <- x$epochs
  nrep <- x$nrep
  flag <- x$flag
  
  # Setup for training
  Y <- factor(metrics$label)
  fraction_0 <- rep(1 - sum(Y == "dead") / length(Y), sum(Y == "dead"))
  fraction_1 <- rep(1 - sum(Y == "live") / length(Y), sum(Y == "live"))
  weights <- numeric(length(Y))
  if (weight) {
    weights[Y == "dead"] <- fraction_0
    weights[Y == "live"] <- fraction_1
    score_weights <- log1p(abs(metrics$score - label.thrs))
    score_weights <- (score_weights - min(score_weights)) / (max(score_weights) - min(score_weights))
    qc_weight <- rep(1, length(Y))
    qc_weight[metrics$qc == "fail"] <- fail.weight
    weights <- weights * score_weights * qc_weight
  } else {
    weights[1:length(weights)] <- 1
  }
  
  # Determine correlation threshold
  if (is.null(cor.threshold)) {
    stop <- 0
    retry <- 0
    while (stop == 0) {
      stats <- as.data.frame(matrix(ncol=6, nrow = cor.steps))
      counter <- 1
      sequence <- 2^seq(log2(cor.min), log2(cor.max), length.out = cor.steps)
      for (cor.threshold in sequence) {
        # Set probabilities
        metrics$prob = 1
        metrics[ metrics$label == "live", "prob"] <- 0
        
        # Extract SVDs
        X = svd$u[,1:npcs]
        rownames(X) <- x$names
        
        # Subset SVDs
        cor.coef <- c()
        for (dim in 1:npcs) { cor.coef <- c(cor.coef, pcaPP::cor.fk(X[,dim], metrics$label)) }
        if (length(which(cor.coef^2 >= cor.threshold)) >= 2) {
          X = X[,which(cor.coef^2 >= cor.threshold)]
          
          # Loop across replicates
          prob.mat <- as.data.frame(matrix(ncol = nrep.cor, nrow = nrow(metrics)))
          for (rep in 1:nrep.cor) {
            # Extract and combine the samples
            idx.dead.fail <- as.numeric(sample(x=as.character(which(metrics$qc == "fail" & metrics$label == "dead")), size=max(min.dead, length(which(metrics$qc == "fail" & metrics$label == "dead"))), replace=TRUE, prob=metrics[ which(metrics$qc == "fail" & metrics$label == "dead"), "prob"]))
            idx.dead.pass <- as.numeric(sample(x=as.character(which(metrics$qc == "pass" & metrics$label == "dead")), size=max(min.dead, length(which(metrics$qc == "pass" & metrics$label == "dead"))), replace=TRUE, prob=metrics[ which(metrics$qc == "pass" & metrics$label == "dead"), "prob"]))
            idx.live.fail <- sample(x=which(metrics$qc == "fail" & metrics$label == "live"), size=min(max.live, length(which(metrics$qc == "fail" & metrics$label == "live"))), replace=TRUE, prob=abs(metrics[ which(metrics$qc == "fail" & metrics$label == "live"), "prob"]-1))
            idx.live.pass <- sample(x=which(metrics$qc == "pass" & metrics$label == "live"), size=min(max.live, length(which(metrics$qc == "pass" & metrics$label == "live"))), replace=TRUE, prob=abs(metrics[ which(metrics$qc == "pass" & metrics$label == "live"), "prob"]-1))
            new.X.fail <- X[c(idx.dead.fail, idx.live.fail),]
            new.X.pass <- X[c(idx.dead.pass, idx.live.pass),]
            new.X <- rbind(new.X.fail, new.X.pass)
            new.Y.fail <- metrics[c(idx.dead.fail, idx.live.fail),"label"]
            new.Y.pass <- metrics[c(idx.dead.pass, idx.live.pass),"label"]
            new.Y <- c(new.Y.fail, new.Y.pass)
            new.Y <- factor(new.Y, levels = c("live","dead"))
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
          ROC <-  suppressMessages(pROC::roc(response = factor(metrics[ metrics$qc == "pass", "label"], levels = c("live","dead")), predictor = as.numeric(metrics[ metrics$qc == "pass", "prob"])))
          cords <- pROC::coords(ROC,"best")
          if (nrow(cords) > 1) {
            metrics$prob <- apply(prob.mat, 1, FUN="mean")
            ROC <-  suppressMessages(pROC::roc(response = factor(metrics[ metrics$qc == "pass", "label"], levels = c("live","dead")), predictor = as.numeric(metrics[ metrics$qc == "pass", "prob"])))
            cords <- pROC::coords(ROC,"best")
            if (nrow(cords) > 1) {
              cords <- pROC::coords(ROC, x = seq(0,1,by=0.001))
              cords <- cords[ cords$specificity == max(cords$specificity),]
              cords <- cords[which.min(cords[,1]),,drop=FALSE]
            }
          }				  
          metrics$prediction <- "live"
          metrics[ metrics$prob > cords$threshold,"prediction"] <- "dead"
          metrics$prediction <- factor(metrics$prediction, levels = c("dead","live"))
          stats[counter,1] <- cor.threshold
          stats[counter,2] <- cords$specificity
          stats[counter,3] <- table(metrics$prediction)[1]
          stats[counter,4] <- table(metrics$prediction, metrics$label)[3]
          stats[counter,5] <- table(metrics$prediction, metrics$label)[4]
          stats[counter,6] <- table(metrics$label)[2]
          counter <- counter + 1
        }
      }
      
      # Evaluate if its necessary to retry
      stats <- stats[ !is.na(stats[,1]),]
      if (nrow(stats[ stats[,3] > 0 & stats[,2] >= 0.99,]) == 0) {
        retry <- retry + 1
        if (retry <= feature.try) {
          ratio <- cor.max / cor.min
          cor.max <- cor.min
          cor.min <- cor.min / ratio
        } else {
          stop <- 1
        }
      } else {
        stop <- 1
      }
    }
    
    # Select best threshold			
    if (nrow(stats[ stats[,3] > 0,]) > 0) {
      stats.subset <- stats[ stats[,3] > 0,]
      if (nrow(stats.subset[ stats.subset[,2] >= 0.99,]) > 0) {
        stats.subset <- stats.subset[ stats.subset[,2] >= 0.99,]
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
        cor.threshold <- stats.subset[ which.max(stats.subset[,2]),1]
      }
    } else {
      cor.threshold <- stats[ which.max(stats[,2]),1]
    }
  }
  
  # Set probabilities
  metrics$prob = 1
  metrics[ metrics$label == "live", "prob"] <- 0
  
  # Setup
  spec.old <- 0
  balance.old <- 0
  relabel.old <- nrow(metrics)
  trigger <- FALSE
  
  # Loop across epochs
  for (eph in 1:epochs) {
    # Get all features
    X = svd$u[,1:npcs]
    rownames(X) <- x$names
    
    # Feature selection
    cor.coef <- c()
    for (dim in 1:npcs) { cor.coef <- c(cor.coef, pcaPP::cor.fk(X[,dim], metrics$label)) }
    if (trigger) { cor.threshold = cor.threshold + cor.min }
    X = X[,which(cor.coef^2 >= cor.threshold)]
    
    # Break if there are no features
    if (dim(X)[2] == 0) {
      break
    }
    
    # Loop across replicates
    prob.mat <- as.data.frame(matrix(ncol = nrep, nrow = nrow(metrics)))
    for (rep in 1:nrep) {
      
      # Extract and combine the samples
      idx.dead.fail <- as.numeric(sample(x=as.character(which(metrics$qc == "fail" & metrics$label == "dead")), size=max(min.dead, length(which(metrics$qc == "fail" & metrics$label == "dead"))), replace=TRUE, prob=metrics[ which(metrics$qc == "fail" & metrics$label == "dead"), "prob"]))
      idx.dead.pass <- as.numeric(sample(x=as.character(which(metrics$qc == "pass" & metrics$label == "dead")), size=max(min.dead, length(which(metrics$qc == "pass" & metrics$label == "dead"))), replace=TRUE, prob=metrics[ which(metrics$qc == "pass" & metrics$label == "dead"), "prob"]))
      idx.live.fail <- sample(x=which(metrics$qc == "fail" & metrics$label == "live"), size=min(max.live, length(which(metrics$qc == "fail" & metrics$label == "live"))), replace=TRUE, prob=abs(metrics[ which(metrics$qc == "fail" & metrics$label == "live"), "prob"]-1))
      idx.live.pass <- sample(x=which(metrics$qc == "pass" & metrics$label == "live"), size=min(max.live, length(which(metrics$qc == "pass" & metrics$label == "live"))), replace=TRUE, prob=abs(metrics[ which(metrics$qc == "pass" & metrics$label == "live"), "prob"]-1))
      new.X.fail <- X[c(idx.dead.fail, idx.live.fail),]
      new.X.pass <- X[c(idx.dead.pass, idx.live.pass),]
      new.X <- rbind(new.X.fail, new.X.pass)
      new.Y.fail <- metrics[c(idx.dead.fail, idx.live.fail),"label"]
      new.Y.pass <- metrics[c(idx.dead.pass, idx.live.pass),"label"]
      new.Y <- c(new.Y.fail, new.Y.pass)
      new.Y <- factor(new.Y, levels = c("live","dead"))
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
    metrics$new_prob <- apply(prob.mat, 1, FUN="median")
    ROC <-  suppressMessages(pROC::roc(response = factor(metrics[ metrics$qc == "pass", "label"], levels = c("live","dead")), predictor = as.numeric(metrics[ metrics$qc == "pass", "new_prob"])))
    cords <- pROC::coords(ROC,"best")
    if (nrow(cords) > 1) {
      metrics$prob <- apply(prob.mat, 1, FUN="mean")
      ROC <-  suppressMessages(pROC::roc(response = factor(metrics[ metrics$qc == "pass", "label"], levels = c("live","dead")), predictor = as.numeric(metrics[ metrics$qc == "pass", "prob"])))
      cords <- pROC::coords(ROC,"best")
      if (nrow(cords) > 1) {
        cords <- pROC::coords(ROC, x = seq(0,1,by=0.001))
        cords <- cords[ cords$specificity == max(cords$specificity),]
        cords <- cords[which.min(cords[,1]),,drop=FALSE]
      }
    }		
    metrics$prediction <- "live"
    metrics[ metrics$new_prob > cords$threshold,"prediction"] <- "dead"
    metrics$prediction <- factor(metrics$prediction, levels = c("dead","live"))
    balance <- nrow(metrics[ metrics$prediction == "dead" & metrics$label == "live",]) / nrow(metrics[ metrics$prediction != metrics$label,])
    relabel <- nrow(metrics[ metrics$prediction != metrics$label,])
    
    # Evaluate how to procede with training
    if (relabel <= nrow(metrics) * 0.002 | relabel >= relabel.old * 2) {
      break
    } else if ((spec.old > cords$specificity | balance.old >= balance * 1.5) & trigger == FALSE) {
      trigger <- TRUE
    } else {
      balance.old <- balance
      spec.old <- cords$specificity 
      metrics$prob <- metrics$new_prob
      metrics$label <- factor(metrics$prediction, levels = c("live","dead"))
      if (min(table(metrics$qc, metrics$label)) == 0) {
        break
        flag <- "Failed"
      }
    }
    relabel.old <- relabel
  }
  
  # Return final labels
  return(list(labels = as.character(metrics$label), flag = as.character(flag)))
})
	  
	  # Extract results
	  labels.df <- as.data.frame(matrix(ncol=rep, nrow=nrow(metrics)))
	  flags <- c()
	  for (iter in 1:length(results)) {
		labels.df[,iter] <- results[[iter]]$labels
		flags <- c(flags, results[[iter]]$flag)
	}
	  
	  # Finalize labels
	  metrics$label <- apply(labels.df,1,FUN = function(x) { ifelse(sum(x == "live") >= n.min, "live", ifelse(sum(x == "dead") >= n.min, "dead", "uncertain"))})
	  metrics$qc <- metrics$original.qc
	  
	  # Finalize flag
	  if (any(flags == "Failed")) {
		flag <- "Failed"
	} else if (any(flags == "Caution")) {
		flag <- "Caution"
	} else {
		flag <- "Succes"
	}
	  # Warn if there are too many uncertain labels
	  uncertain.fraction <- nrow(metrics[ metrics$label == "uncertain" & metrics$qc == "pass",]) / nrow(metrics[ metrics$qc == "pass",])
	  if (uncertain.fraction >= 0.0125 & uncertain.fraction < 0.025) {
	    if (verbose) { message("There is a medium number of barcodes with uncertain classification. Interpret results with caution!") }
		flag <- "Caution"
	  } else if (uncertain.fraction >= 0.025) {
		if (verbose) { message("There is a high number of barcodes with uncertain classification. The model did not converge. Interpret with extreme caution.") }
		flag <- "Failed"
	  }
	}

  # Return
  return(list(labels = labels.df, metrics = metrics, flag = flag))
}
