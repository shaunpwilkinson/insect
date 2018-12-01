#' @noRd
################################################################################
.partition <- function(x, model = NULL, K = 2, allocation = "cluster",
                      refine = "Viterbi", nstart = 20, iterations = 50,
                      kmers = NULL, ksize = NULL, seqweights = "Henikoff",
                      cores = 1, quiet = FALSE, verbose = FALSE, ...){
  ### x is a DNAbin or AAbin object
  # model is a starting model to be trained on each side
  # assumes all seqs are unique
  # allocationis can be integer vector containing the set {1:K} same length as nseq
  DNA <- .isDNA(x)
  dots <- list(...)
  res <- list()
  nseq <- length(x)
  if(K > nseq) K <- nseq
  names(x) <- paste0("S", 1:nseq)
  pointers <- seq_along(x)
  # otud <- FALSE
  if(is.null(seqweights)){
    seqweights <- rep(1, nseq)
  }else if(identical(seqweights, "Henikoff")){
    seqweights <- aphid::weight(x, method = "Henikoff", k = 5)
  }else if(identical(seqweights, "Gerstein")){
    seqweights <- aphid::weight(x, method = "Gerstein", k = 5)
  }else if(length(seqweights) != nseq){
    stop("Invalid sequence weights passed to '.partition'")
  }
  if(nseq == 1) return(NULL)
  tmp <- integer(nseq)
  if(nseq == 2){
    group1 <- c(TRUE, FALSE)
    tmp[1] <- 1
    tmp[2] <- 2
  }else if(nseq == K){
    tmp <- 1:nseq
  }else if(identical(allocation, "cluster") | identical(allocation, "split")){
    if(!quiet & verbose) cat("Clustering sequences into", K, "groups\n")
    if(is.null(kmers) | is.null(ksize)){ # generally will not be true
      # if(!quiet) cat("Counting k-mers\n")
      # kmers <- kmer::kcount(x, k = 5)
      if(!quiet & verbose) cat("Counting", if(is.null(dots$k)) 4 else dots$k, "-mers\n")
      ksize = if(is.null(dots$k)) 4 else dots$k
      suppressMessages(kmers <- kmer::kcount(x, k = ksize))
    }
    if(!quiet & verbose) cat("Assigning sequences to groups ")
    if(identical(allocation, "split")){
      kcounts <- kmer::kcount(x, k = if(!is.null(dots$k)) dots$k else 5)
      infocols <- apply(kcounts, 2, function(v) length(unique(v))) == 2
      if(sum(infocols) > 30 & K == 2){
        if(!quiet & verbose) cat("using k-mer splitting method\n")
        hash1 <- function(v) paste(openssl::md5(as.raw(v == v[1])))
        hashes <- apply(kcounts[, infocols], 2, hash1)
        hfac <- factor(hashes)
        infocol <- match(levels(hfac)[which.max(tabulate(hfac))], hashes)
        tmp <- kcounts[, infocols][, infocol]
        tmp <- unname(tmp == tmp[1]) + 1L ## converts from logical to integer
        rm(kcounts)### free up memory
      }else{
        allocation <- "cluster"
      }
    }
    if(identical(allocation, "cluster")){
      if(!quiet & verbose) cat("using k-means algorithm\n")
      ##########kmers <- kmers/(sapply(x, length) - 4)## k - 1 = 3 ###############
      kmers <- .decodekc(kmers)
      slens <- apply(kmers, 1, sum) # corrected sequence lengths
      kmers <- kmers/(slens + ksize - 1)
      # if(length(x) > 1000){
      #   centroid <- apply(kmers, 2, mean)
      #   errs2 <- (t(t(kmers) - centroid))^2
      #   euclids <- unname(apply(errs2, 1, sum))
      #   dists <- slens/(2 * ksize) * euclids # linearized distances
      #   if(!quiet & verbose) cat("Mean distance from centroid", round(mean(dists), 2), "\n")
      #   if(mean(dists) > 0.1){
      #     if(!quiet & verbose) cat("Original size: ", length(x), "\n")
      #     otus <- .otu(x, k = 4, threshold = 0.99)
      #     if(max(otus) > 20){
      #       pointers <- .point(paste(otus))
      #       centrals <- grepl("\\*$", names(otus)) # central (logical)
      #       cord <- order(pointers[centrals]) # order of central seqs
      #       x <- x[centrals][cord]
      #       nseq <- length(x)
      #       kmers <- kmers[centrals, , drop = FALSE][cord, , drop = FALSE]
      #       #if(is.null(dim(kmers))) kmers <- matrix(kmers, nrow = 1)
      #       #otud <- TRUE
      #       if(!quiet & verbose) cat("Reduced size: ", length(x), "\n")
      #     }else{
      #       if(!quiet & verbose) cat("Too few OTUs to cluster\n")
      #     }
      #   }
      # }
      tmp <- tryCatch(kmeans(kmers, centers = K, nstart = nstart)$cluster,
                      error = function(er) sample(rep(1:K, nseq)[1:nseq]),
                      warning = function(wa) sample(rep(1:K, nseq)[1:nseq]))
      rm(kmers)### free up memory
    }
  }else{
    if(length(allocation) != nseq) stop("Invalid argument given for 'allocation'")
    tmp <- allocation
  }
  res$membership <- integer(0)
  res$init_membership <- tmp
  #res$success <- NA
  res$scores <- numeric(0)   # just so they're in the right order
  if(!quiet & verbose){
    # if(length(tmp) > 50){
    #   cat("Initial membership: ", paste0(head(tmp, 50), collapse = ""), ".......\n")
    # }else{
    #   cat("Initial membership: ", paste0(tmp, collapse = ""), "\n")
    # }
    cat("Group sizes: ")
    for(i in 1:K) cat(i, "=", sum(tmp == i), " ")
    cat("\n")
  }
  scores <- matrix(nrow = K, ncol = nseq)
  pnms <- paste0("model", 1:K) # profile HMM names
  rownames(scores) <- pnms
  colnames(scores) <- names(x) ### these have been changed to S1, S2, ... ?
  ### set up multithread
  if(inherits(cores, "cluster") | identical(cores, 1)){
    stopclustr <- FALSE
  }else{ # create cluster object
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
    # if(cores > navailcores) stop("Number of cores is more than the number available")
    # if(!quiet) cat("Multithreading with", cores, "cores\n")
    if(cores == 1){
      stopclustr <- FALSE
    }else{
      if(!quiet & verbose) cat("Initializing cluster with", cores, "cores\n")
      cores <- parallel::makeCluster(cores)
      stopclustr <- TRUE
    }
  }
  ### Parent model taken care of by 'learn' ?
  if(is.null(model)){
    if(!quiet & verbose) cat("Deriving parent model\n")
    ## this should only happen at top level for tree-learning
    model <- aphid::derivePHMM.list(x, refine = refine, seqweights = seqweights,
                          cores = cores, quiet = !(!quiet & verbose), inserts = "inherited",
                          alignment = TRUE, ... = ...)
  }else if(DNA & nrow(model$E) == 20L){ # transition from AA to DNA
    if(!quiet & verbose) cat("Deriving DNA model\n")
    model <- aphid::derivePHMM.list(x, refine = refine, seqweights = seqweights,
                          cores = cores, quiet = !(!quiet & verbose), inserts = "inherited",
                          alignment = TRUE, ... = ...)
  }
  res[["model0"]] <- model
  for(j in 1:K) res[[pnms[j]]] <- model
  seq_numbers <- integer(K)
  finetune <- FALSE
  md5s <- paste(openssl::md5(as.raw(tmp)))
  fscore <- function(s, model) aphid::forward(model, s, odds = FALSE)$score
  for(i in 1:iterations){
    if(!quiet & verbose) cat("Insect iteration", i, "\n")
    membership <- tmp
    for(j in 1:K) seq_numbers[j] <- sum(membership == j)
    mxcn <- max(seq_numbers)
    for(j in 1:K){
      # if(!quiet) cat("Calculating sequence weights given child model", j, "\n")
      # scale so that weights reflect smallest clade size
      seqweightsj <- seqweights[membership == j]
      seqweightsj <- seqweightsj/mean(seqweightsj) # scale so that mean = 1
      seqweightsj <- seqweightsj * mxcn/seq_numbers[j] ########
      if(!quiet & verbose) cat("Training child model", j, "\n")
      ins <- if(finetune) res[[pnms[j]]]$inserts else model$inserts
      if(is.null(ins)) ins <- TRUE # top level only
      suppressWarnings(
        modelj <- aphid::train(if(finetune) res[[pnms[j]]] else model,
                               x[membership == j], #model
                               method = refine, seqweights = seqweightsj,
                               # pseudocounts = 0.1*length(seqweightsj),#
                               inserts = "inherited",
                               alignment = sum(ins)/length(ins) < 0.5,
                               cores = cores, quiet = !(!quiet & verbose),
                               ... = ...)
      )
      modelj$weights <- NULL
      modelj$mask <- NULL
      modelj$map <- NULL
      modelj$reference <- NULL
      modelj$insertlengths <- NULL
      modelj$name <- NULL
      res[[pnms[j]]] <- modelj
      modelj$alignment <- NULL
      if(!quiet & verbose) cat("Calculating sequence probabilities given child model",
                               j, "\n")
      scores[j, ] <- if(inherits(cores, "cluster")){
        parallel::parSapply(cores, x, fscore, model = modelj)
      }else{
        sapply(x, fscore, model = res[[pnms[j]]])
      }
    }
    tmp <- apply(scores, 2, which.max)
    finetune <- sum(tmp == membership)/nseq > 0.95
    if(finetune & !quiet & verbose) cat("Fine-tune mode active\n")
    if(!quiet & verbose){
      # if(length(tmp) > 50){
      #   cat("Membership: ", paste0(head(tmp, 50), collapse = ""), ".......\n")
      # }else{
      #   cat("Membership: ", paste0(tmp, collapse = ""), "\n")
      # }
      cat("Group sizes ")
      for(i in 1:K) cat(i, ":", sum(tmp == i), " ")
      cat("\n")
    }
    tmpmd5 <- paste(openssl::md5(as.raw(tmp)))
    if(tmpmd5 %in% md5s){
      break
    }else if(length(unique(tmp)) < K){
      if(!quiet & verbose){
        cat("Unsuccessful split into", K, "groups\n")
        cat("Unable to split node\n")
      }
      if(stopclustr) parallel::stopCluster(cores)
      return(NULL)
    }
    md5s <- c(md5s, tmpmd5)
  }
  if(stopclustr) parallel::stopCluster(cores) # not used if called from 'learn'
  res$membership <- membership[pointers]
  res$scores <- scores[, pointers]
  class(res) <- "split"
  return(res)
}
################################################################################




#
# else if(identical(allocation, "cluster2")){
#   if(!quiet) cat("Clustering sequences into", K, "groups\n")
#   qds <- kdistance(x, k = 5, alpha = NULL)
#   tr <- as.dendrogram(hclust(qds, method = "average"))
#   if(K > 2){
#     my_height <- attr(tr, "height")
#     for(i in 1:(K - 2)){
#       which_ls <- which(sapply(tr, is.list))
#       if(length(which_ls) == 1){
#         tr <- c(tr[-which_ls], tr[[which_ls]])
#       }else{
#         branch_lengths <- sapply(tr[which_ls], function(e) my_height - attr(e, "height"))
#         shortest_branch <- which_ls[which.min(branch_lengths)]
#         tr <- c(tr[-shortest_branch], tr[[shortest_branch]])
#       }
#     }
#     attr(tr, "height") <- my_height # necessary?
#     class(tr) <- "dendrogram"
#   }
#   for(i in 1:K){
#     groupi <- names(x) %in% unlist(dendrapply(tr[[i]], function(e) attr(e, "label")))
#     tmp[groupi] <- i
#   }
# }

# else if(identical(allocation, "random")){
#   #randomly allocate to K groups
#   if(!quiet) cat("Randomly allocating sequences to", K, "groups\n")
#   rands <- runif(nseq)
#   quants <- quantile(rands, probs = seq(0, 1, by = 1/K))
#   quants[1] <- 0
#   for(i in 1:K) tmp[rands > quants[i] & rands <= quants[i + 1]] <- i
#   # rmed <- median(rands)
#   # group1 <- rands >= rmed
# }
