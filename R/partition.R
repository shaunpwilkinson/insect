################################################################################
#' Split a group of sequences into subsets with maximum differentiation.
#'
#' Finds the best split of a set of DNA sequences and provides a profile
#'   hidden Markov model for each subset along with the logged full
#'   probability of each sequence given its respective model (calculated
#'   using the forward algorithm).
#'
#' @param x a \code{"DNAbin"} object containing multiple sequences in a list.
#' @param model an optional profile hidden Markov model as the hypothetical
#'   data generating mechanism for the entire sequence set. This is used to
#'   provide the starting parameters to train the top level models.
#' @param needs_training logical indicating whether the profile hidden Markov model
#'   provided above (if applicable) should be trained before attempting to
#'   split the sequence set, using the sequences supplied in \code{x}.
#'   Only necessary if \code{"model"} is not NULL.
#' @param K integer. The number of groups to split the sequences into.
#' @param allocation integer vector giving the initial group membership of the
#'   sequences. The vector length should be the same as that of the sequence set
#'   \code{x}, and all elements should take values between 1 and K.
#'   Alternatively, the character string "cluster" is accepted. In this case
#'   sequences are initially allocated to K groups using K-means clustering.
#'   This is the default option.
#' @param refine the method used to train/refine the models. Valid methods
#'   are "Viterbi" (Viterbi training; default) and "BaumWelch" (a modified
#'   version of the E-M algorithm).
#' @param iterations integer giving the maximum number of training-classification
#'   iterations to be used in the splitting process (see details section).
#'   Note that this is not necessarily the same as the number of Viterbi training
#'   or Baum Welch iterations to be used in model training, which can be set
#'   using the argument \code{"maxiter"} (eventually passed to
#'   \code{\link[aphid]{train}} via the dots argument "...").
#' @param seqweights either NULL (default; all sequences are given an equal
#'   weight of 1), a numeric vector the same length as \code{x} representing
#'   the sequence weights used to derive the model, or a character string giving
#'   the method to derive the weights from the sequences. Currently only the
#'   \code{"Gerstein"} method is supported (the default). For this method, a
#'   tree is first created by k-mer counting (see
#'   \code{\link[phylogram]{kdistance}}),
#'   and sequence weights are derived from the tree using the 'bottom up'
#'   algorithm of Gerstein et al. (1994). The sum of these weights are equal
#'   to the number of sequences in the alignment (so that mean(seqweights) = 1;
#'   Note this does not need to be the case if providing weights as a numeric
#'   vector).
#' @param maxcores integer giving the maximum number of CPUs to be used
#'   when training the models (only applicable if
#'   \code{refine = 'Viterbi'}). Note that the number of cores used may
#'   be less than the number given if the sequence training set is small
#'   or there are fewer cores available.
#' @param quiet logical indicating whether the progress should be printed to
#'   the console.
#' @param ... further arguments to be passed to \code{"aphid::train"} (not
#'   including 'inserts' or 'cores').
#' @return an object of class \code{split}. #############################TODO
#' @details TBA.
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Gerstein M, Sonnhammer ELL, Chothia C (1994) Volume changes in protein evolution.
#'   \emph{Journal of Molecular Biology}, \strong{236}, 1067-1078.
#' @seealso \code{\link{fork}} (parent function)
#' @examples
#'   ## TBA
################################################################################
partition <- function(x, model = NULL, needs_training = FALSE, K = 2,
                      allocation = "cluster", refine = "Viterbi",
                      iterations = 50, seqweights = "Gerstein",
                      maxcores = 1, quiet = FALSE, ...){
  ### x is a DNAbin object
  # model is a starting model to be trained on each side
  # assumes all seqs are unique
  # allocationis can be integer vector containing the set {1:K} same length as nseq
  res <- list()
  nseq <- length(x)
  if(K > nseq) K <- nseq
  names(x) <- paste0("S", 1:nseq) # just ensures all names are unique
  if(is.null(seqweights)){
    seqweights <- rep(1, nseq)
  }else if(identical(seqweights, "Gerstein")){
    seqweights <- aphid::weight(x, method = "Gerstein")
  }else if(length(seqweights) != nseq){
    stop("Invalid sequence weights passed to 'partition'")
  }
  if(nseq == 1) return(NULL)
  tmp <- integer(nseq)
  if(nseq == 2){
    group1 <- c(TRUE, FALSE)
    tmp[1] <- 1
    tmp[2] <- 2
  }else if(nseq == K){
    tmp <- 1:nseq
  }else if(identical(allocation, "cluster")){
    if(!quiet) cat("Clustering sequences into", K, "groups\n")
    x <- lapply(x, function(y) y[y != as.raw(2)])
    freqs <- phylogram::mbed(x, k = 5)
    #tmp <- kmeans(freqs, centers = K)$cluster
    tmp <- tryCatch(kmeans(freqs, centers = K)$cluster,
                    error = function(er) sample(rep(1:K, nseq)[1:nseq]),
                    warning = function(wa) sample(rep(1:K, nseq)[1:nseq]))
  }else{
    if(length(allocation) != nseq) stop("Invalid argument given for 'allocation'")
    tmp <- allocation
  }
  res$membership <- integer(0)
  #res$success <- NA
  res$scores <- numeric(0)   # just so they're in the right order
  if(!quiet) {
    if(length(tmp) > 50){
      cat("Initial membership: ", paste0(head(tmp, 50), collapse = ""), ".......\n")
    }else{
      cat("Initial membership: ", paste0(tmp, collapse = ""), "\n")
    }
    cat("Group sizes ")
    for(i in 1:K) cat(i, ":", sum(tmp == i), " ")
    cat("\n")
  }
  scores <- matrix(nrow = K, ncol = nseq)
  pnms <- paste0("phmm", 1:K) # profile HMM names
  rownames(scores) <- pnms
  colnames(scores) <- names(x)
  if(is.null(model)){
    if(!quiet) cat("Deriving parent model\n")
    nseeds <- ceiling(log(nseq, 2)^2)
    seeds <- sample(1:nseq, size = nseeds)
    #seedsTF <- 1:nseq %in% seeds
    model <- aphid::derivePHMM.list(x, refine = "none", seeds = seeds,
                                    seqweights = seqweights)
    ## just a progressive multiple alignment of seed seqs only
    ## now train the model
    if(!quiet) cat("Training parent model\n")
    optncores <- .optcores(maxcores, nseq)
    para <- optncores > 1
    if(para & !quiet) cat("Multithreading model training over", optncores, "cores\n")
    cores <- if(para) parallel::makeCluster(optncores) else 1
    model <- aphid::train(model, x, method = refine, seqweights = seqweights,
                        cores = cores, quiet = quiet, ... = ...)
    if(para) parallel::stopCluster(cores)
    # model <- aphid::derivePHMM.list(x, refine = refine,
    #                                 seqweights = seqweights, ... = ...)
  }else if(needs_training){
    if(!quiet) cat("Training parent model\n")  # model can change size here
    optncores <- .optcores(maxcores, nseq)
    para <- optncores > 1
    if(para & !quiet) cat("Multithreading model training over", optncores, "cores\n")
    cores <- if(para) parallel::makeCluster(optncores) else 1
    model <- aphid::train(model, x, method = refine, seqweights = seqweights,
                          cores = cores, quiet = quiet, ... = ...)
    if(para) parallel::stopCluster(cores)
  }
  res[["phmm0"]] <- model
  for(j in 1:K) res[[pnms[j]]] <- model
  seq_numbers <- integer(K)
  #seq_proportions <- rep(0, K)
  finetune <- FALSE
  md5s <- paste(openssl::md5(as.raw(tmp)))
  #md5s <- paste(openssl::md5(paste0(tmp, collapse = "")))
  for(i in 1:iterations){
    if(!quiet) cat("HPHMM iteration", i, "\n")
    membership <- tmp
    for(j in 1:K) seq_numbers[j] <- sum(membership == j)
    mcn <- min(seq_numbers)
    #if(mcn < 4) mcn <- 4
    for(j in 1:K){
      # if(!quiet) cat("Calculating sequence weights for child model", j, "\n")
      seqweightsj <- seqweights[membership == j]
      seqweightsj <- seqweightsj/mean(seqweightsj) # scale so that mean = 1
      seqweightsj <- seqweightsj * mcn/seq_numbers[j] # scale so that weights reflect smallest clade size
      #seqweights <- weight(x[membership == j]) * mcn/seq_numbers[j]
      if(!quiet) cat("Training child model", j, "\n")
      # res[[pnms[j]]] <- aphid::train(res[[pnms[j]]], x[membership == j], #model
      #                         method = refine, seqweights = seqweightsj,
      #                         inserts = if(refine == "Viterbi") "inherited" else "map",
      #                         ... = ...)
      optncores <- .optcores(maxcores, nseq = seq_numbers[j])
      para <- optncores > 1
      if(para & !quiet) cat("Multithreading model training over", optncores, "cores\n")
      cores <- if(para) parallel::makeCluster(optncores) else 1
      res[[pnms[j]]] <- aphid::train(if(finetune) res[[pnms[j]]] else model,
                                     x[membership == j], #model
                                     method = refine, seqweights = seqweightsj,
                                     inserts = if(refine == "Viterbi") "inherited" else "map",
                                     cores = cores, quiet = quiet, ... = ...)
      if(para) parallel::stopCluster(cores)
      # res[[pnms[j]]] <- aphid::train(model, x[membership == j], #model
      #                         method = refine, seqweights = seqweightsj,
      #                         inserts = if(refine == "Viterbi") "inherited" else "map",
      #                         ... = ...)
      if(!quiet) cat("Calculating sequence probabilities for child model", j, "\n")
      for(l in 1:nseq){
        scores[j, l] <- aphid::forward(res[[pnms[j]]], x[[l]], odds = FALSE)$score
      }
    }
    tmp <- apply(scores, 2, which.max)
    finetune <- sum(tmp == membership)/nseq > 0.95
    if(finetune & !quiet) cat("Fine-tune mode active\n")
    if(!quiet){
      if(length(tmp) > 50){
        cat("Membership: ", paste0(head(tmp, 50), collapse = ""), ".......\n")
      }else{
        cat("Membership: ", paste0(tmp, collapse = ""), "\n")
      }
      #cat("Membership: ", paste0(head(tmp, 50), collapse = ""), ".......", "\n")
      cat("Group sizes ")
      for(i in 1:K) cat(i, ":", sum(tmp == i), " ")
      cat("\n")
    }
    tmpmd5 <- paste(openssl::md5(as.raw(tmp)))
    # cat("hashes:", md5s, "\n")
    # cat("New hash:", tmpmd5, "\n")
    # cat("In md5s?", tmpmd5 %in% md5s, "\n")
    if(tmpmd5 %in% md5s){
      #cat("matches iteration", match(TRUE, sapply(md5s, identical, tmpmd5)), "\n")
      break
    }else if(length(unique(tmp)) < K){
      if(!quiet){
        cat("Unsuccessful split into", K, "groups\n")
        cat("Unable to split clade\n")
      }
      return(NULL)
    }
    md5s <- c(md5s, tmpmd5)
  }
  res$membership <- membership
  res$scores <- scores
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
