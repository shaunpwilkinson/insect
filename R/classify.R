#' Taxonomic classification for DNA barcode sequences.
#'
#' \code{"classify"} assigns a DNA barcode or other
#'   taxonomically-informative sequence to a node of a classification
#'   tree probabilistically, using a series of nested profile
#'   hidden Markov models.
#'
#' @param x a \code{"DNAbin"} object, either as a vector or a list of vectors.
#' @param tree an object of class \code{"insect"}.
#' @param threshold numeric value between 0 and 1 giving the minimum
#'   Akaike weight for the recursive classification procedure
#'   to continue toward the leaves of the tree.
#' @param decay logical indicating whether the decision to terminate the
#'   classification process should occur based on decaying Akaike weights
#'   (at each node, the Akaike weight of the selected model is multiplied by
#'   the Akaike weight of the selected model at the parent node) or whether
#'   each Akaike weight should be calculated independently of that of the
#'   parent node. Defaults to the latter (FALSE).
#' @param cores integer giving the number of CPUs to parallelize the operation
#'   over. Defaults to 1, and reverts to 1 if x is not a list.
#'   This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   for example by running \code{parallel::stopCluster(cores)}.
#'   The string 'autodetect' is also accepted, in which case the maximum
#'   number of cores to use is one less than the total number of cores available.
#' @return a character string giving the lineage of the input sequence
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @seealso \code{\link{learn}}
#' @examples
#'   ##TBA
################################################################################
classify <- function(x, tree, threshold = 0.9, decay = TRUE, cores = 1){
  classify1 <- function(x, tree, threshold = 0.9, decay = TRUE){
    path <- ""
    akw <- 1
    cakw <- 1
    res <- "" # lineage above the root node
    while(is.list(tree)){
      no_mods <- length(tree)
      sc <- numeric(no_mods) # scores (log probabilities)
      for(i in 1:no_mods){
        sc[i] <- aphid::forward.PHMM(attr(tree[[i]], "model"), x, odds = FALSE)$score
      }
      total_score <- aphid::logsum(sc)
      akwgts <- exp(sc - total_score)
      best_model <- which.max(akwgts)
      newakw <- akwgts[best_model]
      newcakw <- newakw * cakw
      threshold_met <- threshold <= if(decay) newcakw else newakw
      minscore_met <- sc[best_model] >= attr(tree[[best_model]], "minscore")
      if(!(threshold_met & minscore_met)) break
      path <- paste0(path, best_model)
      akw <- newakw
      cakw <- newcakw
      tree <- tree[[best_model]]
      res <- attr(tree, "lineage")
    }
    score <- if(decay) cakw else akw
    res <- paste(c(res, path, paste(score), paste(threshold_met),
                   paste(minscore_met)), collapse = "%")
    return(res)
  }
  unpack <- function(v){# v is a vector of %-delimited strings
    v <- strsplit(v, split = "%")
    out <- sapply(v, function(s) s[1])
    attr(out, "paths") <- sapply(v, function(s) s[2])
    attr(out, "scores") <- as.numeric(sapply(v, function(s) s[3]))
    attr(out, "threshold_met") <- as.logical(sapply(v, function(s) s[4]))
    attr(out, "minscore_met") <- as.logical(sapply(v, function(s) s[5]))
    return(out)
  }
  if(is.list(x)){
    if(inherits(cores, "cluster")){
      res <- parallel::parSapply(cores, x, classify1, tree, threshold, decay)
    }else if(cores == 1){
      res <- sapply(x, classify1, tree, threshold, decay)
    }else{
      navailcores <- parallel::detectCores()
      if(identical(cores, "autodetect")) cores <- navailcores - 1
      if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
      if(cores > navailcores) stop("Insufficient CPUs available")
      if(cores > 1){
        cl <- parallel::makeCluster(cores)
        res <- parallel::parSapply(cl, x, classify1, tree, threshold, decay)
        parallel::stopCluster(cl)
      }else{
        res <- sapply(x, classify1, tree, threshold, decay)
      }
    }
    res <- unpack(res)
  }else{
    res <- unpack(classify1(x, tree, threshold, decay))
  }
  return(res)
}
################################################################################
