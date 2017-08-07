#' Virtual \emph{in situ} hybridization.
#'
#' This function queries a list of DNA sequences with a virtual probe
#'   (either a sequence or a profile hidden Markov model) and returns only
#'   the sequences and regions that are of sufficient similarity based on
#'   log-odds alignment scoring.
#'
#' @param x a list of DNA sequences in \code{DNAbin} format.
#' @param probe a DNA sequence ("DNAbin" object) or profile hidden
#'   Markov model ("PHMM" object) to use as the probe.
#' @param minscore numeric; the minimum specificity (log-odds score
#'   for the optimal alignment) between the query sequence and the probe
#'   for the former to be retained.
#' @param minamplen,maxamplen integers giving the minimum and maximum
#'   acceptable amplicon lengths. Sequences are discarded if the number
#'   of base pairs between the primer-binding sites falls outside of these
#'   limits.
#' @param cores integer giving the number of CPUs to parallelize the operation
#'   over. Defaults to 1, and reverts to 1 if x is not a list.
#'   This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   for example by running \code{parallel::stopCluster(cores)}.
#'   The string 'autodetect' is also accepted, in which case the maximum
#'   number of cores to use is one less than the total number of cores available.
#'   Note that in this case there
#'   may be a tradeoff in terms of speed depending on the number and size
#'   of sequences to be aligned, due to the extra time required to initialize
#'   the cluster.
#' @param quiet logical indicating whether progress should be printed to
#'   the console.
#' @return a list of trimmed sequences, returned as an object of class
#'   \code{DNAbin}.
#' @details TBA
#' @author Shaun Wilkinson
#' @examples
#'   ## TBA
################################################################################
virtualFISH <- function(x, probe, minscore = 0, minamplen = 50,
                        maxamplen = 2000, cores = 1, quiet = FALSE){
  vf1 <- function(s, probe, minscore, minamplen, maxamplen){
    qual <- attr(s, "quality")
    if(!is.null(qual)) stopifnot(length(qual) == length(s))
    vit <- aphid::Viterbi(probe, s, odds = TRUE, type = "semiglobal")
    score <- vit$score
    path <- vit$path
    if(score < minscore) return(NULL)
    if(path[1] == 0L){
      # allow 1 indel
      if(match(1, path) <= 2) match1 <- 0 else return(NULL)
    }else{
      match1 <- match(1, path)
    }
    if(path[length(path)] == 0L){
      # allow 1 indel
      if(match(1, rev(path)) <= 2) match2 <- 0 else return(NULL)
    }else{
      match2 <- match(1, rev(path))
    }
    if(is.na(match1) | is.na(match2)) return(NULL)
    newlength <- length(s) - match1 - match2 + 2
    if(newlength < minamplen | newlength > maxamplen) return(NULL)
    if(match1 > 1){
      s <- s[match1:length(s)]
      if(!is.null(qual)) qual <- qual[match1:length(s)]
    }
    if(match2 > 1){
      s <- s[seq(1, length(s) - match2 + 1)]
      if(!is.null(qual)) qual <- qual[seq(1, length(s) - match2 + 1)]
    }
    attr(s, "score") <- score
    attr(s, "quality") <- qual
    return(s)
  }
  nseq <- length(x)
  tmpattr <- attributes(x)
  whichattr <- which(sapply(tmpattr, length) == length(x))
  if(inherits(cores, "cluster")){
    x <- parallel::parLapply(cores, x, vf1, probe, minscore, minamplen, maxamplen)
  }else if(cores == 1){
    x <- lapply(x, vf1, probe, minscore, minamplen, maxamplen)
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores'")
    if(cores > navailcores) stop("Number of cores is more than available")
    if(cores > 1){
      cl <- parallel::makeCluster(cores)
      x <- parallel::parLapply(cl, x, vf1, probe, minscore, minamplen, maxamplen)
      parallel::stopCluster(cl)
    }else{
      x <- lapply(x, vf1, probe, minscore, minamplen, maxamplen)
    }
  }
  discards <- sapply(x, is.null)
  nseq <- sum(!discards)
  if(nseq > 0){
    if(!quiet) cat("Retained", nseq, "sequences after trimming\n")
    scores <- unlist(lapply(x, function(s) attr(s, "score")), use.names = FALSE)
    x <- x[!discards]
    #x <- lapply(x, as.vector) # removes individual score attrs
    x <- lapply(x, function(s){
      attr(s, "score") <- NULL
      return(s)
    })
    for(i in whichattr) tmpattr[[i]] <- tmpattr[[i]][!discards]
    attributes(x) <- tmpattr
    attr(x, "scores") <- scores
  }else{
    cat("None of the sequences met primer specificity criteria. Returning NULL\n")
    x <- NULL
  }
  return(x)
}
################################################################################
