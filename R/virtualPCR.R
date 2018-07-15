#' Virtual PCR.
#'
#' \code{virtualPCR} queries a list of DNA sequences with virtual primers
#'   (either sequences or profile hidden Markov models) and returns only
#'   the sequences that contain regions of sufficient similarity based on
#'   log-odds alignment scoring.
#'
#' @param x a list of DNA sequences in \code{DNAbin} format.
#' @param up an object of class \code{DNAbin} or \code{PHMM} giving the
#'   forward primer with which to query the sequence list.
#' @param down an optional argument the same type as \code{up} giving the
#'   reverse primer with which to query the sequence list. If NULL only
#'   the forward primer is used.
#' @param rcdown logical indicating whether the reverse primer should be
#'   reverse-complemented prior to aligning with the input sequences. Set
#'   to TRUE only if \code{down} is not NULL, is of class \code{DNAbin}, and
#'   is the reverse complement of the target sequence (e.g. the sequence of
#'   a reverse primer as would be ordered from an oligo supplier).
#' @param trimprimers logical indicating whether the primer-binding sites
#'   should be removed from the sequences in the returned list.
#' @param minfsc numeric, giving the minimum specificity(log-odds score
#'   for the optimal alignment) between the forward primer and a sequence
#'   for that sequence to be retained.
#' @param minrsc numeric, the minimum specificity (log-odds score for
#'   the optimal alignment) between the reverse primer (if provided) and
#'   a sequence for that sequence to be retained.
#' @param minamplen,maxamplen integers giving the minimum and maximum
#'   acceptable amplicon lengths. Sequences are discarded if the number
#'   of base pairs between the primer-binding sites falls outside of these
#'   limits.
#' @param maxNs numeric giving the maximum acceptable proportion
#'   of the ambiguous residue "N" within the output sequences.
#'   Defaults to 0.02.
#' @param partialbind logical indicating whether partial primer matching is
#'   accepted. Defaults to TRUE.
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
#'   of sequences to be processed, due to the extra time required to initialize
#'   the cluster.
#' @param quiet logical indicating whether progress should be printed to
#'   the console.
#' @return a list of trimmed sequences, an object of class
#'   \code{DNAbin}.
#' @author Shaun Wilkinson
#' @examples
#'   ## trim whale sequences using a new set of inner primers
#'   inner_for <- "CGGTTGGGGTGACCTCGGAGTA"
#'   inner_rev <- "GCTGTTATCCCTAGGGTAA"
#'   whales_short <- virtualPCR(whales, up = inner_for, down = inner_rev,
#'                              trimprimers = TRUE)
################################################################################
virtualPCR <- function(x, up, down = NULL, rcdown = TRUE, trimprimers = FALSE,
                       minfsc = 50, minrsc = 50, minamplen = 50,
                       maxamplen = 2000, maxNs = 0.02, partialbind = TRUE, cores = 1,
                       quiet = FALSE){
  ischar <- mode(x) == "character"
  if(ischar) x <- char2dna(x)
  if(mode(up) == "character"){
    up <- char2dna(up)[[1]]
    if(!is.null(down)) down <- char2dna(down)[[1]]
  }
  nseq <- length(x)
  if(nseq == 0) stop("No sequences provided\n")
  if(!quiet) cat("Started with", nseq, "sequences\n")
  if(inherits(cores, "cluster")){
    para <- TRUE
    stopclustr <- FALSE
  }else if(cores == 1){
    para <- FALSE
    stopclustr <- FALSE
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(cores > 1){
      # if(cores > navailcores) stop("Number of cores is more than number available")
      if(!quiet) cat("Multithreading over", cores, "cores\n")
      cores <- parallel::makeCluster(cores)
      para <- TRUE
      stopclustr <- TRUE
    }else{
      para <- FALSE
      stopclustr <- FALSE
    }
  }
  # trim all nucleotides to left of forward primer bind site, including primer if specified
  forfun <- function(s, up, trimprimers, minfsc, partialbind, minamplen){
    vit <- aphid::Viterbi(up, s, type = "semiglobal", odds = TRUE)
    if(vit$score >= minfsc){
      pl <- length(vit$path)
      if(partialbind | vit$path[1] != 0){
        zeroonestarts <- match(0:1, vit$path)
        zeroonestarts <- zeroonestarts[!is.na(zeroonestarts)]
        overhang <- min(zeroonestarts) - 1
        if(overhang > 0) s <- s[-(1:overhang)]
        if(trimprimers){
          zeroonestarts2 <- match(0:1, rev(vit$path))
          zeroonestarts2 <- zeroonestarts2[!is.na(zeroonestarts2)]
          tokeep <- min(zeroonestarts2) - 1
          if(tokeep > 0) s <- rev(rev(s)[1:tokeep]) else return(NULL)
        }
        if(length(s) < minamplen) return(NULL)
        attr(s, "forscore") <- vit$score
        return(s)
      }
    }
    return(NULL)
  }
  tmpattr <- attributes(x)
  whichattr <- which(sapply(tmpattr, length) == length(x))
  if(!quiet) cat("Forward trimming sequences\n")
  x <- if(para){
    parallel::parLapply(cores, x, forfun, up, trimprimers, minfsc, partialbind, minamplen)
  }else{
    lapply(x, forfun, up, trimprimers, minfsc, partialbind, minamplen)
  }
  discards <- sapply(x, is.null)
  nseq <- sum(!discards)
  if(!quiet) cat("Retained", nseq, "sequences after forward trim\n")
  if(nseq > 0){
    forscores <- unlist(lapply(x, function(s) attr(s, "forscore")), use.names = FALSE)
    x <- x[!discards]
    x <- lapply(x, as.vector) # removes individual forscore attrs
    for(i in whichattr) tmpattr[[i]] <- tmpattr[[i]][!discards]
    attributes(x) <- tmpattr
    attr(x, "forscores") <- forscores
  }else{
    cat("None of the sequences met forward primer specificity criteria\n")
    x <- NULL # cant return yet since cluster is still open
  }
  # trim all nucleotides to right of reverse primer bind site, including primer if specified
  if(!is.null(down) & !is.null(x)){
    if(rcdown) down <- ape::complement(down)
    revfun <- function(s, down, trimprimers, minrsc, partialbind, minamplen, maxamplen){
      vit <- aphid::Viterbi(down, s, type = "semiglobal", odds = TRUE)
      if(vit$score >= minrsc){
        pl <- length(vit$path)
        if(partialbind | vit$path[pl] != 0){
          zeroonestarts <- match(0:1, rev(vit$path))
          zeroonestarts <- zeroonestarts[!is.na(zeroonestarts)]
          overhang <- min(zeroonestarts) - 1
          if(overhang > 0) s <- rev(rev(s)[-(1:overhang)])
          if(trimprimers){
            zeroonestarts2 <- match(0:1, vit$path)
            zeroonestarts2 <- zeroonestarts2[!is.na(zeroonestarts2)]
            tokeep <- min(zeroonestarts2) - 1
            if(tokeep > 0) s <- s[1:tokeep] else return(NULL)
          }
          if(length(s) >= minamplen & length(s) <= maxamplen) {
            attr(s, "revscore") <- vit$score
            return(s)
          }
        }
      }
      return(NULL)
    }
    tmpattr <- attributes(x)
    whichattr <- which(sapply(tmpattr, length) == length(x))
    if(!quiet) cat("Reverse trimming sequences\n")
    x <- if(para){
      parallel::parLapply(cores, x, revfun, down, trimprimers, minrsc, partialbind, minamplen, maxamplen)
    }else{
      lapply(x, revfun, down, trimprimers, minrsc, partialbind, minamplen, maxamplen)
    }
    discards <- sapply(x, is.null)
    nseq <- sum(!discards)
    if(nseq > 0){
      revscores <- unlist(lapply(x, function(s) attr(s, "revscore")), use.names = FALSE)
      if(!quiet) cat("Retained", nseq, "sequences after reverse trim\n")
      x <- x[!discards]
      x <- lapply(x, as.vector) # removes individual revscore attrs
      for(i in whichattr) tmpattr[[i]] <- tmpattr[[i]][!discards]
      attributes(x) <- tmpattr
      attr(x, "revscores") <- revscores
    }else{
      if(!quiet) cat("None of the sequences met reverse primer specificity criteria\n")
      x <- NULL
    }
  }
  if(para & stopclustr) parallel::stopCluster(cores)
  if(is.null(x)) return(x)
  if(!quiet) cat("Filtering ambiguous sequences\n")
  discards <- sapply(x, function(s) sum(s == 0xf0)/length(s)) > maxNs
  x <- subset.DNAbin(x, subset = !discards)
  if(!quiet) cat(length(x), "sequences retained after applying ambiguity filter\n")
  if(ischar) x <- dna2char(x)
  if(!quiet) cat("Done\n")
  return(x)
}
################################################################################
