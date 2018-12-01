#' Virtual \emph{in situ} hybridization.
#'
#' This function queries a list of DNA sequences with a virtual probe
#'   (either a sequence or a profile hidden Markov model) and returns only
#'   the sequences and regions that are of sufficient similarity based on
#'   log-odds alignment scoring.
#'
#' @param x a list of DNA sequences in \code{DNAbin} format.
#' @param probe a DNA sequence ("DNAbin" object) or profile hidden
#'   Markov model ("PHMM" object) to be used as the virtual hybridization probe.
#' @param minscore numeric; the minimum specificity (log-odds score
#'   for the optimal alignment) between the query sequence and the probe
#'   for the former to be retained in the output object.
#' @param minamplen,maxamplen integers giving the minimum and maximum
#'   acceptable amplicon lengths.
#' @param up,down optional objects of class \code{DNAbin}
#'   giving the forward and reverse primer sequences with which
#'   to query the sequence list following virtual probe hybridization.
#' @param rcdown logical indicating whether the reverse primer should be
#'   reverse-complemented prior to aligning with the input sequences.
#'   Should be set to TRUE if \code{down}
#'   is the reverse complement of the target sequence (e.g. the sequence of
#'   a reverse primer as would be ordered from an oligo supplier).
#' @param minfsc numeric, giving the minimum specificity(log-odds score
#'   for the optimal alignment) between the forward primer and a sequence
#'   for that sequence to be retained.
#' @param minrsc numeric, the minimum specificity (log-odds score for
#'   the optimal alignment) between the reverse primer (if provided) and
#'   a sequence for that sequence to be retained.
#' @param maxNs numeric giving the maximum acceptable proportion
#'   of the ambiguous residue "N" within the output sequences.
#'   Defaults to 0.02.
#' @param cores integer giving the number of processors for multithreading.
#'   Defaults to 1, and reverts to 1 if \code{x} is not a list.
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
#' @return a list of trimmed sequences, returned as an object of class
#'   \code{DNAbin}.
#' @details This function is generally used when filtering/trimming a
#'   local sequence database,
#'   to mop up any high-scoring sequences with partial/missing primer
#'   bind sites that were not included in the output of the
#'   \code{\link{virtualPCR}}.
#'   For example, this includes sequences that were generated using the same
#'   primer set as used in the virtual PCR, and whose primer binding sites
#'   were trimmed prior to deposition in the sequence database.
#'   Unlike the virtualPCR function, there is
#'   no option to retain the primer-bind sites in the returned sequences.
#' @author Shaun Wilkinson
#' @seealso \code{\link{virtualPCR}}
#' @examples
#'   ## ensure whale sequences are globally alignable
#'   data(whales)
#'   model <- aphid::derivePHMM(whales)
#'   z <- virtualFISH(whales, probe = model)
################################################################################
virtualFISH <- function(x, probe, minscore = 100, minamplen = 50,
                        maxamplen = 500, up = NULL, down = NULL, rcdown = TRUE,
                        minfsc = 60, minrsc = 60, maxNs = 0.02, cores = 1,
                        quiet = FALSE){
  ischar <- mode(x) == "character"
  if(ischar) x <- char2dna(x)
  if(!is.null(up)){
    if(!inherits(up, "DNAbin")){
      if(mode(up) == "character"){
        if(nchar(up[1]) == 1) up <- paste0(up, collapse = "")
        up <- char2dna(up)
      }else{
        if(!inherits(up, "PHMM")) stop("Invalid primer(s)\n")
      }
    }
  }
  if(!is.null(down)){
    if(inherits(down, "DNAbin")){
      if(rcdown) down <- ape::complement(down)
    }else if(mode(down) == "character"){
      if(nchar(down[1]) == 1) down <- paste0(down, collapse = "")
      down <- char2dna(rc(down))
    }else{
      if(!inherits(down, "PHMM")) stop("Invalid primer(s)\n")
    }
  }
  dd1 <- function(s, probe, minscore, minamplen, maxamplen, up, down, minfsc, minrsc){
    s <- s[!s %in% as.raw(c(2, 4))]
    if(!is.null(up)) fpl <- length(up)
    if(!is.null(down)) rpl <- length(down)
    vit <- aphid::Viterbi(probe, s, odds = TRUE, type = "semiglobal")
    if(vit$score < minscore) return(NULL)
    path <- vit$path
    if(path[1] == 0 | path[length(path)] == 0) return(NULL)
    match1 <- match(1, path)
    match2 <- match(1, rev(path))
    if(is.na(match1) | is.na(match2)) return(NULL)
    newlength <- length(s) - match1 - match2 + 2
    if(newlength < minamplen | newlength > maxamplen) return(NULL)
    score <- function(y, z) aphid::Viterbi(y, z, type = "semiglobal", odds = TRUE)$score
    if(is.null(up)) fpl <- 0
    if(is.null(down)) rpl <- 0
    if(match1 > fpl & match2 > rpl){
      if(!is.null(down)){
        rsc <- score(down, s[seq(length(s) - match2 + 2, length(s))])
        if(rsc < minrsc) return(NULL)
      }
      s <- s[seq(1, length(s) - match2 + 1)]
      if(!is.null(up)){
        fsc <- score(up, s[seq(1, match1 - 1)])
        if(fsc < minfsc) return(NULL)
      }
      s <- s[match1:length(s)]
    }else if(match1 <= fpl & match2 <= rpl){ # both sides amputated
      s <- s[match1:length(s)]
      s <- s[seq(1, length(s) - match2 + 1)]
    }else if(match1 <= fpl){ # only left hand amputated
      s <- s[match1:length(s)]
      if(!is.null(down)){
        rsc <- score(down, s[seq(length(s) - match2 + 2, length(s))])
        if(rsc < minrsc) return(NULL)
      }
      s <- s[seq(1, length(s) - match2 + 1)]
    }else if(match2 <= rpl){# only right hand amputated
      s <- s[seq(1, length(s) - match2 + 1)]
      if(!is.null(up)){
        fsc <- score(up, s[seq(1, match1 - 1)])
        if(fsc < minfsc) return(NULL)
      }
      s <- s[match1:length(s)]
    }
    attr(s, "score") <- vit$score
    return(s)
  }
  nseq <- length(x)
  tmpattr <- attributes(x)
  whichattr <- which(sapply(tmpattr, length) == length(x))
  if(inherits(cores, "cluster")){
    x <- parallel::parLapply(cores, x, dd1, probe, minscore, minamplen, maxamplen,
                             up, down, minfsc, minrsc)
  }else if(cores == 1){
    x <- lapply(x, dd1, probe, minscore, minamplen, maxamplen, up, down, minfsc, minrsc)
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores'")
    # if(cores > navailcores) stop("Number of cores is more than available")
    if(cores > 1){
      if(!quiet) cat("Multithreading with", cores, "cores\n")
      cl <- parallel::makeCluster(cores)
      x <- parallel::parLapply(cl, x, dd1, probe, minscore, minamplen, maxamplen,
                               up, down, minfsc, minrsc)
      parallel::stopCluster(cl)
    }else{
      x <- lapply(x, dd1, probe, minscore, minamplen, maxamplen, up, down, minfsc, minrsc)
    }
  }
  discards <- sapply(x, is.null)
  nseq <- sum(!discards)
  if(nseq > 0){
    if(!quiet) cat("Retained", nseq, "sequences after trimming\n")
    scores <- unlist(lapply(x, function(s) attr(s, "score")), use.names = FALSE)
    x <- x[!discards]
    x <- lapply(x, as.vector) # removes individual score attrs
    for(i in whichattr) tmpattr[[i]] <- tmpattr[[i]][!discards]
    attributes(x) <- tmpattr
    attr(x, "scores") <- scores
  }else{
    if(!quiet) cat("None of the sequences met primer/probe specificity criteria. Returning NULL\n")
    x <- NULL
  }
  if(!quiet) cat("Filtering ambiguous sequences\n")
  discards <- sapply(x, function(s) sum(s == 0xf0)/length(s)) > maxNs
  x <- subset.DNAbin(x, subset = !discards)
  if(!quiet) cat(length(x), "sequences retained after applying ambiguity filter\n")
  if(ischar) x <- dna2char(x)
  if(!quiet) cat("Done\n")
  return(x)
}
################################################################################
