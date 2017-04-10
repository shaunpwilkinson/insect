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
#'   a reverse primer).
#' @param trimprimers logical indicating whether the primer-binding sites
#'   should be removed from the sequences in the returned list.
#' @param minfsc numeric, giving the minimum specificity(log-odds score
#'   for the optimal alignment) between the forward primer and a sequence
#'   for that sequence to be retained.
#' @param minrevsc numeric, the minimum specificity (log-odds score for
#'   the optimal alignment) between the reverse primer (if provided) and
#'   a sequence for that sequence to be retained.
#' @param minamplen,maxamplen integers giving the minimum and maximum
#'   acceptable amplicon lengths. Sequences are discarded if the number
#'   of base pairs between the primer-binding sites falls outside of these
#'   limits.
#' @param partialbind logical indicating whether partial primer matching is
#'   accepted.
#' @param reversecheck logical indicating whether the sequences should be
#'   aligned to the primers in both forward and reverse directions.
#'   Increases running time but necessary if some of the sequences in the
#'   input list are suspected to  be oriented in the 3' -> 5' direction
#'   (as is often the case with DNA sequences in GenBank).
#' @param reversethresh numeric, the minimum specificity (log-odds score
#'   for the optimal alignment) between the forward primer and a reverse-
#'   complemented sequence for that sequence to be retained. Only applicable
#'   if \code{reversecheck = TRUE}.
#' @param rm.duplicates logical indicating whether amplicon sequences that
#'   have identical nucleotide composition to one that is already in the
#'   dataset should be removed (TRUE) or retained (FALSE; default).
#' @param quiet logical indicating whether progress should be printed to
#'   the console.
#' @return a list of trimmed sequences, returned as an object of class
#'   \code{DNAbin}.
#' @details TBA
#' @author Shaun Wilkinson
#' @examples
#'   ## TBA
################################################################################
virtualPCR <- function(x, up, down = NULL, rcdown = TRUE, trimprimers = FALSE,
                       minfsc = 70, minrevsc = 70, minamplen = 50,
                       maxamplen = 2000, partialbind = TRUE, reversecheck = TRUE,
                       reversethresh = 80, rm.duplicates = FALSE, quiet = FALSE){

  nseq <- length(x)
  if(nseq == 0) stop("No sequences provided\n")
  if(!quiet) cat("Started with", nseq, "sequences\n")
  if(reversecheck){
    forscoresRC <- numeric(nseq)
    for(i in 1:nseq) {
      forscoresRC[i] <- aphid::Viterbi(up, ape::complement(x[i]),
                                       type = "semiglobal")$score
    }
    numreversed <- sum(forscoresRC > reversethresh)
    if(!quiet) cat("Detected", numreversed, "reversed sequences\n")
    if(numreversed > 0){
      if(!quiet & numreversed > 0) {
        cat("Applying reverse complement function to",
            numreversed, "sequences\n")
      }
      tmpattr <- attributes(x)
      x[forscoresRC > minfsc] <- lapply(x[forscoresRC > minfsc], ape::complement)
      attributes(x) <- tmpattr
      # test forward primer against seqs and select positive hits (eg >70)
      if(!quiet & numreversed > 0) {
        cat(numreversed, "sequences successfully reverse-complemented\n")
      }
    }
  }
  # trim all nucleotides to left of forward primer bind site, including primer if specified
  x1 <- list()
  included <- logical(nseq)
  counter <- 1
  for(i in 1:nseq){
    vit_i <- aphid::Viterbi(up, x[i], type = "semiglobal")
    if(vit_i$score > minfsc){
      pl <- length(vit_i$path)
      if(partialbind | vit_i$path[1] != 0){
        tmp <- x[[i]]
        zeroonestarts <- match(0:1, vit_i$path)
        zeroonestarts <- zeroonestarts[!is.na(zeroonestarts)]
        overhang <- min(zeroonestarts) - 1
        if(overhang > 0) tmp <- tmp[-(1:overhang)]
        if(trimprimers){
          zeroonestarts2 <- match(0:1, rev(vit_i$path))
          zeroonestarts2 <- zeroonestarts2[!is.na(zeroonestarts2)]
          tokeep <- min(zeroonestarts2) - 1
          if(tokeep == 0){
            tmp <- raw(0)
          }else{
            tmp <- rev(rev(tmp)[1:tokeep])
          }
        }
        if(length(tmp) >= minamplen){
          x1[[counter]] <- tmp
          counter <- counter + 1
          included[i] <- TRUE
        }
      }
    }
  }
  nseq <- length(x1)
  if(!quiet) cat("Retained", nseq, "sequences after forward trim\n")
  if(nseq == 0) stop("None of the sequences met forward primer specificity criteria\n")
  names(x1) <- names(x)[included]
  attr(x1, "species") <- attr(x, "species")[included]
  attr(x1, "definition") <- attr(x, "definition")[included]
  attr(x1, "lineage") <- attr(x, "lineage")[included]
  class(x1) <- "DNAbin"
  if(!is.null(down)){
    if(rcdown) down <- ape::complement(down)
    # trim all nucleotides to right of reverse primer bind site, including primer if specified
    x2 <- list()
    included <- logical(nseq)
    counter <- 1
    for(i in 1:nseq){
      vit_i <- aphid::Viterbi(down, x1[i], type = "semiglobal")
      if(vit_i$score > minrevsc){
        pl <- length(vit_i$path)
        if(partialbind | vit_i$path[pl] != 0){
          tmp <- x1[[i]]
          zeroonestarts <- match(0:1, rev(vit_i$path))
          zeroonestarts <- zeroonestarts[!is.na(zeroonestarts)]
          overhang <- min(zeroonestarts) - 1
          if(overhang > 0) tmp <- rev(rev(tmp)[-(1:overhang)])
          if(trimprimers){
            zeroonestarts2 <- match(0:1, vit_i$path)
            zeroonestarts2 <- zeroonestarts2[!is.na(zeroonestarts2)]
            tokeep <- min(zeroonestarts2) - 1
            if(tokeep == 0){
              tmp <- raw(0)
            }else{
              tmp <- tmp[1:tokeep]
            }
          }
          if(length(tmp) >= minamplen & length(tmp) <= maxamplen){
            x2[[counter]] <- tmp
            counter <- counter + 1
            included[i] <- TRUE
          }
        }
      }
    }
    nseq <- length(x2)
    if(!quiet) cat("Retained", nseq, "sequences after reverse trim\n")
    if(nseq == 0) stop("None of the sequences met reverse primer specificity criteria\n")
    names(x2) <- names(x1)[included]
    attr(x2, "species") <- attr(x1, "species")[included]
    attr(x2, "definition") <- attr(x1, "definition")[included]
    attr(x2, "lineage") <- attr(x1, "lineage")[included]
    class(x2) <- "DNAbin"
  }else{
    x2 <- x1
  }
  if(nseq == 0) stop("None of the sequences met primer specificity criteria\n")
  #### remove duplicate sequences
  if(rm.duplicates){
    x2 <- unique.DNAbin(x2)
    if(!quiet) cat("Retained", length(x2), "sequences after duplicate analysis\n")
  }
  if(!quiet) cat("Done\n")
  return(x2)
}
################################################################################
