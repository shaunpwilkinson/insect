#' Concatenate DNAbin objects while preserving attributes.
#'
#' This function joins two or more \code{DNAbin} objects, as well as the
#'   attributes whose lengths match the objects themselves
#'
#' @param ... \code{DNAbin} objects to be concatenated.
#' @return an object of class \code{DNAbin}.
#' @details TBA
#' @author Shaun Wilkinson
#' @seealso \code{\link{subset.DNAbin}}
#' @examples
#'   ## TBA
################################################################################
join <- function(...){
  dots <- list(...)
  nlsts <- length(dots)
  if(nlsts == 0) return(NULL)
  #spec <- defn <- lnge <- character(nseq)
  #res <- vector(mode = "list", length = nseq)
  spec <- attr(dots[[1]], "species")
  defn <- attr(dots[[1]], "definition")
  lnge <- attr(dots[[1]], "lineage")
  res <- if(is.list(dots[[1]])) dots[[1]] else structure(list(dots[[1]]), class = "DNAbin")
  # attrs <- c("species", "definition", "lineage")
  for(i in seq_along(dots)[-1]){
    spec <- c(spec, attr(dots[[i]], "species"))
    defn <- c(defn, attr(dots[[i]], "definition"))
    lnge <- c(lnge, attr(dots[[i]], "lineage"))
    if(!is.list(dots[[i]])) dots[[i]] <- structure(list(dots[[i]]), class = "DNAbin")
    res <- c(res, dots[[i]])
  }
  if(!(length(spec) == length(res) &
       length(defn) == length(res) &
       length(lnge) == length(res))){
    warning("invalid or missing attributes in at lest one of the DNAbin objects")
  }
  attr(res, "species") <- spec
  attr(res, "definition") <- defn
  attr(res, "lineage") <- lnge
  attr(res, "class") <- "DNAbin"
  return(res)
}
################################################################################
# trim DNA
# x is a DNAbin list or vector
# motif can be either a PHMM or a DNAbin vector, must not be rc'd (as ordered)
# resulting sequence includes the motif and anything to the right (if direction
# == "forward) or left (if direction == "reverse")

#' Trim ends of DNA sequences
#'
#' This function uses the Viterbi algorithm to semi-globally align a motif to
#'   a DNA sequence, and removes all nucleotides to the left or right of the
#'   motif.
#'
#' @param x an object of class \code{DNAbin}.
#' @param motif a \code{DNAbin} or \code{PHMM} object.
#' @param direction character string indicating
#'   the direction of the trim. Options are "forward" (trims everything to
#'   the right of the motif), "backward" (trims everything to the left of
#'   the motif) or "both" (retains the motif region only).
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
#' @param ... further arguments to be passed to \code{\link[aphid]{Viterbi}}
#'   (not including 'type').
#' @return an object of class \code{DNAbin}.
#' @details
#'   This functions finds the optimal semiglobal alignment (a.k.a. "glocal"
#'   alignment or global alignment with free end gaps) between a
#'   sequence \code{"x"} and a shorter sequence \code{"motif"}, returning
#'   the motif region of x along with the nucleotides to the left or right
#'   if \code{direction} is set to \code{"reverse"} or \code{"forward"},
#'   respectively.
#' @author Shaun Wilkinson
#' @seealso \code{\link{virtualPCR}}
#' @examples
#'   ## TBA
################################################################################
trim <- function(x, motif, direction = "both", cores = 1, ...){
  trim1 <- function(x, motif, direction){
    vit <- aphid::Viterbi(motif, x, type = "semiglobal", ... = ...)
    if(identical(direction, "forward")){
      ntotrim <- match(c(0, 1), rev(vit$path)) - 1
      ntotrim <- min(ntotrim[!is.na(ntotrim)])
      res <- x[1:(length(x) - ntotrim)]
    }else if(identical(direction, "reverse")){
      ntotrim <- match(c(0, 1), vit$path) - 1
      ntotrim <- min(ntotrim[!is.na(ntotrim)])
      res <- if(ntotrim > 0) x[-(1:ntotrim)] else x
    }else if(identical(direction, "both")){
      ntotrimf <- match(c(0, 1), rev(vit$path)) - 1
      ntotrimf <- min(ntotrimf[!is.na(ntotrimf)])
      last <- length(x) - ntotrimf
      begin <- match(c(0, 1), vit$path)
      begin <- min(begin[!is.na(begin)])
      res <- x[begin:last]
    }else stop("Invalid argument for 'direction'")
    return(res)
  }
  if(is.list(x)){
    if(inherits(cores, "cluster")){
      res <- parallel::parLapply(cores, x, trim1, motif, direction, ...)
    }else if(cores == 1){
      res <- lapply(x, trim1, motif, direction, ...)
    }else{
      #nseq <- length(x)
      navailcores <- parallel::detectCores()
      if(identical(cores, "autodetect")) cores <- navailcores - 1
      if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
      if(cores > navailcores) stop("Number of cores to use is more than number available")
      if(cores > 1){
        cl <- parallel::makeCluster(cores)
        res <- parallel::parLapply(cl, x, trim1, motif, direction, ...)
        parallel::stopCluster(cl)
      }else{
        res <- lapply(x, trim1, motif, direction, ...)
      }
    }
    return(res)
  }else{
    return(trim1(x, motif, direction, ...))
  }
  #res <- if(is.list(x)) lapply(x, trim1, motif, direction) else trim1(x, motif, direction)
  attributes(res) <- attributes(x)
  class(res) <- "DNAbin"
  return(res)
}
################################################################################
# trim <- function(x, motif, direction = "forward", ...){
#   if(!(identical(direction, "forward") | identical(direction, "reverse"))){
#     stop("Direction must be either 'forward' or 'reverse'")
#   }
#   fw <- direction == "forward"
#   trim1 <- function(x, motif, fw){
#     viti <- Viterbi(motif, x, type = "semiglobal", ... = ...)
#     pathi <- if(fw) rev(viti$path) else viti$path
#     ntotrim <- match(c(0, 1), pathi)
#     if(fw) ntotrim <- ntotrim - 1
#     ntotrim <- min(ntotrim[!is.na(ntotrim)])
#     res <- if(fw) x[1:(length(x) - ntotrim)] else x[-(1:(ntotrim - 1))]
#     return(res)
#   }
#   res <- if(is.list(x)) lapply(x, trim1, motif, fw) else trim1(x, motif, fw)
#   class(res) <- "DNAbin"
#   return(res)
# }
