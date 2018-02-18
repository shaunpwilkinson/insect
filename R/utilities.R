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
  DNA <- any(sapply(dots, class) == "DNAbin")
  if(nlsts == 0) return(NULL)
  if(nlsts == 1) return(dots[[1]])
  islist <- sapply(dots, is.list)
  for(i in which(!islist)){
    tmpattr <- attributes(dots[[i]])
    attributes(dots[[i]]) <- NULL
    dots[[i]] <- list(dots[[i]])
    attributes(dots[[i]]) <- tmpattr
  }
  dots <- lapply(dots, unclass)
  findattr <- function(x){
    names(attributes(x))[sapply(attributes(x), length) == length(x)]
  }
  attrlist <- lapply(dots, findattr)
  ual <- unique(unlist(attrlist, use.names = FALSE))
  validattrs <- ual[sapply(ual, function(e) all(sapply(attrlist, function(g) e %in% g)))]
  validattrs <- validattrs[!validattrs %in% c("names", "class")]
  res <- unlist(dots, recursive = FALSE, use.names = TRUE)
  for(i in validattrs){
    attr(res, i) <- unlist(lapply(dots, attr, i), use.names = FALSE)
  }
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
#'   of sequences to be processed, due to the extra time required to initialize
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
  tmpattr <- attributes(x)
  tmpattr <- tmpattr[names(tmpattr) != "quality"]
 #classx <- class(x)
  trim1 <- function(x, motif, direction){
    vit <- aphid::Viterbi(motif, x, type = "semiglobal", ... = ...)
    p <- attr(x, "quality") # trim quality scores too
    if(!is.null(p)) stopifnot(length(x) == length(p))
    if(identical(direction, "forward")){
      ntotrim <- match(c(0, 1), rev(vit$path)) - 1
      ntotrim <- min(ntotrim[!is.na(ntotrim)])
      res <- x[1:(length(x) - ntotrim)]
      if(!is.null(p)) attr(res, "quality") <- p[1:(length(x) - ntotrim)]
    }else if(identical(direction, "reverse")){
      ntotrim <- match(c(0, 1), vit$path) - 1
      ntotrim <- min(ntotrim[!is.na(ntotrim)])
      if(ntotrim > 0){
        res <- x[-(1:ntotrim)]
        if(!is.null(p)) attr(res, "quality") <- p[-(1:ntotrim)]
      }else{
        res <- x
      }
    }else if(identical(direction, "both")){
      ntotrimf <- match(c(0, 1), rev(vit$path)) - 1
      ntotrimf <- min(ntotrimf[!is.na(ntotrimf)])
      last <- length(x) - ntotrimf
      begin <- match(c(0, 1), vit$path)
      begin <- min(begin[!is.na(begin)])
      res <- x[begin:last]
      if(!is.null(p)) attr(res, "quality") <- p[begin:last]
    }else stop("Invalid argument for 'direction'")
    attr(res, "score") <- vit$score
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
      # if(cores > navailcores) stop("Number of cores to use is more than number available")
      if(cores > 1){
        cl <- parallel::makeCluster(cores)
        res <- parallel::parLapply(cl, x, trim1, motif, direction, ...)
        parallel::stopCluster(cl)
      }else{
        res <- lapply(x, trim1, motif, direction, ...)
      }
    }
  }else{
    res <- trim1(x, motif, direction, ...)
    tmpattr$quality <- attr(res, "quality")
    tmpattr$score <- attr(res, "score")
  }
  attributes(res) <- tmpattr
  return(res)
}
################################################################################

