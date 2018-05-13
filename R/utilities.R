#' Concatenate DNAbin objects while preserving attributes.
#'
#' This function joins two or more \code{DNAbin} objects, retaining any
#'   attributes whose lengths match those of the input objects (e.g. "species",
#'   "lineage" and/or "taxID" attributes).
#'
#' @param ... \code{DNAbin} objects to be concatenated.
#' @return an object of class \code{DNAbin}.
#' @author Shaun Wilkinson
#' @seealso \code{\link{subset.DNAbin}}.
#' @examples
#' data(whales)
#' join(whales, whales)
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
  class(res) <- "DNAbin"
  return(res)
}
################################################################################
#' Shave ends from DNA and amino acid sequences
#'
#' This function uses the Viterbi algorithm to semi-globally align a motif to
#'   a DNA or AA sequence, and removes all nucleotides to the left and/or right of the
#'   motif.
#'
#' @param x an object of class \code{DNAbin} or \code{AAbin}.
#' @param motif a \code{DNAbin}, \code{AAbin} or \code{PHMM} object.
#' @param direction character string indicating
#'   the direction of the shave. Options are "forward" (shaves everything to
#'   the right of the motif), "backward" (shaves everything to the left of
#'   the motif) or "both" (retains the motif region only).
#' @param cores integer giving the number of CPUs to parallelize the operation
#'   over. Defaults to 1, and reverts to 1 if \code{x} is not a list.
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
#' @return an object of class \code{DNAbin} or \code{AAbin}
#'   (depending on the input object).
#' @details
#'   This functions finds the optimal semiglobal alignment (a.k.a. "glocal"
#'   alignment or global alignment with free end gaps) between a
#'   sequence \code{"x"} and a shorter sequence \code{"motif"}, returning
#'   the motif region of x along with the nucleotides to the left or right
#'   if \code{direction} is set to \code{"reverse"} or \code{"forward"},
#'   respectively.
#' @author Shaun Wilkinson
#' @seealso \code{\link{virtualPCR}}.
#' @examples
#'   data(whales)
#'   motif = char2dna("AAGTGTAGCATCACTTATTGATCCAAATT")
#'   shave(whales, motif = motif, direction = "both")
################################################################################
shave <- function(x, motif, direction = "both", cores = 1, ...){
  tmpattr <- attributes(x)
  tmpattr <- tmpattr[names(tmpattr) != "quality"]
  shave1 <- function(x, motif, direction){
    vit <- aphid::Viterbi(motif, x, type = "semiglobal", ... = ...)
    p <- attr(x, "quality") # shave quality scores too
    if(!is.null(p)) stopifnot(length(x) == length(p))
    if(identical(direction, "forward")){
      ntoshave <- match(c(0, 1), rev(vit$path)) - 1
      ntoshave <- min(ntoshave[!is.na(ntoshave)])
      res <- x[1:(length(x) - ntoshave)]
      if(!is.null(p)) attr(res, "quality") <- p[1:(length(x) - ntoshave)]
    }else if(identical(direction, "reverse")){
      ntoshave <- match(c(0, 1), vit$path) - 1
      ntoshave <- min(ntoshave[!is.na(ntoshave)])
      if(ntoshave > 0){
        res <- x[-(1:ntoshave)]
        if(!is.null(p)) attr(res, "quality") <- p[-(1:ntoshave)]
      }else{
        res <- x
      }
    }else if(identical(direction, "both")){
      ntoshavef <- match(c(0, 1), rev(vit$path)) - 1
      ntoshavef <- min(ntoshavef[!is.na(ntoshavef)])
      last <- length(x) - ntoshavef
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
      res <- parallel::parLapply(cores, x, shave1, motif, direction, ...)
    }else if(cores == 1){
      res <- lapply(x, shave1, motif, direction, ...)
    }else{
      #nseq <- length(x)
      navailcores <- parallel::detectCores()
      if(identical(cores, "autodetect")) cores <- navailcores - 1
      if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
      # if(cores > navailcores) stop("Number of cores to use is more than number available")
      if(cores > 1){
        cl <- parallel::makeCluster(cores)
        res <- parallel::parLapply(cl, x, shave1, motif, direction, ...)
        parallel::stopCluster(cl)
      }else{
        res <- lapply(x, shave1, motif, direction, ...)
      }
    }
  }else{
    res <- shave1(x, motif, direction, ...)
    tmpattr$quality <- attr(res, "quality")
    tmpattr$score <- attr(res, "score")
  }
  attributes(res) <- tmpattr
  return(res)
}
################################################################################

