#' Dereplicate and rereplicate sequence datasets.
#'
#' These functions are used to extract only the unique sequences from
#'   a set of DNA reads, with the ability to rebuild the original
#'   sequence set at a later time.
#'
#' @param x a list of sequences in \code{DNAbin} or \code{AAbin} format, or a
#'   vector of sequences as concatenated upper-case character strings.
#' @param cores integer giving the number of CPUs to use (defaults to 1).
#'   This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   e.g. by running \code{parallel::stopCluster(cores)}.
#'   The string 'autodetect' is also accepted, in which case the maximum
#'   number of cores to use is one less than the total number of cores
#'   available.
#' @return either a DNAbin/AAbin object, or a vector of concatenated
#'   upper-case character strings, depending on the input object.
#' @author Shaun Wilkinson
#' @examples
#'   data(whales)
#'   tmp <- dereplicate(whales)
#'   whales <- rereplicate(tmp)
#' @name replicate
################################################################################
dereplicate <- function(x, cores = 1){
  hashes <- hash(x, cores = cores)
  pointers <- .point(hashes)
  orignames <- names(x)
  x <- x[!duplicated(pointers)]
  if(.isDNA(x)) for(i in seq_along(x)) attr(x[[i]], "quality") <- NULL
  attr(x, "rerep.names") <- orignames
  attr(x, "rerep.pointers") <- pointers
  return(x)
}
################################################################################
#' @rdname replicate
################################################################################
rereplicate <- function(x){
  if(is.null(attr(x, "rerep.names")) | is.null(attr(x, "rerep.pointers"))){
    stop("x is not rereplicable\n")
  }
  orignames <- attr(x, "rerep.names")
  x <- x[attr(x, "rerep.pointers")]
  names(x) <- orignames
  return(x)
}
################################################################################
