#' Convert sequences to MD5 hashes.
#'
#'  This function converts DNA or amino acid sequences to 128-bit
#'    MD5 hash values for efficient duplicate identification and dereplication.
#'
#' @param x a sequence or list of sequences, either in character string,
#'   character vector, or raw byte format (eg DNAbin or AAbin objects).
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
#' @return a character vector.
#' @details This function uses the \code{md5} function from the openSSL library
#'   (\url{https://www.openssl.org/})
#'   to digest sequences to 128-bit hashes.
#'   These can be compared using base functions such
#'   as \code{duplicated} and \code{unique}, for fast identification and
#'   management of duplicate sequences in large datasets.
#' @author Shaun Wilkinson
#' @references
#'   Ooms J (2017) openssl: toolkit for encryption, signatures and
#'     certificates based on OpenSSL. R package version 0.9.7.
#'     \url{https://CRAN.R-project.org/package=openssl}
#' @examples
#'  data(whales)
#'  hashes <- hash(whales)
#'  sum(duplicated(hashes))
################################################################################
hash <- function(x, cores = 1){
  if(mode(x) == "raw") x <- list(x)
  hash1 <- function(s){
    if(mode(s) != "raw"){
      if(mode(s) == "character"){
        s <- charToRaw(s)
      }else if(mode(s) == "integer"){
        s <-  as.raw(s)
      }else if(mode(s) == "numeric"){
        stop("Can't digest numeric vectors")
      }
    }
    return(paste(openssl::md5(as.vector(s))))
  }
  if(inherits(cores, "cluster")){
    res <- parallel::parSapply(cores, x, hash1)
  }else if(cores == 1){
    res <- sapply(x, hash1)
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' argument")
    if(cores > 1){
      cl <- parallel::makeCluster(cores)
      res <- parallel::parSapply(cl, x, hash1)
      parallel::stopCluster(cl)
    }else{
      res <- sapply(x, hash1)
    }
  }
  return(res)
}
################################################################################
