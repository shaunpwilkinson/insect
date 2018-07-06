#' Convert sequences to MD5 hashes.
#'
#'  This function converts DNA or amino acid sequences to 128-bit
#'    MD5 hash values for efficient duplicate identification and dereplication.
#'
#' @param x a sequence or list of sequences, either in character string,
#'   character vector, or raw byte format (eg DNAbin or AAbin objects).
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
hash <- function(x){
  if(mode(x) == "raw") x <- rawToChar(x)
  if(mode(x) == "list"){
    if(length(x) == 0){
      return(NULL)
    }else{
      if(mode(x[[1]]) == "raw"){
        x <- vapply(x, rawToChar, "", USE.NAMES = TRUE)
      }else if(mode(x[[1]]) == "character"){
        x <- vapply(x, paste0, "", collapse = "", USE.NAMES = TRUE)
      }else{
        stop("Invalid input format\n")
      }
    }
  }
  if(mode(x) == "character"){
    nms <- names(x)
    res <- openssl::md5(x)
    names(res) <- nms
    return(res)
  }else{
    stop("Invalid input format\n")
  }
}
################################################################################
