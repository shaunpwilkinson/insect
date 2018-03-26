#' Reverse complement DNA in character string format.
#'
#' This function reverse complements a DNA sequence or set of DNA
#'   sequences that are in concatenated character string format.
#'
#' @param z a vector of concatenated strings representing DNA sequences
#'   (in upper case).
#' @return a vector of concatenated strings representing DNA sequences
#'   (in upper case).
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples rc("TATTG")
################################################################################
rc <- function(z){
  rc1 <- function(zz){
    s <- strsplit(zz, split = "")[[1]]
    s <- rev(s)
    dchars <- strsplit("ACGTMRWSYKVHDBNI", split = "")[[1]]
    comps <- strsplit("TGCAKYWSRMBDHVNI", split = "")[[1]]
    s <- s[s %in% dchars] # remove spaces etc
    s <- dchars[match(s, comps)]
    s <- paste0(s, collapse = "")
    return(s)
  }
  z <- toupper(z)
  tmpnames <- names(z)
  res <- unname(sapply(z, rc1))
  if(!is.null(attr(z, "quality"))){
    strev <- function(x) sapply(lapply(lapply(unname(x), charToRaw), rev), rawToChar)
    attr(res, "quality") <- unname(sapply(attr(z, "quality"), strev))
  }
  names(res) <- tmpnames
  return(res)
}
################################################################################
