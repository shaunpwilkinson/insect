#' Trim primers.
#'
#' This function trims the primer sequence(s) from a set of DNA sequences.
#'
#' @param x a "DNAbin" object or a vector of concatenated upper-case character strings,
#'   representing the sequences to be trimmed.
#' @param up,down "DNAbin" objects or  vectors of concatenated upper-case character strings,
#'   representing the primer sequences.
#' @return a "DNAbin" object or a vector of concatenated character strings,
#'   depending on the input.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples ##TBA
################################################################################
trim <- function(x, up, down = NULL){ # trims primers and/or indices
  isbin <- .isDNA(x)
  if(isbin){
    x <- dna2char(x)
  }else{
    if(mode(x) != "character") stop("Invalid input\n")
  }
  quals <- attr(x, "quality")
  hasq <- !is.null(quals)
  if(.isDNA(up)) up <- dna2char(up)
  if(.isDNA(down)) down <- dna2char(down)
  up <- disambiguate(up)
  if(!is.null(down)) down <- disambiguate(rc(down))
  strev <- function(x) sapply(lapply(lapply(x, charToRaw), rev), rawToChar)
  isok <- function(z) z[1] > 0 & length(z) == 1
  pattern <- paste0(".*", up, ".+", if(!is.null(down)) paste0(down, ".*") else NULL)
  gx1 <- gregexpr(pattern, x)
  hasp <- sapply(gx1, isok) # has primers?
  x[!hasp] <- rc(x[!hasp]) # complement rejects
  if(hasq) quals[!hasp] <- sapply(quals[!hasp], strev)
  gx2 <-  gregexpr(pattern, x[!hasp]) # length = n rejects
  tmp <- sapply(gx2, isok) # length = n rejects
  gx1[!hasp][tmp] <- gx2[tmp]
  hasp[!hasp] <- tmp
  if(sum(hasp) == 0) return(if(isbin) raw(0) else character(0))
  starts <- unlist(gx1[hasp])
  stops <- sapply(gx1[hasp], function(g) attr(g, "match.length"))
  res <- mapply(substr, x[hasp], starts, stops)
  if(hasq) quals <- unname(mapply(substr, quals[hasp], starts, stops))
  attr(res, "quality") <- quals
  if(isbin) res <- char2dna(res)
  return(res)
}
################################################################################
