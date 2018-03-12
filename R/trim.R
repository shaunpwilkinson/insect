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
#' @details Any sequences not containing the primer(s) in either direction are discarded.
#'   Hence this function can also be used to de-multiplex sequences and remove indices,
#'   though it will generally be faster to do this on the sequencing platform prior to
#'   exporting the FASTQ files.
#' @author Shaun Wilkinson
#' @references TBA
#' @examples ##TBA
################################################################################
trim <- function(x, up, down = NULL){ # trims primers and/or indices
  isbin <- .isDNA(x)
  if(isbin) x <- dna2char(x) else stopifnot(mode(x) == "character")
  quals <- attr(x, "quality")
  hasq <- !is.null(quals)
  if(.isDNA(up)) up <- dna2char(up)
  if(.isDNA(down)) down <- dna2char(down)
  ncup <- nchar(up) # store primer length prior to disambiguation
  ncdn <- nchar(down) # returns empty integer if null
  up <- disambiguate(up)
  if(!is.null(down)) down <- disambiguate(rc(down))
  strev <- function(x) sapply(lapply(lapply(x, charToRaw), rev), rawToChar)
  isok <- function(z) z[1] > 0 & length(z) == 1
  pattern <- paste0(up, ".+", down) ## down can be null
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
  stops <- sapply(gx1[hasp], function(g) attr(g, "match.length")) + starts - 1
  starts <- starts + ncup
  if(!is.null(down)) stops <- stops - ncdn
  discards <- stops - starts < 50
  hasp[discards] <- FALSE
  if(sum(hasp) == 0) return(if(isbin) raw(0) else character(0))
  starts <- starts[!discards]
  stops <- stops[!discards]
  res <- mapply(substr, x[hasp], starts, stops)
  if(hasq) quals <- unname(mapply(substr, quals[hasp], starts, stops))
  attr(res, "quality") <- quals
  if(isbin) res <- char2dna(res)
  return(res)
}
################################################################################
