#' Convert sequences between binary and character string formats.
#'
#' These functions convert DNA and amino acid sequences in
#'   "DNAbin" or "AAbin" format to concatenated character strings,
#'   and vice versa.
#'
#' @param x a "DNAbin" or "AAbin" object.
#' @param z a vector of concatenated strings representing DNA or
#'   amino acid sequences in upper case.
#' @param simplify logical indicating whether length-one "DNAbin"
#'   or "AAbin" objects should be simplified to vectors.
#'   Defaults to FALSE.
#' @return \code{dna2char} and \code{aa2char} return upper case
#'   concatenated character strings.
#'   \code{char2dna} and \code{char2aa} return "DNAbin" and "AAbin" objects,
#'   respectively. These will be lists unless the input object
#'   has length one and simplify = TRUE, in which case the returned object
#'   will be a vector.
#' @author Shaun Wilkinson
#' @details
#'   These functions are used to convert concatenated character strings
#'   (e.g. "TAACGC") to binary format and vice versa.
#'   To convert DNAbin and AAbin objects to non-concatenated
#'   character objects (e.g. \code{c("T", "A", "A", "C", "G", "C")})
#'   refer to the \code{\link[ape]{ape}} package functions
#'   \code{\link[ape]{as.character.DNAbin}} and
#'   \code{\link[ape]{as.character.AAbin}}.
#'   Likewise the \code{\link[ape]{ape}} package functions
#'   \code{\link[ape]{as.DNAbin}} and \code{\link[ape]{as.AAbin}}
#'   are used to convert non-concatenated character
#'   objects to binary format.
#' @references
#'   Paradis E, Claude J, Strimmer K, (2004) APE: analyses of phylogenetics
#'   and evolution in R language. \emph{Bioinformatics} \strong{20}, 289-290.
#'
#'   Paradis E (2007) A bit-level coding scheme for nucleotides. \url{
#'   http://ape-package.ird.fr/misc/BitLevelCodingScheme.html}.
#'
#'   Paradis E (2012) Analysis of Phylogenetics and Evolution with R
#'   (Second Edition). Springer, New York.
#' @examples
#'   char2dna("TAACGC")
#'   char2aa("VGAHAGEY")
#'   dna2char(char2dna("TAACGC"))
#'   aa2char(char2aa("VGAHAGEY"))
#'   char2dna(list(seq1 = "TAACGC", seq2 = "ATTGCG"))
#'   char2aa(list(seq1 = "VGAHAGEY", seq2 = "VNVDEV"))
#' @name conversion
################################################################################
dna2char <- function(x){
  cbytes <- as.raw(c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77, 66, 86, 72, 68, 78, 45, 63))
  indices <- c(136, 24, 72, 40, 96, 144, 192, 48, 80 ,160, 112, 224, 176, 208, 240, 4, 2)
  vec <- raw(240)
  vec[indices] <- cbytes
  if(is.list(x)){
    res <- sapply(x, function(s) rawToChar(vec[as.integer(s)]))
    if(!is.null(attr(x[[1]], "quality"))){
      attr(res, "quality") <- unname(sapply(x, function(s) .qual2char(attr(s, "quality"))))
    }
  }else{
    res <- rawToChar(vec[as.integer(x)])
    if(!is.null(attr(x, "quality"))){
      attr(res, "quality") <- unname(.qual2char(attr(x, "quality")))
    }
  }
  attr(res, "rerep.names") <- attr(x, "rerep.names")
  attr(res, "rerep.pointers") <- attr(x, "rerep.pointers")
  return(res)
}
################################################################################
#' @rdname conversion
################################################################################
aa2char <- function(x){
  res <- if(is.list(x)) sapply(x, rawToChar) else rawToChar(x)
  attr(res, "rerep.names") <- attr(x, "rerep.names")
  attr(res, "rerep.pointers") <- attr(x, "rerep.pointers")
  return(res)
}
################################################################################
#' @rdname conversion
################################################################################
char2dna <- function(z, simplify = FALSE){
  dbytes <- as.raw(c(136, 24, 72, 40, 96, 144, 192, 48, 80 ,160, 112, 224, 176, 208, 240, 240, 4, 2))
  indices <- c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77, 66, 86, 72, 68, 78, 73, 45, 63) # max 89
  vec <- raw(89)
  vec[indices] <- dbytes
  s2d1 <- function(s) vec[as.integer(charToRaw(s))]
  if(length(z) == 1 & simplify){
    res <- s2d1(z)
    attr(res, "quality") <- .char2qual(attr(z, "quality")) # can be null
  }else{
    res <- lapply(z, s2d1)
    if(!is.null(attr(z, "quality"))){
      quals <- lapply(attr(z, "quality"), .char2qual) #list
      for(i in seq_along(res)) attr(res[[i]], "quality") <- quals[[i]]
    }
  }
  attr(res, "rerep.names") <- attr(z, "rerep.names")
  attr(res, "rerep.pointers") <- attr(z, "rerep.pointers")
  class(res) <- "DNAbin"
  return(res)
}
################################################################################
#' @rdname conversion
################################################################################
char2aa <- function(z, simplify = FALSE){
  res <- if(length(z) == 1 & simplify) charToRaw(z) else lapply(z, charToRaw)
  attr(res, "rerep.names") <- attr(z, "rerep.names")
  attr(res, "rerep.pointers") <- attr(z, "rerep.pointers")
  class(res) <- "AAbin"
  return(res)
}
################################################################################
