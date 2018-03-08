#' Convert sequences between binary and string formats.
#'
#' These functions are used to convert DNA and amino acid sequences in
#'   "DNAbin" or "AAbin" format to character strings, and vice versa.
#'
#' @param x a "DNAbin" or "AAbin" object.
#' @param z a vector of concatenated strings representing DNA or amino acid sequences
#'   (upper case).
#' @param simplify logical indicating whether length-one "DNAbin" or "AAbin" objects
#'   should be simplified to vectors. Defaults to FALSE.
#' @return \code{dna2char} and \code{aa2char} return vectors of concatenated strings
#'   (in upper case).
#'   char2dna and char2aa return "DNAbin" and "AAbin" objects, respectively.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples ##TBA
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
    return(res)
  }else{
    res <- rawToChar(vec[as.integer(x)])
    if(!is.null(attr(x, "quality"))){
      attr(res, "quality") <- unname(.qual2char(attr(x, "quality")))
    }
    return(res)
  }
}
################################################################################
#' @rdname conversion
################################################################################
aa2char <- function(x){
  if(is.list(x)) sapply(x, rawToChar) else rawToChar(x)
}
################################################################################
#' @rdname conversion
################################################################################
char2dna <- function(z, simplify = FALSE){
  dbytes <- as.raw(c(136, 24, 72, 40, 96, 144, 192, 48, 80 ,160, 112, 224, 176, 208, 240, 4, 2))
  indices <- c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77, 66, 86, 72, 68, 78, 45, 63) # max 89
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
  class(res) <- "DNAbin"
  return(res)
}
################################################################################
#' @rdname conversion
################################################################################
char2aa <- function(z, simplify = FALSE){
  res <- if(length(z) == 1 & simplify) charToRaw(z) else lapply(z, charToRaw)
  class(res) <- "AAbin"
  return(res)
}
################################################################################
