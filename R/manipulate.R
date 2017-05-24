#' Further bit-level manipulation of DNA sequences.
#'
#' These functions provide additional methods to manipulate objects of class
#'   \code{"DNAbin"}.
#'
#' @param x a \code{DNAbin} object.
#' @param incomparables placeholder, not currently functional.
#' @param point logical indicating whether the duplication indices
#'   should be returned as a \code{"pointers"} attribute of the outputted
#'   logical vector (only applicable for \code{duplicated.DNAbin}).
#'   Note that this can significantly increase the computational
#'   time for larger DNAbin objects.
#' @param attrs logical indicating whether the attributes of the input object
#'   whose length match the object length (or number of rows if it is a matrix)
#'   should be retained and subsetted along with the object.
#'   This is useful if the input object has species, definition and/or lineage
#'   metadata that need to be retained following the duplicate analysis.
#' @param drop logical; indicates whether the input matrix (assuming one is
#'   passed) should be reduced to a vector if subset to a single sequence.
#'   Defaults to FALSE in keeping with the style of the \code{\link{ape}} package.
#' @param subset logical vector giving the elements or rows to be kept.
#' @param ... further agruments to be passed between methods.
#' @return \code{unique.DNAbin} and \code{subset.DNAbin} return a
#'   \code{DNAbin} object. \code{duplicated.DNAbin} returns a logical vector.
#' @details
#' \code{duplicated.DNAbin} returns a logical vector indicating the sequences
#'   that are duplicated, and optionally provides a "pointers" attribute
#'   indicating which of the unique sequences they are duplicates of.
#' \code{unique.DNAbin} takes a list or matrix of DNA sequences and returns
#' the subset of sequences that are unique.
#' @author Shaun Wilkinson
#' @references
#'   Paradis E (2007) A bit-level coding scheme for nucleotides. \url{
#'   http://ape-package.ird.fr/misc/BitLevelCodingScheme.html}.
#' @seealso \code{\link[ape]{DNAbin}}
#' @examples
#'   ## find duplicates in a small subsection of the woodmouse alignment
#'   library(ape)
#'   data(woodmouse)
#'   x <- woodmouse[, 100:150]
#'   duplicates <- duplicated.DNAbin(x, point = TRUE)
#'   duplicates
#'   ## shows that there are 3 unique sequences, No305, No304, and No0910S
#'   ## pointers attribute shows that most are identical to the second one (No304)
#'   ##
#'   ## reduce the alignment to only include the unique sequences
#'   x <- unique.DNAbin(x)
#'   x
#'   ## further subset the alignment to only include the first 2 sequences
#'   x <- subset.DNAbin(x, subset = c(TRUE, TRUE, FALSE))
#' @name manipulate
################################################################################
duplicated.DNAbin <- function(x, incomparables = FALSE, point = TRUE, ...){
  incomparables <- incomparables #TODO
  if(is.list(x)){
    if(point){
      hashes <- sapply(x, function(s) paste(openssl::md5(as.vector(s))))
      dupes <- duplicated(hashes, ... = ...) # logical length x
      pointers <- integer(length(x))
      dupehashes <- hashes[dupes]
      uniquehashes <- hashes[!dupes]
      pointers[!dupes] <- seq_along(uniquehashes)
      pd <- integer(length(dupehashes))
      for(i in unique(dupehashes)) pd[dupehashes == i] <- match(i, uniquehashes)
      pointers[dupes] <- pd
      #pointers[dupes] <- sapply(dupehashes, match, uniquehashes)
      attr(dupes, "pointers") <- pointers
    }else{
      dupes <- duplicated(lapply(x, as.vector)) # just removes attributes
    }
    names(dupes) <- names(x)
  }else if(is.matrix(x)){
    if(point){
      hashes <- apply(x, 1, function(s) paste(openssl::md5(as.vector(s))))
      dupes <- duplicated(hashes, ... = ...)
      pointers <- integer(nrow(x))
      dupehashes <- hashes[dupes]
      uniquehashes <- hashes[!dupes]
      pointers[!dupes] <- seq_along(uniquehashes)
      pd <- integer(length(dupehashes))
      for(i in unique(dupehashes)) pd[dupehashes == i] <- match(i, uniquehashes)
      pointers[dupes] <- pd
      #pointers[dupes] <- sapply(dupehashes, match, uniquehashes)
      attr(dupes, "pointers") <- pointers
    }else{
      dupes <- duplicated(x, MARGIN = 1, ... = ...)
    }
    names(dupes) <- rownames(x)
  }else{
    dupes <- FALSE
    if(point) attr(dupes, "pointers") <- 1
  }
  return(dupes)
}
################################################################################
#' @rdname manipulate
################################################################################
unique.DNAbin <- function(x, incomparables = FALSE, attrs = TRUE,
                          drop = FALSE, ...){
  incomparables <- incomparables
  if(is.list(x)){
    if(attrs){
      tmpattr <- attributes(x)
      if(!is.null(tmpattr)){
        whichattr <- which(sapply(tmpattr, length) == length(x))
      }else attrs <- FALSE
      if(length(whichattr) == 0) attrs <- FALSE
    }
    dupes <- duplicated(lapply(x, as.vector), ... = ...)
    x <- x[!dupes]
    if(attrs) {
      for(i in whichattr) tmpattr[[i]] <- tmpattr[[i]][!dupes]
      attributes(x) <- tmpattr
    }
  }else if(is.matrix(x)){
    if(attrs){
      tmpattr <- attributes(x)
      if(!is.null(tmpattr)){
        whichattr <- which(sapply(tmpattr, length) == nrow(x))
      }else attrs <- FALSE
      if(length(whichattr) == 0) attrs <- FALSE
    }
    dupes <- duplicated(x, MARGIN = 1, ... = ...)
    x <- x[!dupes, , drop = drop]
    if(attrs){
      for(i in whichattr) tmpattr[[i]] <- tmpattr[[i]][!dupes]
      attributes(x) <- tmpattr
    }
  }else{
    if(attr) x <- as.vector(x)
  }
  return(x)
}
################################################################################
#' @rdname manipulate
################################################################################
subset.DNAbin <- function(x, subset, attrs = TRUE, drop = FALSE, ...){
  if(mode(subset) != "logical") stop("Unsupported format for 'subset'")
  if(is.list(x)){
    if(length(subset) != length(x)) stop("'subset' and 'x' must be equal length")
    if(attrs){
      tmpattr <- attributes(x)
      if(!is.null(tmpattr)){
        whichattr <- which(sapply(tmpattr, length) == length(x))
      }else attrs <- FALSE
      if(length(whichattr) == 0) attrs <- FALSE
    }
    x <- x[subset]
    if(attrs){
      for(i in whichattr) tmpattr[[i]] <- tmpattr[[i]][subset]
      attributes(x) <- tmpattr
    }
  }else if(is.matrix(x)){
    if(length(subset) != nrow(x)) stop("'subset' and nrow(x) must be equal length")
    if(attrs){
      tmpattr <- attributes(x)
      if(!is.null(tmpattr)){
        whichattr <- which(sapply(tmpattr, length) == nrow(x))
      }else attrs <- FALSE
      if(length(whichattr) == 0) attrs <- FALSE
    }
    x <- x[subset, , drop = drop]
    if(attrs){
      for(i in whichattr) tmpattr[[i]] <- tmpattr[[i]][subset]
      attributes(x) <- tmpattr
    }
  }
  return(x)
}
################################################################################
