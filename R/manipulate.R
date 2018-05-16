#' Further bit-level manipulation of DNA and amino acid sequences.
#'
#' These functions provide additional methods to manipulate objects of class
#'   \code{"DNAbin"} and \code{"AAbin"} to supplement those available in the
#'   \code{\link[ape]{ape}} package.
#'
#' @param x a \code{"DNAbin"} or \code{"AAbin"} object.
#' @param incomparables placeholder, not currently functional.
#' @param pointers logical indicating whether the re-replication index key
#'   should be returned as a \code{"pointers"} attribute of the output vector
#'   (only applicable for \code{duplicated.DNAbin} and \code{duplicated.AAbin}).
#'   Note that this can increase the computational
#'   time for larger sequence lists.
#' @param attrs logical indicating whether the attributes of the input object
#'   whose length match the object length (or number of rows if it is a matrix)
#'   should be retained and subsetted along with the object.
#'   This is useful if the input object has species, lineage and/or taxon ID
#'   metadata that need to be retained following the duplicate analysis.
#' @param drop logical; indicates whether the input matrix (assuming one is
#'   passed) should be reduced to a vector if subset to a single sequence.
#'   Defaults to FALSE in keeping with the style of the \code{\link{ape}}
#'   package functions.
#' @param subset logical vector giving the elements or rows to be kept.
#' @param ... further arguments to be passed between methods.
#' @return \code{unique} and \code{subset} return a
#'   \code{DNAbin} or \code{AAbin} object. \code{duplicated} returns a logical vector.
#' @author Shaun Wilkinson
#' @references
#'   Paradis E, Claude J, Strimmer K, (2004) APE: analyses of phylogenetics
#'   and evolution in R language. \emph{Bioinformatics} \strong{20}, 289-290.
#'
#'   Paradis E (2007) A bit-level coding scheme for nucleotides. \url{
#'   http://ape-package.ird.fr/misc/BitLevelCodingScheme.html}.
#'
#'   Paradis E (2012) Analysis of Phylogenetics and Evolution with R
#'   (Second Edition). Springer, New York.
#' @seealso \code{\link[ape]{DNAbin}}
#' @examples
#'   data(whales)
#'   duplicates <- duplicated.DNAbin(whales, point = TRUE)
#'   attr(duplicates, "pointers")
#'   ## returned indices show that the last sequence is
#'   ## identical to the second one.
#'   ## subset the reference sequence database to only include unques
#'   whales <- subset.DNAbin(whales, subset = !duplicates)
#'   ## this gives the same result as
#'   unique.DNAbin(whales)
#' @name manipulate
################################################################################
duplicated.DNAbin <- function(x, incomparables = FALSE, pointers = TRUE, ...){
  incomparables <- incomparables #TODO
  if(is.list(x)){
    if(pointers){
      hashes <- sapply(x, function(s) paste(openssl::md5(as.vector(s))))
      dupes <- duplicated(hashes, ... = ...) # logical length x
      pntrs <- .point(hashes)
      attr(dupes, "pointers") <- pntrs
    }else{
      dupes <- duplicated(lapply(x, as.vector)) # as.vector removes attributes
    }
    names(dupes) <- names(x)
  }else if(is.matrix(x)){
    if(pointers){
      hashes <- apply(x, 1, function(s) paste(openssl::md5(as.vector(s))))
      dupes <- duplicated(hashes, ... = ...)
      pntrs <- .point(hashes)
      attr(dupes, "pointers") <- pntrs
    }else{
      dupes <- duplicated(x, MARGIN = 1, ... = ...)
    }
    names(dupes) <- rownames(x)
  }else{
    dupes <- FALSE
    if(pointers) attr(dupes, "pointers") <- 1
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
#' @rdname manipulate
################################################################################
duplicated.AAbin <- function(x, incomparables = FALSE, pointers = TRUE, ...){
  duplicated.DNAbin(x, incomparables, pointers, ...)
}
################################################################################
#' @rdname manipulate
################################################################################
unique.AAbin <- function(x, incomparables = FALSE, attrs = TRUE,
                         drop = FALSE, ...){
  unique.DNAbin(x, incomparables, attrs, drop, ...)
}
################################################################################
#' @rdname manipulate
################################################################################
subset.AAbin <- function(x, subset, attrs = TRUE, drop = FALSE, ...){
  subset.AAbin(x, subset, attrs, drop, ...)
}
################################################################################
