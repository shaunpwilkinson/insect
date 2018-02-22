#' Write sequences to text in FASTA or FASTQ format.
#'
#' These functions take a list of DNA or amino acid sequences in
#'   \code{DNAbin} or \code{AAbin} format
#'   and outputs a text file to a specified directory.
#'
#' @param x a list of sequences in \code{DNAbin} or \code{AAbin} format, or a
#'   vector of sequences as character strings.
#'   For writeFASTQ, only DNAbin objects are accepted, and each element should have
#'   a vector of quality scores of equal length attributed to the sequence.
#'   These vectors are comprised of raw bytes ranging from 00 to 5d
#'   (0 to 93 when converted to integers; this format is more memory
#'   efficient than the ASCII characters used in the FASTQ text files).
#'   See \code{\link{readFASTQ}} for more details.
#' @param file character string giving a valid file path to output the text to.
#'   If file = "" (default setting) the text file is written to the
#'   console.
#' @param ... further options to be passed to \code{cat} (not including \code{"sep"}).
#' @details TBA
#' @author Shaun Wilkinson
#' @references
#'   Illumina help page:
#'   \url{https://help.basespace.illumina.com/articles/descriptive/fastq-files/}
#'
#' @seealso \code{\link{readFASTQ}} for reading FASTQ files into R,
#'   and \code{\link[ape]{write.dna}} in the ape package
#'   for writing DNA to text in FASTA and other formats.
#' @examples
#'   \dontrun{
#'   ## TBA
#'   }
#' @name write
################################################################################
writeFASTQ <- function(x, file = "", ...){
  if(!.isDNA(x)) stop("x must be a 'DNAbin' object")
  reslen <- length(x) * 4
  res <- character(reslen)
  res[seq(1, reslen, by = 2)] <- paste0("@", names(x))
  res[seq(2, reslen, by = 2)]  <- .dna2char(x)
  res[seq(3, reslen, by = 2)]  <- rep("+", length(x))
  res[seq(4, reslen, by = 2)]  <- sapply(lapply(x, attr, "quality"), .qual2char)
  cat(res, file = file, sep = "\n", ... = ...)
}
################################################################################
#' @rdname write
################################################################################
writeFASTA <- function(x, file = "", ...){
  isDNA <- .isDNA(x)
  isAA <- .isAA(x)
  if(!is.null(dim(x))){
    # convert from matrix to list while retaining gaps
    x <- as.list(as.data.frame(t(unclass(x))))
  }
  if(isDNA | isAA){
    if(isDNA){
      tmp <- .dna2char(x)
    }else{
      tmp <- if(is.list(x)) sapply(x, rawToChar) else rawToChar(x)
    }
  }else if(is.list(x)){
    if(length(x[[1]] == 1)){
      tmp <- unlist(x, use.names = TRUE)
    }else{
      tmp <- sapply(x, paste0, collapse = "")
    }
  }
  reslen <- 2 * length(tmp)
  res <- character(reslen)
  res[seq(1, reslen, by = 2)] <- paste0(">", names(tmp))
  res[seq(2, reslen, by = 2)] <- tmp
  cat(res, file = file, sep = "\n", ... = ...)
}
################################################################################
