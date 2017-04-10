#' Write DNA to text in FASTQ format.
#'
#' This function takes a list of DNA sequences in \code{DNAbin} format
#'   with "quality" attributes (e.g. Phred scores) and outputs
#'   a text file in FASTQ format to a specified directory.
#'
#' @param x a list object of class \code{DNAbin} with each element also having
#'   a vector of quality scores of equal length to the sequence. These vectors
#'   are comprised of raw bytes ranging from 00 to 5d (0 to 93 when converted to
#'   integers; this format is more memory efficient than the ASCII characters used
#'   in the FASTQ text files). See \code{\link{readFASTQ}} for more details.
#' @param file character string giving a valid file path to output the text to.
#'   If file = "" (default setting) the text file is written to the
#'   current working directory.
#' @param append logical indicating whether text should be appended below
#'   existing text in the file (TRUE), or whether any existing text should be
#'   overwritten (FALSE; default). Only applicable if the file specified
#'   in \code{file} already exists.
#' @details TBA
#' @author Shaun Wilkinson
#' @references
#'   Illumina help page:
#'   \url{https://help.basespace.illumina.com/articles/descriptive/fastq-files/}
#'
#' @seealso \code{\link{readFASTQ}} for reading FASTQ filed into R,
#'   and \code{\link[ape]{write.dna}} in the ape package
#'   for writing DNA to text in FASTA and other formats.
#' @examples
#'   \dontrun{
#'   ## TBA
#'   }
################################################################################
writeFASTQ <- function(x, file = "", append = FALSE){
  res <- vector(mode = "character", length = length(x) * 4)
  sar <- seq_along(res)
  res[sar %% 4 == 1] <- names(x)
  res[sar %% 4 == 2] <- sapply(x, .dna2char)
  res[sar %% 4 == 3] <- rep("+", length(x))
  res[sar %% 4 == 0] <- sapply(lapply(x, attr, "quality"), .qual2char)
  if(file == ""){
    return(res)
  }else{
    cat(res, file = file, append = append, sep = "\n")
  }
}
