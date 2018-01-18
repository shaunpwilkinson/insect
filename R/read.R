#' Read and filter FASTQ files.
#'
#' \code{readFASTQ} is a text parser/filter that reads files in the FASTQ
#'   format into R, during which it optionally applies various sequence
#'   editing procedures and filters.
#'
#' @param file the name of the FASTQ file from which the sequences are to be read.
#' @param format the format of each element in the returned list. Accepted options are
#'   "string" (each sequence is a concatenated string with a similarly concatenated
#'   quality attribute comprised of metacharacters in the FASTQ coding scheme),
#'   "character" (each sequence is a character vector with one element for each
#'   nucleotide; quality attributes are integer vectors with values
#'   between 0 and 93), and "DNAbin" (each element is a raw vector in DNAbin format)
#' @param ... further arguments to be passed to \code{\link{scan}}.
#' @return Returns a list of DNA sequences, either in \code{DNAbin} format
#'   or as character vectors if \code{DNA} is set to \code{FALSE}. In either
#'   case the sequences in the returned list each have a "quality" attribute,
#'   which is a vector of quality scores between 0 and 93 either in raw bytes
#'   (if DNA = TRUE) or in the same ASCII coding scheme as used in the FASTQ
#'   file format.
#' @details
#'   \subsection{Compatibility:}{The FASTQ convention is somewhat
#'   ambiguous with several slightly different interpretations appearing
#'   in the literature. For now, this function only supports the Illumina
#'   convention for FASTQ files, where each sequence and its associated
#'   metadata occupies four line of the text file as follows : (1) the
#'   run and cluster metadata preceeded by an @ symbol; (2) the sequence
#'   itself in capitals without spaces;
#'   (3) a single "+" symbol; and (4) the Phred quality
#'   scores from 0 to 93 represented as ASCII symbols. For more information
#'   on this convention see the Illumina help page
#'   \href{https://help.basespace.illumina.com/articles/descriptive/fastq-files/}{here}
#'   .}
#'   \subsection{Speed and Memory Requirements:}{This function can
#'   take a while to process larger FASTQ files, a multithreading option
#'   will be available in a future version.}
#' @author Shaun Wilkinson
#' @references
#'   Bokulich NA, Subramanian S, Faith JJ, Gevers D, Gordon JI, Knight R,
#'   Mills DA, Caporaso JG (2013) Quality-filtering vastly improves diversity
#'   estimates from Illumina amplicon sequencing.
#'   \emph{Nat Methods}, \strong{1}, 57-59.
#'
#'   Illumina help page:
#'   \url{https://help.basespace.illumina.com/articles/descriptive/fastq-files/}
#'
#' @seealso \code{\link{writeFASTQ}} for writing DNAbin objects to text
#'   in the FASTQ format, and \code{\link[ape]{read.dna}} in the
#'   \code{\link[ape]{ape}} package for reading DNA in FASTA and other formats
#'   into R.
#' @examples
#'   \dontrun{
#'     ##TBA
#'   }
################################################################################
readFASTQ <- function(file, format = "DNAbin", ...){
  x <- scan(file = file, what = "", sep = "\n", quiet = TRUE, ... = ...)
  x <- x[seq_along(x) %% 4 != 3]
  seqs <- x[seq_along(x) %% 3 == 2]
  names(seqs) <- gsub("^@", "", x[seq_along(x) %% 3 == 1])
  quals <- x[seq_along(x) %% 3 == 0]
  # ape::as.DNAbin(seqs[1])
  if(format == "string"){
    mapply(function(x, y) structure(x, quality = y), seqs, quals, SIMPLIFY = FALSE) # 639.9mb
  }else if(format == "character"){
    seqs2 <- strsplit(seqs, split = "")
    quals2 <- lapply(quals, .char2qual)
    quals2 <- lapply(quals2, as.integer)
    mapply(function(x, y) structure(x, quality = y), seqs2, quals2, SIMPLIFY = FALSE) #2.3gb
  }else if(format == "DNAbin"){
    seqs2 <- lapply(seqs, .char2dna)
    quals2 <- lapply(quals, .char2qual)
    res <- mapply(function(x, y) structure(x, quality = y), seqs2, quals2, SIMPLIFY = FALSE) # 576mb
    class(res) <- "DNAbin"
    return(res)
  }else stop("Invalid format\n")
}

################################################################################
