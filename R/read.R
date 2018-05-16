#' Read FASTA and FASTQ files.
#'
#' Text parsing functions for reading sequences in the FASTA or FASTQ format into R.
#'
#' @param file the name of the FASTA or FASTQ file from which the sequences
#'   are to be read.
#' @param bin logical indicating whether the returned object should be in
#'   binary/raw byte format (i.e. "DNAbin" or "AAbin" objects for
#'   nucleotide and amino acid sequences, respectively).
#'   If FALSE a vector of named character strings is returned.
#' @param residues character string indicating whether the sequences to
#'   be read are composed of nucleotides ("DNA"; default) or amino acids ("AA").
#'   Only required for \code{readFASTA} and if \code{bin = TRUE}.
#' @param alignment logical indicating whether the sequences represent
#'   an alignment to be parsed as a matrix.
#'   Only applies to \code{readFASTA}.
#' @return Either a vector of character strings (if bin = FALSE),
#'   or a list of raw ("DNAbin" or "AAbin") vectors,
#'   with each element having a "quality" attribute.
#' @details
#'   \subsection{Compatibility:}{The FASTQ convention is somewhat
#'   ambiguous with several slightly different interpretations appearing
#'   in the literature. For now, this function supports the Illumina
#'   convention for FASTQ files, where each sequence and its associated
#'   metadata occupies four line of the text file as follows : (1) the
#'   run and cluster metadata preceded by an @ symbol; (2) the sequence
#'   itself in capitals without spaces;
#'   (3) a single "+" symbol; and (4) the Phred quality
#'   scores from 0 to 93 represented as ASCII symbols. For more information
#'   on this convention see the Illumina help page
#'   \href{https://help.basespace.illumina.com/articles/descriptive/fastq-files/}{here}
#'   .}
#'   \subsection{Speed and Memory Requirements:}{
#'   For optimal memory efficiency and compatibility with other functions,
#'   it is recommended to store sequences in raw byte format
#'   as either DNAbin or AAbin objects.
#'   For FASTQ files when bin = TRUE, a vector of quality scores
#'   (also in raw-byte format) is attributed to each sequence.
#'   These can be converted back to numeric quality scores with \code{as.integer}.
#'   For FASTQ files when bin = FALSE the function returns a vector with each
#'   sequence as a concatenated string with a similarly concatenated quality attribute
#'   comprised of the same ASCII metacharacters used in the FASTQ coding scheme.
#'
#'   This function can take a while to process larger FASTQ files,
#'   a multithreading option may be available in a future version.}
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
#' @seealso \code{\link{writeFASTQ}} and \code{\link{writeFASTA}}
#'   for writing sequences to text in the FASTA or FASTQ format.
#'   See also \code{\link[ape]{read.dna}} in the \code{\link[ape]{ape}} package.
#' @examples
#' \donttest{
#'   ## download and extract example FASTQ file to temporary directory
#'   td <- tempdir()
#'   URL <- "https://www.dropbox.com/s/71ixehy8e51etdd/insect_tutorial1_files.zip?dl=1"
#'   dest <- paste0(td, "/insect_tutorial1_files.zip")
#'   download.file(URL, destfile = dest, mode = "wb")
#'   unzip(dest, exdir = td)
#'   x <- readFASTQ(paste0(td, "/COI_sample2.fastq"))
#'  }
#' @name read
################################################################################
readFASTQ <- function(file = file.choose(), bin = TRUE){
  x <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  if(!grepl("^@", x[1])) stop("Not a valid fastq file\n")
  seqs <- toupper(x[seq(2, length(x), by = 4)])
  seqnames <- gsub("^@", "", x[seq(1, length(x), by = 4)])
  quals <- x[seq(4, length(x), by = 4)]
  if(bin){
    seqs2 <- char2dna(seqs)
    quals2 <- lapply(quals, .char2qual) #510 mb total
    res <- mapply(function(x, y) structure(x, quality = y), seqs2, quals2, SIMPLIFY = FALSE)
    names(res) <- seqnames
    class(res) <- "DNAbin"
    return(res)
  }else{
    attr(seqs, "quality") <- quals
    names(seqs) <- seqnames
    return(seqs)
  }
}
################################################################################
#' @rdname read
################################################################################
readFASTA <- function(file = file.choose(), bin = TRUE, residues = "DNA",
                      alignment = FALSE){
  x <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
  if(!grepl("^>", x[1])) stop("Not a valid fasta file\n")
  namelines <- grepl("^>", x)
  f <- cumsum(namelines)
  res <- split(x, f)
  resnames <- sapply(res, function(s) s[1])
  resnames <- gsub("^>", "", resnames)
  # resnames <- gsub("\\|.+", "", resnames)
  res <- toupper(sapply(res, function(s) paste0(s[-1], collapse = "")))
  names(res) <- resnames
  if(bin){
    residues <- toupper(residues)
    DNA <- identical(residues, "DNA")
    if(!DNA){
      if(!identical(residues, "AA")) stop("Invalid 'residues' argument\n")
    }
    res <- if(DNA) char2dna(res) else lapply(res, charToRaw)
    if(alignment){
      if(!all(sapply(res, length) == length(res[[1]]))){
        warning("alignment is TRUE but sequences differ in length\n")
      }
      suppressWarnings(res <- do.call("rbind", res))
      rownames(res) <- resnames
    }
    class(res) <- if(DNA) "DNAbin" else "AAbin"
  }else{
    if(alignment){
      res <- matrix(res, ncol = 1)
      rownames(res) <- resnames
    }
  }
  return(res)
}
################################################################################


# mapply(function(x, y) structure(x, quality = y), seqs, quals, SIMPLIFY = FALSE) # 639.9mb
# }else if(format == "character"){
#   seqs2 <- strsplit(seqs, split = "")
#   quals2 <- lapply(quals, .char2qual)
#   quals2 <- lapply(quals2, as.integer)
#   mapply(function(x, y) structure(x, quality = y), seqs2, quals2, SIMPLIFY = FALSE) #2.3gb
