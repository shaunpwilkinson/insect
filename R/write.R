#' Write sequences to text in FASTA or FASTQ format.
#'
#' These functions take a list of DNA or amino acid sequences in
#'   \code{DNAbin} or \code{AAbin} format
#'   and outputs a text file to a specified directory.
#'
#' @param x a list of sequences in \code{DNAbin} or \code{AAbin} format, or a
#'   vector of sequences as concatenated upper-case character strings.
#'   For writeFASTQ, only DNAbin objects are accepted, and each element should have
#'   a vector of quality scores of equal length attributed to the sequence.
#'   These vectors are comprised of raw bytes ranging from 00 to 5d
#'   (0 to 93 when converted to integers).
#'   See \code{\link{readFASTQ}} for more details.
#' @param file character string giving a valid file path to output the text to.
#'   If file = "" (default setting) the text file is written to the
#'   console.
#' @param compress logical indicating whether the output file should be gzipped.
#' @param append logical indicating whether the output should be appended to the file.
#' @param wrap integer giving the maximum number of characters on each sequence line.
#'   Defaults to NULL (no wrapping).
#' @return NULL (invisibly).
#' @author Shaun Wilkinson
#' @references
#'   Illumina help page:
#'   \url{https://help.basespace.illumina.com/articles/descriptive/fastq-files/}
#'
#' @seealso \code{\link{readFASTQ}} for reading FASTQ files into R,
#'   and \code{\link[ape]{write.dna}} in the ape package
#'   for writing DNA to text in FASTA and other formats.
#' @examples
#' \donttest{
#'   ## download and extract example FASTQ file to temporary directory
#'   td <- tempdir()
#'   URL <- "https://www.dropbox.com/s/71ixehy8e51etdd/insect_tutorial1_files.zip?dl=1"
#'   dest <- paste0(td, "/insect_tutorial1_files.zip")
#'   download.file(URL, destfile = dest, mode = "wb")
#'   unzip(dest, exdir = td)
#'   x <- readFASTQ(paste0(td, "/COI_sample2.fastq"))
#'   ## trim primers from sequences
#'   mlCOIintF <- "GGWACWGGWTGAACWGTWTAYCCYCC"
#'   jgHCO2198 <- "TAIACYTCIGGRTGICCRAARAAYCA"
#'   x <- trim(x, up = mlCOIintF, down = jgHCO2198)
#'   ## quality filter with size selection and singleton removal
#'   x <- qfilter(x, minlength = 250, maxlength = 350)
#'   ## output filtered FASTQ file
#'   writeFASTQ(x, file = paste0(td, "/COI_sample2_filtered.fastq"))
#'   writeFASTA(x, file = paste0(td, "/COI_sample2_filtered.fasta"))
#'  }
#' @name write
################################################################################
writeFASTQ <- function(x, file = "", compress = FALSE, append = FALSE){
  if(.isDNA(x)) x <- dna2char(x)
  if(is.null(attr(x, "quality"))) stop("Sequences are missing quality attributes\n")
  reslen <- length(x) * 4
  res <- character(reslen)
  res[seq(1, reslen, by = 4)] <- if(grepl("^@", names(x)[1])) names(x) else paste0("@", names(x))
  res[seq(2, reslen, by = 4)] <- x
  res[seq(3, reslen, by = 4)] <- rep("+", length(x))
  res[seq(4, reslen, by = 4)] <- attr(x, "quality")
  f <- if(compress) gzfile(file, if(append) "a" else "w") else file(file, if(append) "a" else "w")
  write(res, f, append = append)
  close(f)
  invisible(NULL)
  # cat(res, file = file, sep = "\n", ... = ...)
}
################################################################################
#' @rdname write
################################################################################
writeFASTA <- function(x, file = "", compress = FALSE, append = FALSE, wrap = NULL){
  isDNA <- .isDNA(x)
  isAA <- .isAA(x)
  if(!is.null(dim(x))){ # convert from matrix to list while retaining gaps
    x <- as.list(as.data.frame(t(unclass(x))))
  }
  if(isDNA | isAA){
    tmp <- if(isDNA) dna2char(x) else aa2char(x)
  }else if(is.list(x)){
    if(length(x[[1]] == 1)){
      tmp <- unlist(x, use.names = TRUE)
    }else{
      tmp <- sapply(x, paste0, collapse = "")
    }
  }else{
    tmp <- x
  }
  reslen <- 2 * length(tmp)
  res <- character(reslen)
  res[seq(1, reslen, by = 2)] <- paste0(">", names(tmp))
  res[seq(2, reslen, by = 2)] <- tmp
  if(!is.null(wrap)){
    if(mode(wrap) == "integer") stop("wrap must be an integer")
    if(wrap == 0) stop("wrap must be more than zero")

    wrapfun <- function(s, wrap){ # character vector length 1
      nrows <- ceiling(nchar(s)/wrap)
      if(nrows == 1L) return(s)
      f <- rep(1:nrows, each = wrap)
      f <- f[seq_len(nchar(s))]
      spl <- split(strsplit(s, split = "")[[1]], f)
      spl <- vapply(spl, paste0, "", collapse = "", USE.NAMES = F)
      return(spl)
    }
    res <- as.list(res)
    res[seq(2, reslen, by = 2)] <- lapply(res[seq(2, reslen, by = 2)], wrapfun, wrap = wrap)
    res <- unlist(res, recursive = TRUE, use.names = FALSE)
  }
  #f <- if(compress) gzfile(file, "w") else file(file, "w")
  f <- if(compress) gzfile(file, if(append) "a" else "w") else file(file, if(append) "a" else "w")
  write(res, f, append = append)
  close(f)
  invisible(NULL)
  #cat(res, file = file, sep = "\n", ... = ...)
}
################################################################################
