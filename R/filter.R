#' Quality filtering for amplicon sequences.
#'
#' This function performs several quality checks for FASTQ input files,
#'   removing any sequences that do not conform to the specified
#'   quality standards.
#'   This includes an average quality score assessment, size selection,
#'   singleton removal (or an alternative minimum count) and ambiguous
#'   base-call filtering.
#'
#' @param x a vector of concatenated strings representing DNA sequences
#'   (in upper case) or a DNAbin list object with quality attributes.
#'   This argument will usually be produced by \code{readFASTQ}.
#' @param minqual integer, the minimum average quality score for a
#'   sequence to pass the filter. Defaults to 30.
#' @param maxambigs integer, the maximum number of ambiguities for a sequence
#'   to pass the filter. Defaults to 0.
#' @param mincount integer, the minimum acceptable number of occurrences of a
#'   sequence for it to pass the filter. Defaults to 2 (removes singletons).
#' @param minlength integer, the minimum acceptable sequence length.
#'   Defaults to 50.
#' @param maxlength integer, the maximum acceptable sequence length.
#'   Defaults to 500.
#' @return an object of the same type as the primary input argument
#'   (i.e. a "DNAbin" object if x is a "DNAbin" object, or a vector
#'   of concatenated character strings otherwise).
#' @author Shaun Wilkinson
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
#'   ## filter sequences to remove singletons, low quality & short/long reads
#'   x <- qfilter(x, minlength = 250, maxlength = 350)
#'  }
################################################################################
qfilter <- function(x, minqual = 30, maxambigs = 0, mincount = 2,
               minlength = 50, maxlength = 500){
  if(.isDNA(x)){
    if(!is.list(x)) stop("Sequences in binary format should be in a list\n")
    if(!is.null(minlength)){
      lns <- sapply(x, length)
      keeps <- lns >= minlength
      if(sum(keeps) == 0) return(raw(0))
      x <- subset.DNAbin(x, subset = keeps)
    }
    if(!is.null(maxlength)){
      lns <- sapply(x, length)
      keeps <- lns <= maxlength
      if(sum(keeps) == 0) return(raw(0))
      x <- subset.DNAbin(x, subset = keeps)
    }
    if(!is.null(minqual)){
      if(is.null(attr(x[[1]], "quality"))) stop("Missing quality attributes\n")
      avq <- sapply(x, function(y) mean(as.integer(attr(y, "quality"))))
      keeps <- avq >= minqual
      if(sum(keeps) == 0) return(raw(0))
      x <- subset.DNAbin(x, subset = keeps)
    }
    if(!is.null(maxambigs)){
      nab <- sapply(x, function(y) sum(!(y %in% as.raw(c(136, 24, 72, 40)))))
      keeps <- nab <= maxambigs
      if(sum(keeps) == 0) return(raw(0))
      x <- subset.DNAbin(x, subset = keeps)
    }
    if(!is.null(mincount)){
      hashes <- hash(x)
      hashtab <- table(hashes)
      keeps <- hashes %in% names(hashtab)[hashtab >= mincount]
      if(sum(keeps) == 0) return(raw(0))
      x <- subset.DNAbin(x, subset = keeps)
    }
  }else{
    if(mode(x) != "character" | is.list(x)) stop("Invalid sequence format\n")
    if(!is.null(minlength)){
      lns <- nchar(x)
      keeps <- lns >= minlength
      if(sum(keeps) == 0) return(character(0))
      qs <- attr(x, "quality")
      x <- x[keeps]
      attr(x, "quality") <- qs[keeps]
    }
    if(!is.null(maxlength)){
      lns <- nchar(x)
      keeps <- lns <= maxlength
      if(sum(keeps) == 0) return(character(0))
      qs <- attr(x, "quality")
      x <- x[keeps]
      attr(x, "quality") <- qs[keeps]
    }
    if(!is.null(minqual)){
      qs <- attr(x, "quality")
      if(is.null(qs)) stop("Missing quality attributes\n")
      avq <- sapply(qs, function(y) mean(as.integer(.char2qual(y))))
      keeps <- avq >= minqual
      if(sum(keeps) == 0) return(character(0))
      x <- x[keeps]
      attr(x, "quality") <- qs[keeps]
    }
    if(!is.null(maxambigs)){
      nab <- sapply(gregexpr("[^ACGT]", x), function(y) if(y[1] == -1) 0 else length(y))
      keeps <- nab <= maxambigs
      if(sum(keeps) == 0) return(character(0))
      qs <- attr(x, "quality")
      x <- x[keeps]
      attr(x, "quality") <- qs[keeps]
    }
    if(!is.null(mincount)){
      hashes <- hash(x)
      hashtab <- table(hashes)
      keeps <- hashes %in% names(hashtab)[hashtab >= mincount]
      if(sum(keeps) == 0) return(character(0))
      qs <- attr(x, "quality")
      x <- x[keeps]
      attr(x, "quality") <- qs[keeps]
    }
  }
  return(x)
}
################################################################################
