#' Sequence quality filter.
#'
#' This function checks that the average Phred quality scores of sequences (imported
#'   from a FASTQ file) are above a given threshold, and that the maximum number
#'   of ambiguities is below another given threshold.
#'
#' @param x a vector of concatenated strings representing DNA sequences
#'   (in upper case) or a DNAbin list object with quality attributes.
#'   This argument will usually be the output from \code{readFASTQ}.
#' @param minqual integer, the minimum average quality score for a
#'   sequence to pass the filter. Defaults to 30.
#' @param maxambigs integer, the maximum number of ambiguities for a sequence
#'   to pass the filter. Defaults to 0.
#' @param mincount integer, the minimum acceptable number of occurences of a
#'   sequence for it to pass the filter. Defaults to 2 (removes singletons).
#' @param minlength integer, the minimum acceptable sequence length.
#'   Defaults to 50.
#' @param minlength integer, the maximum acceptable sequence length.
#'   Defaults to 500.
#' @return a possibly shortened object of the same type as the primary
#'   input argument
#'   (i.e. a "DNAbin" object if x is a "DNAbin" object and a vector
#'   of concatenated character strings otherwise).
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples ##TBA
################################################################################
qc <- function(x, minqual = 30, maxambigs = 0, mincount = 2,
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

# a logical vector the same length as the input object, with TRUE for
# sequences passing the filter and FALSE otherwise.
