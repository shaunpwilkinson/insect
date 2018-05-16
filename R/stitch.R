#' Paired-end read stitching.
#'
#' This function aligns forward and reverse reads generated from a
#'   2x amplicon sequencing platform, and produces a consensus sequence with
#'   maximum base-call quality scores attached as "quality" attributes.
#'
#' @param x,y DNAbin objects with quality attributes (see \code{\link{readFASTQ}} to
#'   generate these objects from fastq text files), representing the forward and
#'   reverse reads to be stitched. These objects can be either vectors or lists.
#'   In the latter case, the two objects must be equal length.
#' @param up,down forward and reverse primer sequences (either as concatenated
#'   character strings or "DNAbin" objects).
#'   Either both or neither should be provided (not one or the other).
#' @param mindiff the minimum difference in quality between two different base calls
#'   for the higher quality call to be added to the consensus alignment. In cases where
#'   the quality differences are less than this threshold, the ambiguity code "N" is added
#'   to the consensus sequence. Defaults to 6.
#' @param minoverlap integer giving the minimum number of nucleotides that
#'   should be shared between the forward and reverse sequence alignments.
#'   Defaults to 16.
#' @return a "DNAbin" object or a vector of concatenated character strings,
#'   depending on the input.
#' @author Shaun Wilkinson
#' @seealso \code{\link{readFASTQ}}.
#' @examples
#' \donttest{
#'   ## download and extract example FASTQ file to temporary directory
#'   td <- tempdir()
#'   URL <- "https://www.dropbox.com/s/71ixehy8e51etdd/insect_tutorial1_files.zip?dl=1"
#'   dest <- paste0(td, "/insect_tutorial1_files.zip")
#'   download.file(URL, destfile = dest, mode = "wb")
#'   unzip(dest, exdir = td)
#'   x <- readFASTQ(paste0(td, "/COI_sample1_read1.fastq"), bin = FALSE)
#'   y <- readFASTQ(paste0(td, "/COI_sample1_read2.fastq"), bin = FALSE)
#'   z <- stitch(x, y)
#'   z[1]
#'   attr(z, "quality")[1]
#' }
################################################################################
stitch <- function(x, y, up = NULL, down = NULL, mindiff = 6, minoverlap = 16){
  ## first check that sequences correspond
  stopifnot(!is.null(names(x)) & !is.null(names(y)))
  name11 <- strsplit(names(x)[1], split = "")[[1]]
  name21 <- strsplit(names(y)[1], split = "")[[1]]
  wa <- "x and y may not correspond to the same sequence\n"
  if(length(name11) != length(name21)) warning(wa)
  if(sum(suppressWarnings(name11 != name21)) > 2) warning(wa)
  isbin <- .isDNA(x)
  if(isbin){
    x <- dna2char(x)
    y <- dna2char(y)
  }else{ # should generally be DNAbin
    if(mode(x) != "character" | mode(y) != "character") stop("Invalid input\n")
    if(nchar(x[1]) == 1) stop("Expected concatenated character strings\n")
  }
  Q1 <- attr(x, "quality")
  Q2 <- attr(y, "quality")
  strev <- function(s) sapply(lapply(lapply(unname(s), charToRaw), rev), rawToChar)
  if(!is.null(up)){
    if(is.null(down)) stop("Expected both primers")
    if(.isDNA(up)){
      up <- disambiguate(dna2char(up))
      down <- disambiguate(dna2char(down))
    }else{
      if(mode(up) != "character" | mode(down) != "character") stop("Invalid primer formatting\n")
      up <- disambiguate(up)
      down <- disambiguate(down)
    }
    nchar1 <-  nchar(x)
    nchar2 <- nchar(y)
    lengthok <- nchar1 > 50 & nchar2 > 50
    Q1 <- Q1[lengthok]
    Q2 <- Q2[lengthok]
    x <- x[lengthok]
    y <- y[lengthok]
    g1f <- gregexpr(up, x)
    g2f <- gregexpr(down, y)
    isok <- function(z) z[1] > 0 & length(z) == 1
    forwards <- sapply(g1f, isok) & sapply(g2f, isok) #129
    g1r <- g2r <- vector(mode = "list", length = length(x))
    g1r[!forwards] <- gregexpr(down, x[!forwards]) #shorter
    g2r[!forwards] <- gregexpr(up, y[!forwards])
    reverses <- logical(length(x))
    if(any(!forwards)) reverses[!forwards] <- sapply(g1r[!forwards], isok) & sapply(g2r[!forwards], isok)
    if(sum(forwards) + sum(reverses) == 0) stop("Primer(s) not found in sequences\n")
    tmpR1 <- x
    tmpR2 <- y
    tmpQ1 <- Q1
    tmpQ2 <- Q2
    tmpR1[reverses] <- y[reverses]
    tmpR2[reverses] <- x[reverses]
    tmpQ1[reverses] <- Q2[reverses]
    tmpQ2[reverses] <- Q1[reverses]
    x <- tmpR1
    y <- tmpR2
    Q1 <- tmpQ1
    Q2 <- tmpQ2
    g1f[reverses] <- g1r[reverses]
    g2f[reverses] <- g2r[reverses]
    keeps <- forwards | reverses
    x <- x[keeps]
    y <- y[keeps]
    Q1 <- Q1[keeps]
    Q2 <- Q2[keeps]
    g1 <- g1f[keeps]
    g2 <- g2f[keeps]
    trimfun <- function(s, g) substring(s, first = g + attr(g, "match.length"))
    x <- mapply(trimfun, x, g1)
    y <- mapply(trimfun, y, g2)
    Q1 <- mapply(trimfun, Q1, g1, USE.NAMES = FALSE)
    Q2 <- mapply(trimfun, Q2, g2, USE.NAMES = FALSE)
  }

  y <- rc(y)
  Q2 <- sapply(Q2, strev, USE.NAMES = FALSE)
  # hashes <- hash(paste0(x, y)) # char type
  # pointers <- .point(hashes)
  # dupes <- duplicated(pointers)
  #
  find_alistart1 <- function(r1, r2, minoverlap){
    nchar2 <- nchar(r2)
    r2mers <- substring(r2, first = seq(1, nchar2 - (minoverlap - 1)),
                        last = seq(minoverlap, nchar2))
    res <- NA
    for(i in seq_along(r2mers)){
      g <- gregexpr(r2mers[i], r1)
      if(g[[1]][1] != -1 & length(g[[1]]) == 1){
        tmp <- g[[1]][1] - i + 1
        if(tmp > 0) res <- tmp
        break
      }
    }
    return(res)
  }

  # alistarts <- mapply(find_alistart, x[!dupes], y[!dupes], USE.NAMES = FALSE)[pointers]

  alistarts <- mapply(find_alistart1, x, y, minoverlap, USE.NAMES = FALSE)
  keeps <- !is.na(alistarts)
  x <- x[keeps]
  y <- y[keeps]
  attr(x, "quality") <- Q1[keeps]
  attr(y, "quality") <- Q2[keeps]
  alistarts <- alistarts[keeps]

  x <- char2dna(x)
  y <- char2dna(y)
  stitch1 <- function(r1, r2, alistart, mindiff){
    ## r1 is dnabin vec with qual attr
    ## r2 is same as x but rcd
    ## alistart scalar giving position of nt 1 of y in x
    q1 <- as.integer(attr(r1, "quality"))
    q2 <- as.integer(attr(r2, "quality"))
    alilen <- length(r1) - alistart + 1
    if(alilen > length(r2)) return(NA)
    r1overlap <- c(rep(FALSE, alistart - 1), rep(TRUE, alilen))
    r2overlap <- logical(length(r2))
    r2overlap[1:alilen] <- TRUE
    r1sub <- r1[r1overlap]
    r2sub <- r2[r2overlap]
    q1sub <- q1[r1overlap]
    q2sub <- q2[r2overlap]
    Rsub <- r1sub
    Qsub <- q1sub
    mismatches <- r1sub != r2sub
    if(any(mismatches)){
      # if(sum(mismatches) > 10) return(NA)
      diffs <- q1sub[mismatches] - q2sub[mismatches]
      # if(any(abs(diffs) < mindiff)) return(NA)
      r2better <- diffs < 0
      Rsub[mismatches][r2better] <- r2sub[mismatches][r2better]
      ambigs <- abs(diffs) < mindiff
      Rsub[mismatches][ambigs] <- as.raw(240)
    }
    q2better <- q2sub > Qsub
    Qsub[q2better] <- q2sub[q2better]
    R <- c(r1[!r1overlap], Rsub, r2[!r2overlap])
    Q <- c(q1[!r1overlap], Qsub, q2[!r2overlap])
    attr(R, "quality") <- as.raw(Q)
    return(R)
  }
  res <- mapply(stitch1, x, y, alistarts, mindiff)
  #for(j in seq_along(x)) res <- stitch1(R1s[[j]], R2s[[j]], alistarts[j])
  res <- res[!is.na(res)]
  if(isbin){
    if(length(res) == 0) return(raw(0))
    class(res) <- "DNAbin"
  }else{
    if(length(res) == 0) return(character(0))
    res <- dna2char(res)
  }
  return(res)
}

################################################################################
