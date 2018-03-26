#' Paired-end read stitching.
#'
#' This function aligns forward and reverse reads generated from the Illumina
#'   2x amplicon sequencing platforms, and produces a consensus sequence with
#'   maximum Phred scores attached as "quality" attributes.
#'
#' @param R1,R2 DNAbin objects with quality attributes (see \code{\link{readFASTQ}} to
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
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @seealso \code{\link{readFASTQ}} for reading FASTQ text files into R.
#' @examples
#'   \dontrun{
#'     ##TBA
#'   }
################################################################################
stitch <- function(R1, R2, up = NULL, down = NULL, mindiff = 6, minoverlap = 16){
  ## first check that sequences correspond
  stopifnot(!is.null(names(R1)) & !is.null(names(R2)))
  name11 <- strsplit(names(R1)[1], split = "")[[1]]
  name21 <- strsplit(names(R2)[1], split = "")[[1]]
  wa <- "R1 and R2 may not correspond to the same sequence\n"
  if(length(name11) != length(name21)) warning(wa)
  if(sum(suppressWarnings(name11 != name21)) > 2) warning(wa)
  isbin <- .isDNA(R1)
  if(isbin){
    R1 <- dna2char(R1)
    R2 <- dna2char(R2)
  }else{ # should generally be DNAbin
    if(mode(R1) != "character" | mode(R2) != "character") stop("Invalid input\n")
    if(nchar(R1[1]) == 1) stop("Expected concatenated character strings\n")
  }
  Q1 <- attr(R1, "quality")
  Q2 <- attr(R2, "quality")
  strev <- function(x) sapply(lapply(lapply(unname(x), charToRaw), rev), rawToChar)
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
    nchar1 <-  nchar(R1)
    nchar2 <- nchar(R2)
    lengthok <- nchar1 > 50 & nchar2 > 50
    Q1 <- Q1[lengthok]
    Q2 <- Q2[lengthok]
    R1 <- R1[lengthok]
    R2 <- R2[lengthok]
    g1f <- gregexpr(up, R1)
    g2f <- gregexpr(down, R2)
    isok <- function(z) z[1] > 0 & length(z) == 1
    forwards <- sapply(g1f, isok) & sapply(g2f, isok) #129
    g1r <- g2r <- vector(mode = "list", length = length(R1))
    g1r[!forwards] <- gregexpr(down, R1[!forwards]) #shorter
    g2r[!forwards] <- gregexpr(up, R2[!forwards])
    reverses <- logical(length(R1))
    if(any(!forwards)) reverses[!forwards] <- sapply(g1r[!forwards], isok) & sapply(g2r[!forwards], isok)
    if(sum(forwards) + sum(reverses) == 0) stop("Primer(s) not found in sequences\n")
    tmpR1 <- R1
    tmpR2 <- R2
    tmpQ1 <- Q1
    tmpQ2 <- Q2
    tmpR1[reverses] <- R2[reverses]
    tmpR2[reverses] <- R1[reverses]
    tmpQ1[reverses] <- Q2[reverses]
    tmpQ2[reverses] <- Q1[reverses]
    R1 <- tmpR1
    R2 <- tmpR2
    Q1 <- tmpQ1
    Q2 <- tmpQ2
    g1f[reverses] <- g1r[reverses]
    g2f[reverses] <- g2r[reverses]
    keeps <- forwards | reverses
    R1 <- R1[keeps]
    R2 <- R2[keeps]
    Q1 <- Q1[keeps]
    Q2 <- Q2[keeps]
    g1 <- g1f[keeps]
    g2 <- g2f[keeps]
    trimfun <- function(x, g) substring(x, first = g + attr(g, "match.length"))
    R1 <- mapply(trimfun, R1, g1)
    R2 <- mapply(trimfun, R2, g2)
    Q1 <- mapply(trimfun, Q1, g1, USE.NAMES = FALSE)
    Q2 <- mapply(trimfun, Q2, g2, USE.NAMES = FALSE)
  }

  R2 <- rc(R2)
  Q2 <- sapply(Q2, strev, USE.NAMES = FALSE)
  # hashes <- hash(paste0(R1, R2)) # char type
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

  # alistarts <- mapply(find_alistart, R1[!dupes], R2[!dupes], USE.NAMES = FALSE)[pointers]

  alistarts <- mapply(find_alistart1, R1, R2, minoverlap, USE.NAMES = FALSE)
  keeps <- !is.na(alistarts)
  R1 <- R1[keeps]
  R2 <- R2[keeps]
  attr(R1, "quality") <- Q1[keeps]
  attr(R2, "quality") <- Q2[keeps]
  alistarts <- alistarts[keeps]

  R1 <- char2dna(R1)
  R2 <- char2dna(R2)
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
  res <- mapply(stitch1, R1, R2, alistarts, mindiff)
  #for(j in seq_along(R1)) res <- stitch1(R1s[[j]], R2s[[j]], alistarts[j])
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
