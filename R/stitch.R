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
#'   character strings or "DNAbin" objects). Either both or neither should be provided.
#' @param mindiff the minimum difference in quality between two different base calls
#'   for the higher quality call to be added to the consensus alignment. In cases where
#'   the quality differences are less than this threshold, the ambiguity code "N" is added
#'   to the consensus sequence.
#' @param cores The number of threads to be used to process the operation.
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
stitch <- function(R1, R2, up = NULL, down = NULL, mindiff = 6, cores = 1){
  ## forst check that sequences correspond
  stopifnot(!is.null(names(R1)) & !is.null(names(R2)))
  name11 <- strsplit(names(R1)[1], split = "")[[1]]
  name21 <- strsplit(names(R2)[1], split = "")[[1]]
  wa <- "R1 and R2 may not correspond to the same sequence\n"
  if(length(name11) != length(name21)) warning(wa)
  if(sum(suppressWarnings(name11 != name21)) > 2) warning(wa)
  isbin <- .isDNA(R1)
  if(isbin){
    R1c <- dna2char(R1)
    R2c <- dna2char(R2)
  }else{ # should generally be DNAbin
    if(mode(R1) != "character" | mode(R2) != "character") stop("Invalid input\n")
    if(nchar(R1[1]) < 20) stop("Expected concatenated character strings\n")
    R1c <- R1
    R2c <- R2
    R1 <- char2dna(R1)
    R2 <- char2dna(R2)
  }
  if(!is.null(up)){
    if(is.null(down)) stop("Expected both primers")
    if(.isDNA(up)){
      upc <- disambiguate(dna2char(up))
      downc <- disambiguate(dna2char(down))
    }else{
      if(mode(up) != "character" | mode(down) != "character") stop("Invalid primer formatting\n")
      upc <- disambiguate(up)
      downc <- disambiguate(down)
      up <- char2dna(up) #list w length 1
      down <- char2dna(down)
    }
    lengthok <- nchar(R1c) > 100 & nchar(R2c) > 100
    ## determine which sequences are in which orientation
    g1f <- gregexpr(upc, R1c)
    g2f <- gregexpr(downc, R2c)
    isok <- function(z) z[1] > 0 & length(z) == 1
    forwards <- lengthok & sapply(g1f, isok) & sapply(g2f, isok)
    g1r <- gregexpr(downc, R1c)
    g2r <- gregexpr(upc, R2c)
    reverses <- lengthok & sapply(g1r, isok) & sapply(g2r, isok) & !forwards
    if(sum(forwards) + sum(reverses) == 0) stop("Primer(s) not found in sequences\n")
    ## line all sequences up in a single orientation
    trimprimer <- function(s, b, l){#sequence (DNAbin w/ qual), begin, primer length
      rms <- seq(b, b + l - 1)
      tmpattr <- attr(s, "quality")[-rms]
      s <- s[-rms]
      attr(s, "quality") <- tmpattr
      return(s)
    }
    getl <- function(z) attr(z, "match.length")
    R1out1 <- mapply(trimprimer, R1[forwards], unlist(g1f[forwards]),
                    sapply(g1f[forwards], getl), SIMPLIFY = FALSE)
    R1out2 <- mapply(trimprimer, R2[reverses], unlist(g2r[reverses]),
                     sapply(g2r[reverses], getl), SIMPLIFY = FALSE)
    R1out <- c(R1out1, R1out2)
    R2out1 <- mapply(trimprimer, R2[forwards], unlist(g2f[forwards]),
                    sapply(g2f[forwards], getl), SIMPLIFY = FALSE)
    R2out2 <- mapply(trimprimer, R1[reverses], unlist(g1r[reverses]),
                     sapply(g1r[reverses], getl), SIMPLIFY = FALSE)
    R2out <- c(R2out1, R2out2)
    lengthok <- sapply(R1out, length) > 50 & sapply(R2out, length) > 50
    R1 <- R1out[lengthok]
    R2 <- R2out[lengthok]
  }

  stitch1 <- function(R1, R2, mindiff){
    alig <- aphid::align(R1, ape::complement(R2), type = "semiglobal")
    if(ncol(alig) > length(R1) + length(R2) - 10) return(NULL) ## minimum 10 overlap
    gaps <- alig == as.raw(4)
    globalstart <- max(c(match(FALSE, gaps[1, ]), match(FALSE, gaps[2, ])))
    globalend <- max(c(match(FALSE, rev(gaps[1, ])), match(FALSE, rev(gaps[2, ]))))
    globalend <- ncol(alig) - globalend + 1
    if(globalstart >= globalend) return(NULL)
    quals <- matrix(0, nrow = nrow(alig), ncol = ncol(alig))
    quals[1, !alig[1, ] == as.raw(4)] <- as.integer(attr(R1, "quality"))
    quals[2, !alig[2, ] == as.raw(4)] <- rev(as.integer(attr(R2, "quality")))
    ## impute quality of internal gaps using nearest-neighbors
    quals[, globalstart:globalend][gaps[, globalstart:globalend]] <- NA
    if(any(is.na(quals))) return(NULL)
    # impute <- function(v){} # quick dirty impute function in future?
    # while(any(is.na(quals[1, ]))){
    #   befores <- quals[1, which(is.na(quals[1, ])) - 1]
    #   afters <- quals[1, which(is.na(quals[1, ])) + 1]
    #   binded <- rbind(befores, afters)
    #   quals[1, is.na(quals[1, ])] <- apply(binded, 2, min, na.rm = TRUE)
    # }
    # while(any(is.na(quals[2, ]))){
    #   befores <- quals[2, which(is.na(quals[2, ])) - 1]
    #   afters <- quals[2, which(is.na(quals[2, ])) + 1]
    #   binded <- rbind(befores, afters)
    #   quals[2, is.na(quals[2, ])] <- apply(binded, 2, mean, na.rm = TRUE)
    # }
    maxquals <- apply(quals, 2, max)
    diffs <- quals[1, ] - quals[2, ]
    res <- raw(ncol(alig))
    resq <- integer(ncol(alig))
    bothsame <- apply(alig, 2, function(v) v[1] == v[2])
    res[bothsame] <- alig[1, bothsame]
    resq[bothsame] <- maxquals[bothsame]
    fwdonors <- diffs >= mindiff
    res[fwdonors] <- alig[1, fwdonors]
    resq[fwdonors] <- quals[1, fwdonors]
    rvdonors <- diffs <= -1 * mindiff
    res[rvdonors] <- alig[2, rvdonors]
    resq[rvdonors] <- quals[2, rvdonors]
    unknowns <- res == as.raw(0)
    res[unknowns] <- as.raw(240)
    resq[unknowns] <- maxquals[unknowns]
    gaps <- res == as.raw(4)
    resq <- resq[!gaps]
    res <- res[!gaps]
    attr(res, "quality") <- as.raw(resq)
    return(res)
  }
  if(is.list(R1)){
    stitch2 <- function(x, mindiff) stitch1(x[[1]], x[[2]], mindiff = mindiff) # x is length 2 list
    tmp <- mapply(list, R1, R2, SIMPLIFY = FALSE)
    if(inherits(cores, "cluster")){
      res <- parallel::parLapply(cores, tmp, stitch2, mindiff)
    }else if(cores == 1){
      res <- lapply(tmp, stitch2, mindiff)
    }else{
      navailcores <- parallel::detectCores()
      if(identical(cores, "autodetect")) cores <- navailcores - 1
      if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
      # if(cores > navailcores) stop("Insufficient CPUs available")
      if(cores > 1){
        cl <- parallel::makeCluster(cores)
        res <- parallel::parLapply(cl, tmp, stitch2, mindiff)
        parallel::stopCluster(cl)
      }else{
        res <- lapply(tmp, stitch2, mindiff)
      }
    }
    res <- res[!sapply(res, is.null)]
  }else{
    res <- stitch1(R1, R2)
  }
  if(isbin){
    if(!is.null(res)) class(res) <- "DNAbin"
  }else{
    res <- dna2char(res)
  }
  return(res)
}
################################################################################
