#' Paired-end rread stitching.
#'
#' This function aligns forward and reverse reads generated from the Illumina
#'   2x amplicon sequencing platforms, and produces a consensus sequence with
#'   maximum Phred scores attached as "quality" attributes.
#'
#' @param fw,rv DNAbin objects with quality attributes (see \code{\link{readFASTQ}} to
#'   generate these objects from fastq text files), representing the forward and
#'   reverse reads to be stitched. These objects can be either vectors or lists.
#'   In the latter case, the two objects must be equal length.
#' @param mindiff the minimum difference in quality between two different base calls
#'   for the higher quality call to be added to the consensus alignment. In cases where
#'   the quality differences are less than this threshold, the ambiguity code "N" is added
#'   to the consensus sequence.
#' @param cores The number of threads to be used to process the operation.
#' @param ... further agruments to be passed to \code{\link[aphid]{align}}
#' @return a DNAbin object with quality attributes
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @seealso \code{\link{readFASTQ}} for generating DNAbin objects from FASTQ text files
#' @examples
#'   \dontrun{
#'     ##TBA
#'   }
################################################################################
stitch <- function(fw, rv, mindiff = 6, cores = 1, ...){
  stitch1 <- function(fw, rv){
    alig <- aphid::align(fw, ape::complement(rv), type = "semiglobal", ... = ...)
    if(ncol(alig) > length(fw) + length(rv) - 10) return(NULL) ## minimum 10 overlap
    gaps <- alig == as.raw(4)
    globalstart <- max(c(match(FALSE, gaps[1, ]), match(FALSE, gaps[2, ])))
    globalend <- max(c(match(FALSE, rev(gaps[1, ])), match(FALSE, rev(gaps[2, ]))))
    globalend <- ncol(alig) - globalend + 1
    if(globalstart >= globalend) return(NULL)
    quals <- matrix(0, nrow = nrow(alig), ncol = ncol(alig))
    quals[1, !alig[1, ] == as.raw(4)] <- as.integer(attr(fw, "quality"))
    quals[2, !alig[2, ] == as.raw(4)] <- rev(as.integer(attr(rv, "quality")))
    ## impute quality of internal gaps using nearest-neighbors
    quals[, globalstart:globalend][gaps[, globalstart:globalend]] <- NA
    while(any(is.na(quals[1, ]))){
      befores <- quals[1, which(is.na(quals[1, ])) - 1]
      afters <- quals[1, which(is.na(quals[1, ])) + 1]
      binded <- rbind(befores, afters)
      quals[1, is.na(quals[1, ])] <- apply(binded, 2, min, na.rm = TRUE)
    }
    while(any(is.na(quals[2, ]))){
      befores <- quals[2, which(is.na(quals[2, ])) - 1]
      afters <- quals[2, which(is.na(quals[2, ])) + 1]
      binded <- rbind(befores, afters)
      quals[2, is.na(quals[2, ])] <- apply(binded, 2, mean, na.rm = TRUE)
    }
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
  if(is.list(fw)){
    stitch2 <- function(x) stitch1(x[[1]], x[[2]]) # length 2 list
    tmp <- mapply(list, fw, rv, SIMPLIFY = FALSE)
    if(inherits(cores, "cluster")){
      res <- parallel::parLapply(cores, tmp, stitch2)
    }else if(cores == 1){
      res <- lapply(tmp, stitch2)
    }else{
      navailcores <- parallel::detectCores()
      if(identical(cores, "autodetect")) cores <- navailcores - 1
      if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
      # if(cores > navailcores) stop("Insufficient CPUs available")
      if(cores > 1){
        cl <- parallel::makeCluster(cores)
        res <- parallel::parLapply(cl, tmp, stitch2)
        parallel::stopCluster(cl)
      }else{
        res <- lapply(tmp, stitch2)
      }
    }
    res <- res[!sapply(res, is.null)]
  }else{
    res <- stitch1(fw, rv)
  }
  if(!is.null(res)) class(res) <- "DNAbin"
  return(res)
}
################################################################################
