#' Identify sequences with potentially incorrect lineage metadata
#'
#'  This function evaluates a list of DNA barcode
#'    sequences (a "DNAbin" object with lineage attributes)
#'    and returns the names of the sequences that
#'    may require further checking before being used as training
#'    data for downstream tree-learning operations.
#'
#' @param x a DNAbin object with lineage attributes
#' @param threshold numeric. The minimum sequence weight to
#'    prompt addition to the list of potentially incorrect
#'    sequences. Defaults to 5.
#' @return a character vector giving the names of the potentially incorrect sequences
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
################################################################################
check_metadata <- function(x, threshold = 5){
  lins <- attr(x, "lineage")
  names(lins) <- names(x)

  seqhashes <- unname(sapply(x, function(s) paste(openssl::md5(as.vector(s)))))
  linhashes <- unname(sapply(lins, function(s) paste(openssl::md5(s))))
  comhashes <- paste0(seqhashes, linhashes)
  discards <- duplicated(comhashes)
  xu <- subset.DNAbin(x, subset = !discards) #47664
  linsu <- lins[!discards]
  shu <- seqhashes[!discards]
  lhu <- linhashes[!discards]

  dupes <- duplicated(shu) # logical length x
  seqpointers <- integer(length(shu))
  dupehashes <- shu[dupes]
  uniquehashes <- shu[!dupes]
  seqpointers[!dupes] <- seq_along(uniquehashes)
  pd <- integer(length(dupehashes))
  for(i in unique(dupehashes)) pd[dupehashes == i] <- match(i, uniquehashes)
  seqpointers[dupes] <- pd

  dupes <- duplicated(lhu) # logical length x
  linpointers <- integer(length(lhu))
  dupehashes <- lhu[dupes]
  uniquehashes <- lhu[!dupes]
  linpointers[!dupes] <- seq_along(uniquehashes)
  pd <- integer(length(dupehashes))
  for(i in unique(dupehashes)) pd[dupehashes == i] <- match(i, uniquehashes)
  linpointers[dupes] <- pd


  # duplicates <- duplicated.DNAbin(x, point = TRUE)
  # pointers <- attr(duplicates, "pointers")

  linlist <- split(linsu, seqpointers)
  linlist <- linlist[sapply(linlist, length) > 1]
  getlinlen <- function(l) sum(gregexpr(";", l, fixed = TRUE)[[1]] > 0) + 1
  minlinlens <- sapply(linlist, function(e) min(sapply(e, getlinlen)))
  comlins <- sapply(linlist, insect:::.ancestor)
  comlinlens <- sapply(comlins, getlinlen)
  diflens <- minlinlens - comlinlens
  mismatches <- diflens > 2 & sapply(linlist, length) > threshold
  linlist2 <- linlist[mismatches]
  checkmd1 <- function(l, threshold){
    # l is a vector of semicolon-delimited lineage strings
    d <- as.dist(outer(l, l, FUN = Vectorize(insect:::.ldistance)))
    tree <- as.dendrogram(hclust(d, method = "average"))
    wgts <- aphid::weight(tree)
    dodgyseqs <- names(wgts)[wgts > threshold]
    return(dodgyseqs)
  }
  dodgies <- unlist(lapply(linlist2, checkmd1, threshold), use.names = FALSE)
  # now switch to grouping sequences by lineage


  seqlist <- split(xu, linpointers)
  #hashes <- sapply(x, function(s) paste(openssl::md5(as.vector(s))))
  hashlist <- split(shu, linpointers)
  mismatches2 <- sapply(lapply(hashlist, unique), length) > threshold
  seqlist2 <- seqlist[mismatches2]
  checkmd2 <- function(l, threshold){
    # l is a list of sequences (DNAbin object without lineage metadata)
    d <- phylogram::kdistance(l)
    tree <- as.dendrogram(hclust(d, method = "average"))
    wgts <- aphid::weight(tree)
    dodgyseqs <- names(wgts)[wgts > threshold]
    return(dodgyseqs)
  }
  dodgies2 <- unlist(lapply(seqlist2, checkmd2, threshold), use.names = FALSE)
  return(unique(c(dodgies, dodgies2)))
}
################################################################################
#' Evaluate sequences classification quality.
#' This function is used to check the efficiency of a tree-based sequence
#'   classification given a known lineage string.
#' @param predicted,actual semicolon-delimited lineage strings, where the former
#'   is predicted by a classification tree and the latter is the true lineage.
#' @return A signed integer.
#'   A negative number indicates an underclassification (a correct subset of the actual
#'   lineage), a zero represents a full correct classification,
#'   and a positive number indicates a over- or mis-classification.
check_classification <- function(predicted, actual){
  if(identical(predicted, actual)) return(0)
  splitfun <- function(s) strsplit(s, split = "; ")[[1]]
  actual <- splitfun(gsub("\\.$", "", actual))
  if(identical(predicted, "")) return(length(actual) * -1)
  predicted <- splitfun(gsub("\\.$", "", predicted))
  diflen <- length(predicted) - length(actual)
  if(diflen > 0) actual <- c(actual, rep("", diflen))
  matched <- TRUE
  for(i in seq_along(predicted)){
    if(actual[i] != predicted[i]){
      matched = FALSE
      break
    }
  }
  if(i == length(predicted) & matched){
    return(diflen)
  }else{
    return(length(predicted) - i + 1)
  }
}
################################################################################
