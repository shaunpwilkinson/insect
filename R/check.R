#' Identify sequences with potentially incorrect lineage metadata
#'
#'  This function evaluates a list of DNA barcode
#'    sequences (a "DNAbin" object with lineage attributes)
#'    and returns a table of sequences that
#'    may require further checking before being used as training
#'    data for downstream tree-learning operations.
#'
#' @param x a DNAbin object with lineage attributes.
#' @param db the taxon database.
#' @param ranks character vector of taxonomic ranks to compare, in order of priority.
#' @param threshold numeric between 0 and 1 giving the OTU similarity cutoff.
#' @return a character vector giving the names of the potentially incorrect sequences
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
################################################################################
sift <- function(x, db, ranks = c("phylum", "class", "order"), threshold = 0.97){
  otus <- otu(x, threshold = threshold, k = 5)
  lins <- lapply(attr(x, "taxon"), get_lineage, db)
  sift1 <- function(y){# a vector of named taxa
    y <- y[!duplicated(paste0(gsub("(^.{4}).+", "\\1", names(y)), y))]
    # only compare one seq from each study
    if(length(unique(y)) < 2) return(NULL)
    tab <- sort(table(y), decreasing = TRUE)
    if(tab[1] == tab[2]) return(NULL)
    consensus <- names(tab)[1]
    dodgies <- y != consensus
    n <- length(dodgies)
    res <- data.frame(listed = y[dodgies],
                      suggested = rep(consensus, sum(dodgies)),
                      confidence = sum(!dodgies)/n, n = n)
    return(res)
  }
  out <- data.frame(rank = character(0), listed = character(0),
                    suggested = character(0), confidence = numeric(0), n = integer(0))
  for(r in ranks){
    rnk <- sapply(lins, function(e) e[r])
    names(rnk) <- names(x)
    rnkotus <- as.factor(otus)[!is.na(rnk)]
    rnk <- rnk[!is.na(rnk)]
    rnklist <- split(rnk, rnkotus)
    rnklist <- rnklist[tabulate(rnkotus) > 2]
    rnkout <- lapply(rnklist, sift1)
    rnkout <- rnkout[!sapply(rnkout, is.null)]
    if(length(rnkout) > 0){
      names(rnkout) <- NULL
      rnkout <- do.call("rbind", rnkout)
      rnkout <- cbind(data.frame(rank = rep(r, nrow(rnkout))), rnkout)
      rnkout <- rnkout[!rownames(rnkout) %in% rownames(out), ] #
      rnkout <- rnkout[order(rnkout$confidence, decreasing = TRUE), ]
      out <- rbind(out, rnkout)
    }
  }
  return(out)
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
