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



# check_metadata <- function(x, threshold = 5){
#   lins <- attr(x, "lineage")
#   names(lins) <- names(x)
#   seqhashes <- unname(sapply(x, function(s) paste(openssl::md5(as.vector(s)))))
#   linhashes <- unname(sapply(lins, function(s) paste(openssl::md5(s))))
#   comhashes <- paste0(seqhashes, linhashes)
#   discards <- duplicated(comhashes)
#   xu <- subset.DNAbin(x, subset = !discards) #42092
#   linsu <- lins[!discards]
#   shu <- seqhashes[!discards]
#   lhu <- linhashes[!discards]
#   dupes <- duplicated(shu) # logical length x
#   seqpointers <- integer(length(shu))
#   dupehashes <- shu[dupes]
#   uniquehashes <- shu[!dupes]
#   seqpointers[!dupes] <- seq_along(uniquehashes)
#   pd <- integer(length(dupehashes))
#   for(i in unique(dupehashes)) pd[dupehashes == i] <- match(i, uniquehashes)
#   seqpointers[dupes] <- pd
#   dupes <- duplicated(lhu) # logical length x
#   linpointers <- integer(length(lhu))
#   dupehashes <- lhu[dupes]
#   uniquehashes <- lhu[!dupes]
#   linpointers[!dupes] <- seq_along(uniquehashes)
#   pd <- integer(length(dupehashes))
#   for(i in unique(dupehashes)) pd[dupehashes == i] <- match(i, uniquehashes)
#   linpointers[dupes] <- pd
#   linlist <- split(linsu, seqpointers)
#   linlist <- linlist[sapply(linlist, length) > 1]
#   getlinlen <- function(l) sum(gregexpr(";", l, fixed = TRUE)[[1]] > 0) + 1
#   minlinlens <- sapply(linlist, function(e) min(sapply(e, getlinlen)))
#   comlins <- sapply(linlist, .ancestor)
#   comlinlens <- sapply(comlins, getlinlen)
#   diflens <- minlinlens - comlinlens
#   mismatches <- diflens > 2 & sapply(linlist, length) > threshold
#   linlist2 <- linlist[mismatches]
#   wgt1 <- function(l){
#     # l is a vector of semicolon-delimited lineage strings
#     d <- as.dist(outer(l, l, FUN = Vectorize(.ldistance)))
#     tree <- as.dendrogram(hclust(d, method = "average"))
#     return(aphid::weight(tree))
#   }
#   wgts1 <- lapply(linlist2, wgt1)
#   names(wgts1) <- NULL
#   wgts1 <- unlist(wgts1, use.names = TRUE)
#   # now switch to grouping sequences by lineage
#   seqlist <- split(xu, linpointers)
#   hashlist <- split(shu, linpointers)
#   mismatches2 <- sapply(lapply(hashlist, unique), length) > threshold
#   seqlist2 <- seqlist[mismatches2]
#   wgt2 <- function(l){
#     # l is a list of sequences (DNAbin object without lineage metadata)
#     d <- phylogram::kdistance(l)
#     tree <- as.dendrogram(hclust(d, method = "average"))
#     return(aphid::weight(tree))
#   }
#   wgts2 <- lapply(seqlist2, wgt2)
#   names(wgts2) <- NULL
#   wgts2 <- unlist(wgts2, use.names = TRUE)
#   uniquenames <- unique(c(names(wgts1), names(wgts2)))
#   df1 <- data.frame(accession = names(wgts1), weight = wgts1, row.names = NULL)
#   df2 <- data.frame(accession = names(wgts2), weight = wgts2, row.names = NULL)
#   res <- merge(df1, df2, by = "accession", all.x = TRUE, all.y = TRUE)
#   colnames(res) <- c("accession", "lineage_weight", "sequence_weight")
#   l1 <- res$lineage_weight > threshold
#   l1[is.na(l1)] <- FALSE
#   l2 <- res$sequence_weight > threshold
#   l2[is.na(l2)] <- FALSE
#   res <- res[l1 | l2, ]
#   res$accession <- as.character(res$accession)
#   indices <- match(res$accession, names(x))
#   reshashes <- comhashes[indices]
#   #res$hash <- seqhashes[indices]
#   res$sequence <- sapply(x[indices], .dna2char)
#   for(i in 1:nrow(res)){
#     hashmatches <- which(comhashes == reshashes[i])
#     if(length(hashmatches) > 1){
#       for(j in 2:length(hashmatches)){
#         newrow <- res[i, ]
#         newrow$accession <- names(x)[hashmatches[j]]
#         res <- rbind(res, newrow)
#       }
#     }
#   }
#   res <- res[order(res$lineage_weight, res$sequence_weight, decreasing = TRUE),]
#   return(res)
# }
