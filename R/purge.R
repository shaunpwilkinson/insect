#' Identify and remove erroneous reference sequences.
#'
#'  This function evaluates a DNA reference database (a "DNAbin" object)
#'    and removes any sequences whose taxonomic metadata appear to be inconsistent
#'    with those of their most closely related sequences.
#'
#' @param x a DNAbin list object whose names include taxonomic identification numbers
#'   (see \code{\link{searchGB}} for details).
#' @param db a valid taxonomy database containing the taxonomic identification numbers
#'   included in the "names" attribute of the primary input object (a data.frame object;
#'   see \code{\link{taxonomy}}).
#' @param level character string giving the taxonomic level at which
#'   heterogeneity within a cluster will flag a sequence as potentially erroneous.
#'   This should be a recognized rank within the taxonomy database.
#' @param confidence numeric, the minimum confidence value for a sequence to be purged.
#'   For example, if \code{confidence = 0.8} (the default value) a sequence will only be
#'   purged if its taxonomy differs from at least four other independent sequences
#'   in its cluster.
#' @param cores integer giving the number of processors for multithreading. Defaults to 1.
#'   This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   for example by running \code{parallel::stopCluster(cores)}.
#'   The string 'autodetect' is also accepted, in which case the maximum
#'   number of cores to use is one less than the total number of cores available.
#'   Note that in this case there
#'   may be a tradeoff in terms of speed depending on the number and size
#'   of sequences to be processed, due to the extra time required to initialize
#'   the cluster.
#' @param quiet logical indicating whether progress should be printed to the console.
#' @param ... further arguments to pass to \code{\link[kmer]{otu}} (not including
#'   \code{nstart}).
#' @return a "DNAbin" object.
#' @details This function first clusters the sequence dataset into operational
#'   taxonomic units (OTUs) based on a given genetic similarity threshold
#'   using the \code{\link[kmer]{otu}} function from the \code{\link{kmer}}
#'   package.
#'   Each cluster is then checked for taxonomic homogeneity at a given rank,
#'   and any sequences that appear out of place are removed.
#'   The criteria for sequence removal are that at least two other independent
#'   studies should contradict the taxonomic metadata attributed to the sequence.
#' @author Shaun Wilkinson
#' @examples
#'   data(whales)
#'   data(whale_taxonomy)
#'   whales <- purge(whales, db = whale_taxonomy, level = "species",
#'                   threshold = 0.97, method = "farthest")
################################################################################
purge <- function(x, db, level = "order", confidence = 0.8,
                  cores = 1, quiet = FALSE, ...){
  db$rank <- tolower(db$rank)
  level <- tolower(level)
  stopifnot(level %in% db$rank)
  if(is.null(attr(x, "OTU"))){
    if(!quiet) cat("Clustering OTUs\n")
    otus <- kmer::otu(x, nstart = 20, ... = ...)
  }else{
    if(!quiet) cat("Obtaining OTU membership from input object\n")
    otus <- attr(x, "OTU")
    stopifnot(length(x) == length(otus))
  }
  if(!quiet) cat("Comparing lineage metadata within OTUs\n")
  purge1 <- function(y){
    ## input a vector of named lineages (non delimited) of a certain rank
    ## returns a df
    hashes <- paste0(gsub("(^.{4}).+", "\\1", names(y)), y)
    ## only compare one seq from each study
    yu <- y[!duplicated(hashes)]
    if(length(unique(yu)) < 2) return(NULL)
    tab <- sort(table(yu), decreasing = TRUE)
    if(tab[1] == tab[2]) return(NULL)
    consensus <- names(tab)[1]
    dodgies <- y != consensus
    dodgiesu <- yu != consensus
    nu <- length(dodgiesu)
    res <- data.frame(listed = y[dodgies],
                      suggested = rep(consensus, sum(dodgies)),
                      confidence = sum(!dodgiesu)/nu, nstudies = nu)
    return(res)
  }
  taxIDs <- as.integer(gsub(".+\\|", "", names(x)))
  lins <- get_lineage(taxIDs, db, cores = cores)
  lins <- sapply(lins, function(e) e[level])
  lins[is.na(lins)] <- ""
  names(lins) <- names(x)
  f <- as.factor(otus)# [!is.na(lins)]
  # lins <- lins[!is.na(lins)]
  splitlist <- split(lins, f)
  splitlist <- splitlist[tabulate(f) > 2]
  dodgytab <- lapply(splitlist, purge1)
  dodgytab <- dodgytab[!vapply(dodgytab, is.null, logical(1))]
  if(length(dodgytab) == 0){
    # dodgytab <- data.frame(listed = character(0), suggested = character(0),
    #                   confidence = numeric(0), nstudies = integer(0))
    if(!quiet) cat("No erroneous sequences to remove\n")
    return(x)
  }
  names(dodgytab) <- NULL
  dodgytab <- do.call("rbind", dodgytab)
  dodgytab <- dodgytab[dodgytab$confidence >= confidence, ]
  if(nrow(dodgytab) == 0){
    if(!quiet) cat("No erroneous sequences to remove\n")
    return(x)
  }
  dodgytab <- dodgytab[order(dodgytab$confidence, decreasing = TRUE), ]
  if(!quiet) cat("Removed", nrow(dodgytab), "potentially erroneous sequences: \n")
  if(!quiet) print(dodgytab)
  out <- subset.DNAbin(x, subset = !names(x) %in% rownames(dodgytab))
  return(out)
}
################################################################################

