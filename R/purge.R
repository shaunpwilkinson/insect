#' Identify and remove sequences with potentially erroneous lineage metadata.
#'
#'  This function evaluates a DNA reference database
#'    (a "DNAbin" object with "taxID", "lineage" and/or "species" attributes;
#'    see \code{\link{searchGB}} for details)
#'    and removes any sequences whose taxonomic metadata appear to be inconsistent
#'    with those of the most closely related sequences.
#'
#' @param x a DNAbin list object with "taxID", "lineage" and/or "species" attributes
#'   (see \code{\link{searchGB}} for details).
#' @param db a copy of the NCBI taxonomy database as a data.frame object
#'   (see \code{\link{download_taxon}}).
#' @param level character string giving the taxonomic level at which
#'   heterogeneity within a cluster will flag a sequence as potentially erroneous.
#'   This should be a recognized rank within the NCBI taxonomy database.
#' @param threshold numeric between 0 and 1 giving the OTU similarity cutoff value
#'   with which to cluster the sequences. Defaults to 0.97.
#' @param cores integer giving the number of CPUs to parallelize the operation
#'   over. Defaults to 1.
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
#' @return a "DNAbin" object with "taxID", "lineage", and/or "species" attributes,
#'   depending on the attributes of the input object.
#' @details This function first clusters the sequence dataset into operational
#'   taxonomic units (OTUs) based on a given genetic similarity threshold
#'   using the \code{\link[kmer]{otu}} function from the \code{\link{kmer}}
#'   package.
#'   Each cluster is then checked for lineage homogeneity at a given taxonomic rank,
#'   and any sequences that appear out of place based on the taxon/lineage metadata
#'   of the other OTU members are removed.
#'   The criteria for sequence removal are that at least two other independent
#'   studies should contradict the taxonomic metadata attributed to the sequence.
#' @author Shaun Wilkinson
#' @examples
#'   data(whales)
#'   data(whale_taxa)
#'   whales <- purge(whales, db = whale_taxa, level = "species")
################################################################################
purge <- function(x, db, level = "order", threshold = 0.97, cores = 1,
                  quiet = FALSE){
  level <- tolower(level)
  if(is.null(attr(x, "OTU"))){
    if(!quiet) cat("Clustering OTUs\n")
    otus <- kmer::otu(x, threshold = threshold, k = 5, nstart = 20)
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
  if(is.null(attr(x, "taxID"))){
    if(is.null(attr(x, "lineage"))){
      if(is.null(attr(x, "species"))) {
        stop("taxID, lineage or species attribute required\n")
      }else{
        attr(x, "taxID") <- sapply(attr(x, "species"), get_taxID, db)
      }
    }else{
      attr(x, "taxID") <- sapply(attr(x, "lineage"), get_taxID, db)
    }
  }
  lins <- get_lineage(attr(x, "taxID"), db, cores = cores)
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
  dodgytab <- dodgytab[dodgytab$confidence > 0.66, ]
  if(length(dodgytab) == 0){
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

