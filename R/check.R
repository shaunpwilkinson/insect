#' Identify sequences with potentially erroneous lineage metadata.
#'
#'  This function evaluates a list of DNA barcode sequences
#'    (a "DNAbin" object with "taxID", "lineage" and/or "species" attributes;
#'    see \code{\link{searchGB}} for details)
#'    and returns a table identifying the sequences that
#'    may require further checking before being used as training
#'    data for downstream tree-learning operations.
#'
#' @param x a DNAbin list object with "taxID", "lineage" and/or "species" attributes
#'   (see \code{\link{searchGB}} for details).
#' @param db a copy of the NCBI taxonomy database as a data.frame object
#'   (see \code{\link{download_taxon}}).
#' @param level character string giving the taxonomic level at which
#'   heterogeneity within a cluster will flag a sequence as potentially erroneous.
#'   Should be a recognized rank within the NCBI taxonomy database.
#' @param threshold numeric between 0 and 1 giving the OTU similarity cutoff value
#'   with which to cluster the sequences.
#' @param quiet logical indicating whether progress should be printed to the console.
#' @return a data frame containing the names of the potentially erroneous sequences.
#'   The output object has zero rows if no sequences are flagged.
#' @details This function first clusters the sequence dataset into operational
#'   taxonomic units (OTUs) based on a given genetic distance threshold using the
#'   \code{\link[kmer]{otu}} function in the \code{\link{kmer}} package.
#'   It then proceeds to check each cluster
#'   for lineage homogeneity at a given taxonomic rank,
#'   flagging any records that appear out of place based on the taxa/lineages
#'   of the other OTU members.
#' @author Shaun Wilkinson
#' @examples
#'   data(whales)
#'   data(whale_taxa)
#'   dodgy_seqs <- check(whales, db = whale_taxa, level = "species")
#'   whales <- subset.DNAbin(whales, subset = !names(whales) %in% rownames(dodgy_seqs))
################################################################################
check <- function(x, db, level = "order", threshold = 0.97, quiet = FALSE){
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
  check1 <- function(y){
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
  lins <- get_lineage(attr(x, "taxID"), db)
  lins <- sapply(lins, function(e) e[level])
  names(lins) <- names(x)
  f <- as.factor(otus)[!is.na(lins)]
  lins <- lins[!is.na(lins)]
  splitlist <- split(lins, f)
  splitlist <- splitlist[tabulate(f) > 2]
  out <- lapply(splitlist, check1)
  out <- out[!sapply(out, is.null)]
  if(length(out) == 0){
    out <- data.frame(listed = character(0), suggested = character(0),
                      confidence = numeric(0), nstudies = integer(0))
    if(!quiet) cat("No sequences flagged\n")
    return(out)
  }
  names(out) <- NULL
  out <- do.call("rbind", out)
  out <- out[order(out$confidence, decreasing = TRUE), ]
  if(!quiet) cat("Found", nrow(out), "potentially erroneous sequences\n")
  return(out)
}
################################################################################

