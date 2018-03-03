#' Identify sequences with potentially incorrect lineage metadata.
#'
#'  This function evaluates a list of DNA barcode
#'    sequences (a "DNAbin" object with lineage attributes)
#'    and returns a table of sequences that
#'    may require further checking before being used as training
#'    data for downstream tree-learning operations.
#'
#' @param x a DNAbin object with lineage attributes.
#' @param db the taxon database.
#' @param level taxonomic level at which heterogeneity within an OTUs flags
#'   a sequence as potentially erroneous. Must be a recognized rank
#'   in the NCBI Taxon database.
#' @param threshold numeric between 0 and 1 giving the OTU similarity cutoff value.
#' @param quiet logical indicating whether progress should be printed to the console.
#' @return a data frame containing the names of the potentially erroneous sequences.
#' @details This function first clusters the sequence dataset into operational
#'   taxonomic units (OTUs) based on a given genetic distance threshold,
#'   and then checks for lineage homogeneity within OTUS, flagging any records
#'   that appear out of place based on the lineage metadata of other OTU members.
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
################################################################################
check <- function(x, db, level = "order", threshold = 0.97, quiet = FALSE){
  level <- tolower(level)
  if(!quiet) cat("Clustering OTUs\n")
  otus <- kmer::otu(x, threshold = threshold, k = 5)
  if(!quiet) cat("Comparing lineage metadata within OTUs\n")
  check1 <- function(y){
    ## input a vector of named lineages (non delimited) of a given rank
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
  lins <- lapply(attr(x, "taxID"), get_lineage, db)
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
    if(!quiet) cat("Found no erroneous sequences\n")
    return(out)
  }
  names(out) <- NULL
  out <- do.call("rbind", out)
  out <- out[order(out$confidence, decreasing = TRUE), ]
  if(!quiet) cat("Found", nrow(out), "potentially erroneous sequences\n")
  return(out)
}
################################################################################

