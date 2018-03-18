#' Tabulate classification output.
#'
#' This function presents taxonomic classifications in tabular format
#'   for ease of interpretation, including sequence counts for each sample.
#'   This makes for easy exporting to csv or xlsx files.
#'
#' @param y a character vector of semicolon-delimited lineage strings
#'   or a list of such vectors, as output from the function \code{classify}.
#' @param db a copy of the NCBI taxonomy database.
#'   See \code{\link{download_taxon}} for details.
#' @param aggregated logical indicating whether the output table
#'   should be aggregated by taxon ID, rather than having one row
#'   for each unique sequence. Defaults to FALSE.
#' @param ranks character vector giving the taxonomic ranks to be
#'   included in the output table. Must be a valid rank from the
#'   NCBI taxonomy database (see https://www.ncbi.nlm.nih.gov/taxonomy).
#' @return a data frame.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples ##TBA
################################################################################
tabulize <- function(y, db, aggregated = FALSE,
                     ranks = c("kingdom", "phylum", "class", "order",
                               "family", "genus", "species")){
  ## y is output from classify
  if(!is.list(y)) y <- list(sample1 = y)
  nsites <- length(y)
  if(is.null(names(y))) names(y) <- paste0("sample", seq_along(y))
  origins <- rep(seq_along(y), sapply(y, length)) # original length, tabulatable
  #orignames <- names(y)
  newnames <- unlist(lapply(y, names), use.names = FALSE)
  yul <- unlist(y, use.names = FALSE)
  names(yul) <- newnames
  hshul <- unlist(lapply(y, attr, "hash"), use.names = FALSE) #unlisted hashes
  scrul <- unlist(lapply(y, attr, "score"), use.names = FALSE) # unlisted scores
  pointers <- .point(hshul)
  dupes <- duplicated(hshul)
  yulu <- yul[!dupes]# y unlisted & unique - length 13
  clashes <- hash(yulu)
  clpointers <- .point(clashes)
  taxIDs <- sapply(yulu[!duplicated(clpointers)], get_taxID, db = db)
  if(any(is.na(taxIDs))){
    warning(paste0(c("The following taxa were not found in database ",
                     "and were changed to root taxon:",
                     yul[!duplicated(clashes)][is.na(taxIDs)]), collapse = "\n"))
    taxIDs[is.na(taxIDs)] <- 0 # just change to root for now
  }
  taxvecs <- lapply(taxIDs, get_lineage, db = db)
  taxout <- data.frame(representative = names(yulu),
                       # md5 = hshul[!dupes],
                       taxID = taxIDs[clpointers],
                       taxon = sapply(taxvecs, tail, 1)[clpointers],
                       rank = sapply(taxvecs, function(e) tail(names(e), 1))[clpointers],
                       probability = scrul[!dupes], stringsAsFactors = FALSE)
  for(i in seq_along(ranks)){
    taxout[ranks[i]] <- sapply(taxvecs, function(v) if(is.na(v[ranks[i]])) "" else v[ranks[i]])[clpointers]
  }
  qout <- matrix(0L, nrow = nrow(taxout), ncol = nsites)
  ## get sample names from MiSeq dataset
  colnames(qout) <- names(y) # (list names)
  for(i in 1:nrow(qout)) qout[i, ] <- tabulate(origins[pointers == i], nbins = nsites)
  out <- cbind(taxout, qout)
  out$probability <- round(out$probability, 4)
  if(aggregated){
    aggtax <- aggregate(taxout[!colnames(taxout) %in% c("representative", "md5", "taxID", "probability")],
                        by = taxout["taxID"], FUN = head, 1)
    aggtax$n_unique <- aggregate(taxout$representative, by = taxout["taxID"], FUN = length)[[2]]
    aggquants <- aggregate(qout, by = taxout["taxID"], FUN = sum)[-1]
    out <- cbind(aggtax, aggquants)
  }
  rownames(out) <- NULL
  return(out)
}
################################################################################
