#' Derive full lineage details from taxonomy ID number.
#'
#' This function returns the full lineage of a taxon ID number
#'   using the NCBI taxonomy database.
#'
#' @param taxIDs integer or vector of integers giving the taxon ID number(s).
#' @param db a copy of the NCBI taxonomy database (a data.frame object).
#'   See \code{\link{download_taxon}} for details.
#' @param simplify logical indicating whether a single lineage
#'   derived from a length-one input
#'   should be simplified from a list to a named character vector.
#'   Defaults to TRUE.
#' @param cores integer giving the number of CPUs to parallelize the operation
#'   over (Defaults to 1). This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   for example by running \code{parallel::stopCluster(cores)}.
#'   The string 'autodetect' is also accepted, in which case the maximum
#'   number of cores to use is one less than the total number of cores available.
#'   Note that in this case there
#'   may be a tradeoff in terms of speed depending on the number and size
#'   of sequences to be processed, due to the extra time required to initialize
#'   the cluster.
#' @return the full lineage as a named character vector, or list of named character
#'   vectors if the length of the input object is > 1 or simplify = FALSE.
#'   "names" attributes are taxonomic ranks.
#' @author Shaun Wilkinson
#' @references
#'  Federhen S (2012) The NCBI Taxonomy database.
#'  \emph{Nucleic Acids Research}
#'  \strong{40}, D136-D143. doi:10.1093/nar/gkr1178.
#'
#'  \url{https://www.ncbi.nlm.nih.gov/taxonomy/}
#' @examples
#' data(whales)
#' data(whale_taxa)
#' get_lineage(attr(whales, "taxID")[1], db = whale_taxa)
################################################################################
get_lineage <- function(taxIDs, db, simplify = TRUE, cores = 1){
  gl1 <- function(taxID, db){
    stopifnot(length(taxID) == 1 & mode(taxID) == "numeric")
    res <- resnames <- character(100)
    counter <- 1
    index <- match(taxID, db$tax_id)
    if(is.na(index)){
      # warning(paste("Taxon ID", taxID, "not found in database\n"))
      return(NA)
    }
    repeat{
      if(is.na(index)) break
      if(length(index) > 1) cat(index, "\n")
      res[counter] <- db$name[index]
      resnames[counter] <- db$rank[index]
      index <- db$parent_tax_index[index]
      counter <- counter + 1
    }
    res <- res[1:(counter - 1)]
    names(res) <- resnames[1:(counter - 1)]
    return(rev(res))
  }
  db$parent_tax_index <- match(db$parent_tax_id, db$tax_id)
  ## multithreading
  if(inherits(cores, "cluster")){
    res <- parallel::parLapply(cores, taxIDs, gl1, db)
  }else if(cores == 1){
    res <- lapply(taxIDs, gl1, db)
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' argument")
    if(cores > 1){
      cl <- parallel::makeCluster(cores)
      res <- parallel::parLapply(cl, taxIDs, gl1, db)
      parallel::stopCluster(cl)
    }else{
      res <- lapply(taxIDs, gl1, db)
    }
  }
  if(length(res) == 1 & simplify) res <- res[[1]]
  return(res)
}
################################################################################
#' Derive taxon ID from a lineage string or species name.
#'
#' This function returns the unique taxon ID associated with a given
#'   semicolon-delimited lineage string or taxon,
#'   by looking up the NCBI taxon database.
#'
#' @param lineage A semicolon-delimited lineage string or lineage name.
#' @param db the NCBI taxon database (as a data.frame object).
#'   See download_taxon for details.
#' @param multimatch character, the value to return if the query matches multiple
#'   entries in the database. Accepted values are "NA" (default), and "first"
#'   (the first match).
#' @return The unique taxon database ID (integer).
#' @details This function will return NA if the lineage is not found in the
#'   database or it matches multiple entries.
#' @author Shaun Wilkinson
#' @references
#'  Federhen S (2012) The NCBI Taxonomy database.
#'  \emph{Nucleic Acids Research}
#'  \strong{40}, D136-D143. doi:10.1093/nar/gkr1178.
#'
#'  \url{https://www.ncbi.nlm.nih.gov/taxonomy/}
#' @examples
#' data(whales)
#' data(whale_taxa)
#' get_taxID(attr(whales, "lineage")[1], db = whale_taxa)
################################################################################
get_taxID <- function(lineage, db, multimatch = "NA"){
  if(identical(lineage, "")) lineage <- "root"
  linvec <- rev(strsplit(lineage, split = "; ")[[1]])
  indices <- which(db$name == linvec[1])
  taxids <- db$tax_id[indices]
  if(length(taxids) == 1){
    return(taxids)
  }else if(length(taxids) > 1){
    if(!grepl(";", lineage)){
      # warning(paste("Taxon name matches multiple entries in database,",
      #            "returning first matching entry.\n",
      #            "Tip: try entering semicolon-delimited lineage string\n"))
      if(identical(multimatch, "NA")){
        return(NA)
      }else if(identical(multimatch, "first")){
        return(taxids[1])
      }else return(as.integer(multimatch))
    }else{
      tmp <- lapply(taxids, get_lineage, db = db)
      nmatches <- sapply(tmp, function(l) sum(sapply(l, grepl, lineage)))
      return(taxids[which.max(nmatches)])
    }
  }else{
    #warning("Not found in database\n"))
    return(NA)
  }
}
################################################################################
#' Download NCBI taxonomy database.
#'
#' This function accesses the NCBI API and gets an up-to-date copy of the taxonomy
#'   database.
#'
#' @param synonyms logical indicating whether synonyms should be included.
#'   Note that this increases the size of the returned object by around 10\%.
#' @param quiet logical indicating whether progress should be printed to the console.
#' @return a dataframe with the following elements:
#'   "tax_id", "parent_tax_id", "rank", "name".
#' @details
#'   This function downloads the NCBI taxon database as a data.frame
#'   object with the following columns:
#'   "tax_id", "parent_tax_id", "rank", "name".
#'   As of early 2018 the zip archive to download is approximately
#'   40Mb in size, and the output dataframe object is around
#'   200Mb in memory. Once downloaded, the dataframe can be pruned
#'   for increased speed and memory efficiency using the function
#'   \code{\link{prune_taxon}}.
#' @author Shaun Wilkinson
#' @references
#'  Federhen S (2012) The NCBI Taxonomy database.
#'  \emph{Nucleic Acids Research}
#'  \strong{40}, D136-D143. doi:10.1093/nar/gkr1178.
#'
#'  \url{https://www.ncbi.nlm.nih.gov/taxonomy/}
#' @examples
#' \donttest{
#'   taxonomy <- download_taxon()
#' }
################################################################################
download_taxon <- function(synonyms = FALSE, quiet = FALSE){
  tmp <- tempdir()
  fn <- "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
  download.file(fn, destfile = paste0(tmp, "/tmp.tar.gz"), quiet = quiet)
  if(!quiet) cat("Extracting data\n")
  test <- untar(tarfile = paste0(tmp, "/tmp.tar.gz"), exdir = tmp)
  if(!identical(test, 0L)) stop(cat(test))
  if(!quiet) cat("Building data frame\n")
  x <- scan(file = paste0(tmp, "/nodes.dmp"), what = "", sep = "\n", quiet = TRUE)
  x <- strsplit(x, split = "\t")
  x <- sapply(x, function(s) s[c(1, 3, 5)])
  nodes <- as.data.frame(t(x), stringsAsFactors = FALSE)
  nodes[[1]] <- as.integer(nodes[[1]])
  nodes[[2]] <- as.integer(nodes[[2]])
  colnames(nodes) <- c("tax_id", "parent_tax_id", "rank")
  #if(!quiet) cat("Parsing data frame\n")
  x <- scan(file = paste0(tmp, "/names.dmp"), what = "", sep = "\n", quiet = TRUE)
  if(synonyms){
    syn <- x[grepl("synonym", x)]
    syn <- strsplit(syn, split = "\t")
    syn <- sapply(syn, function(s) s[c(1, 3)])
    syn <- as.data.frame(t(syn), stringsAsFactors = FALSE)
    syn[[1]] <- as.integer(syn[[1]])
    colnames(syn) <- c("tax_id", "name")
  }
  x <- x[grepl("scientific name", x)] # 1637100 elements ~200 Mb, ~5 sec
  x <- strsplit(x, split = "\t")
  x <- sapply(x, function(s) s[c(1, 3)])
  namez <- as.data.frame(t(x), stringsAsFactors = FALSE)
  namez[[1]] <- as.integer(namez[[1]])
  colnames(namez) <- c("tax_id", "name")
  # merge node and names data frames together
  taxa <- merge(nodes, namez, by = "tax_id")
  taxa$parent_tax_id[taxa$tax_id == 1] <- 0
  if(synonyms){
    taxapp <- taxa[match(syn$tax_id, taxa$tax_id), ]
    taxapp$name <- syn$name
    rownames(taxapp) <- NULL
    taxa <- rbind(taxa, taxapp)
  } # attr(taxa, "synonyms") <- syn
  if(!quiet) cat("Done\n")
  return(taxa)
}
################################################################################
#' Prune taxa from taxonomy database.
#'
#' This function prunes the taxon database, removing specified taxa as
#'   desired to improve speed and memory efficiency.
#'
#' @param db a copy of the NCBI taxon database, obtained by running
#'   \code{\link{download_taxon}}.
#' @param taxIDs the names or taxon ID numbers of the taxa to be retained
#'   (or discarded if \code{keep = FALSE}).
#' @param keep logical, indicates whether the specified taxa should be
#'   kept and the rest of the database removed or vice versa. Defaults to TRUE.
#' @return a data.frame with the same column names as
#'   the input object
#'   ("tax_id", "parent_tax_id", "rank", "name").
#' @details TBA
#' @author Shaun Wilkinson
#' @references
#'  Federhen S (2012) The NCBI Taxonomy database.
#'  \emph{Nucleic Acids Research}
#'  \strong{40}, D136-D143. doi:10.1093/nar/gkr1178.
#'
#'  \url{https://www.ncbi.nlm.nih.gov/taxonomy/}
#' @examples
#' ## remove Odontoceti suborder from cetacean taxonomy database
#' data(whale_taxa)
#' prune_taxon(whale_taxa, taxIDs = 9722, keep = FALSE)
################################################################################
prune_taxon <- function(db, taxIDs, keep = TRUE){
  ### taxIDs can be character or taxIDs
  if(mode(taxIDs) == "character"){
    tmp <- db$tax_id[match(taxIDs, db$name)]
    if(any(is.na(tmp))) stop("Not found in db:", taxIDs[is.na(tmp)], "\n")
    taxIDs <- tmp
  }else{
    if(any(is.na(match(taxIDs, db$tax_id)))) {
      stop("Taxon IDs not found in db:",
           paste0(taxIDs[is.na(taxIDs)], collapse = " "),
           "\n")
    }
  }
  taxIDs <- unique(taxIDs)
  if(keep){
    rowstokeep <- integer(0)
    repeat{
      indices <- match(taxIDs, db$tax_id)
      if(any(is.na(indices))) stop("Invalid database\n")
      rowstokeep <- c(rowstokeep, indices)
      taxIDs <- unique(db$parent_tax_id[indices])
      if(all(taxIDs == 0L)) break
      taxIDs <- taxIDs[taxIDs != 0]
    }
    rowstokeep <- sort(unique(rowstokeep))
    db <- db[rowstokeep, ]
  }else{
    rowstodiscard <- match(taxIDs, db$tax_id)
    repeat{
      tmp <- which(db$parent_tax_id %in% taxIDs)
      if(length(tmp) == 0) break
      rowstodiscard <- c(rowstodiscard, tmp)
      taxIDs <- db$tax_id[tmp]
    }
    db <- db[-rowstodiscard, ]
  }
  rownames(db) <- NULL
  return(db)
}
################################################################################
