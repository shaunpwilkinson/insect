#' Get full lineage details from a taxonomic ID number.
#'
#' This function derives the full lineage of a taxon ID number
#'   from a given taxonomy database.
#'
#' @param taxIDs integer or vector of integers giving the taxonomic ID number(s).
#' @param db a taxonomy database (a data.frame object).
#'   See \code{\link{taxonomy}} for details.
#' @param simplify logical indicating whether a single lineage
#'   derived from a length-one input
#'   should be simplified from a list to a named character vector.
#'   Defaults to TRUE.
#' @param numbers logical indicating whether the output string(s) should
#'    be comprised of the taxonomic ID numbers rather than taxon names.
#'    Defaults to FALSE.
#' @param cores integer giving the number of processors for multithreading (Defaults to 1).
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
#' data(whale_taxonomy)
#' taxIDs <- as.integer(gsub(".+\\|", "", names(whales)[1:2]))
#' get_lineage(taxIDs, db = whale_taxonomy)
################################################################################
get_lineage <- function(taxIDs, db, simplify = TRUE, numbers = FALSE, cores = 1){
  db$rank <- as.character(db$rank)
  db$name <- as.character(db$name) # avoid stringsasfactor issues
  pointers <- .point(paste(taxIDs))
  taxIDs <- taxIDs[!duplicated(pointers)]
  if("tax_id" %in% colnames(db)){
    colnames(db)[colnames(db) == "tax_id"] <- "taxID"
    colnames(db)[colnames(db) == "parent_tax_id"] <- "parent_taxID"
    ## for backwards compatibility
  }
  gl1 <- function(taxID, db){
    if(is.na(taxID)) return(NA)
    stopifnot(length(taxID) == 1 & mode(taxID) == "numeric")
    res <- if(numbers) integer(100) else character(100)
    resnames <- character(100)
    counter <- 1
    index <- match(taxID, db$taxID)
    if(is.na(index)){
      # warning(paste("Taxon ID", taxID, "not found in database\n"))
      return(NA)
    }
    repeat{
      if(is.na(index)) break
      # if(length(index) > 1) cat(index, "\n")
      res[counter] <- if(numbers) db$taxID[index] else db$name[index]
      resnames[counter] <- db$rank[index]
      index <- db$parent_tax_index[index]
      counter <- counter + 1
    }
    res <- res[1:(counter - 1)]
    names(res) <- resnames[1:(counter - 1)]
    return(rev(res))
  }
  db$parent_tax_index <- match(db$parent_taxID, db$taxID)
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
  res <- res[pointers]
  if(length(res) == 1 & simplify) res <- res[[1]]
  return(res)
}
################################################################################
#' Get taxon ID from taxonomy database.
#'
#' This function returns the unique ID for a specified taxon name by looking
#'   up a taxonomy database.
#'
#' @param lineage character vector of taxon names or semicolon-delimited
#'   lineage strings.
#' @param db a valid taxonomy database (as a data.frame object).
#'   See \code{\link{taxonomy}} for details.
#' @param multimatch integer giving the value to return if the query
#'   matches multiple entries. Defaults to NA_integer_.
#' @return An integer giving the unique taxon ID,
#'   or NA if the taxon is not found in the database.
#' @author Shaun Wilkinson
#' @references
#'  Federhen S (2012) The NCBI Taxonomy database.
#'  \emph{Nucleic Acids Research}
#'  \strong{40}, D136-D143. doi:10.1093/nar/gkr1178.
#'
#'  \url{https://www.ncbi.nlm.nih.gov/taxonomy/}
#' @examples
#' data(whale_taxonomy)
#' get_taxID("Odontoceti", db = whale_taxonomy)
################################################################################
get_taxID <- function(lineage, db, multimatch = NA){
  pointers <- .point(lineage)
  lineage <- lineage[!duplicated(pointers)]
  db <- db[!duplicated(db), ]
  multimatch <- as.integer(multimatch)
  res <- integer(length(lineage))
  linlist <- strsplit(lineage, split = "; ")
  spps <- vapply(linlist, tail, "", 1)
  duplicates <- duplicated(db$name)
  dupnames <- db$name[duplicates]
  multimatches <- spps %in% dupnames
  if(any(multimatches)){
    db2 <- prune_taxonomy(db, taxIDs = unique(db$taxID[db$name %in% spps[multimatches]]))
    gtid_slow <- function(x, db2, multimatch){ # x is a character vector
      indices <- which(db2$name == tail(x, 1)) # should be length >= 2 vector
      candidates <- get_lineage(db2$taxID[indices], db2, simplify = FALSE) # length 2 list
      candidates <- vapply(candidates, function(cdd, x) all(x %in% cdd), logical(1), x)
      return(if(sum(candidates) == 1) db2$taxID[indices[candidates]] else multimatch)
    }
    res[multimatches] <- vapply(linlist[multimatches], gtid_slow, 0L, db2, multimatch)
    res[!multimatches] <- db$taxID[match(spps[!multimatches], db$name)]
  }else{
    res <- db$taxID[match(spps, db$name)]
  }
  res <- res[pointers]
  return(res)
}
################################################################################
#' Download taxonomy database.
#'
#' This function downloads an up-to-date copy of the taxonomy database.
#'
#' @param db character string specifying which taxonomy database to download.
#'   Currently only "NCBI" is supported.
#' @param synonyms logical indicating whether synonyms should be included.
#'   Note that this increases the size of the returned object by around 10\%.
#' @return a dataframe with the following elements:
#'   "taxID", "parent_taxID", "rank", "name".
#' @details
#'   This function downloads the specified taxonomy database as a data.frame
#'   object with the following columns:
#'   "taxID", "parent_taxID", "rank", "name".
#'   As of early 2018 the zip archive to download the NCBI taxonomy database
#'   is approximately 40Mb in size, and the output dataframe object is around
#'   200Mb in memory. Once downloaded, the dataframe can be pruned
#'   for increased speed and memory efficiency using the function
#'   \code{\link{prune_taxonomy}}.
#' @author Shaun Wilkinson
#' @references
#'  Federhen S (2012) The NCBI Taxonomy database.
#'  \emph{Nucleic Acids Research}
#'  \strong{40}, D136-D143. doi:10.1093/nar/gkr1178.
#'
#'  \url{https://www.ncbi.nlm.nih.gov/taxonomy/}
#' @examples
#' \donttest{
#'   # db <- taxonomy()
#' }
################################################################################
taxonomy <- function(db = "NCBI", synonyms = FALSE){
  if(!identical(db, "NCBI")){
    stop("Only the NCBI taxonomy database is available in this version\n")
  }
  tmp <- tempdir()
  fn <- "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
  download.file(fn, destfile = paste0(tmp, "/tmp.tar.gz"))
  message("Extracting data\n")
  test <- untar(tarfile = paste0(tmp, "/tmp.tar.gz"), exdir = tmp)
  if(!identical(test, 0L)) stop(cat(test))
  message("Building data frame\n")
  x <- scan(file = paste0(tmp, "/nodes.dmp"), what = "", sep = "\n", quiet = TRUE)
  x <- strsplit(x, split = "\t")
  x <- sapply(x, function(s) s[c(1, 3, 5)])
  nodes <- as.data.frame(t(x), stringsAsFactors = FALSE)
  nodes[[1]] <- as.integer(nodes[[1]])
  nodes[[2]] <- as.integer(nodes[[2]])
  colnames(nodes) <- c("taxID", "parent_taxID", "rank")
  x <- scan(file = paste0(tmp, "/names.dmp"), what = "", sep = "\n", quiet = TRUE)
  if(synonyms){
    syn <- x[grepl("synonym", x)]
    syn <- strsplit(syn, split = "\t")
    syn <- sapply(syn, function(s) s[c(1, 3)])
    syn <- as.data.frame(t(syn), stringsAsFactors = FALSE)
    syn[[1]] <- as.integer(syn[[1]])
    colnames(syn) <- c("taxID", "name")
  }
  x <- x[grepl("scientific name", x)] # 1637100 elements ~200 Mb, ~5 sec
  x <- strsplit(x, split = "\t")
  x <- sapply(x, function(s) s[c(1, 3)])
  namez <- as.data.frame(t(x), stringsAsFactors = FALSE)
  namez[[1]] <- as.integer(namez[[1]])
  colnames(namez) <- c("taxID", "name")
  # merge node and names data frames together
  taxa <- merge(nodes, namez, by = "taxID")
  taxa$parent_taxID[taxa$taxID == 1] <- 0
  if(synonyms){
    taxapp <- taxa[match(syn$taxID, taxa$taxID), ]
    taxapp$name <- syn$name
    rownames(taxapp) <- NULL
    taxa <- rbind(taxa, taxapp)
  } # attr(taxa, "synonyms") <- syn
  message("Done\n")
  return(taxa)
}
################################################################################
#' Prune taxonomy database.
#'
#' This function prunes the taxon database, removing specified taxa as
#'   desired to improve speed and memory efficiency.
#'
#' @param db a valid taxonomy database, e.g. obtained by running the
#'   \code{\link{taxonomy}} function.
#' @param taxIDs the names or taxon ID numbers of the taxa to be retained
#'   (or discarded if \code{keep = FALSE}).
#' @param keep logical, indicates whether the specified taxa should be
#'   kept and the rest of the database removed or vice versa. Defaults to TRUE.
#' @return a data.frame with the same column names as
#'   the input object
#'   ("taxID", "parent_taxID", "rank", "name").
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
#' data(whale_taxonomy)
#' prune_taxonomy(whale_taxonomy, taxIDs = 9722, keep = FALSE)
################################################################################
prune_taxonomy <- function(db, taxIDs, keep = TRUE){
  ### taxIDs can be character or taxIDs
  if("tax_id" %in% colnames(db)){
    colnames(db)[colnames(db) == "tax_id"] <- "taxID"
    colnames(db)[colnames(db) == "parent_tax_id"] <- "parent_taxID"
    ## for backwards compatibility
  }
  if(mode(taxIDs) == "character"){
    tmp <- db$taxID[match(taxIDs, db$name)]
    if(any(is.na(tmp))) stop("Not found in db:", taxIDs[is.na(tmp)], "\n")
    taxIDs <- tmp
  }else{
    if(any(is.na(match(taxIDs, db$taxID)))) {
      stop("Taxon IDs not found in db:",
           paste0(taxIDs[is.na(taxIDs)], collapse = " "),
           "\n")
    }
  }
  taxIDs <- unique(taxIDs)
  if(keep){
    rowstokeep <- integer(0)
    repeat{
      indices <- match(taxIDs, db$taxID)
      if(any(is.na(indices))) stop("Invalid database format\n")
      rowstokeep <- c(rowstokeep, indices)
      taxIDs <- unique(db$parent_taxID[indices])
      if(all(taxIDs == 0L)) break
      taxIDs <- taxIDs[taxIDs != 0]
    }
    rowstokeep <- sort(unique(rowstokeep))
    db <- db[rowstokeep, ]
  }else{
    rowstodiscard <- match(taxIDs, db$taxID)
    repeat{
      tmp <- which(db$parent_taxID %in% taxIDs)
      if(length(tmp) == 0) break
      rowstodiscard <- c(rowstodiscard, tmp)
      taxIDs <- db$taxID[tmp]
    }
    db <- db[-rowstodiscard, ]
  }
  rownames(db) <- NULL
  return(db)
}
################################################################################
