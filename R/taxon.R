#' Convert taxon ID to semicolon-delimited lineage string.
#'
#' This function returns the full lineage of a taxon ID number
#'   using the NCBI Taxon database.
#'
#' @param id integer, the taxon ID number.
#' @param db the NCBI taxon database (as a data.frame object).
#'   See download_taxon for details.
#' @return the full lineage as a semicolon-delimited lineage string.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
################################################################################
get_lineage <- function(id, db){
  res <- resnames <- character(100)
  counter <- 1
  index <- match(id, db$tax_id)
  if(is.na(index)){
    warning(paste("Taxon id", id, "not found in database\n"))
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
################################################################################
#' Derive taxon ID from a given lineage.
#'
#' This function returns the unique taxon ID associated with a given
#'   semicolon-delimited lineage string or lineage name,
#'   using the NCBI taxon database.
#'
#' @param lineage A semicolon-delimited lineage string or lineage name.
#' @param db the NCBI taxon database (as a data.frame object).
#'   See download_taxon for details.
#' @return The unique taxon database ID (integer).
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
################################################################################
get_taxID <- function(lineage, db){
  if(identical(lineage, "")) lineage <- "root"
  #linvec <- rev(strsplit(gsub("\\.$", "", lineage), split = "; ")[[1]])
  linvec <- rev(strsplit(lineage, split = "; ")[[1]])
  indices <- which(db$name == linvec[1])
  taxids <- db$tax_id[indices]
  if(length(taxids) == 1){
    return(taxids)
  }else if(length(taxids) > 1){
    if(!grepl(";", lineage)){
      stop(paste("Taxon name matches multiple entries in database,",
                 "try entering semicolon-delimited lineage string"))
    }
    tmp <- lapply(taxids, get_lineage, db = db)
    nmatches <- sapply(tmp, function(l) sum(sapply(l, grepl, lineage)))
    return(taxids[which.max(nmatches)])
  }else{
    #warning("Not found in database\n"))
    return(NA)
  }
}
################################################################################
#' Download taxon database.
#'
#' This function accesses the NCBI API and gets an up-to-date copy of the taxon
#'   database.
#'
#' @param synonyms logical indicating whether a dataframe of synonyms and their
#'   associated taxon IDs should be attributed to the output object.
#'   Note that this increases the size of the returned object by around 10 percent.
#' @return a dataframe with the following elements:
#'   "tax_id", "parent_tax_id", "rank", "name", "parent_tax_index".
#' @details
#'   This function downloads the NCBI taxon database as a data.frame
#'   object with the following columns:
#'   "tax_id", "parent_tax_id", "rank", "name", "parent_tax_index".
#'   As of early 2018 the zip archive to download is approximately
#'   40Mb in size, and the output dataframe object is around
#'   200Mb in memory. Once downloaded, the dataframe can be pruned
#'   for increased speed and memory efficiency using the function
#'   prune_taxon.
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
################################################################################
download_taxon <- function(synonyms = FALSE){
  tmp <- tempdir()
  fn <- "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
  download.file(fn, destfile = paste0(tmp, "/tmp.tar.gz"))
  test <- untar(tarfile = paste0(tmp, "/tmp.tar.gz"), exdir = tmp)
  if(!identical(test, 0L)) stop(cat(test))
  x <- scan(file = paste0(tmp, "/nodes.dmp"), what = "", sep = "\n", quiet = TRUE)
  x <- strsplit(x, split = "\t")
  x <- sapply(x, function(s) s[c(1, 3, 5)])
  nodes <- as.data.frame(t(x), stringsAsFactors = FALSE)
  nodes[[1]] <- as.integer(nodes[[1]])
  nodes[[2]] <- as.integer(nodes[[2]])
  colnames(nodes) <- c("tax_id", "parent_tax_id", "rank")
  x <- scan(file = paste0(tmp, "/names.dmp"), what = "", sep = "\n", quiet = TRUE)
  if(synonyms){
    synonyms <- x[grepl("synonym", x)]
    synonyms <- strsplit(synonyms, split = "\t")
    synonyms <- sapply(synonyms, function(s) s[c(1, 3)])
    synonyms <- as.data.frame(t(synonyms), stringsAsFactors = FALSE)
    synonyms[[1]] <- as.integer(synonyms[[1]])
    colnames(synonyms) <- c("tax_id", "synonym")
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
  taxa$parent_tax_index <- match(taxa$parent_tax_id, taxa$tax_id)
  if(synonyms) attr(taxa, "synonyms") <- synonyms
  return(taxa)
}
################################################################################
#' Purge clades from taxon database.
#'
#' This function prunes the taxon database, removing specified clades as
#'   desired to improve speed and memory efficiency.
#'
#' @param db the taxon database, obtained by running
#'   download_taxon.
#' @param clade the name or taxon ID number of the clade to be removed
#'   (or kept if keep is TRUE).
#' @param keep logical, indicates whether the specified clade should be
#'   kept, and the rest of the database purged. Defaults to FALSE.
#' @return a smaller or equal sized dataframe with the same column names as
#'   the input object
#'   ("tax_id", "parent_tax_id", "rank", "name", "parent_tax_index").
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
################################################################################
prune_taxon <- function(db, clade, keep = FALSE){
  if(mode(clade) == "character") clade <- db$tax_id[match(clade, db$name)]
  if(is.na(clade)) stop("Clade not found in database\n")
  indices <- match(clade, db$tax_id)
  keeps <- db[indices, ]
  repeat{
    indices <- which(db$parent_tax_id %in% clade)
    if(length(indices) == 0) break
    clade <- db$tax_id[indices]
    keeps <- rbind(keeps, db[indices, ])
    db <- db[-indices, ]
  }
  if(keep) keeps else db
}
################################################################################
