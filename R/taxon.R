#' Derive lineage from taxon ID.
#'
#' This function returns the full lineage of a taxon ID number,
#'   given a taxonomic database.
#'
#' @param id integer giving the taxon ID number
#' @param db taxonomic database in matrix form
#' @return the full lineage as a named character vector of taxonomic ranks
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
#'
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
#' Derive taxon ID from lineage string.
#'
#' This function returns the unique taxon ID associated with a given
#'   semicolon-delimited lineage string.
#'
#' @param lineage character. A semicolon-delimited lineage string.
#' @param db taxonomic database in matrix form.
#' @return an integer. The unique taxon database ID.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
################################################################################
get_taxon <- function(lineage, db){
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
#' Download taxon database
#'
#' This function accesses the NCBI API and gets the latest copy of the taxon
#'   database
#'
#' @return a dataframe with the following columns: "tax_id", "parent_tax_id",
#'   "rank", "name", "parent_tax_index"
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @examples
#'   ##TBA
################################################################################
download_taxon <- function(){
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
  return(taxa)
}
################################################################################
