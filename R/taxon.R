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
