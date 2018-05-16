#' Query the NCBI GenBank database.
#'
#' \code{searchGB} searches GenBank using the
#'   Entrez search utilities, and downloads the matching sequences
#'   and/or their accession numbers. Alternatively a vector of
#'   GenBank accession numbers can be passed, in which case the function
#'   simply returns the matching sequences from GenBank.
#'   This function makes a series of API calls so internet connectivity is required.
#'
#' @param query an Entrez search query. For help compiling Entrez queries see
#'   \url{https://www.ncbi.nlm.nih.gov/books/NBK3837/#EntrezHelp.Entrez_Searching_Options}
#'   and \url{https://www.ncbi.nlm.nih.gov/books/NBK49540/}.
#' @param accession an optional vector of GenBank accession numbers to be input
#'   in place of a search query. If both query and accession arguments are
#'   provided the function returns an error.
#'   The function currently has a maximum batch size of around 200 accession numbers.
#' @param sequences logical. Should the sequences be returned or only the
#'   GenBank accession numbers? Note that taxon IDs and species
#'   names are not returned if \code{sequences} is set to FALSE.
#' @param bin logical indicating whether the returned sequences should
#'   be in raw-byte format ("DNAbin" or "AAbin" object type) or as a
#'   vector of named character strings. Defaults to TRUE.
#' @param db the NCBI database from which to download the sequences and/or
#'   accession names. Accepted options are "nucleotide" (default) or "protein".
#' @param taxIDs logical indicating whether the NCBI taxon ID numbers should
#'   be attributed to the output object. Defaults to FALSE.
#' @param species logical indicating whether the species names should
#'   be attributed to the output object. Defaults to FALSE.
#' @param lineages logical indicating whether semicolon-delimited lineage
#'   strings should be attributed to the output object. Defaults to FALSE.
#' @param prompt logical indicating whether to check with the user before
#'   downloading sequences.
#' @param contact an optional character string with the users email address.
#'   This is added to the E-utilities URL and may be used by NCBI to contact
#'   the user if the application causes unintended issues.
#' @param quiet logical indicating whether the progress should be printed
#'   to the console.
#' @return a list of sequences as a \code{DNAbin} or \code{AAbin} object,
#'   or a named vector of character strings (if \code{bin = FALSE}).
#' @details
#'   This function uses the Entrez e-utilities API to search and download
#'   sequences from GenBank.
#'   Occasionally users may encounter an unknown error that is not reproducible
#'   and appears to be related to database records being updated in GenBank.
#'   This can generally be remedied by re-running the function. If problems
#'   persist please feel free to raise an issue on the package bug-reports page at
#'   <http://github.com/shaunpwilkinson/insect/issues>.
#'   In some cases depending on the user security settings the R session
#'   may need to be opened with administrator privileges.
#' @author Shaun Wilkinson
#' @references
#'   NCBI Resource Coordinators (2012) Database resources of the National
#'   Center for Biotechnology Information. \emph{Nucleic Acids Research},
#'    \strong{41} (Database issue): D8–D20.
#' @seealso \code{\link[ape]{read.GenBank}} (ape)
#'   for an alternative means of downloading DNA sequences from GenBank
#'   using accession numbers.
#' @examples
#'   ## Query the GenBank database for Eukaryote mitochondrial 16S DNA sequences
#'   ## between 100 and 300 base pairs in length and last modified between
#'   ## 1999 and 2000.
#'   \donttest{
#'     query <- "Eukaryota[ORGN]+AND+16S[TITL]+AND+100:300[SLEN]+AND+1999:2000[MDAT]"
#'     x <- searchGB(query, species = TRUE, prompt = FALSE)
#'   }
################################################################################
searchGB <- function(query = NULL, accession = NULL, sequences = TRUE,
                      bin = TRUE, db = "nucleotide", taxIDs = FALSE, species = FALSE,
                      lineages = FALSE, prompt = TRUE, contact = NULL,
                     quiet = FALSE){
  if((is.null(query) & is.null(accession)) | (!is.null(query) & !is.null(accession))){
    stop("Either a query or accession number(s) must be provided\n")
  }
  if(!is.null(query)){
    URL1 <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
                   "db=", db,
                   "&term=",
                   query,
                   "&usehistory=y",
                   "&tool=R")
    X <- .scanURL(URL1, retmode = "xml")
    N <- as.integer(xml2::xml_text(xml2::xml_find_first(X, "Count")))
    WebEnv <- xml2::xml_text(xml2::xml_find_first(X, "WebEnv"))
    QueryKey <- xml2::xml_text(xml2::xml_find_first(X, "QueryKey"))
    if(N == 0){
      if(!quiet) warning("No sequences found matching query\n")
      if(sequences & bin) raw(0) else character(0)
    }
  }else if(!is.null(accession)){
    if(!sequences) stop("Accession numbers both provided and required\n")
    N <- length(accession)
    URL1 <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?",
                   "db=", db,
                   "&id=", paste(accession, collapse = ","))
    X <- .scanURL(URL1, retmode = "xml")
    WebEnv <- xml2::xml_text(xml2::xml_find_first(X, "WebEnv"))
    QueryKey <- xml2::xml_text(xml2::xml_find_first(X, "QueryKey"))
    Error <- xml2::xml_text(xml2::xml_find_first(X, "ERROR"))
    if(!is.na(Error)) stop(paste0(Error, "\n"))
  }
  if(prompt){
    decision <- readline(paste(N, if(sequences) "sequences" else "accession numbers",
                               "will be downloaded, continue? (y/n): "))
    if(decision == "n") stop("aborted") else if(!identical(decision, "y")) stop("invalid")
  }
  nrequest <- N%/%500 + as.logical(N%%500) # 1 if remainder exists, 0 otherwise
  accs <- vector(mode = "list", length = nrequest)
  if(sequences){
    seqs <- taxs <- spps <- lins <- accs
    if(!quiet){
      cat("Downloading", N, "sequences from GenBank\n")
      progseq <- seq(from = 0, to = N, length.out = 80)
    }
    for (i in 1:nrequest){
      retstart <- (i - 1) * 500
      URL2 <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
                     "efetch.fcgi?",
                     "db=", db,
                     "&WebEnv=", WebEnv,
                     "&query_key=", QueryKey,
                     "&retstart=", retstart,
                     "&retmax=500",
                     "&rettype=gb",
                     "&retmode=xml",
                     if(!is.null(contact)) paste0("&email=", contact) else NULL)
      tmp <- .scanURL(URL2, retmode = "xml")
      tmp2 <- .extractXML(tmp, species, lineages, taxIDs)
      if(is.null(tmp2)) tmp2 <- .extractXML2(tmp, species, lineages, taxIDs)
      accs[[i]] <- tmp2$accs
      seqs[[i]] <- tmp2$seqs
      if(species) spps[[i]] <- tmp2$spps
      if(lineages) lins[[i]] <- tmp2$lins
      if(taxIDs) taxs[[i]] <- tmp2$taxs
      if(!quiet & retstart > 0){
        nlessthan <- sum(progseq <= retstart)
        cat(paste0(rep("=", nlessthan), collapse = ""))
        if(nlessthan > 0) progseq <- progseq[-(1:nlessthan)]
      }
    }
    res <- unlist(seqs, use.names = FALSE)
    accs <- unlist(accs, use.names = FALSE)
    if(taxIDs) taxs <- as.integer(unlist(taxs, use.names = FALSE))
    if(species) spps <- unlist(spps, use.names = FALSE)
    if(lineages) lins <- unlist(lins, use.names = FALSE)
    if(length(res) == 0) stop("No valid sequences to return\n")
    if(bin){
      res <- if(db == "nucleotide") char2dna(res) else char2aa(res)
    }
    if(taxIDs) attr(res, "taxID") <- taxs
    if(species) attr(res, "species") <- spps
    if(lineages) attr(res, "lineage") <- lins
    names(res) <- accs
  }else{
    if(!quiet) {
      cat("Downloading", N, "accession numbers from GenBank\n")
      progseq <- seq(from = 0, to = N, length.out = 80)
    }
    for (i in 1:nrequest){
      retstart <- (i - 1) * 500
      URL2 <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
                     "db=", db,
                     "&WebEnv=", WebEnv,
                     "&query_key=", QueryKey,
                     "&retstart=", retstart,
                     "&retmax=500",
                     "&rettype=acc",
                     "&retmode=text",
                     "&tool=R",
                     if(!is.null(contact)) paste("&email=", contact, sep = "") else NULL)
      tmp <- .scanURL(URL2, retmode = "text", what = "", sep = "\n", quiet = TRUE)
      accs[[i]] <- gsub("\\.[0123456789]+", "", tmp)
      if(!quiet){
        nlessthan <- sum(progseq <= retstart)
        cat(paste0(rep("=", nlessthan), collapse = ""))
        if(nlessthan > 0) progseq <- progseq[-(1:nlessthan)]
      }
    }
    res <- unlist(accs, use.names = FALSE)
  }
  return(res)
}
################################################################################


