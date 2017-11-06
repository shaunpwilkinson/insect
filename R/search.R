#' Query GenBank and download matching sequences.
#'
#' \code{searchGB} searches the NCBI GenBank database using the Entrez
#'   utilities, and downloads the matching sequences as a
#'   \code{"DNAbin"} object. Note
#'   that an internet connection is required when using this function.
#'   In some cases depending on the user security settings the R session
#'   may need to be opened with administrator privileges.
#'
#' @param query an Entrez query provided as a character string. For help
#'   compiling Entrez queries see
#'   \url{https://www.ncbi.nlm.nih.gov/books/NBK3837/#EntrezHelp.Entrez_Searching_Options}
#'   and \url{https://www.ncbi.nlm.nih.gov/books/NBK49540/}.
#' @param onlyID logical. Should just the GenBank accession numbers be
#'   returned (FALSE) or should the function return the entire sequences
#'   and their metadata (TRUE; default)?
#' @param prompt logical indicating whether to check with the user before
#'   downloading sequences.
#' @param contact an optional character string with the users email address.
#'   This is added to the E-utilities URL and may be used by NCBI to contact
#'   the user if the application causes unintended issues.
#' @param quiet logical indicating whether the progress should be printed
#'   to the console.
#' @return a list of DNA sequences as a \code{DNAbin} object.
#' @details
#'   This function uses the Entrez e-utilities API to search and download
#'   sequences from GenBank.
#'   Occasionally users may encounter an unknown error that is not reproducible
#'   and appears to be related to database records being updated in GenBank.
#'   This can generally be remedied by re-running the function. If problems
#'   persist please feel free to raise an issue on the package bug-reports page at
#'   <http://github.com/shaunpwilkinson/insect/issues>.
#' @author Shaun Wilkinson
#' @references
#'   NCBI Resource Coordinators (2012) Database resources of the National
#'   Center for Biotechnology Information. \emph{Nucleic Acids Research},
#'    \strong{41} (Database issue): D8–D20.
#' @seealso \code{\link{readGB}} and \code{\link[ape]{read.GenBank}} (ape)
#'   for downloading DNA sequences from GenBank using accession numbers.
#' @examples
#'   ## Query the GenBank database for Eukaryote mitochondrial 16S sequences
#'   ## between 100 and 300 base pairs in length and last modified between
#'   ## 1999 and 2000.
#'   \dontrun{
#'     query <- "Eukaryota[ORGN]+AND+16S[TITL]+AND+100:300[SLEN]+AND+1999:2000[MDAT]"
#'     x <- searchGB(query)
#'   }
################################################################################
searchGB <- function(query, onlyID = FALSE, prompt = TRUE,
                     contact = NULL, quiet = FALSE){
  query <- paste0("term=", query, "&", collapse = "")
  URL1 <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
                "esearch.fcgi?",
                "db=nucleotide&",
                query,
                "usehistory=y&",
                "tool=R",
                if(!is.null(contact)) paste("&email=", contact, sep = "") else NULL,
                sep = "")
  X <- scan(file = URL1, what = "", sep = "\n", quiet = TRUE)
  Count <- gsub(".+<Count>([[:digit:]]+)</Count>.+", "\\1", X[3])
  WebEnv <- gsub(".+<WebEnv>(.+)</WebEnv>.+", "\\1", X[3])
  QueryKey <- gsub(".+<QueryKey>(.+)</QueryKey>.+", "\\1", X[3])
  N <- as.numeric(Count)
  if(prompt){
    decision <- readline(paste(N, if(onlyID) "sequence IDs" else "sequences",
                               "will be downloaded, continue? (y/n): "))
    if(decision == "n") stop("aborted") else if(!identical(decision, "y")) stop("invalid")
  }
  nrequest <- N%/%500 + as.logical(N%%500) # 1 if remainder exists, 0 otherwise
  if(onlyID){
    obj <- character(N)
  }else{
    obj <- vector("list", N)
    objnames <- character(N)
    objdef <- character(N)
    objorg <- character(N)
    objtax <- integer(N)
    objlin <- character(N)
    discards <- logical(N)
  }
  counter <- 1
  discat <- ""
  if(!quiet) {
    cat("Downloading", N, "DNA sequences from GenBank\n")
    cat("Showing download progress\n")
    cat(paste0(rep("_", 80), collapse = ""), "\n")
    progseq <- seq(from = 0, to = N, length.out = 80)
  }
  scanURL <- function(x){
    res <- tryCatch(scan(file = x, what = "", sep = "\n", quiet = TRUE),
                    error = function(er) return(NULL),
                    warning = function(wa) return(NULL))
    return(res)
  }
  if(onlyID){
    for (i in 1:nrequest){
      a <- (i - 1) * 500
      b <- 500 * i - 1
      if (i == nrequest) b <- N
      URL2 <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
                    "efetch.fcgi?",
                    "db=nucleotide",
                    "&WebEnv=", WebEnv,
                    "&query_key=", QueryKey,
                    "&retstart=", a,
                    "&retmax=", 500,
                    "&rettype=acc",
                    "&retmode=text",
                    "&tool=R",
                    if(!is.null(contact)) paste("&email=", contact, sep = "") else NULL,
                    sep = "")
      #tmp <- scan(file = URL2, what = "", sep = "\n", quiet = TRUE)
      for(l in 1:10){
        tmp <- scanURL(URL2)
        if(!is.null(tmp)) break else Sys.sleep(5)
      } # this is to stop process aborting if unable to connect
      if(is.null(tmp)) stop("No internet connection")
      tmp <- scan(file = URL2, what = "", sep = "\n", quiet = TRUE)
      tmp <- gsub("\\.[0123456789]+", "", tmp)
      obj[(a + 1):(a + length(tmp))] <- tmp
      counter <- counter + 500
      if(!quiet){
        nlessthan <- sum(progseq <= counter)
        cat(paste0(rep("=", nlessthan), collapse = ""))
        if(nlessthan > 0) progseq <- progseq[-(1:nlessthan)]
      }
    }
  }else{
    for (i in 1:nrequest){
      a <- (i - 1) * 500
      b <- 500 * i - 1
      if (i == nrequest) b <- N
      URL2 <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
                    "efetch.fcgi?",
                    "db=nucleotide",
                    "&WebEnv=", WebEnv,
                    "&query_key=", QueryKey,
                    "&retstart=", a,
                    "&retmax=", 500,
                    "&rettype=gb",
                    "&retmode=text",
                    sep = "")
      for(l in 1:10){
        tmp <- scanURL(URL2)
        if(!is.null(tmp)) break else Sys.sleep(5)
      }
      if(is.null(tmp)) stop("No internet connection")
      #tmp <- scan(file = URL2, what = "", sep = "\n", quiet = TRUE)
      # break long vector into list with 1 element for each sequence
      ends <- which(grepl("^//", tmp))
      starts <- c(1, ends[-length(ends)] + 1)
      # if(length(ends) != (b - a + 1) | length(starts) != (b - a + 1)){
      #   stop("Error in sequence download, check internet connection and try again\n")
      # }
      if(length(starts) != length(ends)) stop("Unknown error occurred, please try again")
      runs <- mapply(":", starts, ends, SIMPLIFY = FALSE)
      tmp2 <- lapply(runs, function(e) tmp[e])
      for(j in seq_along(tmp2)){
        seqok <- TRUE
        deflines <- grep("^DEFINITION", tmp2[[j]])
        if(length(deflines) != 1) seqok <- FALSE
        acclines <- grep("^ACCESSION", tmp2[[j]])
        if(length(acclines) != 1) seqok <- FALSE
        # verlines <- acclines + 1 # can throw error if acclines is char(0)
        orglines <- grep("^  ORGANISM", tmp2[[j]])
        if(length(orglines) != 1) seqok <- FALSE
        reflines <- grep("^REFERENCE", tmp2[[j]])
        if(length(reflines) > 1) reflines <- reflines[1] else if(length(reflines) < 1) seqok <- FALSE
        taxlines <- grep("/db_xref=\"taxon:", tmp2[[j]])
        if(length(taxlines) != 1) seqok <- FALSE
        seqlines <- grep("^ORIGIN", tmp2[[j]])
        if(length(seqlines) == 1) seqlines <- seqlines + 1 else seqok <- FALSE
        endlines <- grep("^//", tmp2[[j]])
        if(length(endlines) == 1) endlines <- endlines - 1 else seqok <- FALSE
        if(seqok){
          seqj <- paste0(tmp2[[j]][seqlines:endlines], collapse = "")
          seqj <- gsub("[ 0123456789]", "", seqj)
          seqj <- strsplit(seqj, split = "")[[1]]
          obj[[counter]] <- unclass(ape::as.DNAbin(seqj))
          # objnames[counter] <- gsub("^ACCESSION *", "", tmp2[[j]][acclines])
          name <- gsub("^VERSION *", "", tmp2[[j]][acclines + 1])
          objnames[counter] <- gsub("\\.[0123456789]+ *", "", name)
          defn <- paste0(tmp2[[j]][deflines:(acclines - 1)], collapse = "")
          defn <- gsub("^DEFINITION *", "", defn)
          defn <- gsub("  +", " ", defn)
          objdef[counter] <- defn
          org <- tmp2[[j]][orglines]
          org <- gsub("  ORGANISM *", "", org)
          objorg[counter] <- org
          tax <- tmp2[[j]][taxlines]
          tax <- as.integer(gsub(".+taxon:([[:digit:]]+).+", "\\1", tax))
          objtax[counter] <- tax
          linseq <- seq(orglines + 1, reflines - 1)
          for(l in linseq){
            if(!(grepl("^ {12}Eukaryota", tmp2[[j]][l]) |
                 grepl("^ {12}Bacteria", tmp2[[j]][l]) |
                 grepl("^ {12}Archaea", tmp2[[j]][l]))){
              linseq <- linseq[-1]
            }else break
          }
          if(length(linseq) > 0){
            # linlines_s <- orglines + 1
            # linlines_e <- linlines_s
            # while(!grepl("\\.", tmp2[[j]][linlines_e])) linlines_e <- linlines_e + 1
            # lineage <- paste0(tmp2[[j]][linlines_s:linlines_e], collapse = "")
            lineage <- paste0(tmp2[[j]][linseq], collapse = "")
            lineage <- gsub("^ +", "", lineage)
            lineage <- gsub("  +", " ", lineage)
            lineage <- gsub("\\(.+\\)", "", lineage)
            lineage <- gsub(" +\\.", "\\.", lineage)
            objlin[counter] <- lineage
          }
        }else{
          discards[counter] <- TRUE
          # if(!quiet) cat(tmp2[[j]][acclines], "is invalid and will be discarded\n")
          discat <- paste0(discat, gsub("ACCESSION   ", "", tmp2[[j]][acclines]), ", ")
        }
        counter <- counter + 1
      }
      if(!quiet){
        nlessthan <- sum(progseq <= counter)
        cat(paste0(rep("=", nlessthan), collapse = ""))
        if(nlessthan > 0) progseq <- progseq[-(1:nlessthan)]
      }
    }
    obj <- obj[!discards]
    names(obj) <- objnames[!discards]
    objlin <- objlin[!discards]
    objorg <- objorg[!discards]
    objlin <- gsub("\\.$", "", objlin)
    # attr(obj, "species") <- objorg[!discards]
    attr(obj, "lineage") <- paste0(objlin, "; ", objorg, ".")
    attr(obj, "definition") <- objdef[!discards]
    attr(obj, "taxon") <- objtax[!discards]
    class(obj) <- "DNAbin"
  }
  if(!quiet){
    if(any(discards)){
      cat("\nThe following sequences are invalid and were discarded: ", gsub(", $", "", discat))
      #for(i in which(discards)) cat(i, " ")
      #cat("\n")
    }
  }
  nolength <- sapply(obj, length) == 0
  noclass <- sapply(obj, class) == "NULL"
  if(any(nolength | noclass)) warning("Some sequences failed to download\n")
  return(obj)
}
################################################################################
#' Read DNA sequences from GenBank.
#'
#' This function downloads DNA sequence data from GenBank via the Entrez API.
#'
#' @param accs character vector or GenBank accession numbers.
#' @param prompt logical indicating whether the user should be prompted with
#'  the total number of sequences to download before commencing.
#' @param contact character string, giving users the option to provide an
#'   emil address to send on the web form.
#' @param quiet logical indicating whether feedback should be printed
#'   to the console.
#' @return a list of DNA sequences as a \code{DNAbin} object.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @seealso \code{\link{searchGB}}
#' @examples
#'   ## TBA
################################################################################
readGB <- function(accs, prompt = FALSE, contact = NULL, quiet = FALSE){
  N <- length(accs)
  if(prompt){
    decision <- readline(paste(N, "sequences will be downloaded, continue? (y/n): "))
    if(decision == "n") stop("aborted") else if(!identical(decision, "y")) stop("invalid")
  }
  nrequest <- N%/%500 + as.logical(N%%500)
  obj <- vector("list", N)
  objnames <- character(N)
  objdef <- character(N)
  objorg <- character(N)
  objtax <- integer(N)
  objlin <- character(N)
  discards <- logical(N)
  counter <- 1
  if(!quiet) {
    cat(paste0(rep("_", 80), collapse = ""), "\n")
    progseq <- seq(from = 0, to = N, length.out = 80)
  }
  for (i in 1:nrequest){
    a <- (i - 1) * 500 + 1
    b <- 500 * i - 1
    if (i == nrequest) b <- N
    URL2 <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
                  "efetch.fcgi?",
                  "db=nucleotide",
                  "&id=",
                  paste(accs[a:b], collapse = ","),
                  "&rettype=gb",
                  "&retmode=text",
                  sep = "")
    tmp <- scan(file = URL2, what = "", sep = "\n", quiet = TRUE)
    # break long vector into list with 1 element for each sequence
    ends <- which(grepl("^//", tmp))
    starts <- c(1, ends[-length(ends)] + 1)
    runs <- mapply(":", starts, ends, SIMPLIFY = FALSE)
    tmp2 <- lapply(runs, function(e) tmp[e])
    for(j in seq_along(tmp2)){
      seqok <- TRUE
      deflines <- grep("^DEFINITION", tmp2[[j]])
      if(length(deflines) != 1) seqok <- FALSE
      acclines <- grep("^ACCESSION", tmp2[[j]])
      if(length(acclines) != 1) seqok <- FALSE
      # verlines <- acclines + 1 # can throw error if acclines is char(0)
      orglines <- grep("^  ORGANISM", tmp2[[j]])
      if(length(orglines) != 1) seqok <- FALSE
      reflines <- grep("^REFERENCE", tmp2[[j]])
      if(length(reflines) > 1) reflines <- reflines[1] else if(length(reflines) < 1) seqok <- FALSE
      taxlines <- grep("/db_xref=\"taxon:", tmp2[[j]])
      if(length(taxlines) != 1) seqok <- FALSE
      seqlines <- grep("^ORIGIN", tmp2[[j]])
      if(length(seqlines) == 1) seqlines <- seqlines + 1 else seqok <- FALSE
      endlines <- grep("^//", tmp2[[j]])
      if(length(endlines) == 1) endlines <- endlines - 1 else seqok <- FALSE
      if(seqok){
        seqj <- paste0(tmp2[[j]][seqlines:endlines], collapse = "")
        seqj <- gsub("[ 0123456789]", "", seqj)
        seqj <- strsplit(seqj, split = "")[[1]]
        obj[[counter]] <- unclass(ape::as.DNAbin(seqj))
        name <- gsub("^VERSION *", "", tmp2[[j]][acclines + 1])
        objnames[counter] <- gsub("\\.[0123456789]+ *", "", name)
        defn <- paste0(tmp2[[j]][deflines:(acclines - 1)], collapse = "")
        defn <- gsub("^DEFINITION *", "", defn)
        defn <- gsub("  +", " ", defn)
        objdef[counter] <- defn
        org <- tmp2[[j]][orglines]
        org <- gsub("  ORGANISM *", "", org)
        objorg[counter] <- org
        tax <- tmp2[[j]][taxlines]
        tax <- as.integer(gsub(".+taxon:([[:digit:]]+).+", "\\1", tax))
        objtax[counter] <- tax
        linseq <- seq(orglines + 1, reflines - 1)
        for(l in linseq){
          if(!(grepl("^ {12}Eukaryota", tmp2[[j]][l]) |
               grepl("^ {12}Bacteria", tmp2[[j]][l]) |
               grepl("^ {12}Archaea", tmp2[[j]][l]))){
            linseq <- linseq[-1]
          }else break
        }
        if(length(linseq) > 0){
          lineage <- paste0(tmp2[[j]][linseq], collapse = "")
          lineage <- gsub("^ +", "", lineage)
          lineage <- gsub("  +", " ", lineage)
          lineage <- gsub("\\(.+\\)", "", lineage)
          lineage <- gsub(" +\\.", "\\.", lineage)
          objlin[counter] <- lineage
        }
      }else{
        discards[counter] <- TRUE
        if(!quiet) cat("Sequence", counter, "is invalid and will be discarded\n")
      }
      counter <- counter + 1
    }
    if(!quiet){
      nlessthan <- sum(progseq <= counter)
      cat(paste0(rep("=", nlessthan), collapse = ""))
      if(nlessthan > 0) progseq <- progseq[-(1:nlessthan)]
    }
  }
  obj <- obj[!discards]
  names(obj) <- objnames[!discards]
  objlin <- objlin[!discards]
  objorg <- objorg[!discards]
  objlin <- gsub("\\.$", "", objlin)
  # attr(obj, "species") <- objorg[!discards]
  attr(obj, "lineage") <- paste0(objlin, "; ", objorg, ".")
  attr(obj, "definition") <- objdef[!discards]
  attr(obj, "taxon") <- objtax[!discards]
  class(obj) <- "DNAbin"
  if(!quiet){
    if(any(discards)){
      cat("\n", sum(discards), "sequences are invalid and were discarded:\n ")
      for(i in which(discards)) cat(i, " ")
      cat("\n")
    }
  }
  nolength <- sapply(obj, length) == 0
  noclass <- sapply(obj, class) == "NULL"
  if(any(nolength | noclass)) warning("Some sequences failed to download\n")
  return(obj)
}
################################################################################
#' Query the BOLD database and download matching sequences.
#'
#' This function searches the barcode of life (BOLD) database using the public
#'   API, and downloads the matching sequences as a \code{"DNAbin"} object.
#'   Note that an internet connection is required when using this function.
#'   In some cases depending on the user security settings the R session
#'   may need to be opened with administrator privileges.
#' @param taxon a recognized scientific taxon name.
#' @param GB logical indicating whether sequences also present in GenBank should
#'   be included in the output object.
#' @param markers character string or vector giving the genetic marker(s)
#'   to limit the result to. Can be any or all of: "COI-5P", "COI-3P", "28S",
#'   "16S", "COII", "CYTB", "atp6", "COXIII", "18S", "ITS", "ITS2", "ITS1", "5.8S".
#' @return a list of DNA sequences as a \code{DNAbin} object.
#' @details
#'   This function calls the BOLD API
#'   (\url{http://www.boldsystems.org/index.php/resources/api?type=webservices})
#'   to search and download DNA barcode sequences from the BOLD database
#'   (Ratnasingham & Hebert 2007).
#' @author Shaun Wilkinson
#' @references
#'   Ratnasingham S & Hebert PDN (2007) BOLD: the barcode of life data
#'   system (www.barcodinglife.org). \emph{Molecular Ecology Notes},
#'   \strong{7}, 355–364.
#' @seealso \code{\link{searchGB}} for downloading DNA sequences from GenBank.
#' @examples
#'   ## Query the BOLD database for Saccharomycetes ITS2 sequences
#'   \dontrun{
#'     x <- searchBOLD("Saccharomycetes", markers = "ITS2")
#'   }
################################################################################
searchBOLD <- function(taxon, GB = TRUE, markers = NULL){
  if(is.null(markers)){
    markers <- c("COI-5P", "COI-3P", "28S", "16S", "COII", "CYTB", "atp6",
                "COXIII", "18S", "ITS", "ITS2", "ITS1", "5.8S")
  }
  mk <- paste0(markers, collapse = "|")
  bq <- "http://boldsystems.org/index.php/API_Public/combined?taxon="
  eq <- "&format=tsv"
  URL <- paste0(bq, taxon, "&marker=", mk, eq)
  # X <- scan(file = URL, what = "", quote = "", sep = "\t", quiet = TRUE)
  X <- httr::GET(URL)
  X <- paste0(rawToChar(X$content, multiple = TRUE), collapse = "")
  Xm <- utils::read.delim(text = X, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # ncols <- grep("marker_codes", X)[1]
  # if(is.na(ncols)) stop("Error 1")
  # if(length(X) %% ncols != 0) stop("Error 2")
  # Xm <- matrix(X, ncol = ncols, byrow = TRUE)
  #colnames(Xm) <- Xm[1, ]
  #Xm <- Xm[-1, ]
  if(nrow(Xm) == 0) return(NULL)
  # spc <- match("species_name", colnames(Xm))
  Xm <- Xm[grepl("[[:alpha:]]", Xm$species_name), ]
  if(nrow(Xm) == 0) return(NULL)
  if(!GB) Xm <- Xm[!grepl("[[:alnum:]]", Xm$genbank_accession), ]
    # gbc <- match("genbank_accession", colnames(Xm))
  if(nrow(Xm) == 0) return(NULL)
  #mac <- match("markercode", colnames(Xm))
  Xm <- Xm[Xm$markercode %in% markers, ] # poss bug in API returns additional markers
  if(nrow(Xm) == 0) return(NULL)
  # dnac <- match("nucleotides", colnames(Xm))
  seqs <- gsub("-", "", Xm$nucleotides)
  seqs <- strsplit(seqs, split = "")
  out <- ape::as.DNAbin(seqs)
  rec <- match("recordID", colnames(Xm))
  recs <- Xm[, rec]
  names(out) <- paste0("BOLD", Xm$recordID)
  lins <- paste(Xm$phylum_name, Xm$class_name, Xm$order_name, Xm$family_name,
                Xm$genus_name, Xm$species_name, sep = "; ")
  attr(out, "lineage") <- paste0(lins, ".")
  return(out)
}
################################################################################
