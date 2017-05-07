#' Query GenBank and download matching sequences.
#'
#' \code{searchGB} searches the NCBI GenBank database using the Entrez
#'   utilities, and downloads the matching sequences either as a
#'   \code{"DNAbin"} object or a list of character vectors. Note
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
#' @param DNA logical indicating whether the returned value should be of
#'   class \code{DNAbin} (TRUE; defult). If \code{FALSE} a list of character
#'   vectors is returned. \code{DNA = TRUE} is recommended for memory and speed
#'   efficiency and compatibility with other functions in the \code{insect}
#'   package.
#' @param prompt logical indicating whether to check with the user before
#'   downloading sequences.
#' @param contact an optional character string with the users email address.
#'   This is added to the E-utilities URL and may be used by NCBI to contact
#'   the user if the application causes unintended issues.
#' @param quiet logical indicating whether the progress should be printed
#'   to the console.
#' @return a list of DNA sequences, either in \code{DNAbin} format or as
#'   character vectors if \code{DNA} is set to \code{FALSE}.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
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
searchGB <- function(query, onlyID = FALSE, DNA = TRUE, prompt = TRUE,
                     contact = NULL, quiet = FALSE){
  # input a string
  # output a DNAbin list object or character equivalent
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
  nrequest <- N%/%500 + as.logical(N%%500)
  if(onlyID){
    obj <- character(N)
  }else{
    obj <- vector("list", N)
    objnames <- character(N)
    objdef <- character(N)
    objorg <- character(N)
    objlineage <- character(N)
    discards <- logical(N)
  }
  counter <- 1
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
      runs <- mapply(":", starts, ends, SIMPLIFY = FALSE)
      tmp2 <- lapply(runs, function(e) tmp[e])
      for(j in seq_along(tmp2)){
        acclines <- which(grepl("^ACCESSION", tmp2[[j]]))
        deflines <- which(grepl("^DEFINITION", tmp2[[j]]))
        orglines <- which(grepl("^  ORGANISM", tmp2[[j]]))
        seqlines <- which(grepl("^ORIGIN", tmp2[[j]])) + 1
        endlines <- which(grepl("^//", tmp2[[j]])) - 1
        if(length(acclines) == 1 & length(deflines) == 1 & length(orglines) == 1 &
           length(seqlines) == 1 & length(endlines) == 1){
          seqj <- paste0(tmp2[[j]][seqlines:endlines], collapse = "")
          seqj <- gsub("[ 0123456789]", "", seqj)
          seqj <- strsplit(seqj, split = "")[[1]]
          obj[[counter]] <- if(DNA) unclass(ape::as.DNAbin(seqj)) else seqj
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
          if(grepl(";", tmp2[[j]][orglines + 1])){
            linlines_s <- orglines + 1
            linlines_e <- linlines_s
            while(!grepl("\\.", tmp2[[j]][linlines_e])) linlines_e <- linlines_e + 1
            lineage <- paste0(tmp2[[j]][linlines_s:linlines_e], collapse = "")
            lineage <- gsub("^ +", "", lineage)
            lineage <- gsub("  +", " ", lineage)
            lineage <- gsub("\\(.+\\)", "", lineage)
            lineage <- gsub(" +\\.", "\\.", lineage)
            objlineage[counter] <- lineage
          }
        }else{
          discards[counter] <- TRUE
          #if(!quiet) cat("Sequence", counter, "is invalid and will be discarded\n")
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
    attr(obj, "species") <- objorg[!discards]
    attr(obj, "definition") <- objdef[!discards]
    attr(obj, "lineage") <- objlineage[!discards]
    if(DNA) class(obj) <- "DNAbin"
  }
  if(!quiet){
    if(any(discards)){
      cat("The following sequences are invalid and were discarded:\n")
      for(i in which(discards)){
        cat(i, "\n")
      }
    }
  }
  return(obj)
}
################################################################################
#' Read DNA sequences from GenBank.
#'
#' This function downloads DNA sequence data from GenBank via the Entrez API.
#'
#' @param accs character vector or GenBank accession numbers.
#' @param DNA logical indicating whether the sequences should be returned
#'   as a \code{DNAbin} object.
#' @param prompt logical indicating whether the user should be prompted with
#'  the total number of sequences to download before commencing.
#' @param contact character string, giving users the option to provide an
#'   emil address to send on the web form.
#' @param quiet logical indicating whether feedback should be printed
#'   to the console.
#' @return an object of class \code{DNAbin} if \code{DNA = TRUE}, or a
#'   list of character vectors otherwise.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @seealso \code{\link{searchGB}}
#' @examples
#'   ## TBA
################################################################################
readGB <- function(accs, DNA = TRUE, prompt = FALSE,
                   contact = NULL, quiet = FALSE){
  # input a character vector of access IDs
  # output a DNAbin list object or character equivalent
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
  objlineage <- character(N)
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
      acclines <- which(grepl("^ACCESSION", tmp2[[j]]))
      deflines <- which(grepl("^DEFINITION", tmp2[[j]]))
      orglines <- which(grepl("^  ORGANISM", tmp2[[j]]))
      seqlines <- which(grepl("^ORIGIN", tmp2[[j]])) + 1
      endlines <- which(grepl("^//", tmp2[[j]])) - 1
      if(length(acclines) == 1 & length(deflines) == 1 & length(orglines) == 1 &
         length(seqlines) == 1 & length(endlines) == 1){
        seqj <- paste0(tmp2[[j]][seqlines:endlines], collapse = "")
        seqj <- gsub("[ 0123456789]", "", seqj)
        seqj <- strsplit(seqj, split = "")[[1]]
        obj[[counter]] <- if(DNA) unclass(ape::as.DNAbin(seqj)) else seqj
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
        if(grepl(";", tmp2[[j]][orglines + 1])){
          linlines_s <- orglines + 1
          linlines_e <- linlines_s
          while(!grepl("\\.", tmp2[[j]][linlines_e])) linlines_e <- linlines_e + 1
          lineage <- paste0(tmp2[[j]][linlines_s:linlines_e], collapse = "")
          lineage <- gsub("^ +", "", lineage)
          lineage <- gsub("  +", " ", lineage)
          lineage <- gsub("\\(.+\\)", "", lineage)
          lineage <- gsub(" +\\.", "\\.", lineage)
          objlineage[counter] <- lineage
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
  attr(obj, "species") <- objorg[!discards]
  attr(obj, "definition") <- objdef[!discards]
  attr(obj, "lineage") <- objlineage[!discards]
  if(DNA) class(obj) <- "DNAbin"
  return(obj)
}
################################################################################
