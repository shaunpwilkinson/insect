#' Read and filter FASTQ files.
#'
#' \code{readFASTQ} is a text parser/filter that reads file in the FASTQ
#'   format into R, during which it optionally applies various sequence
#'   editing procedures and filters. If a valid sample sheet is provided
#'   (\code{sheet} is not NULL), any sequence that does not feature a motif
#'   matching one of the index-primer-primer-index combinations from the
#'   sheet are discarded during the parse.
#'
#' @param path a character string giving the path to the FASTQ text file
#'   to be parsed.
#' @param sheet either a character string giving the path to a valid sample
#'   sheet in csv format, or a similarly valid data frame object already
#'   stored in the R environment. The sheet should have one row for each
#'   unique index-primer-primer-index combination, and 13 cloumns labeled
#'   in the following order: 'Number', (sample number, these will be
#'   appended to the end of the names of the output sequences after the
#'   final colon), 'PCR #' (the PCR number), 'Extract' (a description of
#'   the sample), 'Sample' (the environment from where the sample was
#'   extracted), 'Target' (a name given to the primer pair), 'F Tag',
#'   'R Tag', (unique names given to the forward and reverse indexing tags,
#'   respectively), F Sequence', 'R Sequence' (the sequences of the forward
#'   and reverse index tags, respectively), 'Forward Primer', 'Reverse Primer'
#'   (the sequences of the forward and reverse primers, respectively), and
#'   'P7 Sequence'. All data in the tag, primer and P7 Sequence fields
#'   should be in upper case and should not contain gaps. If
#'   \code{sheet = NULL}, all sequences will be parsed without matching
#'   or trimming index tags and primers, and sequences will be named
#'   according to their names in the input fastq text file (the lines
#'   beginning with "@").
#' @param filter logical indicating whether sequences should undergo length,
#'   quality and ambiguity filtering during the parse. Defaults to TRUE.
#' @param minlength integer giving the minimum acceptable length (the
#'   number of base pairs after the removal of index tags and primers
#'   if applicable) for sequences to be retained. Defaults to 50.
#' @param minqual numeric giving the minimum average Phred quality
#'   score (after trimming index tags and primers) for sequences to
#'   be retained. Defaults to 30.
#' @param maxambigs integer giving the maximum acceptable number of
#'   ambiguous base calls (after trimming index tags and primers) for
#'   sequences to be retained. Defaults to 0.
#' @param mincount integer giving the minimum number of occurrences of a
#'   sequence in the dataset for it to be retained. Defaults to 2
#'   (singletons are discarded).
#' @param DNA logical indicating whether the returned object should be of
#'   class \code{"DNAbin"} as output by the \code{\link[ape]{as.DNAbin}}
#'   function in the \code{\link[ape]{ape}} package (defaults to TRUE;
#'   recommended for memory and speed efficiency). If set to FALSE, the
#'   function returns a list of character vectors containing residues from the
#'   set {ACGTN} with an associated character vector of "quality"
#'   attributes in the same ASCII coding scheme used in the FASTQ file format.
#' @param nlines integer, the maximum number of text lines to be read and
#'   analysed in one batch. This may be reduced if importing large files and
#'   memory is limited. Typically a ten-million-line text file of ~500bp
#'   sequences with associated metadata (quality, run number, cluster
#'   coordinates, etc) will require around 2 GB of memory for the initial
#'   text import, and be reduced to around a quarter or this after the
#'   removal of index-mismatching sequences, short reads and singletons,
#'   and once the sequences and quality scores are compressed to
#'   raw bytes (assuming \code{DNA} is set to \code{TRUE}). Thus reducing
#'   \code{nlines} can prevent memory limits from being exceeded; however
#'   this can increase the computation time.
#' @param quiet logical indicating whether the progress of the parse
#'   operation should be printed to the console.
#' @param ... further arguments to be passed to \code{\link{scan}}.
#' @return Returns a list of DNA sequences, either in \code{DNAbin} format
#'   or as character vectors if \code{DNA} is set to \code{FALSE}. In either
#'   case the sequences in the returned list each have a "quality" attribute,
#'   which is a vector of quality scores between 0 and 93 either in raw bytes
#'   (if DNA = TRUE) or in the same ASCII coding scheme as used in the FASTQ
#'   file format.
#' @details
#'   \subsection{Compatibility:}{The FASTQ convention remains somewhat
#'   ambiguous with several slightly different interpretations appearing
#'   in the literature. For now, this function only supports the Illumina
#'   convention for FASTQ files, where each sequence and its associated
#'   metadata occupies four line of the text file as follows : (1) the
#'   run and cluster metadata preceeded by an @ symbol; (2) the sequence
#'   itself in capitals without spaces and containing only residues from
#'   the set {ACGTN}; (3) a single "+" symbol; and (4) the Phred quality
#'   scores from 0 to 93 represented as ASCII symbols. For more information
#'   on this convention see the Illumina help page
#'   \href{https://help.basespace.illumina.com/articles/descriptive/fastq-files/}{here}
#'   .}
#'   \subsection{Speed and Memory Requirements:}{This function can
#'   take a while to process larger FASTQ files, particularly
#'   if index matching and filtering is required (\code{sheet} is
#'   not \code{NULL} and \code{filter = TRUE}). As an example, a FASTQ text
#'   file with four-million lines (one million 500bp sequences with
#'   associated metadata), with a sample sheet specifying
#'   48 unique index-primer-primer-index combinations, and using the
#'   default fitering parameters took around 10 minutes on a single core
#'   of a Lenovo T530 Thinkpad laptop (Intel i7 2.6GHz w 16GB RAM running
#'   64 bit Ubuntu 16.04). The operation required ~ 2.5 GB of free memory and
#'   the output \code{"DNAbin"} object containing ~ 400,000 sequences and
#'   asociated quality scores occupied ~250 Mb memory.}
#' @author Shaun Wilkinson
#' @references
#'   Bokulich NA, Subramanian S, Faith JJ, Gevers D, Gordon JI, Knight R,
#'   Mills DA, Caporaso JG (2013) Quality-filtering vastly improves diversity
#'   estimates from Illumina amplicon sequencing.
#'   \emph{Nat Methods}, \strong{1}, 57-59.
#'
#'   Illumina help page:
#'   \url{https://help.basespace.illumina.com/articles/descriptive/fastq-files/}
#'
#' @seealso \code{\link{writeFASTQ}} for writing DNAbin objects to text
#'   in the FASTQ format, and \code{\link[ape]{read.dna}} in the
#'   \code{\link[ape]{ape}} package for reading DNA in FASTA and other formats
#'   into R.
#' @examples
#'   \dontrun{
#'     ##TBA
#'   }
################################################################################
readFASTQ <- function(path, sheet = NULL, filter = TRUE, minlength = 50,
                      minqual = 30, maxambigs = 0, mincount = 2,
                      DNA = TRUE, nlines = 1E07, quiet = FALSE, ...){
  maxlen <- as.integer(file.size(path) * 0.002)
  res <- vector(mode = "list", length = maxlen)
  resnames <- vector(mode = "character", length = maxlen)
  tnseq <- 0 # total number of seqs
  read1 <- function(s, qual, pattern, AL, DNA){
    tmp <- s
    if(!(nchar(s) == nchar(qual))) stop("lengths of sequence and quality scores don't match")
    if(!grepl(pattern, s)) stop(qual)
    s <- gsub(pattern, "\\1", s)[[1]]
    s <- if(DNA) .char2dna(s) else strsplit(s, "")[[1]]
    #p2 <- if(AL == 0) "(.*)" else paste0(".{", AL, "}(.{", length(s), "}).*")
    #qual <- gsub(p2, "\\1", qual)[[1]]
    qual <- if(DNA) .char2qual(qual) else strsplit(qual, "")[[1]]
    if(AL > 0) qual <- qual[-(1:AL)]
    qual <- qual[1:length(s)]
    return(structure(s, quality = qual))
  }
  if(is.null(sheet)){
    if(!quiet) cat("No sample sheet supplied, skipping index & primer trim\n")
    if(!quiet) cat("Reading FASTQ file\n")
  }else{
    if(!quiet) cat("Checking format of sample sheet\n")
    if(length(sheet) == 1 & mode(sheet) == "character"){
      sheet <- read.csv(sheet, header = TRUE, stringsAsFactors = FALSE)
    }
    if(!ncol(sheet) == 13) stop("Sample sheet is invalid, please see
                                ?readFASTQ for instructions on correct
                                formatting")
    facs <- sapply(sheet, is.factor)
    sheet[facs] <- lapply(sheet[facs], as.character)
    counts <- integer(nrow(sheet))
    names(counts) <- sheet[, 1]
    if(!quiet) cat("Reading FASTQ file\n")
    if(!quiet) cat("Matching and trimming indices and primers\n")
  }
  a <- 0
  b <- 1
  repeat{
    x <- scan(file = path, what = "", sep = "\n", skip = a, nlines = nlines,
              quiet = TRUE, ... = ...)
    # For 747555 sequences, created 2990220 element vector at 616.6Mb. Took ~ 1 min
    if(identical(x, character(0))) break
    stopifnot(length(x) %% 4 == 0)
    nseq <- length(x)/4
    tnseq <- tnseq + nseq
    # cat("nseq =", nseq, "\n")
    # cat("tnseq =", tnseq, "\n")
    # cat("a =", a, "\n")
    #if(!quiet) cat(nseq, "sequences detected\n")
    x <- x[seq_along(x) %% 4 != 3] # gets rid of +'s
    # Instrument <- gsub("@([[:alnum:]]+):.+", "\\1", x[1])
    # RunID <- gsub("@[[:alnum:]]+:([[:alnum:]]+):.+", "\\1", x[1])
    # FlowID <- gsub("@[[:alnum:]]+:[[:alnum:]]+:([^:]+):.+", "\\1", x[1])
    if(is.null(sheet)){
      seqalongx <- seq_along(x)
      tmp <- mapply(read1, x[seqalongx %% 3 == 2], x[seqalongx %% 3 == 0],
                    "(.*)", 0, DNA, SIMPLIFY = FALSE, USE.NAMES = FALSE)
      #names(tmp) <- x[seqalongx %% 3 == 1]
      nseqi <- nseq
      res[b:(b + nseqi - 1)] <- tmp
      resnames[b:(b + nseqi - 1)] <- x[seqalongx %% 3 == 1]
      # b <- b + nseqi
      #res <- c(res, tmp)
    }else{
      for(i in 1:nrow(sheet)){
        ##TODO should move all this out of the loop to increase speed
        ftag <- strsplit(sheet[i, 8], "")[[1]]
        ftag <- ftag[ftag != " "]
        fprim <- strsplit(sheet[i, 10], "")[[1]]
        fprim <- fprim[fprim != " "]
        AL <- length(ftag) + length(fprim) # total length of index
        index1 <- paste0("^", .disambiguate(sheet[i, 8]), .disambiguate(sheet[i, 10]))
        index2 <- paste0(.disambiguate(.rc(sheet[i, 11])), .rc(sheet[i, 9]))
        pattern <- paste0(index1, "[ACGTN]+", index2, "[ACGTN]*")
        whichquals <- grep(pattern, x[seq_along(x) %% 3 == 2]) * 3
        nseqi <- length(whichquals)
        if(nseqi > 0){
          whichseqs <- whichquals - 1
          whichnames <- whichseqs - 1
          # if(!quiet) cat("Extracted", length(whichseqs), "sequences from",
          #                sheet[i, 3], "with target", sheet[i, 5], "(tags",
          #                sheet[i, 6], "&", sheet[i, 7], ")\n")
          counts[i] <- counts[i] + nseqi
          # AL <- length(strsplit(index1, split = "")[[1]]) - 1
          pattern <- paste0(index1, "([ACGTN]+)", index2, "[ACGTN]*")
          tmp <- mapply(read1, x[whichseqs], x[whichquals], pattern, AL, DNA,
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
          # nameprefix <- paste(sheet[i, 5], sheet[i, 6], sheet[i, 7], sep = ":")
          # #Target:Ftag:Rtag
          # namesuffix <- gsub(".+(:[^:]:[^:]+:[^:]+:[^:]+) .+", "\\1", x[whichnames])
          # #Lane:TileNo:Xcoord:Ycoord
          # names(tmp) <- paste0(nameprefix, namesuffix)
          # names(tmp) <- gsub(":[0123456789 ]+$", paste0(":", sheet[i, 1]), x[whichnames])
          res[b:(b + nseqi - 1)] <- tmp
          resnames[b:(b + nseqi - 1)] <- gsub(":[0123456789 ]+$",
                                              paste0(":", sheet[i, 1]), x[whichnames])
          # b <- b + nseqi
          #res <- c(res, tmp)
          x <- x[-(c(whichnames, whichseqs, whichquals))]
        }else{
          # if(!quiet) cat("Extracted 0 sequences from",
          #                sheet[i, 3], "with target", sheet[i, 5], "(tags",
          #                sheet[i, 6], "and", sheet[i, 7], ")\n")
        }
      }
    }
    ## moved if() stmt from here Sep 4 2017
    a <- a + nlines
    b <- b + nseqi
    if(nseq < nlines/4) break
  }
  if(length(res) > 0){
    res <- res[1:b]
    names(res) <- resnames[1:b]
    if(!quiet){
      if(!is.null(sheet)){
        cat("Retained", b - 1, "out of", tnseq,
            "sequences after index/primer matching\n")
        cat("Distribution of retained sequences across samples:\n")
        for(i in 1:nrow(sheet)) cat("Sample", sheet[i, 1], "=", counts[i], "sequences\n")
      }else{
        cat("Successfully imported", b, "sequences from file\n")
      }
    }
  }else return(NULL)
  if(filter){
    if(!is.null(minlength)){
      if(!quiet) cat("Filtering sequences by length\n")
      keeps <- sapply(res, length) >= minlength
      res <- res[keeps]
      if(!quiet) cat(length(res), "sequences retained after applying length filter\n")
    }
    if(!is.null(minqual)){
      if(!quiet) cat("Filtering sequences by quality\n")
      keeps <- sapply(res, function(s) mean(as.integer(attr(s, "quality"))) >= minqual)
      res <- res[keeps]
      if(!quiet) cat(length(res), "sequences retained after applying quality filter\n")
    }
    if(!is.null(maxambigs)){
      if(!quiet) cat("Filtering ambiguous sequences\n")
      keeps <- sapply(res, function(s) sum(s == as.raw(240)) <= maxambigs)
      res <- res[keeps]
      if(!quiet) cat(length(res), "sequences retained after applying ambiguity filter\n")
    }
    if(!is.null(mincount)){
      if(!quiet) cat("Filtering sequences by abundance\n")
      hashes <- sapply(res, function(x) paste(openssl::md5(as.vector(x))))
      hashtab <- table(as.factor(hashes))
      keeps <- hashes %in% names(hashtab)[hashtab >= mincount]
      res <- res[keeps]
      if(!quiet) cat(length(res), "sequences retained after removing low-abundance reads\n")
    }
  }
  if(DNA){
    if(!quiet) cat("Creating 'DNAbin' object\n")
    class(res) <- "DNAbin"
  }
  # attr(res, "Instrument") <- Instrument
  # attr(res, "RunID") <- RunID
  # attr(res, "FlowID") <- FlowID
  # 357187 seqs retained
  if(!quiet) cat("Done\n")
  return(res)
}

