# Internal 'insect' functions

.dna2char <- function(z){
  cbytes <- as.raw(c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77, 66, 86, 72, 68, 78, 45, 63))
  indices <- c(136, 24, 72, 40, 96, 144, 192, 48, 80 ,160, 112, 224, 176, 208, 240, 4, 2)
  vec <- raw(240)
  vec[indices] <- cbytes
  if(is.list(z)){
    sapply(z, function(s) rawToChar(vec[as.integer(s)]))
  }else{
    rawToChar(vec[as.integer(z)])
  }
}


.char2dna <- function(z){
  dbytes <- as.raw(c(136, 24, 72, 40, 96, 144, 192, 48, 80 ,160, 112, 224, 176, 208, 240, 4, 2))
  indices <- c(65, 84, 71, 67, 83, 87, 82, 89, 75, 77, 66, 86, 72, 68, 78, 45, 63) # max 89
  vec <- raw(89)
  vec[indices] <- dbytes
  s2d1 <- function(s) vec[as.integer(charToRaw(s))]
  res <- if(length(z) > 1) lapply(z, s2d1) else s2d1(z)
  class(res) <- "DNAbin"
  return(res)
}


.qual2char <- function(x){
  qbytes <- as.raw(0:93)
  qchars <- strsplit(paste0("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP",
                            "QRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"),
                     split = "")[[1]]
  return(paste0(qchars[match(x, qbytes)], collapse = ""))
}


.char2qual <- function(x){
  qbytes <- as.raw(0:93)
  qchars <- strsplit(paste0("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP",
                            "QRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"),
                     split = "")[[1]]
  return(qbytes[match(strsplit(x, split = "")[[1]], qchars)])
}


.rc <- function(s){
  # a string of A, C, G, Ts etc
  s <- strsplit(s, split = "")[[1]]
  s <- rev(s)
  dchars <- strsplit("ACGTMRWSYKVHDBN", split = "")[[1]]
  comps <- strsplit("TGCAKYWSRMBDHVN", split = "")[[1]]
  s <- s[s %in% dchars] # remove spaces etc
  s <- dchars[match(s, comps)]
  s <- paste0(s, collapse = "")
  return(s)
}

.disambiguate <- function(s){
  if(is.null(s)) return(NULL)
  s <- gsub("[^ACGTMRWSYKVHDBNI]", "", s)
  if(grepl("[MRWSYKVHDBNI]", s)){
    s <- gsub("M", "[AC]", s)
    s <- gsub("R", "[AG]", s)
    s <- gsub("W", "[AT]", s)
    s <- gsub("S", "[CG]", s)
    s <- gsub("Y", "[CT]", s)
    s <- gsub("K", "[GT]", s)
    s <- gsub("V", "[ACG]", s)
    s <- gsub("H", "[ACT]", s)
    s <- gsub("D", "[AGT]", s)
    s <- gsub("B", "[CGT]", s)
    s <- gsub("N|I", "[ACGT]", s)
  }
  return(s)
}


## tree is a "dendrogram" object (can be a node)
## x is a DNAbin object - should not contain duplicates, must have lineage attrs
.forkr <- function(tree, x, lineages, refine = "Viterbi", nstart = 10,
                   iterations = 50, minK = 2, maxK = 2,
                   minscore = 0.9, probs = 0.05, retry = TRUE, resize = TRUE,
                   maxsize = NULL, kmers = NULL,
                   seqweights = "Gerstein", cores = 1, quiet = FALSE, ...){
  tree <- fork(tree, x, lineages, refine = refine, nstart = nstart,
               iterations = iterations, minK = minK,
               maxK = maxK, minscore = minscore, probs = probs,
               retry = retry, resize = resize, maxsize = maxsize,
               kmers = kmers, seqweights = seqweights, cores = cores,
               quiet = quiet, ... = ...)
  if(is.list(tree)) tree[] <- lapply(tree, .forkr, x = x, lineages = lineages,
                                     refine = refine, nstart = nstart,
                                     iterations = iterations, minK = minK, maxK = maxK,
                                     probs = probs, retry = retry, resize = resize,
                                     maxsize = maxsize, minscore = minscore,
                                     kmers = kmers, seqweights = seqweights,
                                     cores = cores, quiet = quiet, ... = ...)
  return(tree)
}



.reindex <- function(tree){
  #indices <- character(length(attr(tree, "sequences")))
  fun <- function(node){
    if(is.leaf(node)){
      newnode <- attr(node, "sequences")
      names(newnode) <- rep(attr(node, "clade"), length(newnode))
      node <- newnode
    }
    return(node)
  }
  tmp <- dendrapply(tree, fun)
  indices <- unlist(tmp, use.names = TRUE)
  indices <- sort(indices)
  return(names(indices))
}


.digest <- function(x, simplify = TRUE){
  digest1 <- function(s){
    if(mode(s) != "raw"){
      if(mode(s) == "character"){
        s <- sapply(s, charToRaw)
      }else if(mode(s) == "integer"){
        s <- sapply(s, as.raw)
      }else if(mode(s) == "numeric"){
        s <- unlist(sapply(s, function(a) charToRaw(paste(round(a, 4)))), use.names = FALSE)
        #stop("Can't digest numeric vectors")
      }
    }
    return(paste(openssl::md5(as.vector(s))))
  }
  if(is.list(x)){
    if(simplify) return(sapply(x, digest1)) else return(lapply(x, digest1))
  }else{
    return(digest1(x))
  }
}

.ancestor <- function(lineages){
  # input and output both semicolon-delimited character string(s)
  if(all(lineages == lineages[1])) return(lineages[1])
  # lineages <- gsub("\\.$", "", lineages)
  splitfun <- function(s) strsplit(s, split = "; ")[[1]]
  linvecs <- lapply(lineages, splitfun)
  guide <- linvecs[[which.min(sapply(linvecs, length))]]
  a <- 0
  for(l in guide) if(all(sapply(linvecs, function(e) l %in% e))) a <- a + 1
  guide <- if(a > 0) guide[1:a] else character(0)
  lineage <- paste(guide, collapse = "; ")
  # lineage <- paste0(lineage, ".")
  return(lineage)
}

.ldistance <- function(x, y){
  # x and y are lineage strings (semicolon delim)
  # returns a rough distance measure based on relatedness
  getlinlen <- function(l) sum(gregexpr(";", l, fixed = TRUE)[[1]] > 0) + 1
  if(x == y) return(0)
  minlinlen <- min(sapply(c(x, y), getlinlen))
  anc <- .ancestor(c(x, y))
  anlen <- getlinlen(anc)
  return(1 - anlen/minlinlen)
}


################################################################################
# Get rereplication indices.
#
# This function takes a list of sequences containing duplicates, and returns
#   an vector of named indices that can be used to rereplicate the list if it is
#   dereplicated using the \code{unique} function.
#
# @param x a character vector of sequences or sequence hashes.
# @return an integer vector, named if the input object is named.
# @details TBA
# @author Shaun Wilkinson
# @examples ##TBA
################################################################################
.point <- function(h){
  uh <- unique(h)
  pointers <- seq_along(uh)
  names(pointers) <- uh
  unname(pointers[h])
}
################################################################################

.trimprimer <- function(s, nbases){
  qual <- attr(s, "quality")
  s <- s[-(1:nbases)]
  attr(s, "quality") <- qual[-(1:nbases)]
  return(s)
}

.scanURL <- function(x, retmode = "xml", ...){
  scanURL <- function(z, retmode = "xml", ...){
    res <- tryCatch(if(retmode == "xml") xml2::read_xml(z, ... = ...) else scan(file = z, ... = ...),
                    error = function(er) return(NULL),
                    warning = function(wa) return(NULL))
    return(res)
  }
  for(l in 1:10){
    res <- scanURL(x, retmode = retmode, ... = ...)
    if(!is.null(res)) break else Sys.sleep(5)
  }
  if(is.null(res)) stop("Unable to reach URL, please check internet connection\n")
  return(res)
}


.extractXML <- function(x, species = FALSE, lineages = FALSE, taxIDs = FALSE){
  # x is an xml document
  res <- list()
  x <- xml2::xml_children(x)
  res$accs <- xml2::xml_text(xml2::xml_find_all(x, "GBSeq_locus"))
  res$seqs <- toupper(xml2::xml_text(xml2::xml_find_all(x, "GBSeq_sequence")))
  if(species) res$spps <- xml2::xml_text(xml2::xml_find_all(x, "GBSeq_organism"))
  if(lineages) res$lins <- xml2::xml_text(xml2::xml_find_all(x, "GBSeq_taxonomy"))
  if(taxIDs){
    feattab <- xml2::xml_text(xml2::xml_find_all(x, "GBSeq_feature-table"))
    res$taxs <- gsub(".+taxon:([[:digit:]]+).+", "\\1", feattab)
  }
  if(!all(sapply(res, length) == length(res[[1]]))) res <- NULL
  return(res)
}

# this function is much slower but can return NAs
.extractXML2 <- function(x, species = FALSE, lineages = FALSE, taxIDs = FALSE){
  find_accession <- function(e){
    accession <- e$GBSeq_locus[[1]]
    if(is.null(accession)) accession <- NA
    return(accession)
  }
  find_sequence <- function(e){
    seqnc <- e$GBSeq_sequence[[1]]
    if(is.null(seqnc)) seqnc <- NA
    return(toupper(seqnc))
  }
  find_taxID <- function(s){
    tmp <- unlist(s$`GBSeq_feature-table`$GBFeature$GBFeature_quals, use.names = FALSE)
    if(is.null(tmp)){
      taxID <- NA
    }else{
      taxID <- tmp[grepl("^taxon:", tmp)]
      taxID <- as.integer(gsub("taxon:", "", taxID))
    }
    return(taxID)
  }
  find_species <- function(e){
    spp <- e$GBSeq_organism[[1]]
    if(is.null(spp)) spp <- NA
    return(spp)
  }
  find_lineage <- function(e){
    lin <- e$GBSeq_taxonomy[[1]]
    if(is.null(lin)) lin <- NA
    return(lin)
  }
  res <- list()
  tmp <- xml2::as_list(tmp)[[1]]
  accs <- unname(sapply(tmp, find_accession))
  seqs <- unname(sapply(tmp, find_sequence))
  if(species) spps <- unname(sapply(tmp, find_species))
  if(lineages) lins <- unname(sapply(tmp, find_lineage))
  if(taxIDs) taxs <- unname(sapply(tmp, find_taxID))
  discards <- is.na(accs) | is.na(seqs)
  if(species) discards <- discards | is.na(spps)
  if(lineages) discards <- discards | is.na(lins)
  if(taxIDs) discards <- discards | is.na(taxs)
  # if(sum(!discards) == 0) stop("Error 1\n")
  res$accs  <- accs[!discards]
  res$seqs <- seqs[!discards]
  if(species) res$spps <- spps[!discards]
  if(lineages) res$lins <- lins[!discards]
  if(taxIDs) res$taxs <- taxs[!discards]
  return(res)
}

#
#
#
# query <- "Symbiodinium[ORGN]+AND+Internal[TITL]+AND+2014[MDAT]"
#
# URL1 <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
#                "db=nucleotide&term=",
#                query,
#                "&usehistory=y",
#                "&tool=R")
# X <- .scanURL(URL1, retmode = "xml")
# N <- as.integer(xml2::xml_text(xml2::xml_find_first(X, "Count")))
# WebEnv <- xml2::xml_text(xml2::xml_find_first(X, "WebEnv"))
# QueryKey <- xml2::xml_text(xml2::xml_find_first(X, "QueryKey"))
# nrequest <- N%/%500 + as.logical(N%%500) # 1 if remainder exists, 0 otherwise
# accs <- vector(mode = "list", length = nrequest)
#
# i <- 1
# retstart <- (i - 1) * 500
# contact <- NULL
# URL2 <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
#                "efetch.fcgi?",
#                "db=nucleotide",
#                "&WebEnv=", WebEnv,
#                "&query_key=", QueryKey,
#                "&retstart=", retstart,
#                "&retmax=500",
#                "&rettype=gb",
#                "&retmode=xml",
#                if(!is.null(contact)) paste0("&email=", contact) else NULL)
# tmp <- .scanURL(URL2, retmode = "xml")
#
# test <- .extractXML(tmp, species = TRUE, lineages = TRUE, taxIDs = TRUE)
# test[[5]][1:2]
#
# test2 <- .extractXML2(tmp, species = TRUE, lineages = TRUE, taxIDs = TRUE)
# test[[1]][1:2]
