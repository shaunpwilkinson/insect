# Internal 'insect' functions
#' @noRd
.qual2char <- function(x){
  if(is.null(x)) return(NULL) # needed for dna2char
  qbytes <- as.raw(0:93)
  qchars <- strsplit(paste0("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP",
                            "QRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"),
                     split = "")[[1]]
  return(paste0(qchars[match(x, qbytes)], collapse = ""))
}

#' @noRd
.char2qual <- function(x){
  if(is.null(x)) return(NULL)# needed for chardna
  qbytes <- as.raw(0:93)
  qchars <- strsplit(paste0("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP",
                            "QRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"),
                     split = "")[[1]]
  return(qbytes[match(strsplit(x, split = "")[[1]], qchars)])
}


#' @noRd
## tree is a "dendrogram" object (can be a node)
## x is a DNAbin object - should not contain duplicates, must have lineage attrs
.forkr <- function(tree, x, lineages, refine = "Viterbi", nstart = 10,
                   iterations = 50, minK = 2, maxK = 2,
                   minscore = 0.9, probs = 0.05, retry = TRUE, resize = TRUE,
                   maxsize = NULL, kmers = NULL,
                   seqweights = "Gerstein", cores = 1, quiet = FALSE, ...){
  tree <- .fork(tree, x, lineages, refine = refine, nstart = nstart,
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




#' @noRd
.ancestor <- function(lineages){
  # input and output both semicolon-delimited character string(s)
  if(length(lineages) == 1) return(lineages)
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
#' @noRd
.point <- function(h){
  uh <- unique(h)
  pointers <- seq_along(uh)
  names(pointers) <- uh
  unname(pointers[h])
}
################################################################################

#' @noRd
.scanURL <- function(x, retmode = "xml", ...){
  scanURL <- function(z, retmode = "xml", ...){
    errfun <- function(er){
      closeAllConnections()
      return(NULL)
    }
    res <- tryCatch(if(retmode == "xml") xml2::read_xml(z, ... = ...) else
      scan(file = z, ... = ...), error = errfun, warning = errfun)
    return(res)
  }
  for(l in 1:10){
    res <- scanURL(x, retmode = retmode, ... = ...)
    if(!is.null(res)) break else Sys.sleep(5)
  }
  if(is.null(res)) stop("Unable to reach URL, please check connectivity\n")
  return(res)
}

#' @noRd
.extractXML <- function(x, taxIDs = TRUE, species = FALSE, lineages = FALSE){
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

#' @noRd
# this function is much slower but can return NAs
.extractXML2 <- function(x, taxIDs = TRUE, species = FALSE, lineages = FALSE){
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
  tmp <- xml2::as_list(x)[[1]]
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

#' @noRd
.encode1 <- function(xx){ ## a number between 0.000001 and 100
  ## encodes as two bytes to (almost) 4 signif figs
  ## first 13 for the digits (2^13 = 8192) last 3 bits for exponent
  #(1e-6 -> 0, 1e+1 -> 7)
  if(!is.finite(xx)) return(as.raw(c(0, 0)))
  stopifnot(xx >= 1e-06 & xx < 100)
  expo <- floor(log10(xx))
  digi <- round((xx/10^(expo - 3) - 1000) * 0.9102122, 0)
  ## scale factor 8191/8999 = 0.9102122, 3 to get from 1 to 1000
  return(packBits(c(intToBits(digi)[1:13], intToBits(expo + 6)[1:3])))
}

#' @noRd
#'
.decode1 <- function(zz){ ## zz is a 2-byte raw vec
  res <- rawToBits(zz)
  expo <- as.integer(packBits(c(res[14:16], raw(5)))) - 6
  digi <- as.integer(packBits(c(res[1:13], raw(3)))) # 2 ints between 0 and 8191
  digi <- sum(digi* c(1, 256))
  digi <- round(digi/0.9102122, 0) + 1000 # integer between 1000 and 9999
  return(digi * 10^(expo - 3)) # the 3 reduces from 1000 to 1
}

#' @noRd
.encodek <- function(xx){  ## kmer presence/absence matrix
  if(mode(xx) == "raw") return(xx)
  dims <- dim(xx)
  xx <- as.integer(xx != 0)
  dim(xx) <- dims
  rem <- ncol(xx) %% 8 #remainder
  if(rem > 0){
    newcol <- matrix(0L, nrow = nrow(xx), ncol = 8 - rem)
    xx <- cbind(xx, newcol)
  }
  res <- packBits(t(xx))

  dims[2] <- dims[2]/8
  dim(res) <- dims[c(2, 1)]
  res <- t(res)
  return(res)
}

#' @noRd
.decodek <- function(zz){ ## kmer presence/absence matrix
  if(mode(zz) != "raw") return(zz)
  dims <- dim(zz)
  dims[2] <- dims[2] * 8
  res <- rawToBits(t(zz))
  res <- as.integer(res)
  dim(res) <- dims[c(2, 1)]
  res <- t(res)
  return(res)
}

#' @noRd
# Check if object is DNA
.isDNA <- function(x){
  if(inherits(x, "DNAbin")){
    return(TRUE)
  }else if(inherits(x, "AAbin")){
    return(FALSE)
  }else if(mode(x) == "character"){
    return(FALSE)
  }else if(mode(x) == "raw"){
    return(all(x %in% as.raw(c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48,
                               224, 176, 208, 112, 240, 4, 2))))
  }else if(mode(x) == "list"){
    if(length(x) > 0){
      return(all(unlist(x, use.names = FALSE) %in%
                   as.raw(c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48,
                            224, 176, 208, 112, 240, 4, 2))))
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

#' @noRd
# Check if object is amino acid sequence
.isAA <- function(x){
  if(inherits(x, "AAbin")){
    return(TRUE)
  }else if(inherits(x, "DNAbin")){
    return(FALSE)
  }else if(mode(x) == "character"){
    return(FALSE)
  }else if(mode(x) == "raw"){
    return(all(x %in% as.raw(c(65:90, 42, 45, 63))))
  }else if(mode(x) == "list"){
    if(length(x) > 0){
      return(all(unlist(x, use.names = FALSE) %in% as.raw(c(65:90, 42, 45, 63))))
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}
