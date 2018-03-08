# Internal 'insect' functions

.qual2char <- function(x){
  if(is.null(x)) return(NULL) # needed for dna2char
  qbytes <- as.raw(0:93)
  qchars <- strsplit(paste0("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP",
                            "QRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"),
                     split = "")[[1]]
  return(paste0(qchars[match(x, qbytes)], collapse = ""))
}


.char2qual <- function(x){
  if(is.null(x)) return(NULL)# needed for chardna
  qbytes <- as.raw(0:93)
  qchars <- strsplit(paste0("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP",
                            "QRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"),
                     split = "")[[1]]
  return(qbytes[match(strsplit(x, split = "")[[1]], qchars)])
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
