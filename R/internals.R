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
                   maxsize = NULL, kmers = NULL, ksize = NULL, seqweights = NULL,
                   cores = 1, quiet = FALSE, ...){
  tree <- .fork(tree, x, lineages, refine = refine, nstart = nstart,
               iterations = iterations, minK = minK,
               maxK = maxK, minscore = minscore, probs = probs,
               retry = retry, resize = resize, maxsize = maxsize,
               kmers = kmers, ksize = ksize, seqweights = seqweights, cores = cores,
               quiet = quiet, ... = ...)
  if(is.list(tree)) tree[] <- lapply(tree, .forkr, x = x, lineages = lineages,
                                     refine = refine, nstart = nstart,
                                     iterations = iterations, minK = minK, maxK = maxK,
                                     probs = probs, retry = retry, resize = resize,
                                     maxsize = maxsize, minscore = minscore,
                                     kmers = kmers, ksize = ksize, seqweights = seqweights,
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
    if(!is.null(res)) break else Sys.sleep(1)
  }
  if(is.null(res)){
    suppressMessages(tmp <- scan(x, what = "", sep = "\n"))
    tmp <- gsub("\xAE", "\x20", tmp) ## annoying trademark char in some genbank files
    tmpf <- tempfile()
    cat(tmp, file = tmpf, sep = "\n")
    errfun <- function(er){
      closeAllConnections()
      return(NULL)
    }
    res <- tryCatch(if(retmode == "xml") xml2::read_xml(tmpf, ... = ...) else
      scan(file = tmpf, ... = ...), error = errfun, warning = errfun)
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
.encodekc <- function(xx){  ## kmer count matrix
  if(mode(xx) == "raw") return(xx)
  xx <- round(xx)
  if(ncol(xx) %% 2 == 1) xx <- cbind(xx, 0L)
  xx[xx > 15] <- 15
  dimnames(xx) <- NULL
  dims <- dim(xx)
  xx <- as.raw(xx)
  dim(xx) <- dims
  fun <- function(xxx){ # 2-col int mat
    m1 <- matrix(rawToBits(xxx[, 1]), ncol = 8, byrow = TRUE)[, 1:4]
    m2 <- matrix(rawToBits(xxx[, 2]), ncol = 8, byrow = TRUE)[, 1:4]
    return(packBits(t(cbind(m2, m1))))
  }
  res <- matrix(as.raw(0), nrow = dims[1], ncol = dims[2]/2)
  guide <- split(seq(1, dims[2]), f = rep(seq(1, dims[2]/2), each = 2))
  for(i in seq_along(guide)){
    res[, i] <- fun(xx[, guide[[i]]])
  }
  return(res)
}

#' @noRd
.decodekc <- function(zz){ ## kmer count matrix (max count 15)
  if(mode(zz) != "raw") return(zz)
  dims <- dim(zz)
  fun <- function(zzz){
    # takes a raw vector
    # returns a 2 col matrix
    mymat <- matrix(rawToBits(zzz), ncol = 8L, byrow = TRUE)
    m1 <- m2 <- matrix(as.raw(0L), ncol = 8L, nrow = nrow(mymat))
    m1[, 1:4] <- mymat[, 5:8]
    m2[, 1:4] <- mymat[, 1:4]
    m1 <- as.integer(packBits(t(m1)))
    m2 <- as.integer(packBits(t(m2)))
    return(cbind(m1, m2))
  }
  guide <- split(seq(1, dims[2] * 2), f = rep(seq(1, dims[2]), each = 2))
  out <- matrix(0L, nrow = dims[1], ncol = dims[2] * 2)
  for(i in seq_along(guide)) out[, guide[[i]]] <- fun(zz[, i])
  return(out)
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


#' @noRd
# Nearest neighbor OTU search
.otu <- function(z, k = 4, threshold = 0.97){
  if(is.null(names(z))) names(z) <- paste0("S", seq_along(z))
  stopifnot(!any(duplicated(names(z))))
  hashes <- insect::hash(z)
  pointers <- .point(hashes)
  catchnames <- names(z)
  z <- z[!duplicated(hashes)]
  derepnames <- names(z)
  zl <- vapply(z, length, 0L)
  kmers <- kmer::kcount(z, k = k)
  kmers <- kmers/(apply(kmers, 1, sum) + k - 1)
  myseq <- c(2, 5, (2:20)^4)
  myseq <- myseq[myseq < length(z)]
  res <- list()
  counter <- 1L
  takens <- logical(length(z))
  names(takens) <- names(z)
  for(m in myseq){
    if(nrow(kmers) - 1 < m) m <- nrow(kmers) - 1
    neighbors <- RANN::nn2(kmers, k = m)
    idx <- neighbors$nn.idx
    dis <- 1 - ((neighbors$nn.dists)^2 * zl/(2*k))
    keeps <- dis[, m] < threshold
    incursions <- idx[, 1] != seq_along(z) #
    wincs <- which(incursions)
    for(i in wincs){
      ind <- match(i, idx[i, ])
      if(!is.na(ind)){
        tmp <- idx[i, ind]
        idx[i, ind] <- idx[i, 1]
        idx[i, 1] <- tmp
        tmp <- dis[i, ind]
        dis[i, ind] <- dis[i, 1]
        dis[i, 1] <- tmp
      }else{
        idx[i, ] <- c(i, idx[i, seq_len(m - 1)])
        dis[i, ] <- c(0, dis[i, seq_len(m - 1)])
        ## warning("code whero\n")
      }
    }
    if(sum(keeps) == 0){ ## final biggest cluster
      degs <- apply(dis, 1, function(v) sum(v[v >= threshold]))
      res[[counter]] <- names(z)[order(degs, decreasing = TRUE)]
      break
    }
    idx <- idx[keeps, , drop = FALSE]
    dis <- dis[keeps, , drop = FALSE]
    degs <- apply(dis, 1, function(v) sum(v[v >= threshold]))
    names(degs) <- names(z)[keeps]
    degorder <- order(degs, decreasing = TRUE)
    for(i in degorder){ # indices in reduced set
      if(!takens[names(degs)[i]]){
        nmmbs <- sum(dis[i, ] >= threshold)
        members <- idx[i, 1:nmmbs] # indices in semi-full set
        members <- names(z)[members]
        members <- members[!takens[members]] # character matching
        takens[members] <- TRUE # character matching
        res[[counter]] <- members
        counter <- counter + 1
      }
    }
    kmers <- kmers[!keeps, , drop = FALSE]
    z <- z[!keeps]
    zl <- zl[!keeps]
    if(length(z) == 0L) {
      break
    }else if(length(z) == 1L){ # cant do nn search with 1 row
      res[[counter]] <- names(z)
      break
    }
  }
  res <- res[order(sapply(res, length), decreasing = TRUE)]
  out <- rep(seq_along(res), sapply(res, length))
  names(out) <- unlist(res, use.names = FALSE)
  centrals <- vapply(res, "[", "", 1)
  out <- out[derepnames]
  out <- out[pointers]
  centrals <- catchnames %in% centrals
  catchnames[centrals] <- paste0(catchnames[centrals], "*")
  names(out) <- catchnames
  return(out)
}

#' @noRd
# Load balancing for parallel apply functions
## x is vector of load sizes
## ncl is number of clusters
## returns an index vector e.g:
## idx <- .balance(x, ncl)
## x <- x[idx]
## y <- parLapply(cl, x, ...)
## y <- y[order(idx)]
.balance <- function(x, ncl){
  nx <- length(x)
  if(ncl >= nx) return(seq_along(x))
  i <- seq_len(nx)
  fuzz <- min((nx - 1L)/1000, 0.4 * nx/ncl)
  breaks <- seq(1 - fuzz, nx + fuzz, length.out = ncl + 1L)
  tmp <- structure(split(i, cut(i, breaks)), names = NULL)
  tmplens <- vapply(tmp, length, 0L, USE.NAMES = FALSE)
  maxlen <- max(vapply(tmp, length, 0L))
  idx <- vector(mode = "list", length = maxlen)
  for(i in seq_along(idx)){
    idx[[i]] <- sapply(tmp, "[", i)
  }
  idx <- unlist(idx, use.names = FALSE)
  idx <- idx[!is.na(idx)]
  idx <- idx[order(order(x, decreasing = T))]
  return(idx)
}
