# Internal 'insect' functions

.dna2char <- function(x){
  dbytes <- as.raw(c(136, 40, 72, 24, 240))
  dchars <- c("A", "C", "G", "T", "N")
  return(paste0(dchars[match(x, dbytes)], collapse = ""))
}

.qual2char <- function(x){
  qbytes <- as.raw(0:93)
  qchars <- strsplit(paste0("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP",
                            "QRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"),
                     split = "")[[1]]
  return(paste0(qchars[match(x, qbytes)], collapse = ""))
}

.char2dna <- function(x){
  dbytes <- as.raw(c(136, 40, 72, 24, 240)) # A, C, G, T, N
  dchars <- c("A", "C", "G", "T", "N")
  return(dbytes[match(strsplit(x, split = "")[[1]], dchars)])
}

.char2qual <- function(x){
  qbytes <- as.raw(0:93)
  qchars <- strsplit(paste0("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP",
                            "QRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"),
                     split = "")[[1]]
  return(qbytes[match(strsplit(x, split = "")[[1]], qchars)])
}

.rc <- function(s){
  # a string of A, C, G, T, and N's
  s <- strsplit(s, split = "")[[1]]
  s <- rev(s)
  dchars <- strsplit("ACGTMRWSYKVHDBN", split = "")[[1]]
  comps <- strsplit("TGCAKYWSRMBDHVN", split = "")[[1]]
  s <- s[s %in% dchars] # rmove spaces etc
  s <- dchars[match(s, comps)]
  s <- paste0(s, collapse = "")
  return(s)
}

.disambiguate <- function(s){
  s <- gsub("[^ACGTMRWSYKVHDBN]", "", s)
  if(grepl("[MRWSYKVHDBN]", s)){
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
    s <- gsub("N", "[ACGT]", s)
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
  lineages <- gsub("\\.$", "", lineages)
  splitfun <- function(s) strsplit(s, split = "; ")[[1]]
  linvecs <- lapply(lineages, splitfun)
  guide <- linvecs[[which.min(sapply(linvecs, length))]]
  a <- 0
  for(l in guide) if(all(sapply(linvecs, function(e) l %in% e))) a <- a + 1
  guide <- if(a > 0) guide[1:a] else character(0)
  lineage <- paste(guide, collapse = "; ")
  lineage <- paste0(lineage, ".")
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
