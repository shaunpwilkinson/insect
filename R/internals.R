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


.learn1 <- function(tree, x, refine = "Viterbi", iterations = 50, minK = 2, maxK = 2,
                    minscore = 0.9, probs = 0.05, resize = TRUE, seqweights = "Gerstein",
                    quiet = FALSE, ...){
  #tree is a "dendrogram" object (can be a node)
  # x is a DNAbin object - should not contain duplicates
  tree <- fork(tree, x = x, refine = refine, iterations = iterations, minK = minK,
               maxK = maxK, minscore = minscore, probs = probs, resize = resize,
               seqweights = seqweights, quiet = quiet, ... = ...)
  if(is.list(tree)) tree[] <- lapply(tree, .learn1, x = x, refine = refine,
                                     iterations = iterations, minK = minK, maxK = maxK,
                                     probs = probs, resize = resize, minscore = minscore,
                                     seqweights = seqweights, quiet = quiet, ... = ...)
  return(tree)
}
