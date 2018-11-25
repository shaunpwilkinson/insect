#' Expand an existing classification tree.
#'
#' This function is used to grow an existing classification tree, typically
#'   using more relaxed parameter settings than those used when the tree was
#'   created, or if fine-scale control over the tree-learning operation
#'   is required.
#'
#' @param tree an object of class \code{"insect"}.
#' @param clades a vector of character strings giving the binary indices
#'   matching the labels of the nodes that are to be expanded.
#'   Defaults to "0", meaning all subclades are expanded.
#'   See below for further details on clade indexing.
#' @param recursive logical indicating whether the splitting process
#'   should continue recursively until the discrimination criteria
#'   are not met (TRUE; default), or whether a single split should
#'   take place at each of the nodes specified in \code{clades}.
#' @param ... further arguments to be passed on to \code{\link[aphid]{train}}).
#' @inheritParams learn
#' @return an object of class \code{"insect"}.
#' @details
#'   The clade indexing system used here is based on character strings,
#'   where "0" refers to the root node,
#'   "01" is the first child node, "02" is the second child node,
#'   "011" is the first child node of the first child node, etc.
#'   Note that this means each node cannot have more than 9 child nodes.
#' @author Shaun Wilkinson
#' @seealso \code{\link{learn}}.
#' @examples
#' \donttest{
#'   data(whales)
#'   data(whale_taxonomy)
#'   ## split the first node
#'   set.seed(123)
#'   tree <- learn(whales, db = whale_taxonomy, recursive = FALSE)
#'   ## expand only the first clade
#'   tree <- expand(tree, clades = "1")
#'  }
################################################################################
expand <- function(tree, clades = "0", refine = "Viterbi", iterations = 50,
                   nstart = 20, minK = 2, maxK = 2, minscore = 0.9, probs = 0.5,
                   retry = TRUE, resize = TRUE, maxsize = 1000,
                   recursive = TRUE, cores = 1, quiet = FALSE, verbose = FALSE,
                   ...){
  dots <- list(...)

  if(is.null(attr(tree, "numcode"))){
    x <- attr(tree, "trainingset")
    pointers <- attr(x, "rerep.pointers")
    xnames <- attr(x, "rerep.names")
    lineages <- attr(x, "lineages") # same length as full set, includes full strings. can be null
    kmers <- attr(tree, "kmers")
    ksize <- attr(tree, "k")
    if(is.null(kmers) | is.null(ksize)){
      ksize <- if(is.null(dots$k)) 4 else dots$k
      kmers <- .encodekc(kmer::kcount(x, k = ksize))
    }else{
      stopifnot(nrow(kmers) == length(x))
    }
  }else{ #amino classifier
    if(is.null(attr(tree, "frame"))) stop("Error code 0227\n")
    x <- attr(tree, "xaa")
    if(is.null(x)){
      x <- rereplicate(attr(tree, "trainingset"))
      xlengths <- vapply(x, length, 0L, USE.NAMES = FALSE)
      xrems <- xlengths %% 3
      if(length(unique(xrems)) > 1) stop("Error code 3274\n")
      x <- ape::as.character.DNAbin(x)
      x <- lapply(x, seqinr::translate, numcode = attr(tree, "numcode"), frame = attr(tree, "frame"))
      x <- ape::as.AAbin(x)
      keeps <- sapply(x, function(v) !any(v == as.raw(42)))
      if(any(!keeps)) stop("Error code 0442\n")
      hashes <- hash(x)
      xnames <- names(x)
      pointers <- .point(hashes)
      x <- x[!duplicated(hashes)]
    }else{
      pointers <- attr(x, "rerep.pointers")
      xnames <- attr(x, "rerep.names")
    }
    lineages <- attr(attr(tree, "trainingset"), "lineages")
    ksize <- 2
    kmers <- .encodekc(kmer::kcount(x, k = ksize))
  }
  #attr(tree, "kmers") <- NULL ## replaced later

  if(is.null(lineages)){
    taxIDs <- as.integer(gsub(".+\\|", "", xnames))
    lineages <- get_lineage(taxIDs, db = attr(tree, "taxonomy"), numbers = TRUE)
    lineages <- vapply(lineages, paste0, "", collapse = "; ")
  }
  seqweights <- attr(x, "rerep.weights") # length of derep'd set
  if(is.null(seqweights)) seqweights <- aphid::weight(x, k = if(is.null(attr(tree, "numcode"))) 5 else 2)

  #seqweights <- NULL # assigned to .partition 20181024

  ## Establish which parts of the tree to expand
  clades <- gsub("0", "", clades)
  indices <- gsub("([[:digit:]])", "[[\\1]]", clades)
  findnestedleaves <- function(node){
    if(!is.list(node) & is.null(attr(node, "lock"))){
      nestedleaves <<- c(nestedleaves, attr(node, "clade"))
    }
    return(node)
  }
  fnlr <- function(node){ # find nested leaves recursively
    node <- findnestedleaves(node)
    if(is.list(node) & is.null(attr(node, "lock"))) node[] <- lapply(node, fnlr)
    return(node)
  }
  allnestedleaves <- vector(mode = "list", length = length(indices))
  for(i in seq_along(indices)){
    validclade <- TRUE
    toeval <- paste0("validclade <- is.leaf(tree", indices[i],
                     ") | is.list(tree", indices[i], ")")
    eval(parse(text = toeval))
    if(!validclade) stop("Clade ", clades[i], " is out of bounds\n")
    nestedleaves <- character(0)
    eval(parse(text = paste0("tmp <- fnlr(tree", indices[i], ")")))
    rm(tmp)
    allnestedleaves[[i]] <- nestedleaves
  }
  rm(nestedleaves)
  clades <- unlist(allnestedleaves, use.names = FALSE)
  if(length(clades) == 0) return(tree)
  indices <- gsub("([[:digit:]])", "[[\\1]]", clades)
  duplicates <- duplicated(pointers)
  if(any(duplicates)){
    lineages <- sapply(split(lineages, pointers), .ancestor)
    #seqweights <- sapply(split(seqweights, pointers), sum) ### removed 20181014
    rmduplicates <- function(node, whchunq, pointers){
      tmp <- attr(node, "sequences")[attr(node, "sequences") %in% whchunq]
      attr(node, "sequences") <- pointers[tmp]
      return(node)
    }
    tree <- dendrapply(tree, rmduplicates, which(!duplicates), pointers)
  }
  ### set up multithread if required
  if(inherits(cores, "cluster")){
    ncores <- length(cores)
    stopclustr <- FALSE
  }else if(identical(cores, 1)){
    ncores <- 1
    stopclustr <- FALSE
  }else{ # create cluster object
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
    # if(cores > navailcores) stop("Number of cores is more than the number available")
    # if(!quiet) cat("Multithreading over", cores, "cores\n")
    if(cores == 1){
      ncores <- 1
      stopclustr <- FALSE
    }else{
      ncores <- cores
      if(!quiet) cat("Initializing cluster with", ncores, "cores\n")
      cores <- parallel::makeCluster(cores)
      stopclustr <- TRUE
    }
  }

  # if(ncores == 1){
  #   if(!quiet) cat("Counting k-mers\n")
  #   dots <- list(...)
  #   kmers <- kmer::kcount(x, k = if(is.null(dots$k)) 5 else dots$k)
  #   kmers <- kmers/(sapply(x, length) - 4) #k - 1 = 4
  # }else kmers <- NULL
  # kmers <- NULL # need this if using partition to count kmers

  ## prev line commented to prevent k-mer stripping
  ### recursively split nodes
  switchpoint <- max(round((length(x) * 0.8)/ncores), 50L)
  if(switchpoint > 500L) switchpoint <- 500L
  if(switchpoint < ncores * 4) switchpoint <- ncores * 4
  ## switch from basal to terminal node recursion method
  if(ncores > 1 & recursive){
    #if(length(clades) < ncores){
    lockleaves <- function(node, exceptions){
      if(!is.list(node)){
        if(!attr(node, "clade") %in% exceptions) attr(node, "lock") <- TRUE
      }
      return(node)
    }
    ll1 <- function(node, exceptions){
      node <- lockleaves(node, exceptions)
      if(is.list(node)) node[] <- lapply(node, ll1, exceptions)
      return(node)
    }
    tree <- ll1(tree, exceptions = clades)
    findnmembers <- function(node){
      if(!is.list(node)){
        numberofseqs <- length(attr(node, "sequences"))
        names(numberofseqs) <- attr(node, "clade")
        nmembers <<- c(nmembers, numberofseqs)
        eligible <<- c(eligible, is.null(attr(node, "lock")))
      }
      return(node)
    }
    fm1 <- function(node){
      node <- findnmembers(node)
      if(is.list(node)) node[] <- lapply(node, fm1)
      return(node)
    }
    if(!quiet) cat("Recursively splitting basal tree nodes\n")
    repeat{
      nmembers <- integer(0)
      eligible <- logical(0)
      tmp <- fm1(tree)
      rm(tmp)
      nmembers <- nmembers[eligible]
      # if(!any(eligible) | length(nmembers) >= 2 * ncores) break
      if(!any(eligible) | all(nmembers < switchpoint)) break
      whichclade <- names(nmembers)[which.max(nmembers)]
      index <- gsub("([[:digit:]])", "[[\\1]]", whichclade)
      toeval <- paste0("tree", index, "<- .fork(tree",
                       index, ", x, lineages, refine = refine, ",
                       "nstart = nstart, iterations = iterations, minK = minK, maxK = maxK, ",
                       "minscore = minscore, probs = probs, retry = retry, resize = resize, ",
                       "maxsize = maxsize, kmers = kmers, ksize = ksize, seqweights = seqweights, ",
                       "cores = cores, quiet = quiet, verbose = verbose, ... = ...)")
      eval(parse(text = toeval))
      ss <- FALSE # split success; prevents build note due to lack of visible binding
      eval(parse(text = paste0("ss <- is.list(tree", index, ")")))
      # prevent multiple attempts to split the same node
      if(!ss) eval(parse(text = paste0("attr(tree", index, ", 'lock') <-TRUE")))
    }
    clades <- names(nmembers)
    indices <- gsub("([[:digit:]])", "[[\\1]]", clades)
    indices <- indices[.balance(nmembers, ncores)]
    rm(nmembers)
    rm(eligible)
    trees <- vector(mode = "list", length = length(clades)) # = length(indices)
    for(i in seq_along(indices)){
      eval(parse(text = paste0("trees[[", i, "]] <- tree", indices[i])))
    }
    if(!quiet){
      cat("Recursively splitting terminal tree nodes\n")
      if(verbose) cat("Feedback suppressed, this could take a while...\n")
    }
    trees <- parallel::parLapply(cores, trees, .forkr, x, lineages, refine = refine,
                                 nstart = nstart, iterations = iterations,
                                 minK = minK, maxK = maxK, minscore = minscore,
                                 probs = probs, retry = retry, resize = resize,
                                 maxsize = maxsize,
                                 kmers = kmers, ksize = ksize, # large matrix could cause memory probs
                                 seqweights = seqweights,
                                 cores = 1,
                                 quiet = TRUE, ... = ...)
    for(i in seq_along(trees)){
      eval(parse(text = paste0("tree", indices[i], "<- trees[[", i, "]]")))
    }
    rm(trees)
    gc()
  }else{
    indices <- gsub("([[:digit:]])", "[[\\1]]", clades)
    for(i in seq_along(indices)){
      toeval <- paste0("tree", indices[i],
                       if(recursive) "<-.forkr(tree" else "<-.fork(tree",
                       indices[i], ", x, lineages, refine = refine, nstart = nstart, ",
                       "iterations = iterations, minK = minK, maxK = maxK, ",
                       "minscore = minscore, probs = probs, retry = retry, resize = resize, ",
                       "maxsize = maxsize, kmers = kmers, ksize = ksize, seqweights = seqweights, ",
                       "cores = cores, quiet = quiet, verbose = verbose, ... = ...)")
      eval(parse(text = toeval))
    }
  }
  if(stopclustr) parallel::stopCluster(cores)
  ### remove kmers since can be memory hungry, prevent next operations
  #rm(kmers)
  gc()
  ### fix midpoints, members, heights and leaf integers
  ### note changes here also apply to 'learn' function
  if(!quiet) cat("Setting midpoints and members attributes\n")
  #tree <- settreeattr(tree)
  tree <- phylogram::remidpoint(tree)
  class(tree) <- "dendrogram"
  rm_locks <- function(node){
    attr(node, "lock") <- NULL
    return(node)
  }
  tree <- dendrapply(tree, rm_locks)
  reduplicate <- function(node, pointers){
    seqs <- attr(node, "sequences")
    attr(node, "nunique") <- length(seqs)
    attr(node, "sequences") <- which(pointers %in% seqs)
    attr(node, "ntotal") <- length(attr(node, "sequences"))
    return(node)
  }
  if(!quiet) cat("Repatriating duplicate sequences with tree\n")
  tree <- dendrapply(tree, reduplicate, pointers = pointers)
  encodemods <- function(node){
    if(is.leaf(node)) attr(node, "model") <- encodePHMM(attr(node, "model"))
    return(node)
  }
  tree <- dendrapply(tree, encodemods) # more memory efficient
  truncatetaxIDs <- function(node){
    if(is.null(attr(node, "taxID"))){ # in case tree manually expanded
      stopifnot(!is.null(attr(node, "lineage")))
      attr(node, "taxID") <- as.integer(gsub(".+; ", "", attr(node, "lineage")))
      attr(node, "lineage") <- NULL
    }
    return(node)
  }
  tree <- dendrapply(tree, truncatetaxIDs) # more memory efficient
  #if(has_duplicates) x <- fullseqset
  # attributes(x) <- tmpxattr
  if(!quiet) cat("Resetting node heights\n")
  tree <- phylogram::reposition(tree)
  #attr(tree, "kmers") <- kmers # can no longer do this due to temporary AA kmers
  rm(kmers)
  rm(x)
  class(tree) <- c("insect", "dendrogram")
  return(tree)
}
################################################################################

