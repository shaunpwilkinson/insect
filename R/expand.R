#' Expand an existing classification tree.
#'
#' This function is used to grow an existing classification tree, typically
#'   using more relaxed parameter settings than those used when the tree was
#'   created, or if fine-scale control over the tree-learning operation
#'   is required.
#'   Note that the same reference sequence database used to
#'   build the original tree is required for the second argument.
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
#'   ## split the first node
#'   tree <- learn(whales, recursive = FALSE, quiet = FALSE)
#'   ## expand only the first clade
#'   tree <- expand(tree, whales, clades = "1", quiet = TRUE)
#'  }
################################################################################
expand <- function(tree, x, clades = "0", refine = "Viterbi", iterations = 50,
                   nstart = 20, minK = 2, maxK = 2, minscore = 0.9, probs = 0.5,
                   retry = TRUE, resize = TRUE, maxsize = max(sapply(x, length)),
                   recursive = TRUE, cores = 1, quiet = TRUE, ...){
  ## Establish which parts fof the tree to expand
  clades <- gsub("0", "", clades)
  if(!(identical(attr(tree, "sequences"), seq_along(x)))){
    stop("tree is incompatible with sequences\n")
  }
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
  ## following lines are for trees that have been stripped of memory-intensive elements
  if(!quiet) cat("Dereplicating sequences\n")
  hashes <- attr(x, "hashes")
  if(is.null(hashes)) hashes <- hash(x)
  # attr(tree, "hashes") <- NULL
  duplicates <- attr(x, "duplicates")
  if(is.null(duplicates)) duplicates <- duplicated(hashes)
  if(!quiet) cat("Found", sum(!duplicates), "unique sequences\n")
  # attr(tree, "duplicates") <- NULL
  pointers <- attr(x, "pointers")
  if(is.null(pointers)) pointers <- .point(hashes)
  seqweights <- attr(x, "weights")
  if(is.null(seqweights)) seqweights <- aphid::weight(x, k = 5)
  # attr(tree, "weights") <- NULL ## replaced later
  # lineages <- gsub("\\.$", "", attr(x, "lineage"))
  # lineages <- paste0(lineages, "; ", attr(x, "species"))
  lineages <- attr(x, "lineage")
  # lineages <- paste0(lineages, "; ~", attr(x, "species"), "~")
  x <- x[!duplicates]# strip attrs regardless of duplicates
  if(any(duplicates)){
    # fullseqset <- x # needed??
    # x <- x[!duplicates]  ## attributes not needed
    lineages <- sapply(split(lineages, pointers), .ancestor)
    seqweights <- sapply(split(seqweights, pointers), sum)
    ###
    rmduplicates <- function(node, whchunq, pointers){
      tmp <- attr(node, "sequences")[attr(node, "sequences") %in% whchunq]
      attr(node, "sequences") <- pointers[tmp]
      return(node)
    }
    tree <- dendrapply(tree, rmduplicates, which(!duplicates), pointers)
    # if(nrow(distances) == nseq) distances <- distances[!duplicates, ]
    # if(nrow(kmers) == nseq) kmers <- kmers[!duplicates, ]
    # rm attrs, condensed version in final tree
  }# else distances <- distances[ , ]

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
  ## calculate k-mers once but only if run on single core
  ## otherwise uses excessive memory (since this matrix can be > 1GB)
  # kmers <- attr(tree, "kmers")
  #if(is.null(kmers)) kmers <- kmer::mbed(x)
  if(ncores == 1){
    if(!quiet) cat("Counting k-mers\n")
    kmers <- kmer::kcount(x, k = 5)/(sapply(x, length) - 4) #k - 1 = 4
  # }else if(has_duplicates & nrow(kmers) == nseq & ncores = 1){
  #   kmers <- kmers[!duplicates, ]
  }else kmers <- NULL
  # attr(tree, "kmers") <- NULL ## replaced later
  ## prev line commented to prevent k-mer stripping
  ### recursively split nodes
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
    if(!quiet) cat("Recursively splitting basal tree branches\n")
    repeat{
      nmembers <- integer(0)
      eligible <- logical(0)
      tmp <- fm1(tree)
      rm(tmp)
      nmembers <- nmembers[eligible]
      # if(!any(eligible) | length(nmembers) >= 2 * ncores) break
      if(!any(eligible) | all(nmembers < 50)) break
      whichclade <- names(nmembers)[which.max(nmembers)]
      index <- gsub("([[:digit:]])", "[[\\1]]", whichclade)
      toeval <- paste0("tree", index, "<- .fork(tree",
                       index, ", x, lineages, refine = refine, ",
                       "nstart = nstart, iterations = iterations, minK = minK, maxK = maxK, ",
                       "minscore = minscore, probs = probs, retry = retry, resize = resize, ",
                       "maxsize = maxsize, kmers = kmers, seqweights = seqweights, ",
                       "cores = cores, quiet = quiet, ... = ...)")
      eval(parse(text = toeval))
      ss <- FALSE # split success; prevents build note due to lack of visible binding
      eval(parse(text = paste0("ss <- is.list(tree", index, ")")))
      # prevent multiple attempts to split the same node
      if(!ss) eval(parse(text = paste0("attr(tree", index, ", 'lock') <-TRUE")))
    }
    clades <- names(nmembers)
    indices <- gsub("([[:digit:]])", "[[\\1]]", clades)
    rm(nmembers)
    rm(eligible)
    #}
    trees <- vector(mode = "list", length = length(clades)) # = length(indices)
    for(i in seq_along(indices)){
      eval(parse(text = paste0("trees[[", i, "]] <- tree", indices[i])))
    }
    if(!quiet){
      cat("Recursively splitting terminal tree branches\n")
      cat("Feedback suppressed, this could take a while...\n")
    }
    trees <- parallel::parLapply(cores, trees, .forkr, x, lineages, refine = refine,
                                 nstart = nstart, iterations = iterations,
                                 minK = minK, maxK = maxK, minscore = minscore,
                                 probs = probs, retry = retry, resize = resize,
                                 maxsize = maxsize,
                                 kmers = kmers, # large matrix could cause memory probs
                                 seqweights = seqweights, cores = 1,
                                 quiet = TRUE, ... = ...)
    for(i in seq_along(trees)){
      eval(parse(text = paste0("tree", indices[i], "<- trees[[", i, "]]")))
      trees[[i]] <- NA
      gc()
    }
  }else{
    indices <- gsub("([[:digit:]])", "[[\\1]]", clades)
    for(i in seq_along(indices)){
      toeval <- paste0("tree", indices[i],
                       if(recursive) "<-.forkr(tree" else "<-.fork(tree",
                       indices[i], ", x, lineages, refine = refine, nstart = nstart, ",
                       "iterations = iterations, minK = minK, maxK = maxK, ",
                       "minscore = minscore, probs = probs, retry = retry, resize = resize, ",
                       "maxsize = maxsize, kmers = kmers, seqweights = seqweights, ",
                       "cores = cores, quiet = quiet, ... = ...)")
      eval(parse(text = toeval))
    }
  }
  if(stopclustr) parallel::stopCluster(cores)
  ### remove kmers since can be memory hungry, prevent next operations
  rm(kmers)
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
  # rmaligs <- function(node){
  #   if(is.leaf(node)) attr(node, "model")$alignment <- NULL
  #   return(node)
  # }
  # tree <- dendrapply(tree, rmaligs) # more memory efficient
  encodemods <- function(node){
    if(is.leaf(node)) attr(node, "model") <- encodePHMM(attr(node, "model"))
    return(node)
  }
  tree <- dendrapply(tree, encodemods) # more memory efficient
  #if(has_duplicates) x <- fullseqset
  # attributes(x) <- tmpxattr
  if(!quiet) cat("Resetting node heights\n")
  tree <- phylogram::reposition(tree)
  if(!quiet) cat("Making tree ultrametric\n")
  tree <- phylogram::ultrametricize(tree)
  if(!quiet) cat("Done\n")
  class(tree) <- c("insect", "dendrogram")
  return(tree)
}
################################################################################

