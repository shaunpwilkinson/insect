################################################################################
################################################################################
#' Expand an existing classification tree.
#'
#' This function is used to grow an existing classification tree, typically
#'   using more relaxed settings to those used when the tree was created.
#'
#' @param tree an object of class \code{"insect"}.
#' @param clades a vector of character strings giving the binary indices
#'   matching the labels of the nodes that are to be collapsed.
#' @param recursive logical indicating whether the splitting process
#'   should continue recursively until the discrimination criteria
#'   are not met (TRUE; default), or whether a single split should
#'   take place at each of the nodes specified in \code{clades}.
#' @param ... further arguments to be passed to \code{\link{fork}}.
#' @inheritParams learn
#' @return an object of class \code{"insect"}.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @seealso \code{\link{contract}}, \code{\link{learn}}, \code{\link{purge}}
#' @examples
#'   ## TBA
################################################################################
expand <- function(tree, clades = "", refine = "Viterbi", iterations = 50,
                   minK = 2, maxK = 2, minscore = 0.9, probs = 0.05,
                   resize = TRUE, recursive = TRUE, cores = 1,
                   quiet = FALSE, ...){
  x <- attr(tree, "sequences") #full sequence set
  tmpxattr <- attributes(x)
  x <- x[seq_along(x)] # removes attributes which can be memory hungry
  attr(tree, "sequences") <- seq_along(x) # temporary
  ## following lines are for trees that have been stripped of memory-intensive elements
  distances <- attr(tree, "distances")
  if(is.null(distances)) distances <- phylogram::mbed(x)
  attr(tree, "distances") <- NULL ## replaced later
  duplicates <- attr(tree, "duplicates")
  if(is.null(duplicates)) duplicates <- attr(distances, "duplicates")
  if(is.null(duplicates)) duplicates <- duplicated.DNAbin(x, point = TRUE)
  attr(tree, "duplicates") <- NULL ## replaced later
  pointers <- attr(tree, "pointers")
  if(is.null(pointers)) pointers <- attr(distances, "pointers")
  if(is.null(pointers)) pointers <- attr(duplicates, "pointers")
  attr(tree, "pointers") <- NULL ## replaced later
  hashes <- attr(tree, "hashes")
  if(is.null(hashes)) hashes <- attr(distances, "hashes")
  if(is.null(hashes)) hashes <- sapply(x, function(s) paste(openssl::md5(as.vector(s))))
  attr(tree, "hashes") <- NULL ## replaced later
  distances <- distances[!duplicates, ] # rm attrs, condensed version in final tree
  seqweights <- attr(tree, "weights")
  #ok if weights are null
  attr(tree, "weights") <- NULL ## replaced later
  has_duplicates <- any(duplicates)
  if(has_duplicates){
    fullseqset <- x
    #x <- x[!duplicates]
    x <- x[!duplicates]  # attributes not needed
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
    if(cores > navailcores) stop("Number of cores is more than the number available")
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
  ### recursively split nodes
  if(ncores > 1 & recursive){
    if(length(clades) < ncores){
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
      if(!quiet) cat("Recursively partitioning basal tree branches")
      repeat{
        nmembers <- integer(0)
        eligible <- logical(0)
        tmp <- fm1(tree)
        rm(tmp)
        nmembers <- nmembers[eligible]
        if(!any(eligible) | length(nmembers) >= 2 * ncores) break
        whichclade <- names(nmembers)[which.max(nmembers)]
        index <- gsub("([[:digit:]])", "[[\\1]]", whichclade)
        toeval <- paste0("tree", index, "<- fork(tree",
                         index, ", x, refine = refine, ",
                         "iterations = iterations, minK = 2, maxK = 2, ",
                         "minscore = minscore, probs = probs, resize = resize, ",
                         "distances = distances, seqweights = seqweights, ",
                         "cores = cores, quiet = quiet, ... = ...)")
        eval(parse(text = toeval))
        ss <- FALSE # split success; prevents build note due to lack of visible binding
        eval(parse(text = paste0("ss <- is.list(tree", index, ")")))
        # prevent multiple attempts to split the same node
        if(!ss) eval(parse(text = paste0("attr(tree", index, ", lock) <-TRUE")))
      }
      clades <- names(nmembers)
      indices <- gsub("([[:digit:]])", "[[\\1]]", clades)
      rm(nmembers)
      rm(eligible)
    }
    trees <- vector(mode = "list", length = length(clades)) # = length(indices)
    for(i in seq_along(indices)){
      eval(parse(text = paste0("trees[[", i, "]] <- tree", indices[i])))
    }
    if(!quiet) {
      cat("Recursively partitioning terminal tree branches\n")
      cat("Feedback suppressed\n")
      cat("This could take a while...\n")
    }
    trees <- parallel::parLapply(cores, trees, .learn1,
                                 x, refine = refine, iterations = iterations,
                                 minK = minK, maxK = maxK, minscore = minscore,
                                 probs = probs, resize = resize,
                                 distances = distances, # large matrix could cause probs
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
                       if(recursive) "<-.learn1(tree" else "<-fork(tree",
                       indices[i], ", x, refine = refine, ",
                       "iterations = iterations, minK = minK, maxK = maxK, ",
                       "minscore = minscore, probs = probs, resize = resize, ",
                       "distances = distances, seqweights = seqweights, ",
                       "cores = cores, quiet = quiet, ... = ...)")
      eval(parse(text = toeval))
    }
  }
  if(stopclustr) parallel::stopCluster(cores)
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
    akws <- attr(node, "Akweights")
    scrs <- attr(node, "scores")
    attr(node, "nunique") <- length(seqs)
    #newseqs <- which(pointers %in% seqs)
    newseqs <- newakws <- newscrs <- vector(mode = "list", length = length(seqs))
    for(i in seq_along(seqs)){
      newseqs[[i]] <- which(pointers == seqs[i])
      if(!is.null(akws)) newakws[[i]] <- rep(akws[i], length(newseqs[[i]]))
      if(!is.null(scrs)) newscrs[[i]] <- rep(scrs[i], length(newseqs[[i]]))
      # newakweights[newseqs == seqs[i]] <- akweights[i]
    }
    attr(node, "sequences") <- unlist(newseqs, use.names = FALSE)
    if(!is.null(akws)) attr(node, "Akweights") <- unlist(akws, use.names = FALSE)
    if(!is.null(scrs)) attr(node, "scores") <- unlist(scrs, use.names = FALSE)
    attr(node, "ntotal") <- length(newseqs)
    return(node)
  }
  if(!quiet) cat("Repatriating duplicate sequences with tree\n")
  tree <- dendrapply(tree, reduplicate, pointers = pointers)
  if(has_duplicates) x <- fullseqset
  attributes(x) <- tmpxattr
  if(!quiet) cat("Resetting node heights\n")
  tree <- phylogram::reposition(tree)
  if(!quiet) cat("Making tree ultrametric\n")
  tree <- phylogram::ultrametricize(tree)
  if(!quiet) cat("Labelling nodes\n")
  lineages <- gsub("\\.", "", attr(x, "lineage"))
  lineages <- paste0(lineages, "; ", attr(x, "species"))
  attachlins <- function(node, lineages){
    splitfun <- function(s) strsplit(s, split = "; ")[[1]]
    linvecs <- lapply(lineages[attr(node, "sequences")], splitfun)
    guide <- linvecs[[which.min(sapply(linvecs, length))]]
    counter <- 0
    for(l in guide){
      if(all(sapply(linvecs, function(e) l %in% e))) counter <- counter + 1
    }
    guide <- if(counter > 0) guide[1:counter] else character(0)
    lineage <- paste(guide, collapse = "; ")
    attr(node, "lineage") <- lineage
    attr(node, "label") <- paste0(guide[length(guide)], " (",
                                  attr(node, "nunique"), ",",
                                  attr(node, "ntotal"), ")")
    return(node)
  }
  tree <- dendrapply(tree, attachlins, lineages)
  attr(tree, "sequences") <- x # must happen after attaching lineages
  attr(tree, "distances") <- distances
  attr(tree, "duplicates") <- duplicates
  attr(tree, "pointers") <- pointers
  attr(tree, "weights") <- seqweights
  attr(tree, "hashes") <- hashes # length is length(x)
  if(!quiet) cat("Done\n")
  class(tree) <- c("insect", "dendrogram")
  return(tree)
}


#   attr(tree, "sequences") <- if(has_duplicates) fullseqset else x
#   label <- function(node, x){ # node id dendro, x is DNAbin
#     if(is.leaf(node)){
#       if(length(attr(node, "sequences")) > 1){
#         attr(node, "label") <- paste0(names(x)[attr(node, "sequences")[1]],
#                                       "...(", attr(node, "nunique"), ",",
#                                       attr(node, "ntotal"), ")")
#       }else{
#         attr(node, "label") <- names(x)[attr(node, "sequences")]
#       }
#     }
#     return(node)
#   }
#   if(!quiet) cat("Labeling leaf nodes\n")
#   tree <- dendrapply(tree, label, x = attr(tree, "sequences"))
#   if(!quiet) cat("Done\n")
#   return(tree)
# }
################################################################################
################################################################################
#' Contract a classification tree based on pattern matching.
#'
#' This function is used to collapse over-extended nodes.
#'
#' @param tree an object of class \code{"insect"}.
#' @param clades a vector of character strings giving the binary indices
#'   matching the labels of the nodes that are to be collapsed.
#' @param quiet logical indicating whether feedback should be printed
#'   to the console.
#' @return an object of class \code{"insect"}.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @seealso \code{\link{expand}}, \code{\link{purge}}
#' @examples
#'   ## TBA
################################################################################
contract <- function(tree, clades = "", quiet = FALSE){
  # tree is a hphmm
  # clades is a character vector eg c("22111", "22112")
  if(identical(clades, "")){
    if(!is.leaf(tree)){
      tmpattr <- attributes(tree)
      tree <- 1
      attributes(tree) <- tmpattr
      attr(tree, "height") <- 0
      attr(tree, "members") <- 1L
      attr(tree, "midpoint") <- 0
      attr(tree, "leaf") <- TRUE
    }
  }else{
    for(i in seq_along(clades)){
      clade_i <- as.numeric(strsplit(clades[i], split = "")[[1]])
      subtree <- tree
      for(j in clade_i) subtree <- subtree[[j]]
      if(is.list(subtree)){
        tmpattr <- attributes(subtree)
        subtree <- 1
        attributes(subtree) <- tmpattr
        attr(subtree, "height") <- 0
        attr(subtree, "members") <- 1L
        attr(subtree, "midpoint") <- NULL
        attr(subtree, "leaf") <- TRUE
        toeval <- paste0(paste0("tree[[",
                                paste0(clade_i, collapse = "]][["),
                                "]]"), " <- subtree")
        eval(parse(text = toeval))
      }else stop("Clade", i, "is a leaf or non-existent index\n")
    }
    if(!quiet) cat("Setting midpoints and members attributes\n")
    tree <- phylogram::remidpoint(tree)
    class(tree) <- "dendrogram"
    if(!quiet) cat("Resetting node heights\n")
    tree <- phylogram::reposition(tree)
    if(!quiet) cat("Making tree ultrametric\n")
    tree <- phylogram::ultrametricize(tree)
  }
  return(tree)
}
################################################################################
################################################################################
#' Find and collapse over-extended nodes
#'
#' This function tests each node of a classification tree
#'   for over-extension by checking if all sequences belonging
#'   to the node have identical lineage metadata. Over-extended
#'   nodes are collapsed and the simplified tree returned.
#'
#' @param tree an object of class \code{"insect"}.
#' @param quiet logical indicating whether feedback should be printed
#'   to the console.
#' @return an object of class \code{"insect"}.
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @seealso \code{\link{expand}}, \code{\link{contract}}
#' @examples
#'   ## TBA
################################################################################
purge <- function(tree, quiet = FALSE){
  purge1 <- function(node, lineages){
    if(is.list(node)){
      nodelins <- lineages[attr(node, "sequences")]
      if(all(nodelins == nodelins[1])){
        for(i in seq_along(node)) node[[i]] <- contract(node[[i]])
        #node <- contract(node)
      }
    }
    return(node)
  }

  purge2 <- function(node, lineages){
    node <- purge1(node, lineages)
    if(is.list(node)) node[] <- lapply(node, purge2, lineages)
    return(node)
  }
  x <- attr(tree, "sequences")
  attr(tree, "sequences") <- seq_along(x)
  if(!quiet) cat("Collapsing over-extended nodes\n")
  # lineages <- attr(x, "lineage")
  lineages <- gsub("\\.", "", attr(x, "lineage"))
  lineages <- paste0(lineages, "; ", attr(x, "species"))
  if(!(length(x) == length(lineages))) stop("Error 1")
  tree <- purge2(tree, lineages)
  if(all(lineages == lineages[1])) tree <- contract(tree)
  if(!quiet) cat("Setting midpoints and members attributes\n")
  tree <- phylogram::remidpoint(tree)
  #cat(class(tree), "\n") ##############
  class(tree) <- "dendrogram"
  if(!quiet) cat("Resetting node heights\n")
  tree <- phylogram::reposition(tree)
  if(!quiet) cat("Making tree ultrametric\n")
  tree <- phylogram::ultrametricize(tree)
  attr(tree, "sequences") <- x
  return(tree)
}
################################################################################


# purge1 <- function(node, lineages){
#   if(is.list(node)){
#     for(i in seq_along(node)){
#       if(is.list(node[[i]])){
#         nodelinsi <- lineages[attr(node[[i]], "sequences")]
#         if(all(nodelinsi == nodelinsi[1])) node[[i]] <- contract(node[[i]])
#       }
#     }
#   }
#   return(node)
# }

# purge <- function(tree, quiet = FALSE){
  # purge1 <- function(node, lineages){
  #   if(is.list(node)){
  #     nodelins <- lineages[attr(node, "sequences")]
  #     if(all(nodelins == nodelins[1])){
  #       #for(i in seq_along(node)) node[[i]] <- contract(node[[i]])
  #       node <- contract(node)
  #     }
  #   }
  #   return(node)
  # }
#   purge2 <- function(node, lineages){
#     node <- purge1(node, lineages)
#     if(is.list(node)) node[] <- lapply(node, purge2, lineages)
#     return(node)
#   }
#   x <- attr(tree, "sequences")
#   attr(tree, "sequences") <- seq_along(x)
#   if(!quiet) cat("Collapsing over-extended nodes\n")
#   lineages <- attr(x, "lineage")
#   tree <- purge2(tree, lineages)
#   if(!quiet) cat("Setting midpoints and members attributes\n")
#   tree <- phylogram::remidpoint(tree)
#   cat(class(tree), "\n") ##############
#   #class(tree) <- "dendrogram"
#   if(!quiet) cat("Resetting node heights\n")
#   tree <- phylogram::reposition(tree)
#   if(!quiet) cat("Making tree ultrametric\n")
#   tree <- phylogram::ultrametricize(tree)
#   attr(tree, "sequences") <- x
#   return(tree)
# }







# ### This was last right before returning tree in contract
# x <- attr(tree, "sequences") #full sequence set
# has_duplicates <- any(attr(tree, "duplicates"))
# if(has_duplicates) x <- x[!attr(tree, "duplicates")]
# label <- function(node, x){ # node is dendro, x is DNAbin
#   if(is.leaf(node)){
#     if(length(attr(node, "sequences")) > 1){
#       attr(node, "label") <- paste0(names(x)[attr(node, "sequences")[1]],
#                                     "...(", attr(node, "nunique"), ",",
#                                     attr(node, "ntotal"), ")")
#     }else{
#       attr(node, "label") <- names(x)[attr(node, "sequences")]
#     }
#   }
#   return(node)
# }
# if(!quiet) cat("Labeling leaf nodes\n")
# tree <- dendrapply(tree, label, x)
