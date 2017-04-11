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
                   resize = TRUE, recursive = TRUE, quiet = FALSE, ...){
  x <- attr(tree, "sequences") #full sequence set
  attr(tree, "sequences") <- seq_along(x)
  seqweights <- attr(tree, "weights")
  duplicates <- attr(tree, "duplicates")
  has_duplicates <- any(duplicates)
  if(has_duplicates){
    fullseqset <- x #includes attributes
    #x <- x[!duplicates]
    x <- subset(x, subset = !duplicates)  #includes attributes
    rmduplicates <- function(node, whchunq, pointers){
      tmp <- attr(node, "sequences")[attr(node, "sequences") %in% whchunq]
      attr(node, "sequences") <- pointers[tmp]
      return(node)
    }
    tree <- dendrapply(tree, rmduplicates, whchunq = which(!duplicates),
                       pointers = attr(duplicates, "pointers"))
  }
  if(identical(clades, "")){
    if(recursive){
      tree <- .learn1(tree, x, refine = refine, iterations = iterations,
                      minK = minK, maxK = maxK, minscore = minscore, probs = probs,
                      resize = resize, seqweights = seqweights,
                      quiet = quiet, ... = ...)
    }else{
      tree <- fork(tree, x, refine = refine, iterations = iterations,
                   minK = minK, maxK = maxK, minscore = minscore, probs = probs,
                   resize = resize, seqweights = seqweights,
                   quiet = quiet, ... = ...)
    }
  }else{
    for(i in seq_along(clades)){
      clade_i <- as.numeric(strsplit(clades[i], split = "")[[1]])
      index <- paste0(paste0("tree[[", paste0(clade_i, collapse = "]][["), "]]"))
      toeval <- paste0(index, if(recursive) "<-.learn1(" else "<-fork(",
                       index, ", x, refine = refine, ",
                       "iterations = iterations, minK = minK, maxK = maxK, ",
                       "minscore = minscore, probs = probs, resize = resize, ",
                       "seqweights = seqweights, quiet = quiet, ... = ...)")
      eval(parse(text = toeval))
    }
  }
  ### fix midpoints, members, heights and leaf integers
  ### note changes here also apply to 'learn' function
  if(!quiet) cat("Setting midpoints and members attributes\n")
  #tree <- settreeattr(tree)
  tree <- phylogram::remidpoint(tree)
  class(tree) <- "dendrogram"
  add_duplicates <- function(node, pointers){
    seqs <- attr(node, "sequences")
    akws <- attr(node, "Akweights")
    attr(node, "nunique") <- length(seqs)
    #newseqs <- which(pointers %in% seqs)
    newseqs <- vector(mode = "list", length = length(seqs))
    newakws <- vector(mode = "list", length = length(seqs))
    for(i in seq_along(seqs)){
      newseqs[[i]] <- which(pointers == seqs[i])
      if(!is.null(akws)) newakws[[i]] <- rep(akws[i], length(newseqs[[i]]))
      # newakweights[newseqs == seqs[i]] <- akweights[i]
    }
    attr(node, "sequences") <- unlist(newseqs, use.names = FALSE)
    if(!is.null(akws)) attr(node, "Akweights") <- unlist(akws, use.names = FALSE)
    attr(node, "ntotal") <- length(newseqs)
    return(node)
  }
  if(!quiet) cat("Repatriating duplicate sequences with tree\n")
  tree <- dendrapply(tree, add_duplicates,
                     pointers = attr(duplicates, "pointers"))
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
  attr(tree, "duplicates") <- duplicates
  attr(tree, "pointers") <- attr(duplicates, "pointers")
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
    tmpattr <- attributes(tree)
    tree <- 1
    attributes(tree) <- tmpattr
    attr(tree, "height") <- 0
    attr(tree, "members") <- 1L
    attr(tree, "midpoint") <- 0
    attr(tree, "leaf") <- TRUE
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
        node <- contract(node)
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
  lineages <- attr(x, "lineage")
  tree <- purge2(tree, lineages)
  if(!quiet) cat("Setting midpoints and members attributes\n")
  tree <- phylogram::remidpoint(tree)
  cat(class(tree), "\n") ##############
  #class(tree) <- "dendrogram"
  if(!quiet) cat("Resetting node heights\n")
  tree <- phylogram::reposition(tree)
  if(!quiet) cat("Making tree ultrametric\n")
  tree <- phylogram::ultrametricize(tree)
  attr(tree, "sequences") <- x
  return(tree)
}
################################################################################












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
