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
#' @seealso \code{\link{contract}}, \code{\link{learn}}
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
  if(!quiet) cat("Setting midpoints and members attributes\n")
  #tree <- settreeattr(tree)
  tree <- phylogram::remidpoint(tree)
  class(tree) <- "dendrogram"
  add_duplicates <- function(node, pointers){
    attr(node, "nunique") <- length(attr(node, "sequences"))
    attr(node, "sequences") <- which(pointers %in% attr(node, "sequences"))
    attr(node, "ntotal") <- length(attr(node, "sequences"))
    return(node)
  }
  if(!quiet) cat("Repatriating duplicate sequences with tree\n")
  tree <- dendrapply(tree, add_duplicates,
                     pointers = attr(duplicates, "pointers"))
  # fix heights
  if(!quiet) cat("Resetting node heights\n")
  tree <- phylogram::reposition(tree)
  # make ultrametric
  if(!quiet) cat("Making tree ultrametric\n")
  tree <- phylogram::ultrametricize(tree)
  attr(tree, "sequences") <- if(has_duplicates) fullseqset else x
  label <- function(node, x){ # node id dendro, x is DNAbin
    if(is.leaf(node)){
      if(length(attr(node, "sequences")) > 1){
        attr(node, "label") <- paste0(names(x)[attr(node, "sequences")[1]],
                                      "...(", attr(node, "nunique"), ",",
                                      attr(node, "ntotal"), ")")
      }else{
        attr(node, "label") <- names(x)[attr(node, "sequences")]
      }
    }
    return(node)
  }
  if(!quiet) cat("Labeling leaf nodes\n")
  tree <- dendrapply(tree, label, x = attr(tree, "sequences"))
  if(!quiet) cat("Done\n")
  return(tree)
}
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
#' @seealso \code{\link{expand}}
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
    attr(subtree, "midpoint") <- 0
    attr(subtree, "leaf") <- TRUE
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
        # attr(subtree, "label") <- paste0(names(attr(tree, "sequences"))[attr(subtree, "sequences")[1]],
        #                                  "...(", length(attr(subtree, "sequences")), ")")
        toeval <- paste0(paste0("tree[[", paste0(clade_i, collapse = "]][["), "]]"), " <- subtree")
        eval(parse(text = toeval))
      }else stop("Clade", i, "is a leaf\n")
    }
    if(!quiet) cat("Setting midpoints and members attributes\n")
    tree <- phylogram::remidpoint(tree)
    class(tree) <- "dendrogram"
    # fix heights
    if(!quiet) cat("Resetting node heights\n")
    tree <- phylogram::reposition(tree)
    # make ultrametric
    if(!quiet) cat("Making tree ultrametric\n")
    tree <- phylogram::ultrametricize(tree)
    label <- function(node, x){ # node is dendro, x is DNAbin
      if(is.leaf(node)){
        if(length(attr(node, "sequences")) > 1){
          attr(node, "label") <- paste0(names(x)[attr(node, "sequences")[1]],
                                        "...(", attr(node, "nunique"), ",",
                                        attr(node, "ntotal"), ")")
        }else{
          attr(node, "label") <- names(x)[attr(node, "sequences")]
        }
      }
      return(node)
    }
    if(!quiet) cat("Labeling leaf nodes\n")
    x <- attr(tree, "sequences") #full sequence set
    has_duplicates <- any(attr(tree, "duplicates"))
    if(has_duplicates) x <- x[!attr(tree, "duplicates")]
    tree <- dendrapply(tree, label, x)
  }
  return(tree)
}
################################################################################
