#' Allocate sequences for cross validation by identity.
#'
#' This function takes a reference sequence database and allocates each sequence to either
#'   a query set (a.k.a. test set) or a training set, in order to cross validate a
#'   supervised taxon classifier. The method is based on that of Edgar (2018), but uses
#'   recursive divisive clustering and retains all sequences rather than discarding those that
#'   violate the top-hit identity constraint.
#'
#' @param x a set of reference sequences. Can be a
#'   "DNAbin" object or a named vector of upper-case DNA character strings.
#' @param threshold numeric between 0 and 1 giving the identity threshold for sequence allocation.
#' @param allocate character giving the method to use to allocate eligible sequences to the query set.
#'   Options are "max" (default) which chooses the largest node from each pair in order to
#'   maximize the size of the query set,
#'   or "sample", which randomly chooses one node from each eligible pair.
#' @param ... further arguments to pass to "kmeans"
#' @return a logical vector the same length as the input object, indicating which sequences
#'   should be allocated to the query set
#' @author Shaun Wilkinson
#' @references
#'   Edgar RC (2018) Accuracy of taxonomy prediction for 16S rRNA and fungal ITS sequences.
#'   PeerJ 6:e4652. DOI 10.7717/peerj.4652
#' @examples
#' \donttest{
#'   data(whales)
#'   allocateCVI(whales)
#' }
################################################################################
allocateCVI <- function(x, threshold = 0.9, allocate = "max", ...){
  if(mode(x) == "character") x <- char2dna(x, simplify = FALSE)
  cluster <- function(x, threshold = 0.9, ...){ ## outputs a dendrogram
    if(any(duplicated(names(x)))) stop("Sequence names must be unique\n")
    hashes <- hash(x)
    pointers <- .point(hashes)
    x <- x[!duplicated(hashes)]
    xlengths <- vapply(x, length, 0L)
    kmers <- kmer::kcount(x, k = 4)
    kmers <- kmers/(apply(kmers, 1, sum) + 4 - 1)
    tree <- 1
    attr(tree, "leaf") <- TRUE
    attr(tree, "sequences") <- seq_along(x)
    attr(tree, "height") <- 0
    otun <- function(node, ...){
      if(!is.list(node)){ ## fork leaves only
        if(length(attr(node, "sequences")) > 1){
          kmers <- kmers[attr(node, "sequences"), , drop = FALSE] # up 1 envir
          lens <- xlengths[attr(node, "sequences")] # up 1 envir
          centroid <- apply(kmers, 2, mean)
          dists2centroid <- apply((t(t(kmers) - centroid))^2, 1, sum)
          cseq <- which.min(dists2centroid) # index
          attr(node, "central") <- cseq
          dists2central <- apply((t(kmers) - kmers[cseq, ])^2, 2, sum) * lens/(2*4)
          if(max(dists2central) < 1 - threshold) return(node)
          km <- if(nrow(kmers) > 2){
            tryCatch(kmeans(kmers, centers = 2, ... = ...),
                     error = function(er) return(NULL),
                     warning = function(wa) return(NULL))
          }else{
            list(cluster = 1:2)
          }
          if(is.null(km)) return(node)
          tmpattr <- attributes(node)
          node <- vector(mode = "list", length = 2)
          attributes(node) <- tmpattr
          attr(node, "leaf") <- NULL
          for(i in 1:2){
            node[[i]] <- 1
            attr(node[[i]], "height") <- attr(node, "height") - 1
            attr(node[[i]], "leaf") <- TRUE
            attr(node[[i]], "sequences") <- attr(node, "sequences")[km$cluster == i]
          }
        }else{
          attr(node, "central") <- 1L
        }
      }
      return(node)
    }
    otur <- function(tree, ...){
      tree <- otun(tree, ...)
      if(is.list(tree)) tree[] <- lapply(tree, otur, ...)
      return(tree)
    }
    tree <- otur(tree, ... = ...)
    class(tree) <- "dendrogram"
    reduplicate <- function(node, pointers){
      cseq <- attr(node, "sequences")[attr(node, "central")] # index in derepd set
      cseq <- match(cseq, pointers) # index in rerepd set
      attr(node, "sequences") <- which(pointers %in% attr(node, "sequences")) #indices in rerepd set
      attr(node, "central") <- match(cseq, attr(node, "sequences")) #index in nodeseq vect
      stopifnot(!is.na(attr(node, "central")))
      return(node)
    }
    tree <- dendrapply(tree, reduplicate, pointers)
    tree <- phylogram::remidpoint(tree)
    tree <- phylogram::reposition(tree)
    return(tree)
  }
  tree <- cluster(x, threshold = threshold, ... = ...) # cluster within predefined threshold
  res <- logical(length(x))# F or T = training or test
  assign_test <- function(node, allocate){
    if(all(vapply(node, is.leaf, logical(1)))){# only apply to level 1 nodes
      nodesizes <- vapply(node, function(l) length(attr(l, "sequences")), 0L)
      if(allocate == "max"){
        whichnode <- which.max(nodesizes)
      }else if(allocate == "sample"){
        whichnode <- sample(1:2, size = 1)
      }
      res[attr(node[[whichnode]], "sequences")] <<- TRUE
    }
    return(node)
  }
  tmp <- dendrapply(tree, assign_test, allocate)
  rm(tmp)
  gc()
  return(res)
}
################################################################################
