#' Split nodes of a classification tree.
#'
#' This function is used to split a leaf node of a classification tree,
#'   given a set of sequences.
#'
#' @param node an object of class \code{"insect"}.
#' @param x an object of class \code{"DNAbin"}.
#' @param kmers an optional matrix of k-mer frequencies used for
#'   initial assignment of sequences to groups via k-means clustering.
#'   Rows should sum to one to account for differences in sequence length.
#'   Defaults to NULL, in which case k-mers are automatically
#'   counted (using k = 4) and normalized to sequence length.
#' @inheritParams learn
#' @return an object of class \code{"insect"}.
#' @details Note that seqweights argument should have the same length as x.
#' @author Shaun Wilkinson
#' @references
#'   Blackshields G, Sievers F, Shi W, Wilm A, Higgins DG (2010) Sequence embedding
#'   for fast construction of guide trees for multiple sequence alignment.
#'   \emph{Algorithms for Molecular Biology}, \strong{5}, 21.
#' @seealso \code{\link{learn}} (parent function),
#'   \code{\link{partition}} (child function).
#' @examples
#'   ## TBA
################################################################################
fork <- function(node, x, refine = "Viterbi", nstart = 10, iterations = 50,
                 minK = 2, maxK = 2, minscore = 0.9, probs = 0.05,
                 resize = TRUE, maxsize = NULL, kmers = NULL,
                 seqweights = "Gerstein", cores = 1, quiet = FALSE, ...){
  if(!is.list(node) & is.null(attr(node, "lock"))){ # fork leaves only
    seqs <- x[attr(node, "sequences")]
    nseq <- length(seqs)
    if(minK == 1 | nseq < minK){
      if(!is.null(attr(node, "phmm"))) attr(node, "phmm")$alignment <- NULL
      return(node)
    }
    if(nseq < maxK) maxK <- nseq
    if(is.null(seqweights)){
      seqweights <- rep(1, nseq)
    }else if(identical(seqweights, "Gerstein")){
      seqweights <- aphid::weight(x[attr(node, "sequences")], method = "Gerstein")
    }else if(length(seqweights) == length(x)){
      seqweights <- seqweights[attr(node, "sequences")]
      ### scale weights to average 1
      seqweights <- seqweights/mean(seqweights)
    }else stop("Invalid seqweights argument")
    if(!is.null(kmers)){
      if(nrow(kmers) == length(x)){
        kmers <- kmers[attr(node, "sequences"), ]
      }else{
        stopifnot(nrow(kmers) == length(attr(node, "sequences")))
      }
    }

    # if(is.null(kmers)){
    #   # distances <- phylogram::mbed(x[attr(node, "sequences")])
    #   if(!quiet) cat("Counting k-mers\n")
    #   kmers <- phylogram::kcount(seqs, k = 5)/(sapply(seqs, length) - 4)#k-1=3
    # # }else if(nrow(distances) == length(x)){
    #   # distances <- distances[attr(node, "sequences"), ]
    # }else if(nrow(kmers) == length(x)){
    #   kmers <- kmers[attr(node, "sequences"), ]
    # }else stop("Invalid kmers argument")

    ### set up multithread
    if(inherits(cores, "cluster") | identical(cores, 1)){
      stopclustr <- FALSE
    }else{ # create cluster object
      navailcores <- parallel::detectCores()
      if(identical(cores, "autodetect")) cores <- navailcores - 1
      if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
      if(cores > navailcores) stop("Number of cores is more than number available")
      # if(!quiet) cat("Multithreading over", cores, "cores\n")
      if(cores == 1){
        stopclustr <- FALSE
      }else{
        if(!quiet) cat("Initializing cluster with", cores, "cores\n")
        cores <- parallel::makeCluster(cores)
        stopclustr <- TRUE
      }
    }
    ### Split clade
    if(!quiet) cat("\nAttempting to split node", attr(node, "clade"), "\n")
    if(is.null(attr(node, "phmm"))){
      # should only happen at top level and if a PHMM is not provided
      mod <- NULL
    }else{
      mod <- attr(node, "phmm")
      attr(node, "phmm")$alignment <- NULL ## save on memory
      ins <- mod$inserts
      toosparse <- if(is.null(ins)) FALSE else sum(ins)/length(ins) > 0.5
      if(resize & toosparse & !quiet) cat("Skipping resize step\n")
      if(resize & attr(node, "clade") != "" & !toosparse){
        #don't need to retrain top level model
        # if mod was truncated then no need to resize
        if(!quiet) cat("Resizing parent model\n")
        ### model is allowed to change size here
        alig <- mod$alignment
        if(is.null(alig)) alig <- aphid::align(seqs, model = mod, cores = cores)
        mod <- aphid::derivePHMM.default(alig, seqweights = seqweights,
                                         inserts = if(nseq < 1000) "map" else "threshold",
                                         maxsize = maxsize)
        rm(alig)
        gc()
        if(!quiet) cat("New model size :", mod$size, "\n")
      }
    }
    split_node <- FALSE
    nclades <- minK
    allocation <- "cluster"
    repeat{
      seqsplit <- partition(seqs, model = mod, refine = refine, K = nclades,
                            allocation = allocation, nstart = nstart,
                            iterations = iterations, kmers = kmers,
                            seqweights = seqweights, cores = cores, quiet = quiet,
                            ... = ...)
      if(is.null(seqsplit)){
        if(!quiet) cat("Sequence splitting failed, returning unsplit node\n")
        if(stopclustr) parallel::stopCluster(cores)
        if(!is.null(attr(node, "phmm"))) attr(node, "phmm")$alignment <- NULL
        return(node)
      }
      if(is.null(attr(node, "phmm"))){ # should only be TRUE at top level
        if(!quiet) cat("Assigning top-level model\n")
        attr(node, "phmm") <- seqsplit$phmm0
      }
      if(is.null(attr(node, "scores"))){ # should only be TRUE at top level
        if(!quiet) cat("Calculating top-level scores\n")
        fscore <- function(s, model) aphid::forward(model, s, odds = FALSE)$score
        attr(node, "scores") <- if(inherits(cores, "cluster")){
          parallel::parSapply(cores, seqs, fscore, model = seqsplit$phmm0)
        }else{
          sapply(seqs, fscore, model = seqsplit$phmm0)
        }
      }
      membership <- seqsplit$membership
      scores <- seqsplit$scores
      total_scores <- apply(scores, 2, aphid::logsum)
      akwgts <- t(exp(t(scores) - total_scores))
      performances <- numeric(nseq)
      for(i in 1:nseq) performances[i] <- akwgts[membership[i], i]
      if(!quiet) cat("Akaike weights:", performances, "\n")
      minperfs <- numeric(nclades)
      for(i in 1:nclades){
        minperfs[i] <- quantile(performances[membership == i], probs = probs)
        if(!quiet){
          cat("Group", i, "size =", sum(membership == i), "\n")
          cat(sum(performances[membership == i] > minscore), "of",
              sum(membership == i), "correctly predicted with Akaike weight >",
              minscore, "\n")
          cat(sum(performances[membership == i] > 0.5), "of",
              sum(membership == i), "correctly predicted with Akaike weight > 0.5\n")
          cat("Lower", probs, "quantile of Akaike weights:", minperfs[i], "\n")
          cat("Minimum Akaike weight:", min(performances[membership == i]), "\n")
        }
      }
      if(all(minperfs > minscore)){
        split_node <- TRUE
        break
      }else if(nclades == maxK){
        if(!quiet) cat("Minimum performance criteria not reached, unable to split clade\n")
        split_node <- FALSE
        break
      }else{
        if(nclades == 2 & identical(allocation, "cluster")){
          allocation <- "split"
          if(!quiet){
            cat("Minimum performance criteria not reached\n")
            cat("Trying alternative initial grouping method\n")
          }
        }else{
          nclades <- nclades + 1
          allocation <- "cluster"
          if(!quiet){
            cat("Minimum performance criteria not reached\n")
            cat("Attempting", nclades, "way split\n")
          }
        }
      }
    }
    if(stopclustr) parallel::stopCluster(cores) # not used if called from 'learn'
    # Akaike weights should all be close to 1
    # now decorate the node (if disc ability is > a certain threshold?)
    # splitfun <- function(s) strsplit(s, split = ";")[[1]]
    if(split_node){# placeholder for discriminant thresholding
      # change from leaf to inner node
      if(!quiet) cat("Splitting node", attr(node, "clade"), "\n")
      tmpattr <- attributes(node)
      node <- vector(mode = "list", length = nclades)
      attributes(node) <- tmpattr
      attr(node, "leaf") <- NULL
      # attr(node, "Akaike") <- akwgts #(nclades x nseq matrix of Akaike weights for the new sub-models)
      for(i in 1:nclades){
        node[[i]] <- 1
        attr(node[[i]], "height") <- attr(node, "height") - 1
        attr(node[[i]], "leaf") <- TRUE
        #attr(node[[i]], "label") <- paste0(attr(node, "label"), i)
        attr(node[[i]], "clade") <- paste0(attr(node, "clade"), i)
        attr(node[[i]], "sequences") <- attr(node, "sequences")[membership == i]
        attr(node[[i]], "scores") <- scores[i, membership == i]
        attr(node[[i]], "Akweights") <- performances[membership == i]
        attr(node[[i]], "phmm") <- seqsplit[[paste0("phmm", i)]]
      }
    }else{
      attr(node, "phmm")$alignment <- NULL
    }
  }
  return(node)
}
################################################################################
