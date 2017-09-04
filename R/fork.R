#' Split nodes of a classification tree.
#'
#' This function is used to split a leaf node of a classification tree,
#'   given a set of sequences.
#'
#' @param node an object of class \code{"insect"}.
#' @param x an object of class \code{"DNAbin"}.
#' @param lineages a character vector (the same length as x) of
#'   semi-colon delimited strings providing the lineage metadata
#'   for each sequence.
#' @param kmers an optional matrix of k-mer frequencies used for
#'   initial assignment of sequences to groups via k-means clustering.
#'   Rows should sum to one to account for differences in sequence length.
#'   Defaults to NULL, in which case k-mers are automatically
#'   counted (using k = 5) and normalized to sequence length.
#' @param seqweights either NULL (all sequences are given weights
#'   of 1), a numeric vector the same length as \code{x} representing
#'   the sequence weights used to derive the model, or a character string giving
#'   the method to derive the weights from the sequences. Currently only the
#'   \code{"Gerstein"} method is supported (default). For this method, a
#'   tree is first created by k-mer counting (see \code{\link[phylogram]{topdown}}),
#'   and sequence weights are then derived from the tree using the 'bottom up'
#'   algorithm of Gerstein et al (1994).
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
fork <- function(node, x, lineages, refine = "Viterbi", nstart = 10,
                 iterations = 50, minK = 2, maxK = 2, minscore = 0.9,
                 probs = 0.05, retry = TRUE, resize = TRUE, maxsize = NULL,
                 kmers = NULL, seqweights = "Gerstein", cores = 1,
                 quiet = FALSE, ...){
  #if(!is.list(node) & is.null(attr(node, "lock")) & is.null(attr(node, "onespp"))){
  indices <- attr(node, "sequences")
  nseq <- length(indices)
  lineages <- lineages[indices]
  if(!is.list(node) & is.null(attr(node, "lock")) & !all(lineages == lineages[1])){
    # fork leaves only
    # seqs <- x[attr(node, "sequences")]
    #lins <- lineages[attr(node, "sequences")]
    if(minK == 1 | nseq < minK){
      if(!is.null(attr(node, "model"))) attr(node, "model")$alignment <- NULL
      return(node)
    }
    if(nseq < maxK) maxK <- nseq
    if(is.null(seqweights)){
      seqweights <- rep(1, nseq)
    }else if(identical(seqweights, "Gerstein")){
      seqweights <- aphid::weight(x[indices], method = "Gerstein")
    }else if(length(seqweights) == length(x)){
      seqweights <- seqweights[indices]
      seqweights <- seqweights/mean(seqweights) ##scale weights to average 1
    }else stop("Invalid seqweights argument")
    if(is.null(kmers)){# generally will be NULL due to large size
      if(!quiet) cat("\nCounting kmers\n")
      kmers <- phylogram::kcount(x[indices], k = 5)
    }else{
      if(nrow(kmers) == length(x)) kmers <- kmers[indices, ]
    }
    stopifnot(nrow(kmers) == length(indices))
    ## set up multithread
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

    ## Split clade
    if(!quiet) cat("Attempting to split node", attr(node, "clade"), "\n")
    if(is.null(attr(node, "model"))){
      # should only happen at top level and if a PHMM is not provided
      mod <- NULL
    }else{
      mod <- attr(node, "model")
      attr(node, "model")$alignment <- NULL ## save on memory
      ins <- mod$inserts
      toosparse <- if(is.null(ins)) FALSE else sum(ins)/length(ins) > 0.5
      if(resize & toosparse & !quiet) cat("Skipping resize step\n")
      if(resize & attr(node, "clade") != "" & !toosparse){
        ## don't need to retrain top level model
        ## if mod was truncated then no need to resize
        if(!quiet) cat("Resizing parent model\n")
        alig <- mod$alignment
        if(is.null(alig)) alig <- aphid::align(x[indices], model = mod, cores = cores)
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
    alternative <- FALSE
    repeat{
      seqsplit <- partition(x[indices], model = mod, refine = refine, K = nclades,
                            allocation = "cluster", nstart = nstart,
                            iterations = iterations, kmers = kmers,
                            seqweights = seqweights, cores = cores, quiet = quiet,
                            ... = ...)
      if(is.null(seqsplit)){
        if(!quiet) cat("Sequence splitting failed, returning unsplit node\n")
        if(stopclustr) parallel::stopCluster(cores)
        #if(!is.null(attr(node, "model"))) attr(node, "model")$alignment <- NULL
        ## already done in previous block
        return(node)
      }
      # if(is.null(attr(node, "model"))){
      #   # should only be TRUE at top level- and even then not usually
      #   if(!quiet) cat("Assigning top-level model\n")
      #   attr(node, "model") <- seqsplit$model0
      # }
      # if(is.null(attr(node, "scores"))){ # should only be TRUE at top level
      #   if(!quiet) cat("Calculating top-level scores\n")
      #   fscore <- function(s, model) aphid::forward(model, s, odds = FALSE)$score
      #   attr(node, "scores") <- if(inherits(cores, "cluster")){
      #     parallel::parSapply(cores, x[indices], fscore, model = seqsplit$model0)
      #   }else{
      #     sapply(x[indices], fscore, model = seqsplit$model0)
      #   }
      # }
      membership <- seqsplit$membership
      scores <- seqsplit$scores
      total_scores <- apply(scores, 2, aphid::logsum)
      akwgts <- t(exp(t(scores) - total_scores))
      performances <- numeric(nseq)
      for(i in 1:nseq) performances[i] <- akwgts[membership[i], i]
      minperfs <- numeric(nclades)
      for(i in 1:nclades){
        minperfs[i] <- quantile(performances[membership == i], probs = probs)
        if(!quiet){
          if(i == 1) cat("Akaike weights:", performances, "\n")
          cat("Group", i, "size =", sum(membership == i), "\n")
          cat(sum(performances[membership == i] > minscore), "of",
              sum(membership == i), "correctly predicted with Akaike weight >",
              minscore, "\n")
          # cat(sum(performances[membership == i] > 0.5), "of",
          #     sum(membership == i), "correctly predicted with Akaike weight > 0.5\n")
          cat("Lower", probs, "quantile of Akaike weights:", minperfs[i], "\n")
          cat("Minimum Akaike weight:", min(performances[membership == i]), "\n")
        }
      }
      if(nclades == 2 & retry & min(minperfs) < 0.999){
        infocols <- apply(kmers, 2, function(v) length(unique(v))) == 2
        if(sum(infocols) > 50){
          if(!quiet) cat("Comparing result with alternative grouping method\n")
          hash <- function(v) paste(openssl::md5(as.raw(v == v[1])))
          hashes <- apply(kmers[, infocols], 2, hash)
          hfac <- factor(hashes)
          infocol <- match(levels(hfac)[which.max(tabulate(hfac))], hashes)
          alloc2 <- kmers[, infocols][, infocol]
          ## convert from logical to integer
          alloc2 <- unname(alloc2 == alloc2[1]) + 1
          if(!all(alloc2 == membership) & !all(alloc2 == seqsplit$init_membership)){
            seqsplit2 <- partition(x[indices], model = mod, refine = refine, K = 2,
                                   allocation = alloc2, nstart = nstart,
                                   iterations = iterations, kmers = kmers,
                                   seqweights = seqweights, cores = cores,
                                   quiet = quiet, ... = ...)
            if(!is.null(seqsplit2)){
              membership2 <- seqsplit2$membership
              if(!all(membership2 == membership)){
                scores2 <- seqsplit2$scores
                total_scores2 <- apply(scores2, 2, aphid::logsum)
                akwgts2 <- t(exp(t(scores2) - total_scores2))
                performances2 <- numeric(nseq)
                for(i in 1:nseq) performances2[i] <- akwgts2[membership2[i], i]
                minperfs2 <- numeric(nclades)
                for(i in 1:nclades){
                  minperfs2[i] <- quantile(performances2[membership2 == i], probs = probs)
                  if(!quiet){
                    if(i == 1) cat("Akaike weights:", performances2, "\n")
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
                if(quantile(performances2, probs = probs) > quantile(performances, probs = probs)){
                  if(!quiet) cat("Improvement found using alternative grouping\n")
                  seqsplit <- seqsplit2
                  membership <- membership2
                  scores <- scores2
                  total_scores <- total_scores2
                  akwgts <- akwgts2
                  performances <- performances2
                  minperfs <- minperfs2
                  alternative <- TRUE
                }else if(!quiet) cat("No improvement found using alternative grouping\n")
              }else if(!quiet) cat("Alternative grouping method gave same clusters\n")
            }else if(!quiet) cat("Alternative grouping gave null result\n")
          }else if(!quiet) cat("Alternative grouping method duplicated k-means\n")
        }else if(!quiet) cat("Insufficient informative columns for alternative grouping\n")
      }
      if(all(minperfs > minscore)){
        split_node <- TRUE
        break
      }else if(nclades == maxK){
        if(!quiet) cat("Minimum performance criteria not reached, unable to split clade\n")
        split_node <- FALSE
        break
      }else{
        nclades <- nclades + 1
        if(!quiet){
          cat("Minimum performance criteria not reached\n")
          cat("Attempting", nclades, "way split\n")
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
      attr(node, "alternative") <- alternative
      # onespecies <- all(lineages[indices] == lineages[indices[1]])
      for(i in 1:nclades){
        node[[i]] <- 1
        attr(node[[i]], "height") <- attr(node, "height") - 1
        attr(node[[i]], "leaf") <- TRUE
        #attr(node[[i]], "label") <- paste0(attr(node, "label"), i)
        attr(node[[i]], "clade") <- paste0(attr(node, "clade"), i)
        attr(node[[i]], "sequences") <- indices[membership == i]
        attr(node[[i]], "minscore") <- min(scores[i, membership == i]) - 1
        # above nominal amount just to ensure equality
        # attr(node[[i]], "Akweights") <- performances[membership == i]
        attr(node[[i]], "model") <- seqsplit[[paste0("model", i)]]
        # if(onespecies){
        #   attr(node[[i]], "onespp") <- TRUE
        #   attr(node[[i]], "lineage") <- lineages[indices[1]]
        # }else{
        #   attr(node[[i]], "lineage") <- .ancestor(lineages[indices[membership == i]])
        # }
        attr(node[[i]], "lineage") <- .ancestor(lineages[membership == i])
        ## prevents unnecessary recursion, but should have at least one split
        ## within each species
      }
    } # else{
    #   attr(node, "model")$alignment <- NULL
    # }
    attr(node, "model")$alignment <- NULL
  }
  return(node)
}
################################################################################
