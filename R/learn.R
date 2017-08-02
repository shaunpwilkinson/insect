#' Informatic sequence classification tree learning.
#'
#' This function creates classification trees using recursive
#'   partitioning and nested profile hidden Markov models.
#'
#' @param x an object of class\code{"DNAbin"} representing a list of
#'   DNA sequences to be used as the traning data for the tree-learning process.
#'   All sequences should be from the same genetic region of interest
#'   and be globally alignable (i.e. without unjustified end-gaps).
#' @param model an optional object of class \code{"PHMM"} to form the model
#'   at the root node of the classification tree. Used to train (optimize
#'   parameters for) subsequent nested models to be positioned at successive
#'   sub-nodes. If NULL, the root model is derived \emph{de-novo} from the
#'   sequence list prior to the recursive partitioning process.
#' @param refine character string giving the iterative model refinement
#'   method to be used in the partitioning process. Valid options are
#'   \code{"Viterbi"} (Viterbi training; the default option) and
#'   \code{"BaumWelch"} (a modified version of the Expectation-Maximization
#'   algorithm).
#' @param iterations integer giving the maximum number of training-classification
#'   iterations to be used in the splitting process (see details section).
#'   Note that this is not necessarily the same as the number of Viterbi training
#'   or Baum Welch iterations to be used in model training, which can be set
#'   using the argument \code{"maxiter"} (eventually passed to
#'   \code{\link[aphid]{train}} via the dots argument "...").
#' @param nstart integer. The number of random starting sets to be chosen
#'   for initial k-means assignment of sequences to groups.
#' @param minK integer. The minimum number of furications allowed at each inner
#'   node of the tree. Defaults to 2 (all inner nodes are bifuricating).
#' @param maxK integer. The maximum number of furications allowed at each inner
#'   node of the tree. Defaults to 2 (all inner nodes are bifuricating).
#' @param minscore numeric between 0 and 1. The minimum acceptable value
#'   for the \emph{n}th percentile of Akaike weights (where \emph{n} is
#'   the value given in \code{"probs"}, for the node to be split and the
#'   recursion process to continue.
#'   At any given node, if the \emph{n}th percentile of Akaike weights
#'   falls below this threshold, the recursion process for the node will
#'   terminate. As an extreme example, if \code{minscore = 0.9} and
#'   \code{probs = 0.05} (the default settings), and after generating two
#'   candidate PHMMs to occupy the candidate subnodes the lower 5th percentile
#'   of Akaike weights is 0.89, the splitting process will
#'   terminate and the function will simply return the unsplit root node.
#' @param probs numeric between 0 and 1. The percentile of Akaike weights
#'   to test against the minimum score threshold given in \code{"minscore"}.
#' @param resize logical indicating whether the models should be free to
#'   change size during the training process or if the number of modules
#'   should be fixed. Defaults to TRUE. Only applicable if
#'   \code{refine = "Viterbi"}.
#' @param maxsize integer giving the upper bound on the number of modules
#'   in the PHMMs. If NULL (default) no maximum size is enforced.
#' @param seqweights an optional numeric vector the same length as the
#'   sequence list giving the weights to use when training the
#'   PHMMs. Alternatively its length can be the same as the total number
#'   of unique sequences as specified in \code{duplicates}.
#'   The character string "Gerstein" is also valid, which calls
#'   the \code{\link[aphid]{weight}} function in the
#'   \code{\link[aphid]{aphid}} package to derive the weights using the
#'   algorithm of Gerstein et al. 1994.
#'   If NULL, all sequences are given weights of 1.
#' @param recursive logical indicating whether the splitting process
#'   should continue recursively until the discrimination criteria
#'   are not met (TRUE; default), or whether a single split should
#'   take place at the root node.
#' @param cores integer giving the number of CPUs to use
#'   when training the models (only applicable if
#'   \code{refine = 'Viterbi'}). Defaults to 1.
#'   This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   e.g. by running \code{parallel::stopCluster(cores)}.
#'   The string 'autodetect' is also accepted, in which case the maximum
#'   number of cores to use is one less than the total number of cores
#'   available.
#' @param quiet logical indicating whether feedback should be printed
#'   to the console. Note that the output can be rather verbose for
#'   larger trees.
#' @param ... further arguments to be passed to \code{\link{partition}}
#'   (which may then be passed on to \code{\link[aphid]{train}}).
#' @return an object of class \code{"insect"}. This object is a
#'   \code{"dendrogram"} with several additional attributes including:
#'   "sequences" the sequences attribute of the root node is the
#'   original object \code{x}, while the sequences attributes of all
#'   subnodes are indices pointing to which sequences in \code{x}
#'   belong to the node (this saves on memory).
#'   "weights" (root node only) a numeric vector giving the weights of the
#'   \emph{unique} sequences in \code{x} (as indicated by
#'   \code{duplicates}).
#'   "duplicates" (root node only) a logical vector the same length as
#'   \code{x} indicating which sequences are duplicated.
#'   "pointers" (root node only) integer vector the same length as
#'   \code{x} indicating which sequence in \code{x} each sequence
#'   is a duplicate of. If there are no duplicates this is simply
#'   \code{1:length(x)}.
#'
#' @details
#'   TBA.
#' @author Shaun Wilkinson
#' @references
#'   Blackshields G, Sievers F, Shi W, Wilm A, Higgins DG (2010) Sequence embedding
#'   for fast construction of guide trees for multiple sequence alignment.
#'   \emph{Algorithms for Molecular Biology}, \strong{5}, 21.
#'
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Gerstein M, Sonnhammer ELL, Chothia C (1994) Volume changes in protein evolution.
#'   \emph{Journal of Molecular Biology}, \strong{236}, 1067-1078.
#'
#'   Juang B-H, Rabiner LR (1990) The segmental K-means
#'   algorithm for estimating parameters of hidden Markov models.
#'   \emph{IEEE Transactions on Acoustics, Speech, and Signal Processing},
#'   \strong{38}, 1639-1641.
#'
#' @seealso \code{\link{partition}}
#' @examples
#'   ## TBA
################################################################################
learn <- function(x, model = NULL, refine = "Viterbi", iterations = 50,
                  nstart = 10, minK = 2, maxK = 2, minscore = 0.9, probs = 0.1,
                  resize = TRUE, maxsize = NULL, seqweights = "Gerstein",
                  recursive = TRUE, cores = 1, quiet = FALSE, ...){
  # x is a "DNAbin" object
  tmpxattr <- attributes(x)
  nseq <- length(x)
  x <- x[] # removes attributes
  if(is.null(seqweights)) seqweights <- rep(1, length(x))
  # First initialize the tree
  tree <- 1
  attr(tree, "clade") <- ""
  attr(tree, "leaf") <- TRUE
  #distances <- phylogram::mbed(x)
  #duplicates <- attr(distances, "duplicates")
  #pointers <- attr(distances, "pointers")
  #hashes <- attr(distances, "hashes")
  if(!quiet) cat("Finding duplicates\n")
  hashes <- .digest(x, simplify = TRUE)
  duplicates <- duplicated(hashes) # logical length x
  pointers <- integer(length(x))
  dhashes <- hashes[duplicates]
  uhashes <- hashes[!duplicates]
  pointers[!duplicates] <- seq_along(uhashes)
  pd <- integer(length(dhashes))
  for(i in unique(dhashes)) pd[dhashes == i] <- match(i, uhashes)
  pointers[duplicates] <- pd
  #distances <- distances[!duplicates, ]
  nuseq <- sum(!duplicates)
  has_duplicates <- any(duplicates)
  if(has_duplicates){
    fullseqset <- x #excludes attributes
    # x <- x[!duplicates]
    x <- x[!duplicates] #subset.DNAbin(x, subset = !duplicates)  #exc attributes
    if(length(seqweights) == nseq) seqweights <- seqweights[!duplicates]
  }
  if(identical(seqweights, "Gerstein")){
    if(!quiet) cat("Deriving sequence weights\n")
    seqweights <- aphid::weight(x, "Gerstein")
  }
  attr(tree, "sequences") <- seq_along(x) # tmp-eventually replaced by DNAbin
  stopifnot(length(seqweights) == length(x))
  ### integer vector of indices pointing to x arg
  attr(tree, "phmm") <- model
  attr(tree, "height") <- 0
  ### set up multithread
  if(inherits(cores, "cluster")){
    ncores <- length(cores)
    stopclustr <- FALSE
  }else if(identical(cores, 1)){
    ncores <- 1
    stopclustr <- FALSE
  }else{ # create cluster object
    navailcores <- parallel::detectCores() # relatively costly ~ 0.06s
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
  if(ncores == 1){
    if(!quiet) cat("Counting k-mers\n")
    ## kmers are actually frequencies not counts
    ## could offer option to specify k here eventually but 4 is ok for now
    kmers <- phylogram::kcount(x, k = 5)/(sapply(x, length) - 4) #k-1=4
  }else kmers <- NULL
  if(recursive){
    if(!quiet) cat("Learning tree\n")
    if(ncores == 1){
      tree <- .learn1(tree, x = x, refine = refine, iterations = iterations,
                      nstart = nstart, minK = minK, maxK = maxK, minscore = minscore,
                      probs = probs, resize = resize, maxsize = maxsize,
                      kmers = kmers, seqweights = seqweights, cores = cores,
                      quiet = quiet, ... = ...)
    }else{
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
      if(!quiet) cat("Recursively partitioning basal tree branches\n")
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
        toeval <- paste0("tree", index, "<- fork(tree",
                         index, ", x, refine = refine, ",
                         "iterations = iterations, nstart = nstart, minK = minK, maxK = maxK, ",
                         "minscore = minscore, probs = probs, resize = resize, ",
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
      trees <- vector(mode = "list", length = length(clades))
      for(i in seq_along(indices)){
        eval(parse(text = paste0("trees[[", i, "]] <- tree", indices[i])))
      }
      if(!quiet) {
        cat("Recursively partitioning terminal tree branches\n")
        cat("Feedback suppressed, this could take a while...\n")
      }
      trees <- parallel::parLapply(cores, trees, .learn1,
                                   x, refine = refine, iterations = iterations, nstart = nstart,
                                   minK = minK, maxK = maxK, minscore = minscore,
                                   probs = probs, resize = resize, maxsize = maxsize,
                                   kmers = kmers, seqweights = seqweights,
                                   cores = 1, quiet = TRUE, ... = ...)
      for(i in seq_along(trees)){
        eval(parse(text = paste0("tree", indices[i], "<- trees[[", i, "]]")))
        trees[[i]] <- NA
        gc()
      }
    }
  }else{
    tree <- fork(tree, x = x, refine = refine, iterations = iterations,
                 minK = minK, maxK = maxK, minscore = minscore,
                 probs = probs, resize = resize, maxsize = maxsize,
                 kmers = kmers, seqweights = seqweights, cores = cores,
                 quiet = quiet, ... = ...)
  }
  if(stopclustr) parallel::stopCluster(cores)
  ### remove kmers since can be memory hungry, prevent next operations
  rm(kmers)
  gc()
  ### fix midpoints, members, heights and leaf integers
  ### note changes here also apply to 'expand' function
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
    attr(node, "ntotal") <- length(attr(node, "sequences"))
    return(node)
  }
  if(!quiet) cat("Repatriating duplicate sequences with tree\n")
  tree <- dendrapply(tree, reduplicate, pointers)
  if(has_duplicates) x <- fullseqset
  attributes(x) <- tmpxattr
  if(!quiet) cat("Resetting node heights\n")
  tree <- phylogram::reposition(tree)
  if(!quiet) cat("Making tree ultrametric\n")
  tree <- phylogram::ultrametricize(tree)
  if(!quiet) cat("Labelling nodes\n")
  lineages <- gsub("\\.", "", attr(x, "lineage"))
  lineages <- paste0(lineages, "; ~", attr(x, "species"), "~")
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
  attr(tree, "duplicates") <- duplicates # length is length(x)
  attr(tree, "pointers") <- pointers # length is length(x)
  #attr(tree, "kmers") <- kmers # number of rows is number of unique seqs
  attr(tree, "weights") <- seqweights # length is number of unique seqs
  attr(tree, "hashes") <- hashes # length is length(x)
  #attr(tree, "indices") <- .reindex(tree)
  if(!quiet) cat("Done\n")
  class(tree) <- c("insect", "dendrogram")
  return(tree)
}
################################################################################

# lockleaves <- function(node, exceptions){
#   if(!is.list(node)){
#     if(!attr(node, "clade") %in% exceptions) attr(node, "lock") <- TRUE
#   }
#   return(node)
# }
# ll1 <- function(node, exceptions){
#   node <- lockleaves(node, exceptions)
#   if(is.list(node)) node[] <- lapply(node, ll1, exceptions)
#   return(node)
# }
# for(i in seq_along(trees)){
#   trees[[i]] <- ll1(tree, exceptions = names(nmembers)[i])
# }

# for(i in seq_along(trees)){
#   whichclade <- names(nmembers)[i]
#   if(whichclade == ""){
#     index <- ""
#   }else{
#     whichcladei <- as.numeric(strsplit(whichclade, split = "")[[1]])
#     index <- paste0(paste0("[[", paste0(whichcladei, collapse = "]][["), "]]"))
#   }
#   toeval <- paste0("tree", index, "<- trees[[i]]", index)
#   eval(parse(text = toeval))
#   trees[[i]] <- NA
#   gc()
# }
# rm(nmembers)
# rm(eligible)


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
# #tree <- dendrapply(tree, label, x = if(has_duplicates) fullseqset else x)
# tree <- dendrapply(tree, label, x)
# attachlins <- function(node, x){
#   splitfun <- function(s) strsplit(s, split = ";")[[1]]
#   lineages <- lapply(attr(x, "lineage")[attr(node, "sequences")], splitfun)
#   linlengths <- sapply(lineages, length)
#   whichminlen <- which.min(linlengths)
#   minlen <- linlengths[whichminlen]
#   minlin <- lineages[[whichminlen]]
#   lineage <- ""
#   for(l in 1:minlen){
#     inminlin <- sapply(lineages, function(e) minlin[l] %in% e)
#     if(all(inminlin)){
#       lineage <- paste0(lineage, lineages[[whichminlen]][l], ";")
#     }
#   }
#   # lineage <- gsub(";$", "\\.", lineage)
#   lineage <- gsub(";$", "", lineage)
#   attr(node, "lineage") <- lineage
#   return(node)
# }
#


# if(newK == 1){
#   res$membership <- rep(1L, nseq)
#   res$success <- FALSE
#   res$scores <- scores[tmp[1], ]
#   res$phmm1 <- model
#   return(res)
# }else if(newK != K){
#   ####
#   if(!quiet) cat("Reducing number of splits from", K, "to", newK, "\n")
#   res <- split.DNA2(x, model = model, needs_training = FALSE, allocation = allocation,
#              K = newK, refine = refine, iterations = iterations,
#              quiet = quiet, ... = ...)
#   return(res)
# }


# #
#
# hphmm <- function(tree, refine = "Viterbi", parents = FALSE, quiet = FALSE){
#   #tree leaves have sequences attributes
#   passphmms <- function(tree, refine = "Viterbi", parents = FALSE, quiet = FALSE){
#     if(is.list(tree)){
#       has.seqs <- sapply(tree, function(e) !is.null(attr(e, "sequences")))
#       needs.phmm <- sapply(tree, function(e) is.null(attr(e, "phmm")))
#       if(all(has.seqs & needs.phmm)){
#         seq.numbers <- sapply(tree, function(e) length(attr(e, "sequences")))
#         mcn <- min(seq.numbers)
#         if(mcn < 4) mcn <- 4 # otherwise signal gets drowned by pseudocounts
#         attr(tree, "sequences") <- structure(list(), class = "DNAbin")
#         attr(tree, "label") <- ""
#         for(i in seq_along(tree)){
#           # also need to retain atributes species, definition, etc TODO
#           attr(tree, "sequences") <- c(attr(tree, "sequences"), attr(tree[[i]], "sequences"))
#           attr(tree, "label") <- paste(attr(tree, "label"), attr(tree[[i]], "label"), sep = ". ")
#         }
#         attr(tree, "label") <- gsub("^. ", "", attr(tree, "label"))
#         combined.phmm <- derive.PHMM(attr(tree, "sequences"), refine = "Viterbi", quiet = quiet)
#         for(i in seq_along(tree)){
#           if(seq.numbers[i] > 2){
#             qds <- phylogram::kdistance(attr(tree[[i]], "sequences"), k = 5, alpha = NULL)
#             guidetree <- as.dendrogram(hclust(qds, method = "average"))
#             seqweights <- aphid::weight(guidetree, method = "Gerstein")[names(attr(tree[[i]], "sequences"))]
#             seqweights <- seqweights * mcn/seq.numbers[i] #sum(seqweights)
#           }else if(seq.numbers[i] == 2){
#             seqweights <- c(2, 2)
#           }else seqweights <- 4#1
#           if(!quiet) cat("Training model for clade ", attr(tree[[i]], "label"), "\n")
#           attr(tree[[i]], "phmm") <- train(combined.phmm, attr(tree[[i]], "sequences"), method = refine,
#                                            seqweights = seqweights, inserts = "none", deltaLL = 1E-03, quiet = quiet)
#         }
#         if(parents){ #slower but also allows parent model comparison in classify.DNA
#           seq.number <- length(attr(tree, "sequences"))
#           if(seq.number > 2){
#             qds <- phylogram::kdistance(attr(tree, "sequences"), k = 5, alpha = NULL)
#             guidetree <- as.dendrogram(hclust(qds, method = "average"))
#             seqweights <- aphid::weight(guidetree, method = "Gerstein")[names(attr(tree, "sequences"))]
#             seqweights <- seqweights * mcn/seq.number
#           }else if(seq.number == 2){
#             seqweights <- c(2, 2)
#           }
#           attr(tree, "parent") <- train(combined.phmm, attr(tree, "sequences"), method = refine,
#                                         seqweights = seqweights, inserts = "none", deltaLL = 1E-03, quiet = quiet)
#         }
#       }
#     }
#     return(tree)
#   }
#   hphmm1 <- function(tree, refine = "Viterbi", parents = FALSE, quiet = FALSE){
#     #tree leaves have sequences attributes
#     if(is.list(tree)){
#       tree[] <- lapply(tree, passphmms, refine = refine, parents = parents, quiet = quiet)
#       tree[] <- lapply(tree, hphmm1, refine = refine, parents = parents, quiet = quiet)
#     }
#     return(tree)
#   }
#   repeat{
#     stopifnot(is.list(tree))
#     has.seqs <- sapply(tree, function(e) !is.null(attr(e, "sequences")))
#     needs.phmm <- sapply(tree, function(e) is.null(attr(e, "phmm")))
#     if(all(has.seqs & needs.phmm)){
#       tree <- passphmms(tree, refine = refine, parents = parents, quiet = quiet)
#       return(tree)
#     }else{
#       tree <- hphmm1(tree, refine = refine, parents = parents, quiet = quiet)
#     }
#   }
# }


#
# fork <- function(tree, sequences, refine = "Viterbi", iterations = 10, maxK = 5,
#                  minscore = 0.9, probs = 0, resize = TRUE, splinter = FALSE,
#                  quiet = FALSE, ...){
#   # tree is a dendrogram obj, sequences is a DNAbin object
#   # resize, should models be allowed to change size at each lower level?
#   if(!is.list(tree)){ # fork leaves only
#     seqs <- sequences[attr(tree, "sequences")]
#     nseq <- length(seqs)
#     if(nseq == 1) return(tree)
#     if(nseq < maxK) maxK <- nseq
#     if(!quiet) cat("\nAttempting to split clade", attr(tree, "clade"), "\n")
#     if(is.null(attr(tree, "phmm"))){ # should only happen at top level
#       mod <- NULL
#     }else{
#       mod <- attr(tree, "phmm")
#       if(resize){
#         if(!quiet) cat("Retraining parent model\n")
#         # model is allowed to change size here
#         mod <- train(mod, seqs, method = "Viterbi", ... = ...)
#         if(!quiet) cat("New model size :", mod$size, "\n")
#         if(refine == "BaumWelch") mod <- train(mod, seqs, method = "BaumWelch", ... = ...)
#       }
#     }
#     split_node <- FALSE
#     nclades <- 2
#     seqsplit <- partition(seqs, model = mod, needs_training = FALSE, refine = refine,
#                           K = nclades, iterations = iterations, quiet = quiet, ... = ...)
#     if(is.null(seqsplit)) return(tree)
#     if(is.null(attr(tree, "phmm"))){ # should only happen at top level
#       if(!quiet) cat("Assigning top-level model\n")
#       attr(tree, "phmm") <- seqsplit$phmm0
#     }
#     if(is.null(attr(tree, "scores"))){ # should only happen at top level
#       if(!quiet) cat("Calculating top-level scores\n")
#       scores <- numeric(nseq)
#       for(i in 1:nseq){
#         scores[i] <- forward(seqsplit$phmm0, seqs[[i]], odds = FALSE, ... = ...)$score
#       }
#       attr(tree, "scores") <- scores
#     }
#     membership <- seqsplit$membership
#     scores <- seqsplit$scores
#     total_scores <- apply(scores, 2, aphid::logsum)
#     akwgts <- t(exp(t(scores) - total_scores))
#     performances <- numeric(nseq)
#     for(i in 1:nseq) performances[i] <- akwgts[membership[i], i]
#     if(!quiet) cat("Akaike weights:", performances, "\n")
#     minperfs <- numeric(nclades)
#     # for(i in 1:nclades) minperfs[i] <- min(performances[membership == i])
#     for(i in 1:nclades){
#       minperfs[i] <- quantile(performances[membership == i], probs = probs)
#       if(!quiet){
#         cat("Group", i, "\n")
#         cat("Size: ", sum(membership == i), "\n")
#         cat(sum(performances[membership == i] > minscore), "of",
#             sum(membership == i), "correctly predicted with Akaike weight > ", minscore, "\n")
#         cat(sum(performances[membership == i] > 0.5), "of",
#             sum(membership == i), "correctly predicted with Akaike weight > 0.5\n")
#         cat("Lower", probs, "quantile of Akaike weights:", minperfs[i], "\n")
#         cat("Minimum Akaike weight:", min(performances[membership == i]), "\n")
#       }
#     }
#     if(any(minperfs < minscore)){
#       tmpnclades <- nclades + 1
#       tmperfs <- performances
#       tmpmbrs <- membership
#       if(!quiet & tmpnclades <= maxK){
#         cat("Minimum performance threshold not reached, trying multi-way splits\n")
#       }
#       while(tmpnclades <= maxK){
#         if(splinter){
#           whichmintmperfs <- which.min(tmperfs)
#           tmpmbrs[whichmintmperfs] <- tmpnclades # forms a new clade
#           tmperfs[whichmintmperfs] <- 1
#         }
#         tmpsplit <- partition(seqs, model = mod, needs_training = FALSE, refine = refine,
#                               allocation = if(splinter) tmpmbrs else "cluster",
#                               K = tmpnclades, iterations = iterations, quiet = quiet, ... = ...)
#         if(!is.null(tmpsplit)){
#           tmpmbrs <- tmpsplit$membership
#           tpmscrs <- tmpsplit$scores
#           tmptots <- apply(tpmscrs, 2, aphid::logsum)
#           tmpwgts <- t(exp(t(tpmscrs) - tmptots))
#           tmperfs <- numeric(nseq)
#           for(i in 1:nseq) tmperfs[i] <- tmpwgts[tmpmbrs[i], i]
#           if(!quiet) cat("Candidate Akaike weights:", tmperfs, "\n")
#           tmpmpfs <- numeric(tmpnclades)
#           #for(i in 1:tmpnclades) tmpmpfs[i] <- min(tmperfs[tmpmbrs == i])
#           for(i in 1:tmpnclades) tmpmpfs[i] <- quantile(tmperfs[tmpmbrs == i], probs = probs)
#           if(!quiet) cat("Candidate lower quantile Akaike weights by clade:", tmpmpfs, "\n")
#           if(min(tmpmpfs) > min(minperfs)){
#             if(!quiet){
#               cat(tmpnclades, "way split was superior to predecessor\n")
#               cat("testing performance against threshold\n")
#             }
#             nclades <- tmpnclades
#             seqsplit <- tmpsplit
#             membership <- tmpmbrs
#             scores <- tpmscrs
#             total_scores <- tmptots
#             akwgts <- tmpwgts
#             performances <- tmperfs
#             minperfs <- tmpmpfs
#             if(all(minperfs > minscore)){
#               if(!quiet) cat("Performance threshold reached,", nclades, "way split successful\n")
#               split_node <- TRUE
#               break
#             }else{
#               if(!quiet){
#                 cat("Minimum performance threshold not reached\n")
#                 if(tmpnclades == maxK){
#                   cat("Unable to split clade\n")
#                 }else cat("Splitting further\n")
#               }
#             }
#           }else{
#             if(!quiet){
#               if(tmpnclades == maxK){
#                 cat("Higher-level split was unsuccessful, unable to split clade\n")
#               }else{
#                 cat("Higher-level split was suboptimal, splitting further\n")
#               }
#             }
#             # if(!quiet) cat("Higher-level splitting was unsuccessful, unable to split clade\n")
#             # break
#           }
#         }
#         tmpnclades <- tmpnclades + 1
#       }
#     }else split_node <- TRUE
#     # Akaike weights should all be close to 1
#     # now decorate the tree (if disc ability is > a certain threshold?)
#     if(split_node){# placeholder for discriminant thresholding
#       # change from leaf to inner node
#       if(!quiet) cat("Creating new node\n")
#       tmpattr <- attributes(tree)
#       tree <- vector(mode = "list", length = nclades)
#       attributes(tree) <- tmpattr
#       attr(tree, "leaf") <- NULL
#       for(i in 1:nclades){
#         tree[[i]] <- 1
#         attr(tree[[i]], "height") <- attr(tree, "height") - 1
#         attr(tree[[i]], "leaf") <- TRUE
#         #attr(tree[[i]], "label") <- paste0(attr(tree, "label"), i)
#         attr(tree[[i]], "clade") <- paste0(attr(tree, "clade"), i)
#         attr(tree[[i]], "sequences") <- attr(tree, "sequences")[membership == i]
#         attr(tree[[i]], "scores") <- scores[i, membership == i]
#         attr(tree[[i]], "phmm") <- seqsplit[[paste0("phmm", i)]]
#       }
#     }
#   }
#   return(tree)
# }
