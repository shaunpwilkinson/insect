#' Tree-based sequence classification.
#'
#' \code{"classify"} probabilistically assigns a DNA barcode or other
#'   taxonomically-informative sequence to a node of a classification
#'   tree by recursively testing the sequence against a series of nested profile
#'   hidden Markov models. This enables an Akaike weight confidence
#'   score to be output along with the taxon ID.
#'
#' @param x a DNA sequence or set of sequences. Acceptable input formats are
#'   a single sequence (either as a character string or a "DNAbin" vector),
#'   a vector of character strings, a DNAbin list,
#'   a list of character string vectors (one vector for each sample),
#'   or a list of DNAbin list objects (one DNAbin object for each sample).
#'   The latter two options enable easier tabulation of sequence
#'   counts by site/sample by passing the output on to the
#'   \code{\link{tabulize}} function.
#' @param tree an object of class \code{"insect"}
#'   (see \code{\link{learn}} for details about this object class).
#' @param threshold numeric between 0 and 1 giving the minimum
#'   Akaike weight for the recursive classification procedure
#'   to continue toward the leaves of the tree. Defaults to 0.9.
#' @param decay logical indicating whether the decision to terminate the
#'   classification process should occur based on decaying Akaike weights
#'   (at each node, the Akaike weight of the selected model is multiplied by
#'   the Akaike weight of the selected model at the parent node) or whether
#'   each Akaike weight should be calculated independently of that of the
#'   parent node. Defaults to the latter (FALSE).
#' @param ping logical indicating whether an exact-match search should
#'   be carried out before applying the classification algorithm.
#'   If TRUE (the default value) and the query sequence is identical to
#'   at least one of the training sequences used to learn the tree,
#'   the common ancestor of the matching training sequences is returned
#'   with an associated score value of 1.
#'   The output lineage string will generally specify the taxonomic ID to species level
#'   but may be to genus/family, etc in cases where the genetic marker lacks resolution.
#' @param cores integer giving the number of CPUs to parallelize the operation
#'   over. Defaults to 1, and reverts to 1 if x is not a list.
#'   This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   for example by running \code{parallel::stopCluster(cores)}.
#'   The string 'autodetect' is also accepted, in which case the maximum
#'   number of cores to use is one less than the total number of cores available.
#' @param scores logical indicating whether the Akaike weight confidence
#'   values should be attributed to the output object as a numeric vector
#'   called "score". Defaults to TRUE.
#' @param paths logical indicating whether the path of each sequence
#'   through the classification tree should be attributed to the output object
#'   as a character vector called "path". Defaults to FALSE.
#' @param seqlengths logical indicating whether the original sequence
#'   lengths should be attributed to the returned object.
#' @param threshs,minscores,minlengths,maxlengths logicals, indicating
#'   which classification test results should be attributed to the output object
#'   (for advanced use).
#' @return a semicolon-delimited character string giving the lineage of the
#'   input sequence, usually with additional attributes.
#' @details
#'   This function requires a pre-computed classification tree
#'   of class "insect", which is simply a dendrogram object with additional attributes
#'   (see \code{\link{learn}} for details).
#'   Query sequences obtained from the same primer set used to construct
#'   the tree are classified to produce taxonomic
#'   IDs with an associated degree of confidence.
#'   The classification algorithm works as follows:
#'   starting from the root node of the tree,
#'   the log-likelihood of the query sequence
#'   (the log-probability of the sequence given a particular model)
#'   is computed for each of the models occupying the two child nodes using the
#'   forward algorithm (see Durbin et al. (1998)).
#'   The competing likelihood values are then compared by computing
#'   their Akaike weights (Johnson and Omland, 2004).
#'   If one model is overwhelmingly more likely to have produced
#'   the sequence than the other,
#'   that child node is chosen and the classification is updated
#'   to reflect the taxonomic ID stored at the node.
#'   This classification procedure is repeated, continuing down the
#'   tree until either an inconclusive result is returned by a
#'   model comparison test (i.e. the Akaike weight is lower than
#'   a pre-defined threshold, usually 0.9),
#'   or a terminal leaf node is reached,
#'   at which point a species-level classification is generally returned.
#'   The function outputs a taxonomic ID as a semicolon-delimited lineage string,
#'   as well as the Akaike weight of the model at the final node and several
#'   other results depending on the argument values.
#'   Note that the default behavior is for the Akaike weight to ‘decay’
#'   as it moves down the tree, by computing the cumulative product of
#'   all preceding Akaike weight values.
#'   This is perhaps an overly conservative approach
#'   but it minimizes the chance of generating spurious sequence classifications.
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Johnson JB, Omland KS (2004) Model selection in ecology and evolution.
#'   \emph{Trends in Ecology and Evolution}. \strong{19}, 101-108.
#'
#' @seealso \code{\link{learn}} \code{\link{tabulize}}
#' @examples
#' \donttest{
#'   data(whales)
#'   ## use sequences 2-19 to learn the tree
#'   ## note that training data must retain lineage attribute
#'   training_data <- subset.DNAbin(whales, subset = seq_along(whales) > 1)
#'   ## append species names to semicolon-delimited lineage strings for added resolution
#'   attr(whales, "lineage") <- paste0(attr(whales, "lineage"), "; ", attr(whales, "species"))
#'   ## learn the tree
#'   set.seed(999)
#'   tree <- learn(training_data, cores = 2, quiet = FALSE, maxiter = 5)
#'   ## find predicted lineage for sequence #1
#'   classify(whales[[1]], tree)
#'   ## compare with actual lineage
#'   attr(whales, "lineage")[1]
#' }
################################################################################
classify <- function(x, tree, threshold = 0.9, decay = TRUE, ping = TRUE,
                     cores = 1, scores = TRUE, paths = FALSE, seqlengths = FALSE,
                     threshs = FALSE, minscores = FALSE, minlengths = FALSE,
                     maxlengths = FALSE){
  if(is.null(attr(tree, "key"))){
    cat("Note: ping is TRUE but tree has no hash key. Exact matching not possible\n")
  }
  depth <- function(l) ifelse(is.list(l), 1L + max(sapply(l, depth)), 0L)
  ## thanks flodal for depth function
  xdepth <- depth(x)
  stopifnot(xdepth < 3)
  isbin <- if(xdepth == 2) .isDNA(x[[1]]) else .isDNA(x)
  if(xdepth == 2 & !isbin) stop("List has depth 2 but elements are not of class 'DNAbin'")
  lol <- if(isbin) xdepth == 2 else xdepth == 1
  if(lol){
    splinds <- rep(seq_along(x), sapply(x, length)) # split indices
    #for(i in seq_along(x)) names(x[[i]]) <- paste0(names(x[[i]]), ":", i)
    orignames <- names(x)
    newnames <- unlist(lapply(x, names), use.names = FALSE)
    x <- unlist(x, use.names = FALSE, recursive = FALSE)
    names(x) <- newnames
    if(isbin) class(x) <- "DNAbin"
  }
  if(!isbin) x <- char2dna(x)
  hashes <- hash(x, cores = cores)
  pointers <- .point(hashes)
  dupes <- duplicated(pointers)
  needs_rerep <- any(dupes)
  if(needs_rerep){
    rerepnames <- names(x)
    x <- x[!dupes]
  }
  classify1 <- function(x, tree, threshold = 0.9, decay = TRUE, ping = TRUE){
    if(ping){
      xhash <- hash(x)
      xmatch <- attr(tree, "key")[xhash][1]
      if(!is.na(xmatch)) return(paste(c(xmatch, 0, 1, rep(TRUE, 4)), collapse = "%"))
    }
    path <- ""
    akw <- 1
    cakw <- 1
    res <- "" # lineage above the root node
    ## test at root node first - not needed?
    # rootsc <- aphid::forward.PHMM(decodePHMM(attr(tree, "model")), x, odds = FALSE)$score
    # minscore_met <- if(is.null(attr(tree, "minscore"))) TRUE else rootsc > (attr(tree, "minscore") - 4)
    # ## trees made before 27 Mar 18 don't have top level minscore attr
    # minlength_met <- length(x) >= attr(tree, "minlength") - 3
    # maxlength_met <- length(x) <= attr(tree, "maxlength") + 3
    # threshold_met <- TRUE
    # score <- if(decay) cakw else akw
    # if(!(threshold_met & minscore_met & minlength_met & maxlength_met)){
    #   res <- paste(c(res, path, paste(score), paste(threshold_met),
    #                  paste(minscore_met), paste(minlength_met),
    #                  paste(maxlength_met)), collapse = "%")
    #   return(res)
    # }
    while(is.list(tree)){
      no_mods <- length(tree)
      sc <- numeric(no_mods) # scores (log probabilities)
      for(i in 1:no_mods){
        modi <- decodePHMM(attr(tree[[i]], "model"))
        sc[i] <- aphid::forward.PHMM(modi, x, odds = FALSE)$score
        # takes a similar amount of time
      }
      total_score <- aphid::logsum(sc)
      akwgts <- exp(sc - total_score)
      best_model <- which.max(akwgts)
      newakw <- akwgts[best_model]
      newcakw <- newakw * cakw
      threshold_met <- threshold <= if(decay) newcakw else newakw
      minscore_met <- sc[best_model] >= attr(tree[[best_model]], "minscore") - 4
      # 4 is approx asymtote for single bp change as n training seqs -> inf
      minlength_met <- length(x) >= attr(tree[[best_model]], "minlength") - 3
      maxlength_met <- length(x) <= attr(tree[[best_model]], "maxlength") + 3
      if(!(threshold_met & minscore_met & minlength_met & maxlength_met)) break
      path <- paste0(path, best_model)
      akw <- newakw
      cakw <- newcakw
      tree <- tree[[best_model]]
      res <- attr(tree, "lineage")
    }
    score <- if(decay) cakw else akw
    res <- paste(c(res, path, paste(score), paste(threshold_met),
                   paste(minscore_met), paste(minlength_met),
                   paste(maxlength_met)), collapse = "%")
    return(res)
  }
  unpack <- function(v){# v is a vector of %-delimited strings
    v <- strsplit(v, split = "%")
    out <- sapply(v, function(s) s[1])
    attr(out, "path") <- unname(sapply(v, function(s) s[2]))
    attr(out, "score") <- as.numeric(sapply(v, function(s) s[3]))
    attr(out, "threshold_met") <- as.logical(sapply(v, function(s) s[4]))
    attr(out, "minscore_met") <- as.logical(sapply(v, function(s) s[5]))
    attr(out, "minlength_met") <- as.logical(sapply(v, function(s) s[6]))
    attr(out, "maxlength_met") <- as.logical(sapply(v, function(s) s[7]))
    return(out)
  }
  if(is.list(x)){
    if(inherits(cores, "cluster")){
      res <- parallel::parSapply(cores, x, classify1, tree, threshold, decay, ping)
    }else if(cores == 1){
      res <- sapply(x, classify1, tree, threshold, decay, ping)
    }else{
      navailcores <- parallel::detectCores()
      if(identical(cores, "autodetect")) cores <- navailcores - 1
      if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
      # if(cores > navailcores) stop("Insufficient CPUs available")
      if(cores > 1){
        cl <- parallel::makeCluster(cores)
        res <- parallel::parSapply(cl, x, classify1, tree, threshold, decay, ping)
        parallel::stopCluster(cl)
      }else{
        res <- sapply(x, classify1, tree, threshold, decay, ping)
      }
    }
    res <- unpack(res)
  }else{
    res <- unpack(classify1(x, tree, threshold, decay, ping))
  }
  if(!scores) attr(res, "score") <- NULL
  if(!paths) attr(res, "path") <- NULL
  if(!threshs) attr(res, "threshold_met") <- NULL
  if(!minscores) attr(res, "minscore_met") <- NULL
  if(!minlengths) attr(res, "minlength_met") <- NULL
  if(!maxlengths) attr(res, "maxlength_met") <- NULL
  ## override classifications for exact matches (optional) ###
  if(needs_rerep){
    tmpattr <- attributes(res)
    res <- res[pointers]
    for(i in seq_along(tmpattr)) tmpattr[[i]] <- tmpattr[[i]][pointers]
    if(!all(sapply(tmpattr, length) == length(res))) stop("Error 1\n")
    ## careful here in case of future addition of class attribute
    tmpattr$names <- rerepnames
    attributes(res) <- tmpattr
    # attr(res, "score") <- tmpattr$score[pointers]
    # attr(res, "path") <- tmpattr$path[pointers]
    # attr(res, "threshold_met") <- tmpattr$threshold_met[pointers]
    # attr(res, "minscore_met") <- tmpattr$minscore_met[pointers]
    # attr(res, "minlength_met") <- tmpattr$minlength_met[pointers]
    # attr(res, "maxlength_met") <- tmpattr$maxlength_met[pointers]
  }
  attr(res, "hash") <- unname(hashes) # needed for tabulize
  if(seqlengths) attr(res, "length") <- unname(sapply(x, length)[pointers])
  if(lol){
    tmpattr <- attributes(res)
    res <- split(res, f = splinds)
    for(i in seq_along(tmpattr)[names(tmpattr) != "names"]){
      splattr <- split(tmpattr[[i]], f = splinds)
      for(j in seq_along(res)){
        attr(res[[j]], names(tmpattr)[i]) <- splattr[[j]]
      }
    }
    names(res) <- orignames
  }
  return(res)
}
################################################################################
