#' Taxonomic classification for DNA barcode sequences.
#'
#' \code{"classify"} probabilistically assigns a DNA barcode or other
#'   taxonomically-informative sequence to a node of a classification
#'   tree using a series of nested profile hidden Markov models.
#'
#' @param x a \code{"DNAbin"} object, either as a vector or a list of vectors.
#' @param tree an object of class \code{"insect"}.
#' @param threshold numeric value between 0 and 1 giving the minimum
#'   Akaike weight for the recursive classification procedure
#'   to continue toward the leaves of the tree.
#' @param cores integer giving the number of CPUs to parallelize the operation
#'   over. Defaults to 1, and reverts to 1 if x is not a list.
#'   This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   for example by running \code{parallel::stopCluster(cores)}.
#'   The string 'autodetect' is also accepted, in which case the maximum
#'   number of cores to use is one less than the total number of cores available.
#' @return a character string giving the lineage of the input sequence
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @seealso \code{\link{learn}}
#' @examples
#'   ##TBA
################################################################################
classify <- function(x, tree, threshold = 0.999, cores = 1){
  classify1 <- function(x, tree, threshold = 0.9){
    path <- integer(100)
    Akweights <- numeric(100)
    scores = numeric(100) # log probs of seq given best model
    cAkweights <- numeric(100) # cumulative product of akaike weights
    counter <- 1
    seqhash <- paste(openssl::md5(as.vector(x)))
    matches <- which(attr(tree, "hashes") == seqhash)
    if(length(matches) == 0) matches <- NA
    res <- "" # lineage above the root node
    reachedleaf <- TRUE
    while(is.list(tree)){
      no_mods <- length(tree)
      sc <- numeric(no_mods)
      for(i in 1:no_mods) {
        sc[i] <- aphid::forward.PHMM(attr(tree[[i]], "phmm"), x, odds = FALSE)$score
      }
      total_score <- aphid::logsum(sc)
      akwgts <- exp(sc - total_score)
      best_model <- which.max(akwgts)
      path[counter] <- best_model
      Akweights[counter] <- akwgts[best_model]
      cAkweights[counter] <- akwgts[best_model] * if(counter == 1) 1 else cAkweights[counter - 1]
      scores[counter] <- sc[best_model]
      #threshold_met <- akwgts[best_model] >= threshold
      threshold_met <- cAkweights[counter] >= threshold
      minscore_met <- sc[best_model] >= min(attr(tree[[best_model]], "scores"))
      if(!(threshold_met)){# & minscore_met)){# | best_model == no_mods + 1){
        #path <- path[1:counter]
        # path vector should be 2 shorter than scores & weights
        path <- if(counter < 3) integer(0) else path[1:(counter - 2)] ###TODO what if counter = 1?
        Akweights <- Akweights[1:counter]
        cAkweights <- cAkweights[1:counter]
        scores <- scores[1:counter]
        reachedleaf <- FALSE
        break
      }
      res <- attr(tree, "lineage")
      tree <- tree[[best_model]]
      counter <- counter + 1
    }
    # if(res[counter, 1] == 0) res <- res[1:(counter - 1), ]
    # if(path[counter] == 0){# i.e. reached leaf node
    if(reachedleaf){# i.e. reached leaf node
      if(matches[1] %in% attr(tree, "sequences")){ # i.e. there is an exact match
        res <- attr(tree, "lineage") # can confidently drop 1 level
        # path <- path[1:(counter - 1)]
        path <- if(counter < 2) integer(0) else path[1:(counter - 1)]
        ## just a backup,  counter would therefore have to be > 1
        ## unlesss tree is a single leaf
      }else{ # reached leaf but no exact match
        # path <- path[1:(counter - 1)]
        path <- if(counter < 3) integer(0) else path[1:(counter - 2)]
        reachedleaf <- FALSE
      }
      if(counter > 1){
        Akweights <- Akweights[1:(counter - 1)]
        cAkweights <- cAkweights[1:(counter - 1)]
        scores <- scores[1:(counter - 1)]
      }else{
        Akweights <- cAkweights <- scores <- numeric(0)
      }
    }
    attr(res, "path") <- path
    attr(res, "scores") <- scores
    attr(res, "weights") <- Akweights
    attr(res, "cweights") <- cAkweights
    attr(res, "threshold") <- threshold_met
    attr(res, "minscore") <- minscore_met
    attr(res, "matches") <- matches
    attr(res, "reachedleaf") <- reachedleaf
    return(res)
  }
  if(is.list(x)){
    if(inherits(cores, "cluster")){
      res <- parallel::parLapply(cores, x, classify1, tree, threshold)
    }else if(cores == 1){
      res <- lapply(x, classify1, tree, threshold)
    }else{
      #nseq <- length(x)
      navailcores <- parallel::detectCores()
      if(identical(cores, "autodetect")) cores <- navailcores - 1
      if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
      if(cores > navailcores) stop("Number of cores to use is more than number available")
      if(cores > 1){
        cl <- parallel::makeCluster(cores)
        res <- parallel::parLapply(cl, x, classify1, tree, threshold)
        parallel::stopCluster(cl)
      }else{
        res <- lapply(x, classify1, tree, threshold)
      }
    }
    return(res)
  }else{
    return(classify1(x, tree, threshold))
  }
}



# res[counter, 1] <- best_model
# res[counter, 2] <- akwgts[best_model]
# res[counter, 3] <- sc[best_model]

# seqs <- seqs[attr(tree, "sequences")]
# # an integer vector of indices pointing to the DNAbin obj (saves on memory)
# # attributes are stripped during subsetting
# attr(seqs, "species") <- seqattr$species[attr(tree, "sequences")]
# attr(seqs, "description") <- seqattr$description[attr(tree, "sequences")]
# attr(seqs, "lineage") <- seqattr$lineage[attr(tree, "sequences")]
# if(length(unique(attr(seqs, "lineage"))) == 1){
#   lineage <- attr(seqs, "lineage")[1]
# }else{
#   lineages <- lapply(attr(seqs, "lineage"), function(s) strsplit(s, split = ";")[[1]])
#   linlengths <- sapply(lineages, length)
#   whichminlen <- which.min(linlengths)
#   minlen <- linlengths[whichminlen]
#   lineage <- ""
#   for(i in 1:minlen){
#     if(all(sapply(lineages, function(e) lineages[[whichminlen]][i] %in% e))){
#       lineage <- paste0(lineage, lineages[[whichminlen]][i], ";")
#     }
#   }
#   lineage <- gsub(";$", "\\.", lineage)
# }
# out <- list()
# out$lineage <- lineage
# out$sequences <- seqs
# out$results <- res
# out$threshold_met <- threshold_met
# out$minscore_met <- minscore_met
# class(out) <- "classification"
