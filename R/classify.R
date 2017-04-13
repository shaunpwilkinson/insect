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
#' @return a character string giving the lineage of the input sequence
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @seealso \code{\link{learn}}
#' @examples
#'   ##TBA
################################################################################
classify <- function(x, tree, threshold = 0.9){
  classify1 <- function(x, tree, threshold = 0.9){
    # res <- data.frame(Path = integer(100), Akaike_weight = numeric(100),
    #                   Score = numeric(100), stringsAsFactors = FALSE)
    path <- integer(100)
    Akweights <- numeric(100)
    scores = numeric(100)
    counter <- 1
    seqhash <- paste(openssl::md5(as.vector(x)))
    matches <- which(attr(tree, "hashes") == seqhash)
    if(length(matches) == 0) matches <- NA
    res <- "" # lineage above the root node
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
      scores[counter] <- sc[best_model]
      threshold_met <- akwgts[best_model] >= threshold
      minscore_met <- sc[best_model] >= min(attr(tree[[best_model]], "scores"))
      if(!threshold_met){# | !minscore_met){# | best_model == no_mods + 1){
        path <- path[1:counter]
        Akweights <- Akweights[1:counter]
        scores <- scores[1:counter]
        # res <- res[1:counter, ]
        break
      }
      res <- attr(tree, "lineage")
      tree <- tree[[best_model]]
      counter <- counter + 1
    }
    # if(res[counter, 1] == 0) res <- res[1:(counter - 1), ]
    if(path[counter] == 0){
      path <- path[1:(counter - 1)]
      Akweights <- Akweights[1:(counter - 1)]
      scores <- scores[1:(counter - 1)]
    }
    if(counter > 1){
      if(matches[1] %in% attr(tree, "sequences")){
        res <- attr(tree, "lineage") # exact match -> can confidently drop 1 level
      }
    }
    attr(res, "path") <- path
    attr(res, "scores") <- scores
    attr(res, "weights") <- Akweights
    attr(res, "threshold") <- threshold_met
    attr(res, "minscore") <- minscore_met
    attr(res, "matches") <- matches
    return(res)
  }
  if(is.list(x)){
    lapply(x, classify1, tree, threshold)
  }else{
    classify1(x, tree, threshold)
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
