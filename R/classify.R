#' Tree-based sequence classification.
#'
#' \code{"classify"} assigns taxon IDs to DNA sequences using an
#'   informatic sequence classification tree.
#'
#' @param x a DNA sequence or set of sequences. Can be a
#'   "DNAbin" object or a named vector of upper-case character strings.
#' @param tree an object of class \code{"insect"}
#'   (see \code{\link{learn}} for details).
#' @param threshold numeric between 0 and 1 giving the minimum
#'   Akaike weight for the recursive classification procedure
#'   to continue toward the leaves of the tree. Defaults to 0.9.
#' @param decay logical indicating whether the decision to terminate the
#'   classification process should be made based on decaying Akaike weights
#'   (at each node, the Akaike weight of the selected model is multiplied by
#'   the Akaike weight of the selected model at the parent node) or whether
#'   each Akaike weight should be calculated independently of that of the
#'   parent node. Defaults to TRUE (the former).
#' @param ping logical indicating whether an exact-match search should
#'   be carried out before applying the classification algorithm.
#'   If TRUE (the default) and the query sequence is identical to
#'   at least one of the training sequences used to learn the tree,
#'   the common ancestor of the matching training sequences is returned
#'   with an score of NA.
#'   The output lineage string will generally specify the taxonomic ID to species level
#'   but may be to genus/family, etc for low resolution genetic markers.
#' @param ranks character vector giving the taxonomic ranks to be
#'   included in the output table. Must be a valid rank from the
#'   taxonomy database attributed to the classification tree
#'   (\code{attr(tree, "taxonomy")}). Set to NULL to exclude taxonomic ranks
#'   from the output table.
#' @param cores integer giving the number of CPUs to parallelize the operation
#'   over (defaults to 1).
#'   This argument may alternatively be a 'cluster' object,
#'   in which case it is the user's responsibility to close the socket
#'   connection at the conclusion of the operation,
#'   for example by running \code{parallel::stopCluster(cores)}.
#'   The string 'autodetect' is also accepted, in which case the maximum
#'   number of cores to use is one less than the total number of cores available.
#' @return a data.frame.
#' @details
#'   This function requires a pre-computed classification tree
#'   of class "insect", which is a dendrogram object with additional attributes
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
#'   a pre-defined threshold, e.g. 0.9),
#'   or a terminal leaf node is reached,
#'   at which point a species-level classification is generally returned.
#'   The function outputs a table with one row for each input sequence
#'   Output table fields include "name" (the unique sequence identifier),
#'   "taxID" (the taxonomic identification number from the taxonomy database),
#'   "taxon" (the name of the taxon),
#'   "rank" (the rank of the taxon, e.g. species, genus family, etc),
#'   and "score" (the Akaike weight from the model selection procedure).
#'   Note that the default behavior is for the Akaike weight to ‘decay’
#'   as it moves down the tree, by computing the cumulative product of
#'   all preceding Akaike weight values.
#'   This is perhaps an overly conservative approach
#'   but it minimizes the chance of taxon ID errors.
#'   The output table also includes the higher taxonomic ranks specified in the
#'   \code{ranks} argument, and two additional columns called "path"
#'   (the path of the sequence through the classification tree)
#'   and "reason" outlining why the recursive classification procedure was
#'   terminated. Reason codes are as follows:
#'   \itemize{
#'     \item 0 reached leaf node
#'     \item 1 failed to meet minimum score threshold at inner node
#'     \item 2 failed to meet minimum score of training sequences at inner node
#'     \item 3 sequence length shorter than minimum length of training sequences at inner node
#'     \item 4 sequence length exceeded maximum length of training sequences at inner node
#'   }.
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#'
#'   Johnson JB, Omland KS (2004) Model selection in ecology and evolution.
#'   \emph{Trends in Ecology and Evolution}. \strong{19}, 101-108.
#'
#' @seealso \code{\link{learn}}
#' @examples
#' \donttest{
#'   data(whales)
#'   data(whale_taxonomy)
#'   ## use all sequences except first one to train the classifier
#'   set.seed(999)
#'   tree <- learn(whales[-1], db = whale_taxonomy, maxiter = 5, cores = 2)
#'   ## find predicted lineage for first sequence
#'   classify(whales[1], tree)
#'   ## compare with actual lineage
#'   taxID <- as.integer(gsub(".+\\|", "", names(whales)[1]))
#'   get_lineage(taxID, whale_taxonomy)
#' }
################################################################################
classify <- function(x, tree, threshold = 0.9, decay = TRUE, ping = TRUE,
                     ranks = c("phylum", "class", "order", "family", "genus", "species"),
                     cores = 1){
  if(is.null(names(x))) names(x) <- paste0("S", seq_along(x))
  if(mode(x) == "character") x <- char2dna(x, simplify = FALSE)
  if(!is.list(x)){
    if(!mode(x) == "raw") stop("Unrecognized format for x\n")
    nam <- "S1"# deparse(substitute(x))
    x <- list(x)
    names(x) <- nam
    class(x) <- "DNAbin"
  }
  origins <- if(all(grepl("_", names(x)))){
    gsub("_.+", "", names(x))
  }else{
    rep("sample1", length(x))
  }
  usm <- unique(origins)
  origins <- factor(origins, levels = usm)
  ## very important to specify ordering!
  hashes <- hash(x)
  duplicates <- duplicated(hashes)
  x <- x[!duplicates]
  uhashes <- hashes[!duplicates]
  indices <- split(match(hashes, uhashes), f = origins)
  qout <- matrix(NA_integer_, nrow = sum(!duplicates), ncol = length(usm))
  colnames(qout) <- usm
  for(i in seq_along(usm)) qout[, i] <- tabulate(indices[[i]], nbins = length(x))
  if(!inherits(tree, "dendrogram")) stop("Unrecognized tree format\n")
  db <- attr(tree, "taxonomy")
  if(is.null(db)) stop("tree is missing taxonomy database\n")
  attr(tree, "taxonomy") <- NULL # reduce memory usage for multithreading
  if(ping){
    if(is.null(attr(tree, "key"))){
      warning(paste0("ping is TRUE but tree has no hash key. ",
                     "Exact matching not possible\n"))
    }
  }
  # decode second and third tier models for increased speed
  for(i in seq_along(tree)){
    attr(tree[[i]], "model") <- decodePHMM(attr(tree[[i]], "model"))
    if(is.list(tree[[i]])){
      for(j in seq_along(tree[[i]])){
        attr(tree[[i]][[j]], "model") <- decodePHMM(attr(tree[[i]][[j]], "model"))
      }
    }
  }
  classify1 <- function(x, tree, threshold = 0.9, decay = TRUE, ping = TRUE){
    ## takes a single named raw vector
    ## outputs a 1-row dataframe
    xhash <- hash(x)
    if(ping){
      xmatch <- attr(tree, "key")[xhash][1]
      if(!is.na(xmatch)){
        out <- data.frame(taxID = xmatch,
                          score = NA_real_,
                          path = NA_character_,
                          reason = NA_integer_)
        return(out)
      }
    }
    path <- ""
    akw <- 1
    cakw <- 1
    tax <- 0L # lineage above the root node
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
      minlength_met <- length(x) >= attr(tree[[best_model]], "minlength") - 2
      maxlength_met <- length(x) <= attr(tree[[best_model]], "maxlength") + 2
      if(!(threshold_met & minscore_met & minlength_met & maxlength_met)) break
      path <- paste0(path, best_model)
      akw <- newakw
      cakw <- newcakw
      tree <- tree[[best_model]]
      tax <- attr(tree, "taxID")
    }
    score <- round(if(decay) cakw else akw, 4)
    reason <- if(!threshold_met){
      1L
    }else if(!minscore_met){
      2L
    }else if(!minlength_met){
      3L
    }else if(!maxlength_met){
      4L
    }else{
      0L
    }
    out <- data.frame(taxID = tax,
                      score = score,
                      path = path,
                      reason = reason)
    return(out)
  }
  if(inherits(cores, "cluster")){
    res <- parallel::parLapply(cores, x, classify1, tree, threshold, decay, ping)
  }else if(cores == 1){
    res <- lapply(x, classify1, tree, threshold, decay, ping)
  }else{
    navailcores <- parallel::detectCores()
    if(identical(cores, "autodetect")) cores <- navailcores - 1
    if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
    if(cores > 1){
      cl <- parallel::makeCluster(cores)
      res <- parallel::parLapply(cl, x, classify1, tree, threshold, decay, ping)
      parallel::stopCluster(cl)
    }else{
      res <- lapply(x, classify1, tree, threshold, decay, ping)
    }
  }
  res <- do.call("rbind", res) ## changed output to dataframe 20180617
  lineages <- get_lineage(res$taxID, db, cores = cores,
                          simplify = FALSE, numbers = FALSE)
  tmp <- sapply(lineages, tail, 1)
  lhcols <- data.frame(representative = names(x),
                     taxID = res$taxID,
                     taxon = tmp,
                     rank = names(tmp),
                     stringsAsFactors = FALSE)
  if(!is.null(ranks)){
    rnkmat <- matrix(NA_character_, nrow = length(x), ncol = length(ranks))
    rnkmat <- as.data.frame(rnkmat, stringsAsFactors = FALSE)
    colnames(rnkmat) <- ranks
    fun <- function(l, r) if(is.na(l[r])) "" else l[r]
    for(i in seq_along(ranks)){
      rnkmat[ranks[i]] <- vapply(lineages, fun, "", ranks[i])
    }
    lhcols <- cbind(lhcols, rnkmat)
  }
  lhcols <- cbind(lhcols, qout, res[c("score", "reason", "path")])
  lhcols$path <- as.character(lhcols$path)
  rownames(lhcols) <- NULL
  return(lhcols)
}
################################################################################
