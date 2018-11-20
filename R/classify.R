#' Tree-based sequence classification.
#'
#' \code{"classify"} assigns taxon IDs to DNA sequences using an
#'   informatic sequence classification tree.
#'
#' @param x a sequence or set of sequences. Can be a
#'   "DNAbin" or "AAbin" object or a named vector of
#'   upper-case DNA character strings.
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
#' @param ping logical or numeric (between 0 and 1) indicating whether
#'   a nearest neighbor search should
#'   be carried out, and if so,
#'   what the minimum distance to the nearest neighbor
#'   should be for the the recursive classification algorithm to be skipped.
#'   If TRUE and the query sequence is identical to
#'   at least one of the training sequences used to learn the tree,
#'   the common ancestor of the matching training sequences is returned
#'   with an score of NA.
#'   If a value between 0 and 1 is provided, the common ancestor of the
#'   training sequences with similarity greater than or equal to 'ping'
#'   is returned, again with a score of NA.
#'   If \code{ping} is set to 0 or FALSE, the recursive classification
#'   algorithm is applied to all sequences, regardless of proximity to
#'   those in the training set.
#'   For high values (e.g. \code{ping >= 0.99}) the output will generally
#'   specify the taxonomic ID to species level,
#'   but is often to genus/family/etc level for low resolution genetic markers.
#' @param mincount integer, the minimum number of training sequences belonging to a
#'   selected child node for the classification to progress.
#' @param ranks character vector giving the taxonomic ranks to be
#'   included in the output table. Must be a valid rank from the
#'   taxonomy database attributed to the classification tree
#'   (\code{attr(tree, "taxonomy")}). Set to NULL to exclude taxonomic ranks
#'   from the output table.
#' @param tabulize logical indicating whether sequence counts should be
#'   attached to the output table. If TRUE, the output table will have one
#'   row for each unique sequence, and columns will include counts for
#'   each sample (where samples names precede sequence identifiers in the input
#'   object; see details below).
#' @param metadata logical indicating whether to include additional columns
#'   containing the paths, individual node scores and reasons for termination.
#'   Defaults to FALSE. Included for advanced use and debugging.
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
#'   This minimizes the chance of type I taxon ID errors (overclassifications and misclassifications).
#'   The output table also includes the higher taxonomic ranks specified in the
#'   \code{ranks} argument, and if \code{metadata = TRUE} additional columns
#'   are included called "path"
#'   (the path of the sequence through the classification tree), "scores" (the
#'   scores at each node through the tree, UTF-8-encoded),
#'   and "reason" outlining why the recursive classification procedure was
#'   terminated:
#'   \itemize{
#'     \item 0 reached leaf node
#'     \item 1 failed to meet minimum score threshold at inner node
#'     \item 2 failed to meet minimum score of training sequences at inner node
#'     \item 3 sequence length shorter than minimum length of training sequences at inner node
#'     \item 4 sequence length exceeded maximum length of training sequences at inner node
#'     \item 5 nearest neighbor in training set does not belong to selected node (obsolete)
#'     \item 6 node is supported by too few sequences
#'     \item 7 reserved
#'     \item 8 sequence could not be translated (amino acids only)
#'     \item 9 translated sequence contains stop codon(s) (amino acids only)
#'   }
#'   Additional columns detailing the nearest neighbor search include "NNtaxID", "NNtaxon",
#'   "NNrank", and "NNdistance".
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
                     mincount = 3L, ranks = c("kingdom", "phylum", "class",
                                              "order", "family", "genus", "species"),
                     tabulize = FALSE, metadata = FALSE, cores = 1){
  #if(is.null(attr(tree, "nullscore"))) attr(tree, "nullscore") <- -1E08 ###
  if(!is.null(attr(tree, "training_data"))) attr(tree, "training_data") <- NULL
  urnks <- sort(unique(attr(tree, "taxonomy")$rank))
  if(all(grepl("^rank[0123456789]+", urnks))) ranks <- urnks[-1] ##RDP format
  if(is.null(names(x))) names(x) <- paste0("S", seq_along(x))
  if(mode(x) == "character") x <- char2dna(x, simplify = FALSE)
  if(!is.list(x)){
    if(!mode(x) == "raw") stop("Unrecognized format for x\n")
    nam <- "S1"# deparse(substitute(x))
    x <- list(x)
    names(x) <- nam
    class(x) <- "DNAbin"
  }
  DNA <- inherits(x, "DNAbin")
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
  catchnames <- names(x)
  x <- x[!duplicates]
  uhashes <- hashes[!duplicates]
  if(tabulize){
    indices <- split(match(hashes, uhashes), f = origins)
    qout <- matrix(NA_integer_, nrow = sum(!duplicates), ncol = length(usm))
    colnames(qout) <- usm
    for(i in seq_along(usm)) qout[, i] <- tabulate(indices[[i]], nbins = length(x))
  }else{
    pointers <- .point(hashes)
  }
  if(!inherits(tree, "dendrogram")) stop("Unrecognized tree format\n")
  db <- attr(tree, "taxonomy")
  if(is.null(db)) stop("tree is missing taxonomy database\n")
  attr(tree, "taxonomy") <- NULL # reduce memory usage for multithreading
  key <- attr(tree, "key")
  attr(tree, "key") <- NULL
  z <- attr(tree, "trainingset")
  attr(tree, "trainingset") <- NULL
  m <- decodePHMM(attr(tree, "model"))
  ksize <- attr(tree, "k")
  td <- 1 - ping #threshold distance
  ###doNN <- FALSE
  if(is.null(attr(tree, "kmers"))){ ## for backward compatibility
    warning("Tree is missing kmer count matrix, can't complete nearest-neighbor search\n")
    if(ping > 0 & ping < 1){
      ping <- 0
    }else if(ping == 1){
      if(is.null(key)){
        warning("Tree is missing hash key, setting 'ping' to FALSE")
      }else{
        xhash <- hash(x)
        xmatch <- unname(key[xhash])
        for(i in seq_along(x)){
          if(is.na(xmatch[i])){
            attr(x[[i]], "NN") <- NA_integer_ # sequence index
            attr(x[[i]], "NNhit") <- FALSE # -> do full classification procedure
          }else{
            attr(x[[i]], "NN") <- unname(key[xmatch[i]])# taxid
            attr(x[[i]], "NNhit") <- TRUE # -> skip full classification procedure
          }
        }
      }
    }else if(ping < 1){
      if(ping > 0) stop("Can't neighbor-search without kmer count matrix\n")
      for(i in seq_along(x)){
        attr(x[[i]], "NN") <- NA_integer_ # sequence index
        attr(x[[i]], "NNhit") <- FALSE # -> do full classification procedure
      }
    }else stop("Invalid 'ping' argument\n")
  }else{

    if((ping > 0 & ping < 1) | (ping == 1 & is.null(key))){
      zk <- .decodekc(attr(tree, "kmers"))
      zk <- zk/(attr(tree, "seqlengths") - ksize + 1L)
      attr(tree, "kmers") <- NULL
      xl <- vapply(x, length, 0L, USE.NAMES = FALSE)
      xk <- round(kmer::kcount(x, k = ksize))
      xk <- xk/(xl - ksize + 1L)
      neighbors <- RANN::nn2(zk, query = xk, k = min(50, nrow(zk) - 1L))
      neighbors$nn.dists <-  xl/(2 * ksize) * neighbors$nn.dists^2 ## linearize with JC69
      nnidxs <- match(neighbors$nn.idx[, 1], attr(tree, "pointers")) #NN indices in full set

      ###doNN <- TRUE # do nearest neighbor search
      ##########
      # for(i in seq_along(x)){
      #   attr(x[[i]], "idxs") <- neighbors$nn.idx[i, ]
      #   attr(x[[i]], "dists") <- neighbors$nn.dists[i, ]
      #   attr(x[[i]], "NN") <- nnidxs[i]
      # }
      ###td2 <- 2 * td
      for(i in seq_along(x)){
        if(neighbors$nn.dists[i, 1] <= td){
          idxs <- neighbors$nn.idx[i, neighbors$nn.dists[i,] <= td]
          ###alig <- aphid::align(c(x[i], z[idxs]), model = m, quiet = TRUE)
          ###dis <- if(DNA) ape::dist.dna(alig) else ape::dist.aa(alig)/100
          ###dis <- unname(dis[seq_along(idxs)])
          ###idxs <- idxs[dis <= td]
          #taxs <- as.integer(gsub(".+\\|", "", seqnames[idxs]))
          taxs <- unname(key[idxs]) ## length(key) = nrow(kmers)
          if(length(taxs) == 0L){
            attr(x[[i]], "NN") <-  nnidxs[i] #nnidxs[i] # sequence index
            attr(x[[i]], "NNhit") <- FALSE # -> do full classification procedure
          }else if(length(taxs) == 1L){
            attr(x[[i]], "NN") <- taxs # taxid
            attr(x[[i]], "NNhit") <- TRUE # -> skip full classification procedure
            ###neighbors$nn.dists[i, 1] <- dis[dis <= td]
            ###neighbors$nn.idx[i, 1] <- idxs[1]
          }else{
            lins <- get_lineage(taxs, db, numbers = TRUE)
            lins <- vapply(lins, paste0, "", collapse = "; ")
            anc <- .ancestor(lins)
            attr(x[[i]], "NN") <- as.integer(gsub(".+; ", "", anc)) # taxid
            attr(x[[i]], "NNhit") <- TRUE # -> skip full classification procedure
            ###wmd <- which.min(dis[dis <= td])
            ###neighbors$nn.dists[i, 1] <- min(dis)
            ###neighbors$nn.idx[i, 1] <- idxs[wmd]
          }
        }else{
          attr(x[[i]], "NN") <-  nnidxs[i] #nnidxs[i] # sequence index
          attr(x[[i]], "NNhit") <- FALSE # -> do full classification procedure
        }
      }

      #####
    }else if(ping == 1){ # faster hash matching
      xhash <- unname(hash(x))
      xmatch <- unname(key[xhash])
      for(i in seq_along(x)){
        if(is.na(xmatch[i])){
          ######cat(neighbors$nn.dists[i, 1], "\n")
          attr(x[[i]], "NN") <-  NA_integer_ # sequence index
          attr(x[[i]], "NNhit") <- FALSE # -> do full classification procedure
        }else{
          attr(x[[i]], "NN") <- xmatch[i]# taxid
          attr(x[[i]], "NNhit") <- TRUE # -> skip full classification procedure
        }
      }
      ## hash key
    }else{
      ###doNN <- FALSE
      for(i in seq_along(x)){
        attr(x[[i]], "NN") <- NA_integer_  #  nnidxs[i]# sequence index
        attr(x[[i]], "NNhit") <- FALSE # -> do full classification procedure
      }
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
  origx <- x
  if(!is.null(attr(tree, "numcode"))){
    if(is.null(attr(tree, "frame"))) stop("Amino classifier missing frame attribute\n")
    if(is.null(attr(tree, "remainder"))) stop("Amino classifier missing remainder attribute\n")
    for(i in seq_along(x)){
      if(!attr(x[[i]], "NNhit")){ # leave DNA hits intact (in DNAbin format -for mixed hashing)
        tmpattr <- attributes(x[[i]])
        if(length(x[[i]]) %% 3 != attr(tree, "remainder")){
          attr(x[[i]], "CODE8") <- TRUE ## invalid length
        }else{
          x[[i]] <- ape::as.character.DNAbin(x[[i]])
          x[[i]] <- seqinr::translate(x[[i]], numcode = attr(tree, "numcode"), frame = attr(tree, "frame"))
          x[[i]] <- ape::as.AAbin(x[[i]])
          attributes(x[[i]]) <- tmpattr
          if(any(x[[i]] == as.raw(42))) attr(x[[i]], "CODE9") <- TRUE ## contains stop codon(s)
        }
      }
    }
  }

  classify1 <- function(x, tree, threshold = 0.9, decay = TRUE, ping = TRUE, mincount = 3L){
    ## takes a single named raw vector with NN, NNhit, and possibly other attrs
    ## outputs a 1-row dataframe
    if(attr(x, "NNhit")){
      out <- data.frame(taxID = attr(x, "NN"),
                        score = NA_real_,
                        path = NA_character_,
                        scores = NA_character_,
                        reason = -1L,
                        stringsAsFactors = FALSE)
      return(out)
    }
    if(!is.null(attr(x, "CODE8")) | !is.null(attr(x, "CODE9"))){
      out <- data.frame(taxID = 1L,
                        score = NA_real_,
                        path = NA_character_,
                        scores = NA_character_,
                        reason = if(!is.null(attr(x, "CODE8"))) 8L else 9L,
                        stringsAsFactors = FALSE)
      return(out)
    }
    if(is.null(attr(x, "path"))){
      path <- ""
    }else{
      path <- attr(x, "path")
      indices <- gsub("([[:digit:]])", "[[\\1]]", path)
      eval(parse(text = paste0("tree <- tree", indices)))
    }
    scores <- if(is.null(attr(x, "scores"))) "" else attr(x, "scores")
    akw <- 1
    cakw <- if(is.null(attr(x, "cakw"))) 1 else attr(x, "cakw")
    tax <- if(is.null(attr(x, "tax"))) 1L else attr(x, "tax")
    #tax <- 1L # root (cant use 0L due to get_lineage call below)
    threshold_met <- minscore_met <- minlength_met <- maxlength_met <- neighbor_check <- mincount_met <- TRUE
    while(is.list(tree)){
      no_mods <- length(tree)
      sc <- numeric(no_mods)##### + 1L)##### # scores (log probabilities)
      for(i in seq_len(no_mods)){
        modi <- decodePHMM(attr(tree[[i]], "model"))
        sc[i] <- aphid::forward.PHMM(modi, x, odds = FALSE)$score
      }
      #####sc[no_mods + 1L] <- attr(tree, "nullscore") ######
      ##### cat(sc, "\n") ###############
      total_score <- aphid::logsum(sc)
      akwgts <- exp(sc - total_score)
      best_model <- which.max(akwgts)
      newakw <- akwgts[best_model]
      newcakw <- newakw * cakw
      threshold_met <- threshold <= if(decay) newcakw else newakw
      #####if(best_model == no_mods + 1L) threshold_met <- FALSE ###### null model selected
      # 4 is approx asymtote for single bp change as n training seqs -> inf
      if(threshold_met){
        minscore_met <- sc[best_model] >= attr(tree[[best_model]], "minscore")
        minlength_met <- length(x) >= attr(tree[[best_model]], "minlength")
        maxlength_met <- length(x) <= attr(tree[[best_model]], "maxlength")
        neighbor_check <- TRUE
        # neighbor_check <- if(is.na(attr(x, "NN"))){
        #   TRUE
        # }else{
        #   attr(x, "NN") %in% attr(tree[[best_model]], "sequences")
        # }
        mincount_met <- attr(tree[[best_model]], "ntotal") >= mincount
      }else{
        minscore_met <- minlength_met <- maxlength_met <- neighbor_check <- mincount_met <- FALSE
      }
      if(!(threshold_met & minscore_met & minlength_met & maxlength_met & neighbor_check & mincount_met)) break
      path <- paste0(path, best_model)
      scores <- paste0(scores, intToUtf8(round(newakw * 100)))
      akw <- newakw
      cakw <- newcakw
      tree <- tree[[best_model]]
      #####attr(tree, "nullscore") <- aphid::logsum(sc[seq_len(no_mods)][-best_model]) #####
      tax <- attr(tree, "taxID")
      if(!is.null(attr(tree, "aaleaf"))){
        out <- data.frame(taxID = tax,
                          score = round(if(decay) cakw else akw, 4),
                          path = path,
                          scores = scores,
                          reason = 100L,
                          stringsAsFactors = FALSE)
        return(out)
      }
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
    }else if (!neighbor_check){
      5L
    }else if (!mincount_met){
      6L
    }else{
      0L
    }
    out <- data.frame(taxID = tax,
                      score = score,
                      path = path,
                      scores = scores,
                      reason = reason,
                      stringsAsFactors = FALSE)
    return(out)
  }
  classify2 <- function(cores, x, tree, threshold = 0.9, decay = TRUE, ping = TRUE, mincount = 3L){
    ## x is (usually) a DNAbin or AAbin list
    if(inherits(cores, "cluster")){
      ###if(doNN) x <- parallel::parLapply(cores, x, classify0, z, m, td, db, key, DNA)
      res <- parallel::parLapply(cores, x, classify1, tree, threshold, decay, ping, mincount)
    }else if(cores == 1){
      ###if(doNN) x <- lapply(x, classify0, z, m, td, db, key, DNA)
      res <- lapply(x, classify1, tree, threshold, decay, ping, mincount)
    }else{
      navailcores <- parallel::detectCores()
      if(identical(cores, "autodetect")) cores <- navailcores - 1
      if(!(mode(cores) %in% c("numeric", "integer"))) stop("Invalid 'cores' object")
      if(cores > 1){
        cl <- parallel::makeCluster(cores)
        ###if(doNN) x <- parallel::parLapply(cl, x, classify0, z, m, td, db, key, DNA)
        res <- parallel::parLapply(cl, x, classify1, tree, threshold, decay, ping, mincount)
        parallel::stopCluster(cl)
      }else{
        ###if(doNN) x <- lapply(x, classify0, z, m, td, db, key, DNA)
        res <- lapply(x, classify1, tree, threshold, decay, ping, mincount)
      }
    }
    res <- do.call("rbind", res)
    return(res)
  }

  ## next lines are to speed up aa sequence classifications
  hashes2 <- hash(x) # now only finding duplicate aa seqs
  pointers2 <- .point(hashes2)
  #xnames <- names(x) ## need to think about NN attrs?
  gc()

  res <- classify2(cores, x[!duplicated(hashes2)], tree, threshold, decay, ping, mincount)
  res <- res[pointers2, ] # rereplicating AA dupes

  if(any(res$reason == 100L)){ #only happens for hybrid trees
    aaleaf <- which(res$reason == 100L)
    for(i in aaleaf){
      tmpattr <- attributes(x[[i]])
      x[[i]] <- origx[[i]] # change from AA to DNA
      attributes(x[[i]]) <- tmpattr
      attr(x[[i]], "path") <- res$path[i]
      attr(x[[i]], "scores") <- res$scores[i]
      attr(x[[i]], "cakw") <- res$cakw[i]
      attr(x[[i]], "tax") <- res$tax[i]
    }
    res[aaleaf, ] <- classify2(cores, x[aaleaf], tree, threshold, decay, ping, mincount)
  }
  ## do DNA classifs hree

  #res <- do.call("rbind", res) ## changed output to dataframe 20180617
  lineages <- get_lineage(res$taxID, db, cores = cores,
                          simplify = FALSE, numbers = FALSE)
  tmp <- sapply(lineages, tail, 1)
  lhcols <- data.frame(representative = names(x),
                     taxID = res$taxID,
                     taxon = tmp,
                     rank = names(tmp),
                     score = res$score,
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
  if(metadata){
    lhcols <- cbind(lhcols, res[c("path", "scores", "reason")])
    if(exists("neighbors")){
      taxids <- unname(key[neighbors$nn.idx[, 1]])
      dbinds <- match(taxids, db$taxID)
      if(any(is.na(dbinds))) stop("Error code 9283\n")
      nns <- cbind(taxids, db$name[dbinds], db$rank[dbinds], round(neighbors$nn.dists[, 1], 4))
      colnames(nns) <- c("NNtaxID", "NNtaxon", "NNrank","NNdistance")
      lhcols <- cbind(lhcols, nns)
    }
  }
  if(tabulize){
    lhcols <- cbind(lhcols, qout)
  }else{
    lhcols <- lhcols[pointers, ]
    lhcols$representative <- catchnames
  }
  rownames(lhcols) <- NULL
  return(lhcols)
}
################################################################################


# classify0 <- function(xx, z, m, td, db, key, DNA){ ## parallelizable
#   stopifnot(!is.null(attr(xx, "dists")))
#   stopifnot(!is.null(attr(xx, "idxs")))
#   stopifnot(!is.null(attr(xx, "NN")))
#   if(attr(xx, "dists")[1] <= 2 * td){
#     idxs <- attr(xx, "idxs")[attr(xx, "dists") <= 2 * td]
#     alig <- aphid::align(c(list(xx), z[idxs]), model = m, quiet = TRUE)
#     dis <- if(DNA) ape::dist.dna(alig) else ape::dist.aa(alig)/100
#     dis <- unname(dis[seq_along(idxs)])
#     idxs <- idxs[dis <= td]
#     #taxs <- as.integer(gsub(".+\\|", "", seqnames[idxs]))
#     taxs <- unname(key[idxs]) ## length(key) = nrow(kmers)
#     if(length(taxs) == 0L){
#       ##attr(xx, "NN") <- attr(xx, "nnidx") #nnidxs[i] # sequence index, no rplcmnt needed
#       attr(xx, "NNhit") <- FALSE # -> do full classification routine
#       attr(xx, "dists") <- attr(xx, "dists")[1]
#       attr(xx, "idxs") <- attr(xx, "idxs")[1]
#     }else if(length(taxs) == 1L){
#       attr(xx, "NN") <- taxs # taxid
#       attr(xx, "NNhit") <- TRUE # -> skip full classification routine
#       attr(xx, "dists") <- dis[dis <= td]
#       attr(xx, "idxs") <- idxs[1]
#     }else{
#       lins <- get_lineage(taxs, db, numbers = TRUE)
#       lins <- vapply(lins, paste0, "", collapse = "; ")
#       anc <- .ancestor(lins)
#       attr(xx, "NN") <- as.integer(gsub(".+; ", "", anc)) # taxid
#       attr(xx, "NNhit") <- TRUE # -> skip full classification routine
#       attr(xx, "dists") <- min(dis[dis <= td])
#       attr(xx, "idxs") <- idxs[which.min(dis[dis <= td])]
#     }
#   }else{
#     #attr(x[[i]], "NN") <-  nnidxs[i] #nnidxs[i] # sequence index
#     attr(xx, "NNhit") <- FALSE # -> do full classification routine
#     attr(xx, "dists") <- attr(xx, "dists")[1]
#     attr(xx, "idxs") <- attr(xx, "idxs")[1]
#   }
#   return(xx)
# }
