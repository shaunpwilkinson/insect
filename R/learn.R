#' Informatic sequence classification tree learning.
#'
#' This function learns a classification tree from a reference sequence database
#'   using a recursive partitioning procedure.
#'
#' @param x a reference database of class\code{"DNAbin"} representing a list of
#'   DNA sequences to be used as the training data.
#'   All sequences should be from the same genetic region of interest
#'   and be globally alignable (i.e. without unjustified end-gaps).
#'   The sequences must have "names" attributes, either in RDP format
#'   (containing semicolon-delimited lineage strings),
#'   or that include taxonomic ID numbers corresponding with those in the taxonomy
#'   database \code{db} (separated from the sequence ID by a "|" character).
#'   For example: "AF296347|30962", "AF296346|8022", "AF296345|8017", etc.
#'   See \code{\link{searchGB}} for more details on creating the reference
#'   sequence database and \code{\link{taxonomy}} for the associated heirarchical
#'   taxonomic database.
#' @param db a heirarchical taxonomy database in the form of a data.frame.
#'   Cannot be NULL unless training data is in RDP format
#'   (containing semicolon delimited lineage strings).
#'   The object should have
#'   four columns, labeled "taxID", "parent_taxID", "rank" and "name".
#'   The first two should be numeric, and all ID numbers in the
#'   "parent_taxID" column should link to those in the "taxID" column.
#'   This excludes the first row,
#'   which should have \code{parent_taxID = 0} and \code{name = "root"}.
#'   See \code{\link{taxonomy}} for more details.
#' @param model an optional object of class \code{"PHMM"} providing the
#'   starting parameters. Used to train (optimize parameters for)
#'   subsequent nested models to be positioned at successive
#'   sub-nodes. If NULL, the root model is derived from the
#'   sequence list prior to the recursive partitioning process.
#' @param refine character string giving the iterative model refinement
#'   method to be used in the partitioning process. Valid options are
#'   \code{"Viterbi"} (Viterbi training; the default option) and
#'   \code{"BaumWelch"} (a modified version of the Expectation-Maximization
#'   algorithm).
#' @param iterations integer giving the maximum number of training-classification
#'   iterations to be used in the splitting process.
#'   Note that this is not necessarily the same as the number of Viterbi training
#'   or Baum Welch iterations to be used in model training, which can be set
#'   using the argument \code{"maxiter"} (eventually passed to
#'   \code{\link[aphid]{train}} via the dots argument "...").
#' @param nstart integer. The number of random starting sets to be chosen
#'   for initial k-means assignment of sequences to groups. Defaults to 20.
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
#'   terminate. As an example, if \code{minscore = 0.9} and
#'   \code{probs = 0.5} (the default settings), and after generating two
#'   candidate PHMMs to occupy the candidate subnodes the median
#'   of Akaike weights is 0.89, the splitting process will
#'   terminate and the function will simply return the unsplit root node.
#' @param probs numeric between 0 and 1. The percentile of Akaike weights
#'   to test against the minimum score threshold given in \code{"minscore"}.
#' @param retry logical indicating whether failure to split a node based on
#'   the criteria outlined in 'minscore' and 'probs' should prompt a second
#'   attempt with different initial groupings. These groupings are based on
#'   maximum kmer frequencies rather than k-means division, which can give
#'   suboptimal groupings when the cluster sizes are different (due to
#'   the up-weighting of larger clusters in the k-means algorithm).
#' @param resize logical indicating whether the models should be free to
#'   change size during the training process or if the number of modules
#'   should be fixed. Defaults to TRUE. Only applicable if
#'   \code{refine = "Viterbi"}.
#' @param maxsize integer giving the upper bound on the number of modules
#'   in the PHMMs. If NULL, no maximum size is enforced.
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
#'   to the console.
#' @param verbose logical indicating whether extra feedback should be
#'   printed to the console, including progress at each split.
#' @param numcode,frame passed to \code{\link[seqinr]{translate}}.
#'   Set to NULL (default) unless an amino acid sequence classifier is required.
#' @param ... further arguments to be passed on to \code{\link[aphid]{train}}).
#' @return an object of class \code{"insect"}.
#' @details The "insect" object type is a dendrogram
#'   with several additional attributes stored at each node.
#'   These include:
#'   "clade" the index of the node (see further details below);
#'   "sequences" the indices of the sequences in the reference
#'   database used to create the object;
#'   "taxID" the taxonomic identifier of the lowest common taxon
#'   of the sequences belonging to the node (linking to \code{"db"});
#'   "minscore" the lowest likelihood among the training sequences given
#'   the profile HMM stored at the node;
#'   "minlength" the minimum length of the sequences belonging to the node;
#'   "maxlength" the maximum length of the sequences belonging to the node;
#'   "model" the profile HMM derived from the sequence subset belonging to the node;
#'   "nunique" the number of unique sequences belonging to the node;
#'   "ntotal" the total number of sequences belonging to the node (including duplicates);
#'   "key" the hash key used for exact sequence matching
#'   (bypasses the classification procedure if an exact match is found; root node only);
#'   "taxonomy" the taxonomy database containing the taxon ID numbers (root node only).
#'
#'   The clade indexing system used here is based on character strings,
#'   where "0" refers to the root node,
#'   "01" is the first child node, "02" is the second child node,
#'   "011" is the first child node of the first child node, etc.
#'   The leading zero may be omitted for brevity.
#'   Note that each inner node can not have more than 9 child nodes.
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
learn <- function(x, db = NULL, model = NULL, refine = "Viterbi", iterations = 50,
                  nstart = 20, minK = 2, maxK = 2, minscore = 0.9, probs = 0.5,
                  retry = TRUE, resize = TRUE, maxsize = 1000,
                  recursive = TRUE, cores = 1, quiet = FALSE, verbose = FALSE,
                  numcode = NULL, frame = NULL, ...){

  if(!quiet) cat("Training classifier\n")
  if(mode(x) == "character") x <- char2dna(x)
  dots <- list(...)
  if(is.null(dots$k)){
    ksize <- if(inherits(x, "DNAbin")) 4 else 2
  }else{
    if(inherits(x, "DNAbin")){
      if(dots$k > 6) stop("Maximum kmer size of 6 for DNA\n")
      ksize <- dots$k
    }else{
      if(dots$k > 3) stop("Maximum kmer size of 3 for amino acid sequences\n")
      ksize <- dots$k
    }
  }
  if(is.null(db)){
    if(!all(grepl(";", names(x)))){
      stop("Training data must be in RDP format if no taxonomy database is provided\n")
    }else{
      if(!quiet) cat("Building heirarchical taxonomy database\n")
    }
    if(all(grepl("\\\t", names(x)))){
      xnames <- gsub("\\\t.+", "", names(x))
      xlins <- gsub(".+\\\t", "", names(x))
    }else{
      xnames <- paste0("S", seq_along(x))
      xlins <- names(x)
    }
    if(!grepl("^[Rr]oot;", xlins[1])){
      xlins <- paste0("Root;", xlins)
    }
    xlins <- gsub(";$", "", xlins)
    xlins <- gsub("; ", ";", xlins)
    xs <- strsplit(xlins, split = ";") #split
    linlens <- vapply(xs, length, 0L)
    linlen <- max(linlens)
    discards <- linlens < linlen
    if(any(discards)){
      warning("Incomplete lineages detected. Removing", sum(discards), "training sequences")
      x <- x[!discards]
      if(length(x) < 10) stop("Formatting error\n")
    }
    xsa <- lapply(xs, function(l) paste0(l, "__", seq_along(l) - 1)) #appended
    taxa <- as.data.frame(xsa, stringsAsFactors = FALSE)
    taxa <- t(as.matrix(taxa))
    rownames(taxa) <- NULL
    newmat <- matrix(NA_integer_, nrow = nrow(taxa), ncol = ncol(taxa))
    newmat[, 1] <- 1L # taxID for Root
    for(i in seq(2, linlen)){
      taxa[, i] <- paste0(newmat[, i - 1], ";", taxa[, i])
      newmat[, i] <- .point(taxa[, i]) + max(newmat[, i - 1])
      taxa[, i] <- paste0(newmat[, i], ";", taxa[, i])
    }
    db <- unique(as.vector(taxa))
    db[1] <- paste0("1;0;", db[1]) #root
    db <- gsub("__([[:digit:]]+)$", paste0(";rank", "\\1"), db)
    db <- strsplit(db, split = ";")
    db <- t(as.matrix(as.data.frame(db, stringsAsFactors = FALSE)))
    rownames(db) <- NULL
    db <- db[, c(1, 2, 4, 3)]
    colnames(db) <- c("taxID", "parent_taxID", "rank", "name")
    db <- as.data.frame(db, stringsAsFactors = FALSE)
    db$taxID <- as.integer(db$taxID)
    db$parent_taxID <- as.integer(db$parent_taxID)
    taxIDs <- newmat[, linlen]
    lineages <- apply(newmat, 1, paste0, collapse = "; ")
    names(x) <- paste0(xnames, "|", taxIDs)
  }else{
    if(!quiet) cat("Converting taxon IDs to full lineage strings\n")
    if(!grepl("\\|", names(x)[1])){
      stop("Names of input sequences must include taxonomic ID numbers\n")
    }
    taxIDs <- as.integer(gsub(".+\\|", "", names(x)))
    lineages <- get_lineage(taxIDs, db = db, cores = cores, numbers = TRUE)
    lineages <- vapply(lineages, paste0, "", collapse = "; ")
    db <- prune_taxonomy(db, taxIDs = taxIDs, keep = TRUE)
  }

  if(!quiet) cat("Initializing tree\n")
  tree <- 1
  attr(tree, "k") <- ksize #new
  attr(tree, "leaf") <- TRUE
  attr(tree, "height") <- 0
  attr(tree, "midpoint") <- 0
  attr(tree, "members") <- 1
  class(tree) <- "dendrogram"
  attr(tree, "taxonomy") <- db
  attr(tree, "clade") <- ""
  amino <- FALSE
  if(!is.null(numcode)){ #make amino acid classifier
    if(is.null(frame)) stop("Argument provided for numcode but not frame\n")
    amino <- TRUE
    xlengths <- vapply(x, length, 0L, USE.NAMES = FALSE)
    xrems <- xlengths %% 3
    remtab <- sort(table(xrems), decreasing = TRUE)
    remainder <- as.integer(names(remtab)[1]) # expected remainder
    discards <- xrems != remainder
    if(any(discards)) warning(sum(discards), " sequences are not translatable and will be removed from the trainingset")
    x <- x[!discards]
    lineages <- lineages[!discards]
    if(length(x) < 10) stop("Too few training sequences remain\n")
    if(!quiet) cat("Translating sequences\n")
    xaa <- ape::as.character.DNAbin(x)
    xaa <- lapply(xaa, seqinr::translate, numcode = numcode, frame = frame)
    discards <- sapply(xaa, function(s) !all(s %in% LETTERS[-c(2, 10, 15, 21, 24, 26)]))
    nd <- sum(discards)
    if(nd > 0) warning(nd, " sequences contain stop codons/ambiguities and will be removed from the trainingset")
    stopifnot(nd < length(xaa))
    xaa <- ape::as.AAbin(xaa[!discards])
    #keeps <- sapply(xaa, function(v) !any(v == as.raw(42)))
    x <- x[!discards]
    lineages <- lineages[!discards]
    if(length(x) < 10) stop("Too few training sequences remaining\n")
    attr(tree, "numcode") <- numcode
    attr(tree, "frame") <- frame
    attr(tree, "remainder") <- remainder
    xaa <- dereplicate(xaa)
    attr(tree, "xaa") <- xaa
  }
  attr(tree, "sequences") <- seq_along(x)
  #attr(tree, "seqnames") <- names(x) # cant just name sequences attr due to derep-rerep
  if(!quiet) cat("Dereplicating sequences")
  tset <- dereplicate(x)
  if(!quiet) cat(", found", length(tset),"unique sequences\n")
  if(!quiet) cat("Calculating sequence weights\n")
  attr(tset, "rerep.weights") <- aphid::weight(tset, method = "Henikoff", k = 5)
  #if(!quiet) cat("Clustering OTUs\n")
  #otus <- .otu(tset, threshold = 0.97) # for increased partitioning speed at top levels
  #if(!quiet) cat("Found", max(otus), "OTUs\n")
  #names(tset) <- paste0("OTU", otus) ## S1, S2*, S3 etc, with asterisks showing central seqs
  #centralseqs <- grepl("\\*$", names(otus))
  #names(tset)[centralseqs] <- paste0(names(tset)[centralseqs], "*")
  attr(tset, "lineages") <- lineages #full length, passed to 'expand' to save time
  attr(tree, "lineage") <- .ancestor(lineages) # eventually removed during expansion
  attr(tree, "minscore") <- -1E06 # nominal
  xlengths <- vapply(tset, length, 0L, USE.NAMES = FALSE)
  attr(tree, "seqlengths") <- xlengths #used to normalize kmer matrix in `classify``

  # if(is.null(attr(x, "hashes"))) attr(x, "hashes") <- hash(x)
  # if(is.null(attr(x, "duplicates"))) attr(x, "duplicates") <- duplicated(attr(x, "hashes"))
  # if(is.null(attr(x, "pointers"))) attr(x, "pointers") <- .point(attr(x, "hashes"))
  attr(tree, "pointers") <- attr(tset, "rerep.pointers") #new

  if(!quiet) cat("Making hash key for exact sequence matching\n")
  hashes <- hash(x)
  #duplicates <- duplicated(hashes)
  ancestors <- split(lineages, f = factor(hashes, levels = unique(hashes)))
  anclens <- vapply(ancestors, length, 0L, USE.NAMES = FALSE)
  ancestors[anclens > 1] <- lapply(ancestors[anclens > 1], .ancestor)
  tmpnames <- names(ancestors)
  ancestors <- as.integer(gsub(".+; ", "", ancestors))
  names(ancestors) <- tmpnames
  attr(tree, "key") <- ancestors # for exact matching - order matches kmer matrix row wise

  ## same length as hash key, nrow(kmers)
  if(!quiet) cat("Counting ", ksize, "-mers\n", sep = "")
  attr(tree, "kmers") <- .encodekc(kmer::kcount(tset, k = ksize)) #new
  # attr(tree, "kmers") <- kmer::kcount(x[!attr(x, "duplicates")],
  #                                              k = if(is.null(dots$k)) 5 else dots$k)

  if(amino){
    x <- xaa
    xlengths <- vapply(x, length, 0L, USE.NAMES = FALSE)
  }else{
    x <- tset
  }
  attr(tree, "minlength") <- min(xlengths)
  attr(tree, "maxlength") <- max(xlengths)

  if(is.null(model)){
    if(!quiet) cat("Dereplicating sequences\n")
    if(!quiet) cat("Deriving top level model\n")
    if(length(x) > 100){
      # sample model training to avoid excessive memory usage
      samp <- sample(seq_along(x), size = 2 * log(length(x), 2)^2)
      x <- x[samp]
    }
    suppressWarnings(
      model <- aphid::derivePHMM(x, refine = refine,
                                 #seqweights = attr(tset, "derep.weights")[samp],
                                 maxsize = maxsize,
                                 inserts = "inherited", alignment = FALSE,
                                 quiet = TRUE, cores = cores, maxiter = 20,
                                 limit = 0.98, k = if(inherits(x, "DNAbin")) 4 else 2)
    )

    ## strip memory intensive elements but not alignment yet
    model$weights <- NULL
    model$mask <- NULL
    model$map <- NULL
    model$reference <- NULL
    model$insertlengths <- NULL
    model$name <- NULL
  }else if(mode(model) == "raw"){
    model <- decodePHMM(model)
  }
  attr(tree, "model") <- model
  attr(tree, "trainingset") <- tset # includes full numbered lineages
  rm(tset)
  rm(x)
  if(exists("xaa")) rm(xaa)
  if(!quiet) cat("Recursively splitting nodes\n")
  tree <- expand(tree, clades = "", refine = refine, iterations = iterations,
                 nstart = nstart, minK = minK, maxK = maxK, minscore = minscore,
                 probs = probs, retry = retry, resize = resize, maxsize = maxsize,
                 recursive = recursive, cores = cores, quiet = quiet,
                 verbose = verbose, ... = ...)

  if(!is.null(numcode) & recursive){
    aaleaf <- function(node){
      if(is.leaf(node)) attr(node, "aaleaf") <- TRUE
      return(node)
    }
    tree <- dendrapply(tree, aaleaf)
    tmpnc <- attr(tree, "numcode")
    attr(tree, "numcode") <- NULL
    if(!quiet) cat("Transitioning from AA to DNA models\n")
    tmpf <- tempfile(fileext = ".rds")
    if(!quiet & verbose) cat("Saving AA classifier as ", tmpf, "\n")
    saveRDS(tree, file = tmpf)
    tree <- expand(tree, clades = "", refine = refine, iterations = iterations,
                   nstart = nstart, minK = minK, maxK = maxK, minscore = minscore,
                   probs = probs, retry = retry, resize = resize, maxsize = maxsize,
                   recursive = recursive, cores = cores, quiet = quiet,
                   verbose = verbose, ... = ...)
    if(!quiet & verbose) cat("AA classifier saved as ", tmpf, "\n")
    attr(tree, "numcode") <- tmpnc
  }
  attr(attr(tree, "trainingset"), "lineages") <- NULL # no longer required
  attr(tree, "xaa") <- NULL # no longer required
  if(!quiet) cat("Done\n")
  return(tree)
}
################################################################################


