#' Informatic sequence classification tree learning.
#'
#' This function learns a classification tree from a reference sequence database
#'   using a recursive partitioning procedure.
#'
#' @param x an object of class\code{"DNAbin"} representing a list of
#'   DNA sequences to be used as the training data for the tree-learning process.
#'   All sequences should be from the same genetic region of interest
#'   and be globally alignable (i.e. without unjustified end-gaps).
#'   The sequences must have "names" attributes that include
#'   taxonomic ID numbers corresponding with those in the taxonomy
#'   database \code{db} (separated from the sequence ID by a "|" character).
#'   For example: "AF296347|30962", "AF296346|8022", "AF296345|8017", etc.
#'   See \code{\link{searchGB}} for more details on creating the reference
#'   sequence database.
#' @param db a valid taxonomy database (as a data.frame object).
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
#'   in the PHMMs. If NULL (default) no maximum size is enforced.
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
#'   to the console. Note that the output can be verbose.
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
learn <- function(x, db, model = NULL, refine = "Viterbi", iterations = 50,
                  nstart = 20, minK = 2, maxK = 2, minscore = 0.9, probs = 0.5,
                  retry = TRUE, resize = TRUE, maxsize = max(sapply(x, length)),
                  recursive = TRUE, cores = 1, quiet = TRUE, ...){
  if(mode(x) == "character") x <- char2dna(x)
  if(!grepl("\\|", names(x)[1])){
    stop("Names of input sequences must include taxonomic ID numbers\n")
  }
  taxIDs <- as.integer(gsub(".+\\|", "", names(x)))
  lineages <- get_lineage(taxIDs, db = db, cores = cores, numbers = TRUE)
  attr(x, "lineage") <- vapply(lineages, paste0, "", collapse = "; ")
  ## Initialize the tree as a dendrogram object
  tree <- 1
  attr(tree, "leaf") <- TRUE
  attr(tree, "height") <- 0
  attr(tree, "midpoint") <- 0
  attr(tree, "members") <- 1
  class(tree) <- "dendrogram"
  attr(tree, "taxonomy") <- prune_taxonomy(db, taxIDs = taxIDs, keep = TRUE)
  attr(tree, "clade") <- ""
  attr(tree, "sequences") <- seq_along(x)
  attr(tree, "lineage") <- .ancestor(attr(x, "lineage"))
  attr(tree, "minscore") <- -1E06 # nominal
  xlengths <- sapply(x, length)
  attr(tree, "minlength") <- min(xlengths)
  attr(tree, "maxlength") <- max(xlengths)
  if(is.null(attr(x, "hashes"))) attr(x, "hashes") <- hash(x)
  if(is.null(attr(x, "duplicates"))) attr(x, "duplicates") <- duplicated(attr(x, "hashes"))
  if(is.null(attr(x, "pointers"))) attr(x, "pointers") <- .point(attr(x, "hashes"))
  if(is.null(attr(x, "weights"))){
    if(!quiet) cat("Deriving sequence weights\n")
    attr(x, "weights") <- aphid::weight(x, k = 5)
  }
  ancestors <- character(max(attr(x, "pointers")))
  names(ancestors) <- attr(x, "hashes")[!attr(x, "duplicates")]
  for(i in seq_along(ancestors)){
    ancestors[i] <- .ancestor(attr(x, "lineage")[attr(x, "hashes") == names(ancestors)[i]])
  }
  tmpnames <- names(ancestors)
  ancestors <- as.integer(gsub(".+; ", "", ancestors))
  names(ancestors) <- tmpnames
  attr(tree, "key") <- ancestors# for exact matching
  if(is.null(model)){
    if(!quiet) cat("Dereplicating sequences\n")
    xu <- x[!attr(x, "duplicates")]
    xuw <- sapply(split(attr(x, "weights"), attr(x, "pointers")), sum)
    if(!quiet) cat("Deriving top level model\n")
    if(length(xu) > 1000){
      # progressive model training to avoid excessive memory usage
      samp <- sample(seq_along(xu), size = 1000)
      model <- aphid::derivePHMM(xu[samp], refine = refine,
                                 seqweights = xuw[samp], maxsize = maxsize,
                                 inserts = "inherited", alignment = FALSE,
                                 quiet = TRUE, cores = cores, maxiter = 20)
      model <- aphid::train(model, xu, method = refine, seqweights = xuw,
                            inserts = "inherited", alignment = FALSE,
                            cores = cores, quiet = TRUE, maxiter = 20)
    }else{
      model <- aphid::derivePHMM(xu, refine = refine, seqweights = xuw,
                                 inserts = "inherited", alignment = FALSE,
                                 quiet = TRUE, cores = cores, maxiter = 20)
    }
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
  tree <- expand(tree, x, clades = "", refine = refine, iterations = iterations,
                 nstart = nstart, minK = minK, maxK = maxK, minscore = minscore,
                 probs = probs, retry = retry, resize = resize, maxsize = maxsize,
                 recursive = recursive, cores = cores, quiet = quiet, ... = ...)
  return(tree)
}
################################################################################


