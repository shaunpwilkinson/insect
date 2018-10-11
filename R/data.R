#' Cetacean 16S rDNA sequences.
#'
#' A dataset containing 19 mitochondrial 16S rDNA sequences from 18 cetacean species,
#'   downloaded from GenBank on 27 March 2018.
#'
#' @format A "DNAbin" list object containing 19 cetacean mitochondrial sequences
#' in raw-byte format,
#' averaging 130 nucleotides in length. The sequences are named with the GenBank
#' accession numbers followed by a "|" symbol, followed by their NCBI taxonomy ID
#' numbers.
#' The sequences were downloaded using the
#' \code{\link{searchGB}} function on 17 June 2018
#' (query term: "cetacea[ORGN]+AND+16S+rRNA[GENE]"), and trimmed using the
#' \code{\link{virtualPCR}} function with the primers 16Smam1 and 16Smam2
#' (CGGTTGGGGTGACCTCGGA and GCTGTTATCCCTAGGGTAACT, respectively; Taylor 1996).
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/genbank/}
#' @references
#'   Taylor PG (1996) Reproducibility of ancient DNA sequences from extinct
#'   Pleistocene fauna. \emph{Molecular Biology and Evolution}, \strong{13}, 283-285.
#'
"whales"
################################################################################
#' Cetacean section of NCBI taxonomy database.
#'
#' A copy of the NCBI taxonomy reference database, subsetted to include only the
#'   cetacean taxa in the \code{\link{whales}} dataset.
#'
#' @format A data.frame object with 72 rows and four columns, labeled as follows:
#' \describe{
#'   \item{taxID}{the NCBI unique taxon identifier (integer).}
#'   \item{parent_taxID}{the NCBI unique taxon identifier
#'     of the immediate parent taxon (integer).}
#'   \item{rank}{The taxonomic rank (i.e. species, genus, etc; character).}
#'   \item{name}{The scientific name of the taxon (character).}
#' }
#' The database was accessed from
#' \url{ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz} on 17 June 2018
#' using the \code{\link{taxonomy}} function, and pruned using
#' \code{\link{prune_taxonomy}} with
#' \code{taxIDs = as.integer(gsub(".+\\|", "", names(whales)))}.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/taxonomy/}
#'
"whale_taxonomy"
################################################################################
#' Table of marine COI amplicon sequence variants from American Samoa
#'
#'  This matrix contains counts of COI amplicon sequence variants (ASV)
#'  extracted from autonomous reef monitoring structures (ARMS) in
#'  Ofu, American Samoa. Unpublished data courtesy of Molly Timmers (NOAA).
#'
#' @format a 2 x 16 integer matrix containing abundance counts of 16 ASVs from
#'   two sites. This table contains the first 16 rows of the 'seqtab.nochim' output from
#'   the DADA2 pipeline (\url{https://benjjneb.github.io/dada2/tutorial.html}).
#'   The column names contain the DNA sequences of the ASVs, and row names
#'   correspond with site codes.
#'
"samoa"
################################################################################

