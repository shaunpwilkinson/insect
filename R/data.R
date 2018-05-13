#' Cetacean 16S rDNA sequences.
#'
#' A dataset containing 19 mitochondrial 16S rDNA sequences from 18 cetacean species,
#'   downloaded from GenBank on 27 March 2018.
#'
#' @format A "DNAbin" list object containing 19 sequences in raw-byte format,
#' averaging 130 nucleotides in length. The object also contains the following attributes,
#' each of which is a vector the same length as the sequence list:
#' \describe{
#'   \item{"names"}{GenBank accession numbers of the original sequences.}
#'   \item{"taxID"}{NCBI taxonomy database identifiers.}
#'   \item{"lineage"}{semicolon-delimited lineage strings.}
#'   \item{"species"}{species names.}
#' }
#' The GenBank accession numbers given in the "names" attribute correspond to
#' complete mitochondrial genome sequences that were downloaded using the
#' \code{\link{searchGB}} function on 27 march 2018
#' (query term: "cetacea[ORGN]+AND+16S+rRNA[GENE]+AND+1970:2017[MDAT]").
#' The sequences were then trimmed using the
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
#'   \item{tax_id}{the NCBI unique taxon identifier (integer).}
#'   \item{parent_tax_id}{the NCBI unique taxon identifier
#'     of the immediate parent taxon (integer).}
#'   \item{rank}{The taxonomic rank (i.e. species, genus, etc; character).}
#'   \item{name}{The scientific name of the taxon (character).}
#' }
#' The database was accessed from
#' \url{ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz} on 27 March 2018
#' using the \code{\link{download_taxon}} function, and pruned using
#' \code{\link{prune_taxon}} with \code{taxIDs = attr(whales, "taxID")}.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/taxonomy/}
#'
"whale_taxa"
################################################################################
