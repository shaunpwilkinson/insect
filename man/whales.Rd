% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{whales}
\alias{whales}
\title{Cetacean 16S rDNA sequences.}
\format{
A "DNAbin" list object containing 19 cetacean mitochondrial sequences
in raw-byte format,
averaging 130 nucleotides in length. The sequences are named with the GenBank
accession numbers followed by a "|" symbol, followed by their NCBI taxonomy ID
numbers.
The sequences were downloaded using the
\code{\link{searchGB}} function on 17 June 2018
(query term: "cetacea[ORGN]+AND+16S+rRNA[GENE]"), and trimmed using the
\code{\link{virtualPCR}} function with the primers 16Smam1 and 16Smam2
(CGGTTGGGGTGACCTCGGA and GCTGTTATCCCTAGGGTAACT, respectively; Taylor 1996).
}
\source{
\url{https://www.ncbi.nlm.nih.gov/genbank/}
}
\usage{
whales
}
\description{
A dataset containing 19 mitochondrial 16S rDNA sequences from 18 cetacean species,
  downloaded from GenBank on 27 March 2018.
}
\references{
Taylor PG (1996) Reproducibility of ancient DNA sequences from extinct
  Pleistocene fauna. \emph{Molecular Biology and Evolution}, \strong{13}, 283-285.
}
\keyword{datasets}
