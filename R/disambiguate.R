#' Convert oligonucleotide sequences into regular expressions.
#'
#' This function is used to convert an oligonucleotide sequence into a regular
#'   expression that can be used to query a sequence dataset.
#'   This is generally used for finding and removing primer,
#'   adapter and/or index sequences.
#'
#' @param z a concatenated string representing a DNA oligonucleotide sequence,
#'   possibly with IUPAC ambiguity codes (all in upper case).
#' @return a regular expression.
#' @author Shaun Wilkinson
#' @examples disambiguate("GGWACWGGWTGAACWGTWTAYCCYCC")
################################################################################
disambiguate <- function(z){
  if(mode(z) != "character") stop ("Expected a concatenated character string\n")
  z <- gsub(" ", "", z)
  z <- toupper(z)
  if(grepl("[^ACGTMRWSYKVHDBNI]", z)) stop("Sequence contains invalid residues\n")
  if(grepl("[MRWSYKVHDBNI]", z)){
    z <- gsub("M", "[AC]", z)
    z <- gsub("R", "[AG]", z)
    z <- gsub("W", "[AT]", z)
    z <- gsub("S", "[CG]", z)
    z <- gsub("Y", "[CT]", z)
    z <- gsub("K", "[GT]", z)
    z <- gsub("V", "[ACG]", z)
    z <- gsub("H", "[ACT]", z)
    z <- gsub("D", "[AGT]", z)
    z <- gsub("B", "[CGT]", z)
    z <- gsub("N|I", "[ACGT]", z)
  }
  return(z)
}
################################################################################
