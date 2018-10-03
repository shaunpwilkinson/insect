#' Encode and decode profile HMMs in raw byte format.
#'
#' These functions are used to compress and decompress profile
#'   hidden Markov models for DNA to improve memory efficiency.
#'
#' @param x an object of class "PHMM"
#' @param z a raw vector in the encodePHMM schema.
#' @return encodePHMM returns a raw vector. \code{decodePHMM} returns
#'   an object of class "PHMM" (see Durbin et al (1998) and
#'   the \code{\link[aphid]{aphid}}
#'   package for more details
#'   on profile hidden Markov models).
#' @details Profile HMMs used in tree-based classification
#'   usually include many parameters, and hence large trees with
#'   many PHMMs can occupy a lot of memory. Hence a basic encoding
#'   system was devised to store the emission and transition probabilities
#'   in raw-byte format to three (nearly four) decimal places.
#'   This does not seem to significantly affect the accuracy of likelihood scoring,
#'   and has a moderate impact on classification speed, but can
#'   reduce the memory allocation requirements for large trees by up to
#'   95 percent.
#' @author Shaun Wilkinson
#' @references
#'   Durbin R, Eddy SR, Krogh A, Mitchison G (1998) Biological
#'   sequence analysis: probabilistic models of proteins and nucleic acids.
#'   Cambridge University Press, Cambridge, United Kingdom.
#' @examples
#' \donttest{
#'   ## generate a simple classification tree with two child nodes
#'   data(whales)
#'   data(whale_taxonomy)
#'   tree <- learn(whales, db = whale_taxonomy, recursive = FALSE)
#'   ## extract the omnibus profile HMM from the root node
#'   PHMM0 <- decodePHMM(attr(tree, "model"))
#'   ## extract the profile HMM from the first child node
#'   PHMM1 <- decodePHMM(attr(tree[[1]], "model"))
#'  }
#' @name encoding
################################################################################
encodePHMM <- function(x){
  if(mode(x) == "raw") return(x)
  AA <- nrow(x$E) == 20
  if(!AA & nrow(x$E) != 4) stop("Invalid model for DNA or AA sequences \n")

  logibits <- raw(32)
  if(any(is.finite(x$A[3, ]))) logibits[32] <- as.raw(1) # DI trans enabled?
  if(any(is.finite(x$A[7, ]))) logibits[31] <- as.raw(1) # ID trans enabled?
  if(AA) logibits[30] <- as.raw(1)
  logibytes <- packBits(logibits)
  sizebytes <- packBits(intToBits(x$size)) # 4 bytes for model size = 4bil max size
  A <- x$A[-(c(2, 3, 6, 7, 9)), ]
  A <- A[is.finite(A)] * -1 # now a vector
  A[A < 1e-06] <- 1e-06 ## minimum value (max prob remember to fix final DM trans to 1)
  Abytes <- as.vector(sapply(A, .encode1))
  E <- if(AA) x$E[1:19, ] * -1 else x$E[1:3, ] * -1 # drop T row as can be calculated
  E[E < 1e-06] <- 1e-06
  Ebytes <- as.vector(sapply(E, .encode1))
  qa <- x$qa[-c(3, 7)] * -1 # can't drop any
  qabytes <- as.vector(sapply(qa, .encode1))
  qe <- x$qe * -1
  qebytes <- as.vector(sapply(qe, .encode1))
  z <- c(logibytes, sizebytes, Abytes, Ebytes, qabytes, qebytes)
  return(z)
}
################################################################################
#' @rdname encoding
################################################################################
decodePHMM <- function(z){
  if(mode(z) != "raw") return(z)
  logibits <- rawToBits(z[1:4])
  AA <- as.logical(logibits[30])
  ID <- as.logical(logibits[31])
  DI <- as.logical(logibits[32])
  zsize <-  sum(as.integer(z[5:8]) * c(1, 256, 65536, 16777216))
  alength <- (((zsize + 1) * 4) - 3) * 2
  ## 1 extra col at front, 4 rows excl DI, ID etc, 3 -Infs (DD start, DD MD end), 2 bytes per entry
  astart <- 9
  aend <- alength + 9 - 1
  A <- z[astart:aend]
  A <- apply(matrix(A, nrow = 2), 2, .decode1) * -1
  A <- c(-Inf, A) ## append first DD transition
  lcis <- seq(length(A) - 1, length(A)) # last column indices
  lastcol <- c(-Inf, -Inf, A[lcis]) # final DD, MD, MM and IM transitions
  A <- A[-lcis]
  A <- c(A, lastcol)
  A <- matrix(A, nrow = 4) # just DD, MD, MM and IM at the moment
  DMrow <- log(1 - exp(A[1, ]))
  DMrow[1] <- -Inf
  DMrow[length(DMrow)] <- 0
  DIrow <- rep(-Inf, zsize + 1)
  MIrow <- (1 - exp(A[2, ])) - exp(A[3, ])
  MIrow[MIrow < 1e-06] <- 1e-06
  MIrow <- log(MIrow)
  IDrow <- rep(-Inf, zsize + 1)
  IIrow <- 1 - exp(A[4, ])
  IIrow[IIrow < 1e-06] <- 1e-06
  IIrow <- log(IIrow)
  A <- rbind(A[1, ], DMrow, DIrow, A[2:3, ], MIrow, IDrow, A[4, ], IIrow)
  dimnames(A) <- list(type = c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II"),
                      module = paste(0:zsize))
  estart <- aend + 1
  eend <- length(z) - if(AA) 54 else 22 #(4 x 2 for qe, 7 x 2 for qa)
  E <- z[estart:eend]
  E <- apply(matrix(E, nrow = 2), 2, .decode1) * -1
  E <- matrix(E, nrow = if(AA) 19 else 3)
  Trow <- 1 - apply(exp(E), 2, sum)
  #Trow <- apply(exp(E), 2, function(v) 1 - v[1] - v[2] - v[3])
  Trow[Trow < 1e-06] <- 1e-06
  E <- rbind(E, log(Trow))
  colnames(E) <- paste(1:zsize)
  rownames(E) <- if(AA) LETTERS[-c(2, 10, 15, 21, 24, 26)] else c("A", "C", "G", "T")
  qastart <- eend + 1
  qaend <- qastart + 13 #(7 * 2 - 1)
  tmp <- z[qastart:qaend]
  tmp <- apply(matrix(tmp, nrow = 2), 2, .decode1) * -1
  qa <- structure(rep(-Inf, 9), names = c("DD", "DM", "DI", "MD", "MM", "MI", "ID", "IM", "II"))
  qa[c(1, 2, 4, 5, 6, 8, 9)] <- tmp
  qestart <- qaend + 1
  qeend <- length(z)
  qe <- z[qestart:qeend]
  qe <- apply(matrix(qe, nrow = 2), 2, .decode1) * -1
  names(qe) <- if(AA) LETTERS[-c(2, 10, 15, 21, 24, 26)] else c("A", "C", "G", "T")
  res <- list(size = zsize, A = A, E = E, qa = qa, qe = qe)
  class(res) <- "PHMM"
  return(res)
}
################################################################################
