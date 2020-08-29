#' Demultiplex merged FASTQ
#'
#' This function is used to demultiplex FASTQ files
#'   containing sequence reads with index and primer sequences still attached.
#'   The function trims the tags and primers, and exports
#'   two FASTQ files for each forward-reverse index combination.
#'
#' @param R1 character string giving the path to the first FASTQ file
#' @param R2 character string giving the path to the second FASTQ file.
#'   Set to NULL for single-end reads.
#' @param tags named character vector specifying the unique forward:reverse
#'   tag combinations used in the run.
#'   Each tag combination should be entered as an upper case character string delimited by a colon
#'   (e.g. "ATCGACAC:ATGCACTG") and named according to the unique sample ID.
#' @param up,down upper case character strings giving the forward and reverse primer sequences.
#' @param destdir character string giving the path to the directory where
#'   the new FASTQ files should be written.
#' @param adapter1 the forward flowcell adapter sequence to check and trim (set NULL to ignore).
#'   For standard Illumina MiSeq forward adapter set to "AATGATACGGCGACCACCGAGATCTACAC" (paired end sequencing only).
#' @param adapter2 the reverse flowcell adapter sequence to check and trim (set NULL to ignore).
#'   For standard Illumina MiSeq reverse adapter set to "CAAGCAGAAGACGGCATACGAGAT" (single or paired end sequencing).
#' @return NULL (invisibly)
#' @author Shaun Wilkinson
################################################################################
demultiplex <- function(R1, R2 = NULL, tags, up, down, destdir = "demux", adapter1 = NULL, adapter2 = NULL){
  if(is.null(names(tags))) names(tags) <- paste0("sample_", seq_along(tags))
  if(!dir.exists(destdir)) dir.create(destdir)
  #for(i in names(tags)) file.create(paste0(destdir, "/", i, ".fastq"))
  message("Writing FASTQ files")
  ftags <- unique(sub(":.*", "", tags))
  rtags <- unique(sub(".*:", "", tags))
  ftaglengths <- nchar(ftags)
  rtaglengths <- nchar(rtags)
  ## progress bar
  totallines <- file.size(R1)/180 # approx
  nreadins <- ceiling(totallines/1E07)
  nbars <- max(round((80/nreadins)/length(ftags)), 1)
  progbar <- paste0(rep("=", nbars), collapse = "")
  skip <- 0
  breakme <- FALSE
  counter <- 0
  ncup <- nchar(up)
  ncdn <- nchar(down)
  up <- disambiguate(up)
  rcdown <- disambiguate(rc(down))
  down <- disambiguate(down)
  #fregexs <- paste0("^", ftags, up)
  fregexs <- paste0(ftags, up)
  if(is.null(R2)){ # single-end reads
    rregexs <- paste0(rcdown, rc(rtags))
    repeat{
      counter <- counter + 1
      tmp <- scan(R1, nlines = 1E07, skip = skip, what = "", sep = "\n", quiet = TRUE)
      nlines <- length(tmp)
      if(nlines == 0) break
      if(nlines < 1E07) breakme <- TRUE
      skip <- skip + 1E07
      seqs <- tmp[seq(2, nlines, by = 4)]
      names(seqs) <- tmp[seq(1, nlines, by = 4)]
      quals <- tmp[seq(4, nlines, by = 4)]
      if(!is.null(adapter2)){
        oks <- grepl(rc(adapter2), seqs)
        seqs <- seqs[oks]
        quals <- quals[oks]
      }
      oks <- grepl(paste0(".+", up, ".+", rcdown, ".+"), seqs)
      seqs <- seqs[oks]
      quals <- quals[oks]
      for(i in seq_along(ftags)){
        cat(progbar)
        hitsi <- grepl(fregexs[i], seqs) #logical
        if(sum(hitsi) > 0){
          seqstart <- ftaglengths[i] + ncup + 1
          seqsi <- substr(seqs[hitsi], seqstart, 10000)
          qualsi <- substr(quals[hitsi], seqstart, 10000)
          for(j in seq_along(rtags)){
            whichsample <- match(paste0(ftags[i], ":", rtags[j]), tags)
            if(!is.na(whichsample) & length(seqsi) > 0){
              gregs <- gregexpr(rregexs[j], seqsi)
              gregs <- vapply(gregs, head, 0L, 1L) #note cant use unlist due to var list length
              hitsij <- gregs > 0
              nhits <- sum(hitsij)
              if(nhits > 0){
                if(!all(gregs == 1L)){ ## will all be 1 if no r tags used
                  seqsij <- substr(seqsi[hitsij], rep(1, nhits), gregs[hitsij] - 1)
                  qualsij <- substr(qualsi[hitsij], rep(1, nhits), gregs[hitsij] - 1)
                }else{
                  seqsij <- seqsi[hitsij]
                  qualsij <- qualsi[hitsij]
                }
                reslen <- length(seqsij) * 4
                res <- character(reslen)
                res[seq(1, reslen, by = 4)] <- names(seqsij)
                res[seq(2, reslen, by = 4)] <- seqsij
                res[seq(3, reslen, by = 4)] <- rep("+", length(seqsij))
                res[seq(4, reslen, by = 4)] <- qualsij
                fpath <- paste0(destdir, "/", names(tags)[whichsample], ".fastq")
                cat(res, file = fpath, sep = "\n", append = TRUE)
                seqsi <- seqsi[!hitsij] # makes it a bit faster
                qualsi <- qualsi[!hitsij]
              }
            }
          }
          seqs <- seqs[!hitsi] # makes it a bit faster
          quals <- quals[!hitsi]
        }
      }
      if(breakme) break
    }
  }else{# paired-end reads
    rregexs <- paste0("^", rtags, down)
    repeat{
      counter <- counter + 1
      tmp1 <- scan(R1, nlines = 1E07, skip = skip, what = "", sep = "\n", quiet = TRUE)
      tmp2 <- scan(R2, nlines = 1E07, skip = skip, what = "", sep = "\n", quiet = TRUE)
      nlines <- length(tmp1)
      if(nlines == 0) break # unlikely but could happen
      if(nlines < 1E07) breakme <- TRUE
      skip <- skip + 1E07
      seqs1 <- tmp1[seq(2, nlines, by = 4)]
      names(seqs1) <- tmp1[seq(1, nlines, by = 4)]
      quals1 <- tmp1[seq(4, nlines, by = 4)]
      seqs2 <- tmp2[seq(2, nlines, by = 4)]
      names(seqs2) <- tmp2[seq(1, nlines, by = 4)]
      quals2 <- tmp2[seq(4, nlines, by = 4)]

      if(!is.null(adapter2)){
        oks <- grepl(rc(adapter2), seqs1)
        seqs1 <- seqs1[oks]
        seqs2 <- seqs2[oks]
        quals1 <- quals1[oks]
        quals2 <- quals2[oks]
      }

      if(!is.null(adapter1)){
        oks <- grepl(rc(adapter1), seqs2)
        seqs1 <- seqs1[oks]
        seqs2 <- seqs2[oks]
        quals1 <- quals1[oks]
        quals2 <- quals2[oks]
      }

      oks <- grepl(paste0(".+", up, ".+"), seqs1) & grepl(paste0(".+", down, ".+"), seqs2)
      seqs1 <- seqs1[oks]
      seqs2 <- seqs2[oks]
      quals1 <- quals1[oks]
      quals2 <- quals2[oks]
      for(i in seq_along(ftags)){
        cat(progbar)
        hitsi <- grepl(fregexs[i], seqs1)
        if(sum(hitsi) > 0){
          seqstart1i <- ftaglengths[i] + ncup + 1
          seqs1i <- substr(seqs1[hitsi], seqstart1i, 10000)
          quals1i <- substr(quals1[hitsi], seqstart1i, 10000)
          seqs2i <- seqs2[hitsi]
          quals2i <- quals2[hitsi]
          for(j in seq_along(rtags)){
            seqstart2j <- rtaglengths[j] + ncdn + 1
            whichsample <- match(paste0(ftags[i], ":", rtags[j]), tags)
            if(!is.na(whichsample) & length(seqs1i) > 0){
              hitsij <- grepl(rregexs[j], seqs2i)
              nhits <- sum(hitsij)
              if(nhits > 0){
                seqs1ij <- seqs1i[hitsij]
                quals1ij <- quals1i[hitsij]
                seqs2ij <- substr(seqs2i[hitsij], seqstart2j, 10000)
                quals2ij <- substr(quals2i[hitsij], seqstart2j, 10000)
                reslen <- length(seqs1ij) * 4
                res1 <- res2 <- character(reslen)
                res1[seq(1, reslen, by = 4)] <- names(seqs1ij)
                res1[seq(2, reslen, by = 4)] <- seqs1ij
                res1[seq(3, reslen, by = 4)] <- rep("+", length(seqs1ij))
                res1[seq(4, reslen, by = 4)] <- quals1ij
                res2[seq(1, reslen, by = 4)] <- names(seqs2ij)
                res2[seq(2, reslen, by = 4)] <- seqs2ij
                res2[seq(3, reslen, by = 4)] <- rep("+", length(seqs2ij))
                res2[seq(4, reslen, by = 4)] <- quals2ij
                path1 <- paste0(destdir, "/", names(tags)[whichsample], "_R1_001.fastq")
                path2 <- paste0(destdir, "/", names(tags)[whichsample], "_R2_001.fastq")
                cat(res1, file = path1, sep = "\n", append = TRUE)
                cat(res2, file = path2, sep = "\n", append = TRUE)
                seqs1i <- seqs1i[!hitsij] # makes it a bit faster
                seqs2i <- seqs2i[!hitsij]
                quals1i <- quals1i[!hitsij]
                quals2i <- quals2i[!hitsij]
              }
            }
          }
          seqs1 <- seqs1[!hitsi]
          seqs2 <- seqs2[!hitsi]
          quals1 <- quals1[!hitsi]
          quals2 <- quals2[!hitsi]
        }
      }
      if(breakme) break
    }
  }
  message("\nDone")
  invisible(NULL)
}
################################################################################
