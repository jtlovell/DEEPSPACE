#' @title prep_fastas
#' @description
#' \code{prep_fastas} prep_fastas
#' @name prep_clwFastas
#'
#' @param queryFasta character string coercible to a file path specifying the
#' location of the query (to-be-windowed) fasta file
#' @param targetFasta character string coercible to a file path specifying the
#' location of the target fasta file
#' @param queryID character string to name the query (windowed) genome
#' @param targetID character string to name the target genome
#' @param windowSize integer, the basepair size of the window in the query
#' genome
#' @param stepSize integer, the basepair distance between adjacet window start
#' sites. Defaults to NULL and is calculated from the genome size / nWindows
#' @param nWindows integer, the number of windows to chop the query genome.
#' @param outDir character string coercible to a file path specifying the
#' output directory
#' @param stripChrname character string coercible to a regular expression passed
#' to gsub that re-names the sequences in both the query and target fastas
#' @param minChrLen integer, minimum length of a sequence to be considered
#' @param overwrite logical, should existing files be overwritten?
#'
#' @details coming soon
#'
#' @import data.table
#' @import Biostrings
#' @import GenomicRanges
#' @export
prep_fastas <- function(queryFasta,
                        targetFasta,
                        outQueryFasta,
                        outTargetFasta,
                        windowSize,
                        stepSize,
                        nWindows,
                        stripChrname,
                        minChrLen){

  ##############################################################################
  # 1. check that both fasta files exist and are DNA sequences
  queryFasta <- check_isDNAFasta(queryFasta)
  targetFasta <- check_isDNAFasta(targetFasta)


  ##############################################################################
  # 3. check parameters
  # -- regex to strip text from chromosome names
  stripChrname <- as.character(stripChrname[1])
  if(is.null(stripChrname) || is.na(stripChrname))
    stop("stripChrname must be a character string coercible to a regex\n")

  # -- window size
  windowSize <- as.integer(windowSize)
  if(!is.integer(windowSize))
    stop("windowSize could not be coerced to an integer\n")
  if(windowSize < 50)
    stop("small windows is not what DEEPSPACE is designed for. Try gepard\n")

  # -- min chr size
  minChrLen <- as.integer(minChrLen)
  if(!is.integer(minChrLen))
    stop("minChrLen could not be coerced to an integer\n")
  if(minChrLen < 0)
    minChrLen <- 0

  # -- overwrite
  overwrite <- is.logical(overwrite[1])
  if(!is.logical(overwrite) || is.na(overwrite))
    stop("overwrite must be TRUE/FALSE\n")

  ##############################################################################
  # 4. Prep fastas
  # -- read in the fastas
  if(!file.exists(outQueryFasta) || overwrite){
    qss <- readDNAStringSet(queryFasta)
    qss <- qss[width(qss) >= minChrLen]
    names(qss) <- gsub(stripChrname, "", names(qss))
    if(any(duplicated(names(qss))))
      stop("query chromosome names are not unique after applying stripChrname\n")
    if(any(names(qss) == ""))
      stop("one query chromosome name was entirely stripped away after applying stripChrname\n")


    # -- determine the step size
    if(is.null(stepSize)){
      nWindows <- as.integer(nWindows)
      if(!is.integer(nWindows))
        stop("nWindows could not be coerced to an integer and stepSize was not provided\n")
      nbp <- sum(width(qss))
      stepSize <- round(nbp / nWindows)
      if(stepSize < 1)
        stepSize <- 1
    }

    # -- check step size
    stepSize <- as.integer(stepSize)
    if(!is.integer(stepSize))
      stop("stepSize could not be coerced to an integer\n")

    # -- window the fasta file
    qssw <- window_ss(
      ss = qss,
      windowSize = windowSize,
      windowStep = stepSize)

    # -- write the windows query genome
    writeXStringSet(qssw, filepath = outQueryFasta)
  }else{
    fl <- fasta.seqlengths(outQueryFasta)
    nWindows <- length(fl)
    tmp <- readLines(outQueryFasta, 1)
    tmp <- strsplit(tmp, "XstartX")[[1]][2]
    stepSize <- diff(as.numeric(strsplit(tmp, "_XendX")[[1]]))
    windowSize <- fl[1]
  }

  # -- subset and copy over target fasta if it doesn't exist
  if(!file.exists(outTargetFasta) || overwrite){
    tss <- readDNAStringSet(targetFasta)
    tss <- tss[width(tss) >= minChrLen]
    names(tss) <- gsub(stripChrname, "", names(tss))
    if(any(duplicated(names(tss))))
      stop("target chromosome names are not unique after applying stripChrname\n")
    if(any(names(tss) == ""))
      stop("one target chromosome name was entirely stripped away after applying stripChrname\n")
    writeXStringSet(tss, filepath = outTargetFasta)
  }

  pafFile <- file.path(outDir, sprintf(
    "%s_windows__vs__%s.paf",
    queryID, targetID))

  ##############################################################################
  # 5. Return parameters
  outParam <- data.table(
    queryFasta = outQueryFasta,
    targetFasta = outTargetFasta,
    minChrLen = minChrLen,
    windowSize = windowSize,
    stepSize = stepSize,
    nWindows = nWindows,
    pafFile = pafFile)
  return(outParam)
}

