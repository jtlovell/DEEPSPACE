#' @title minimap2_windows
#' @description
#' \code{minimap2_windows} minimap2_windows
#' @name minimap2_windows
#'
#' @param queryFasta character string coercible to a file path specifying the
#' location of the query (to-be-windowed) fasta file
#' @param targetFasta character string coercible to a file path specifying the
#' location of the target fasta file
#' @param pafFile character string coercible to a file path specifying the
#' location of the output (paf) file
#' @param windowSize integer, the basepair size of the window in the query
#' genome
#' @param filterTopReps integer, passed to minimap2 -f. Setting to 0 can
#' significantly slow operations but will allow visualization of all repeats
#' @param maxTopHits integer, passed to minimap2 -N
#' @param minSecRat numeric, passed to minimap2 -p
#' @param mismatchPenalty numeric, passed to minimap2 -B
#' @param kmerSize integer, passed to minimap2 -k
#' @param kmerStep integer, passed to minimap2 -w
#' @param nCores integer, passed to minimap2 -t
#' @param overwrite logical, should existing files be overwritten?
#' @param basicParams character string of additional parameters passed to
#' minimap2.
#' @param verbose logical, should updates be printed to the console?
#'
#' @details coming soon
#'
#' @import data.table
#' @import Biostrings
#' @import GenomicRanges
#' @export
minimap2_windows <- function(queryFasta,
                             targetFasta,
                             pafFile,
                             windowSize,

                             filterTopReps = 0.0001,
                             maxTopHits = 50,
                             minSecRat = 0.9,
                             mismatchPenalty = 4,
                             kmerSize = 25,
                             kmerStep = 20,
                             nCores = 1,

                             overwrite = FALSE,
                             basicParams = "-A1 -U50,500 --no-long-join -r0,0 --rmq=no -E2,1 -n1 -m0 --frag=yes",
                             minimap2call = "minimap2",
                             verbose = TRUE){

  ##############################################################################
  ##############################################################################
  # ad hoc functions
  ##############################################################################
  ##############################################################################
  check_basicParams <- function(x){
    x <- as.character(x[1])
    if(!is.character(x))
      stop("basicParams must be coercible to a character vector\n")
    if(grepl(" -B", x))
      stop("cannot specify -B in basicParams, use mismatchPenalty\n")
    if(grepl(" -f", x))
      stop("cannot specify -f in basicParams, use filterTopReps\n")
    if(grepl(" -p", x))
      stop("cannot specify -p in basicParams, use minSecRat\n")
    if(grepl(" -N", x))
      stop("cannot specify -N in basicParams, use maxTopHits\n")
    if(grepl(" -k", x))
      stop("cannot specify -k in basicParams, use kmerSize\n")
    if(grepl(" -w", x))
      stop("cannot specify -w in basicParams, use kmerStep\n")
    if(grepl(" -t", x))
      stop("cannot specify -t in basicParams, use nCores\n")
    return(x)
  }

  ##############################################################################
  make_userParams <- function(mismatchPenalty,
                              filterTopReps,
                              minSecRat,
                              maxTopHits,
                              kmerSize,
                              kmerStep,
                              nCores){
    mismatchPenalty <- as.integer(mismatchPenalty[1])
    if(!is.integer(mismatchPenalty))
      stop("mismatchPenalty must be coercible to an integer\n")
    filterTopReps <- as.numeric(filterTopReps[1])
    if(!is.numeric(mismatchPenalty))
      stop("mismatchPenalty must be coercible to an numeric\n")
    minSecRat <- as.numeric(minSecRat[1])
    if(!is.numeric(minSecRat))
      stop("minSecRat must be coercible to an numeric\n")
    maxTopHits <- as.integer(maxTopHits[1])
    if(!is.integer(maxTopHits))
      stop("maxTopHits must be coercible to an integer\n")
    kmerSize <- as.integer(kmerSize[1])
    if(!is.integer(kmerSize))
      stop("kmerSize must be coercible to an integer\n")
    kmerStep <- as.integer(kmerStep[1])
    if(!is.integer(kmerStep))
      stop("kmerStep must be coercible to an integer\n")
    nCores <- as.integer(nCores[1])
    if(!is.integer(nCores))
      stop("nCores must be coercible to an integer\n")

    userParams <- sprintf(
      "-B%s -f%s -p%s -N%s -k%s -w%s -t%s",
      format(mismatchPenalty, scientific = FALSE),
      format(as.numeric(filterTopReps), scientific = FALSE),
      format(as.numeric(minSecRat), scientific = FALSE),
      format(as.integer(maxTopHits), scientific = FALSE),
      format(as.integer(kmerSize), scientific = FALSE),
      format(as.integer(kmerStep), scientific = FALSE),
      format(as.integer(nCores), scientific = FALSE))
    return(userParams)
  }
  ##############################################################################
  ##############################################################################

  ##############################################################################
  # 1. Check input parameters
  # -- user defined parameters
  userParams <- make_userParams(
    mismatchPenalty = mismatchPenalty,
    filterTopReps = filterTopReps,
    minSecRat = minSecRat,
    maxTopHits = maxTopHits,
    kmerSize = kmerSize,
    kmerStep = kmerStep,
    nCores = nCores)

  # -- basic parameters
  basicParams <- check_basicParams(basicParams)

  # -- parameters that are taken from the windows
  windowSize <- as.integer(windowSize[1])
  if(!is.integer(windowSize))
    stop("windowSize must be coercible to an integer\n")
  calcParams <- sprintf(
    "-g%s -F%s",
    format(as.integer(ceiling(windowSize/2)), scientific = FALSE),
    format(as.integer(windowSize * 2), scientific = FALSE))

  ##############################################################################
  # 2. call minimap2
  # -- check that output directory exists
  if(!dir.exists(dirname(pafFile)))
    stop("path to pafFile is not valid\n")

  queryFasta <- check_isDNAFasta(queryFasta)
  targetFasta <- check_isDNAFasta(targetFasta)

  if(overwrite | !file.exists(pafFile)){
    ers <- system2(
      command = minimap2call,
      args = c(userParams,
               calcParams,
               basicParams,
               sprintf("-o %s", gsub(".gz$", "",pafFile)),
               targetFasta,
               queryFasta),
      stdout = TRUE,
      stderr = TRUE)
    if(verbose)
      cat(ers, sep = "\n")

    system2("gzip", sprintf("-f %s", gsub(".gz$", "", pafFile)))
  }

  return(data.table(
    queryFasta = queryFasta,
    targetFasta = targetFasta,
    params = paste(userParams, calcParams, basicParams, collapse = " "),
    pafFile = pafFile))
}
