#' @title windowed minimap2 alignments
#' @description
#' \code{mm2_windows} takes two DNAStringSets and aligns the query (windowed)
#' against the full reference genome. Returns a data.table representation of
#' the paf file.
#' @name mm2_windows
#' @param query DNAStringSet containing the fasta to be windowed
#' @param target DNAStringSet containing the fasta to serve as the target
#' @param wd character string coercible to a file.path where the temporary
#' files should be written. Taken as the current working directory by default.
#' @param speedPreset character string giving the speed preset to use. If
#' 'custom', ignore.
#' @param xPreset character string giving the divergence preset to use. If
#' 'custom', ignore.
#' @param windowSize integer, the basepair size of the window in the query genome
#' @param stepSize integer, the basepair size of the step between windows
#' @param maxPropNs maximum proportion of Ns for a window to be included
#' @param filterTopReps integer, passed to minimap2 -f. Setting to 0 can
#' significantly slow operations but will allow visualization of all repeats
#' @param maxTopHits integer, passed to minimap2 -N
#' @param minSecRat numeric, passed to minimap2 -p
#' @param mismatchPenalty numeric, passed to minimap2 -B
#' @param gapOpen integer, passed to minimap2 -O
#' @param gapExt integer, passed to minimap2 -e
#' @param nWindows integer, number of windows in the analysis (overrides step)
#' @param batchSize integer, passed to minimap2 -I, currently ignored
#' @param mm2verbose logical, should the full mm2 output be printed?
#' @param kmerLength integer, passed to minimap2 -k
#' @param nCores integer, passed to minimap2 -t
#' @param overwrite logical, should existing files be overwritten?
#' @param basicParams character string of additional parameters passed to
#' minimap2.
#' @param minimap2call character coercible to a file path pointing to the
#' minimap2 exectutable. If minimap2 is the path, the default 'minimap2' will
#' suffice
#' @param overwrite logical, if tmp files exist, should they be overwritten?
#' @param verbose logical, should updates be printed to the console?
#'
#' @details coming soon
#'
#' @import data.table
#' @import Biostrings
#' @import GenomicRanges
#' @export
mm2_windows <- function(query,
                        target,
                        speedPreset = c("custom", "all", "default", "fast", "draft"),
                        xPreset = c("custom", "5", "10", "20"),
                        wd = getwd(),
                        windowSize = NULL,
                        stepSize = NULL,
                        nWindows = NULL,
                        filterTopReps = 0.0001,
                        maxTopHits = 100,
                        minSecRat = 0.5,
                        mismatchPenalty = 4,
                        maxPropNs = .1,
                        gapOpen = "6,26",
                        gapExt = "2,1",
                        batchSize = "4G", # not currently used
                        kmerLength = 19,
                        nCores = 1,
                        basicParams = "-A1 -U50,500 --no-long-join --rmq=no -n1 -m0 --frag=yes",
                        minimap2call = "minimap2",
                        overwrite = FALSE,
                        mm2verbose = FALSE,
                        verbose = TRUE){

  ##############################################################################
  ##############################################################################
  # ad hoc functions
  ##############################################################################
  ##############################################################################
  calc_kmerLengthPreset <- function(x, xPreset){
    if(xPreset == "5")
      kmerLength <- kmerLength + 2
    if(xPreset == "20")
      kmerLength <- kmerLength - 2
    return(kmerLength)
  }

  calc_widthPreset <- function(x, xPreset){
    if(xPreset %in% c("custom", "10"))
      return(x)
    if(xPreset == "5")
      windowSize <- ceiling(x / 2)
    if(xPreset == "10")
      windowSize <- 500
    if(xPreset == "20")
      windowSize <- floor(x * 2)
    return(windowSize)
  }
  check_basicParams <- function(x){
    x <- as.character(x[1])
    if(!is.character(x))
      stop("basicParams must be coercible to a character vector\n")
    if(grepl(" -B", x))
      stop("cannot specify -B in basicParams, use mismatchPenalty\n")
    if(grepl(" -r", x))
      stop("cannot specify -r in basicParams, use gapChain\n")
    if(grepl(" -O", x))
      stop("cannot specify -O in basicParams, use gapOpen\n")
    if(grepl(" -E", x))
      stop("cannot specify -E in basicParams, use gapExt\n")
    if(grepl(" -f", x))
      stop("cannot specify -f in basicParams, use filterTopReps\n")
    if(grepl(" -p", x))
      stop("cannot specify -p in basicParams, use minSecRat\n")
    if(grepl(" -N", x))
      stop("cannot specify -N in basicParams, use maxTopHits\n")
    if(grepl(" -k", x))
      stop("cannot specify -k in basicParams, use kmerLength\n")
    if(grepl(" -w", x))
      stop("cannot specify -w in basicParams, use minWinSize\n")
    if(grepl(" -t", x))
      stop("cannot specify -t in basicParams, use nCores\n")
    return(x)
  }

  ##############################################################################
  make_userParams <- function(gapChain,
                              gapOpen,
                              gapExt,
                              mismatchPenalty,
                              filterTopReps,
                              minSecRat,
                              maxTopHits,
                              kmerLength,
                              minWinSize,
                              chainStop,
                              insertLength,
                              nCores){
    gapChain <- as.character(gapChain[1])
    if(!is.character(gapChain))
      stop("gapChain must be coercible to a character\n")
    gapOpen <- as.character(gapOpen[1])
    if(!is.character(gapOpen))
      stop("gapOpen must be coercible to an character\n")
    gapExt <- as.character(gapExt[1])
    if(!is.character(gapExt))
      stop("gapExt must be coercible to an character\n")
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
    kmerLength <- as.integer(kmerLength[1])
    if(!is.integer(kmerLength))
      stop("kmerLength must be coercible to an integer\n")
    minWinSize <- as.integer(minWinSize[1])
    if(!is.integer(minWinSize))
      stop("minWinSize must be coercible to an integer\n")
    chainStop <- as.integer(chainStop[1])
    if(!is.integer(chainStop))
      stop("chainStop must be coercible to an integer\n")
    insertLength <- as.integer(insertLength[1])
    if(!is.integer(insertLength))
      stop("insertLength must be coercible to an integer\n")
    nCores <- as.integer(nCores[1])
    if(!is.integer(nCores))
      stop("nCores must be coercible to an integer\n")
    userParams <- sprintf(
      "-r%s -O%s -e%s -B%s -f%s -p%s -N%s -k%s -w%s -g%s -F%s -t%s",
      gapChain, gapOpen, gapExt,
      format(mismatchPenalty, scientific = FALSE),
      format(as.numeric(filterTopReps), scientific = FALSE),
      format(as.numeric(minSecRat), scientific = FALSE),
      format(as.integer(maxTopHits), scientific = FALSE),
      format(as.integer(kmerLength), scientific = FALSE),
      format(as.integer(minWinSize), scientific = FALSE),
      format(as.integer(chainStop), scientific = FALSE),
      format(as.integer(insertLength), scientific = FALSE),
      format(as.integer(nCores), scientific = FALSE))
    return(userParams)
  }
  ##############################################################################
  ##############################################################################
  # end Ad Hoc functions
  ##############################################################################
  ##############################################################################

  ##############################################################################
  # 1. Get the presets in order
  ##############################################################################
  # -- check the args
  speedPreset <- match.arg(
    as.character(speedPreset),
    choices = c("custom", "all", "default", "fast", "draft"))
  xPreset <- match.arg(
    as.character(xPreset),
    choices = c("custom", "5", "10", "20"))

  # -- width
  windowSize <- as.integer(windowSize[1])
  if(length(windowSize) == 0)
    windowSize <- calc_widthPreset(x = 500, xPreset = xPreset)

  # -- speed preset defaults
  if(speedPreset == "draft"){
    maxTopHits <- 2; filterTopReps <- .01;minSecRat <- .95
    maxPropNs <- 0; kmerLength <- 25; nWindows <- 100e3
  }
  if(speedPreset == "fast"){
    maxTopHits <- 4; filterTopReps <- .001; minSecRat <- .75
    maxPropNs <- .1; kmerLength <- 23; nWindows <- 200e3

  }
  if(speedPreset == "default"){
    maxTopHits <- 100; filterTopReps <- .0001; minSecRat <- .5
    maxPropNs <- .1; kmerLength <- 21; nWindows <- 500e3
  }
  if(speedPreset == "all"){
    maxTopHits <- 10000; filterTopReps <- 0; minSecRat <- 0.1
    maxPropNs <- .1; kmerLength <- 19; nWindows <- 1e6
  }

  # -- if custom, ensure all parameters look ok, if not, adjust kmerLength by
  # -- evoultionary distance
  if(speedPreset == "custom"){
    maxTopHits <- as.integer(maxTopHits[1])
    filterTopReps <- as.numeric(filterTopReps[1])
    minSecRat <- as.numeric(minSecRat[1])
    maxPropNs <- as.numeric(maxPropNs[1])
    kmerLength <- as.integer(kmerLength[1])
    windowSize <- as.integer(windowSize[1])
    if(length(c(kmerLength, maxPropNs, minSecRat,
                filterTopReps, maxTopHits, windowSize)) != 6)
      stop("something is wrong with the mm2 parameters ...
           they should be single value numeric or integer\n")
  }else{
    kmerLength <- calc_kmerLengthPreset(kmerLength, xPreset)
    if(length(windowSize) == 0)
      windowSize <- calc_widthPreset(x = 500, xPreset = xPreset)
  }

  # -- if nWindows is provided, re-calculate step size
  nWindows <- as.integer(nWindows)[1]
  nWindows <- nWindows[!is.na(nWindows)]
  if(length(nWindows) == 1){
    stepSize <- floor(sum(width(query)) / nWindows)
  }else{
    stepSize <- windowSize * 2
  }

  ##############################################################################
  # 2 -- get all parameters formatted for mm2
  ##############################################################################

  # -- check mm2 install
  check_minimap2Install(minimap2call)

  # -- adjust parameters based on those specified by user
  gapChain <- sprintf("%s,%s", windowSize, windowSize * 10)
  insertLength <- windowSize * 10
  chainStop <- windowSize * 2
  minWinSize <- floor((2 / 3) * kmerLength)

  # -- convert user defined parameters to character string of parameters
  userParams <- make_userParams(
    gapChain = gapChain,
    gapOpen = gapOpen,
    gapExt = gapExt,
    mismatchPenalty = mismatchPenalty,
    filterTopReps = filterTopReps,
    minSecRat = minSecRat,
    maxTopHits = maxTopHits,
    kmerLength = kmerLength,
    minWinSize = minWinSize,
    chainStop = chainStop,
    insertLength = insertLength,
    nCores = nCores)

  # -- check basic parameters, given as a character string
  basicParams <- check_basicParams(basicParams)

  # -- check input sequences
  if(!is(query, "DNAStringSet"))
    stop("query must be a DNAStringSet object\n")
  if(!is(target, "DNAStringSet"))
    stop("target must be a DNAStringSet object\n")
  if(verbose)
    cat(sprintf(
      "Preparing input files:\n\tTarget: %sMb in %s sequences\n",
      round(sum(width(target))/1e6, 1), length(target)))

  # -- check window size parameters
  windowSize <- as.integer(windowSize[1])
  if(!is.integer(windowSize))
    stop("windowSize must be coercible to an integer\n")
  stepSize <- as.integer(stepSize[1])
  if(!is.integer(stepSize))
    stop("stepSize must be coercible to an integer\n")

  # -- check other parameters
  verbose <- verbose[1]
  if(!is.logical(verbose))
    verbose <- TRUE

  ##############################################################################
  # 3. Sliding window on the query
  ##############################################################################

  # -- do the sliding window on the query
  queryGr <- window_it(
    x = query,
    windowSize = windowSize,
    stepSize = stepSize)
  queryGr$index <- 1:length(queryGr)
  querySS <- query[queryGr]
  names(querySS) <- queryGr$index
  if(verbose)
    cat(sprintf(
      "\tQuery: %sMb in %s %sbp windows (%sbp step)\n",
      round(sum(width(querySS))/1e6, 1), length(querySS), windowSize, stepSize))

  # -- check the number of Ns and drop any windows with enough
  ncnt <- vcountPattern("N", querySS) / width(querySS)
  fewEnufNs <- ncnt <= maxPropNs
  querySS <- querySS[fewEnufNs]

  ##############################################################################
  # 4. Final checks on input/output directory structure
  ##############################################################################

  # -- check that the working directory exists
  wd <- path.expand(wd[1])
  if(!dir.exists(wd) && dir.exists(dirname(wd))){
    dir.create(wd)
  }else{
    if(!dir.exists(wd))
      stop("working directory (wd) does not exist")
  }

  # -- check that the tmp files do not exist, or if they do if overwrite
  overwrite <- overwrite[1]
  if(!is.logical(overwrite))
    overwrite <- FALSE
  qf <- file.path(wd, "qtmp.fa")
  tf <- file.path(wd, "ttmp.fa")
  pf <- file.path(wd, "tmp.paf")
  if(file.exists(qf) && !overwrite)
    stop("query tmp file (wd/qtmp.fa) exists and !overwrite ...
         either set overwrite = TRUE or remove this file manually\n")
  if(file.exists(tf) && !overwrite)
    stop("target tmp file (wd/ttmp.fa) exists and !overwrite ...
         either set overwrite = TRUE or remove this file manually\n")
  if(file.exists(qf) && !overwrite)
    stop("tmp paf output file (wd/tmp.paf) exists and !overwrite ...
         either set overwrite = TRUE or remove this file manually\n")

  ##############################################################################
  # 5 Run minimap2
  ##############################################################################

  # -- ensure the files get cleaned up on exit
  on.exit(
    expr = c(unlink(qf), unlink(tf), unlink(pf)),
    add = T)

  # -- write the stringsets to file
  writeXStringSet(x = querySS, filepath = qf)
  writeXStringSet(x = target, filepath = tf)

  # -- build the minimap2 call
  mm2args <- c(userParams, basicParams, sprintf("-o %s", pf), tf, qf)

  # -- call minimap2
  if(verbose)
    cat(sprintf("running minimap2 as:\n\t %s %s %s\n",
                minimap2call, userParams, basicParams))
  ers <- system2(
    command = minimap2call,
    args = mm2args,
    stdout = TRUE,
    stderr = TRUE)

  if(mm2verbose)
    cat(ers, sep = "\n")

  ##############################################################################
  # 6. Read in the paf and return
  ##############################################################################

  # -- read in the paf file
  p <- fread_paf(pf)

  # -- if empty, warn
  if(nrow(p) == 0){
    warning("mm2 run complete, but output paf file is empty ...
            something is wrong\n")
  }else{
    if(verbose)
      cat("Done! Produced a paf file with %s entries", nrow(p))
    # -- get vectors of starts/ends etc. from granges
    qnv <- as.character(seqnames(queryGr))
    qlv <- width(query); names(qlv) <- names(query)
    qsv <- start(queryGr)
    qev <- end(queryGr)
    gri <- as.character(queryGr$index)

    # -- get matching coordinates between the window number and granges
    mtch <- match(p$qname, gri)

    # -- add data back to the paf
    qname <- qlen <- NULL
    p[,`:=`(windowID = qname,
            qname = qnv[mtch],
            qstart = qsv[mtch],
            qend = qev[mtch])]
    p[,qlen := qlv[qname]]
  }

  return(p)
}
