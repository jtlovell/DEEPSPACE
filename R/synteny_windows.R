#' @title Subset to collinear windowed hits
#' @description
#' \code{synteny_windows} takes a data.table representation of a paf file,
#' where the hits are from evenly spaced 'windowed' alignments and uses
#' MCScanX and dbscan to subset to collinear hits and block coordinates. Most of
#' the functioning under-the-hood are modified directly from GENESPACE.
#' @name synteny_windows
#'
#' @param paf data.table with the 12 standard paf columns and a 13th column
#' title "windowID", containing the unique window identifier
#' @param width integer, the basepair size of the window in the query genome
#' @param step integer, the basepair size of the step between windows
#' @param topHits2keep maximum proportion of Ns for a window to be included
#' @param minMapq integer, passed to minimap2 -f. Setting to 0 can
#' significantly slow operations but will allow visualization of all repeats
#' @param minBlkSize integer, minimum block size
#' @param minSplitSize integer, minimum block size after splitting (maxGapInBlk)
#' @param initNgaps initial number of gaps, passed to MCScanX
#' @param finNgaps final number of gaps, passed to MCScanX
#' @param tmpDir file.path where the temporary files should be written
#' @param MCScanX_hCall file.path to the executable
#' @param maxGapInBlk integer, the maximum size of a gap within a block
#' @param verbose logical, should updates be printed to the console?
#'
#'    ‘’ ‘’ ‘’ ‘’
#'
#' @details coming soon
#'
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @import GenomicRanges
#' @export
synteny_windows <- function(paf,
                            tmpDir = getwd(),
                            MCScanX_hCall,
                            step = NULL,
                            width = NULL,
                            topHits2keep = 1,
                            minMapq = 12,
                            minBlkSize = 5,
                            minSplitSize = 2,
                            initNgaps = minBlkSize * 10,
                            finNgaps = minBlkSize,
                            maxGapInBlk = NULL,
                            verbose = TRUE){
  ##############################################################################
  ##############################################################################
  # -- ad hoc functions
  ##############################################################################
  convert_pos2ord <- function(start, end, group, to = 0){
    x <- (as.numeric(start) + as.numeric(end))/2
    to <- as.numeric(to)
    if(to > 0){
      x <- round_toInteger(x, to)
    }
    ord <- frankv(
      list(group, x),
      ties.method = "dense")
    return(ord)
  }

  project_anchorsOnPaf <- function(anchorPaf, paf, round2, buffer){
    hits <- data.table(anchorPaf)
    raw <- data.table(paf)

    hits <- hits[,colnames(raw), with = F]
    hits[,`:=`(
      qord = convert_pos2ord(
        start = qstart, end = qend, group = qname, to = round2),
      tord = convert_pos2ord(
        start = tstart, end = tend, group = tname, to = round2),
      isAnchor = TRUE)]

    qord <- tord <- qname <- tname <- index <- NULL
    setkey(hits, qord, tord)
    chrs <- with(hits, unique(paste(qname, tname)))

    tmp <- subset(raw, paste(qname, tname) %in% chrs)
    tmp <- subset(tmp, !index %in% hits$index)
    tmp[,isAnchor := FALSE]
    hits <- rbind(hits, tmp, fill = T)

    qname <- qstart <- qend <- qord <- NULL
    setkey(hits, qname, qstart, qend)
    hits[,qiord := interp_xy(x = (qstart + qend)/2, y = qord), by = "qname"]

    tname <- tstart <- tend <- tord <- NULL
    setkey(hits, tname, tstart, tend)
    hits[,tiord := interp_xy(x = (tstart + tend)/2, y = tord), by = "tname"]

    qiord <- tiord <- qrndi <- trndi <- isAnchor <- inBuffer <- NULL
    hits[,`:=`(qrndi = round_toInteger(qiord, buffer * 2),
               trndi = round_toInteger(tiord, buffer * 2))]
    rnds <- with(subset(hits, isAnchor), paste(qrndi, trndi))
    hits <- subset(hits, paste(qrndi, trndi) %in% rnds)

    hits[,inBuffer := pull_proxHits(
      x = qiord,
      y = tiord,
      isAnchor = isAnchor,
      radius = buffer),
      by = c("qname", "tname")]

    hits <- subset(hits, inBuffer)
    hits <- hits[,colnames(paf), with = F]
    return(hits)
  }
  ##############################################################################

  ##############################################################################
  pull_proxHits <- function(x, y, isAnchor, radius){

    # -- parameter checking
    if(length(x) != length(y))
      stop("x and y must be vectors of the same length\n")
    if(sum(isAnchor) == 0)
      return(rep(FALSE, length(x)))
    if(sum(!isAnchor) == 0)
      return(rep(TRUE, length(x)))
    if(length(x) < 2)
      stop("x and y must have lengths > 1\n")
    x <- as.numeric(x)
    if(any(is.na(x)))
      stop("x values must all be numeric or integer\n")
    y <- as.numeric(y)
    if(any(is.na(y)))
      stop("y values must all be numeric or integer\n")

    isAnchor <- as.logical(isAnchor)
    if(any(is.na(isAnchor)))
      stop("values given to isAnchor (logical vector of isAnchor observations
           are anchors) must all be coercible to logical\n")

    radius <- as.numeric(radius[1])
    if(is.na(radius) || is.null(radius))
      stop("radius must be a single numeric value > 0\n")
    radius <- as.numeric(radius[1])
    if(radius <= 0)
      stop("radius must be a single numeric value > 0\n")

    # -- function to do it
    # -- get fixed radius nearest neighbors
    if(all(isAnchor) || all(!isAnchor)){
      return(isAnchor)
    }else{
      xy <- data.frame(x = x, y = y, isAnchor = isAnchor, inBuffer = isAnchor)
      isAnchor <- NULL
      fn <- frNN(
        x = subset(xy, isAnchor)[c("x", "y")],
        query = subset(xy, !isAnchor)[c("x", "y")],
        eps = radius)

      # -- subset to positions with any nearest neighbors in radius
      hasAnch <- fn$id[sapply(fn$id, length) > 0]
      wh <- as.numeric(names(hasAnch))
      # -- set these observations to true and return
      if(length(wh) > 0)
        xy$inBuffer[wh] <- TRUE
      return(xy$inBuffer)
    }
  }
  ##############################################################################
  ##############################################################################
  # -- main function starts here

  # -- get windows width and steps from the paf if set to NULL
  width <- as.integer(width[1])
  step <- as.integer(step[1])

  if(is.na(width) || is.na(step) || length(width) < 1 || length(step) < 1){
    qstart <- qend <- st <- en <- NULL
    qcoord <- paf[,list(st = min(qstart), en = max(qend)),
                  by = c("qname", "windowID")]
    width <- with(qcoord, max(en - st)) + 1
    setkey(qcoord, qname, st, en)
    step <- min(qcoord[,list(di = min(diff(st))), by = "qname"]$di)
  }

  maxGapInBlk <- as.integer(maxGapInBlk[1])
  if(is.na(maxGapInBlk) || length(maxGapInBlk) < 1){
    maxGapInBlk <- (width + step) * 10
  }

  fuzzy <- (step * 2) + (width * 2)
  if(!"windowID" %in% names(paf))
    stop("paf object must have window IDs")

  # -- pull hit set following rules of top hits and mapq
  windowID <- nmatch <- mapq <- alen <- tname <- tstart <- tend <- wrank <- NULL
  setorder(paf, windowID, -nmatch, -mapq, -alen)
  paf[,wrank := 1:.N, by = windowID]

  setkey(paf, tname, tstart, tend)
  paf[,index := 1:.N]

  cull <- subset(paf, wrank <= topHits2keep)
  cull <- subset(cull, mapq >= minMapq)

  if(verbose)
    cat(sprintf("\tPulled %sk seed hits (topN = %s, min mapq = %s)\n",
                round(nrow(cull)/1e3, 1), topHits2keep, minMapq))

  # -- run MCScanX
  cull1 <- mcscanx_paf(
    paf = cull,
    nGaps = initNgaps,
    blkSize = minBlkSize,
    MCScanX_hCall = MCScanX_hCall,
    tmpDir = tmpDir)
  cull1 <- subset(cull1, !is.na(mcscanxBlkID))
  if(verbose)
    cat(sprintf("\tFound %sk initial anchor hits\n",
                round(nrow(cull1)/1e3, 1)))

  # -- use initial anchors as coordinate system
  mapq <- index <- NULL
  cull2a <- project_anchorsOnPaf(
    anchorPaf = cull1,
    paf = subset(paf, mapq >= minMapq),
    round2 = step + width,
    buffer = minBlkSize * 2)
  cull2b <- project_anchorsOnPaf(
    anchorPaf = cull1,
    paf = subset(paf, mapq < minMapq),
    round2 = step + width,
    buffer = 1)
  cull2b <- subset(cull2b, !index %in% cull2a$index)
  cull2 <- rbind(cull2a, cull2b)
  if(verbose)
    cat(sprintf("\tPulled %sk hits proximate to initial anchors\n",
                round(nrow(cull2)/1e3, 1)))

  # -- re-run mcscanx
  mcscanxBlkID <- NULL
  cull3 <- mcscanx_paf(
    paf = cull2,
    nGaps = initNgaps,
    blkSize = minBlkSize,
    MCScanX_hCall = MCScanX_hCall,
    tmpDir = tmpDir)
  cull3 <- subset(cull3, !is.na(mcscanxBlkID))

  # -- convert to hit order and cluster
  qstart <- qend <- qname <- qpos <- tname <- tstart <- tend <- tpos <- NULL
  cull3[,tpos := convert_pos2ord(
    start = tstart, end = tend, group = tname, to = fuzzy)]
  cull3[,qpos := convert_pos2ord(
    start = qstart, end = qend, group = qname, to = fuzzy)]

  regID <- NULL
  cull3[,regID := dbscan(frNN(
    x = cbind(tpos, qpos),
    eps = initNgaps),
    minPts = minBlkSize)$cluster,
    by = c("qname", "tname")]
  cull3 <- subset(cull3, regID != 0)
  cull3[,regID := paste(qname, tname, regID)]
  cull3[,mcscanxBlkID := NULL]

  qstart <- qend <- qname <- qpos <- tname <- tstart <- tend <- tpos <- NULL
  cull3[,tpos := convert_pos2ord(
    start = tpos, end = tpos, group = regID, to = 2)]
  cull3[,qpos := convert_pos2ord(
    start = qpos, end = qpos, group = regID, to = 2)]

  blkID <- NULL
  cull3[,blkID := dbscan(frNN(
    x = cbind(tpos, qpos),
    eps = finNgaps),
    minPts = minBlkSize)$cluster,
    by = c("regID")]
  cull3 <- subset(cull3, blkID != 0)
  cull3[,blkID := paste(blkID, regID)]
  if(verbose)
    cat(sprintf("\tRe-culled to %sk proximate hits ",
                round(nrow(cull3)/1e3, 1)))

  mcscanxBlkID <- NULL
  cull4 <- rbindlist(lapply(split(cull3, by = "blkID"), function(x)
    mcscanx_paf(
      paf = x,
      nGaps = finNgaps,
      blkSize = minBlkSize,
      MCScanX_hCall = MCScanX_hCall,
      tmpDir = tmpDir)))
  cull4 <- subset(cull4, !is.na(mcscanxBlkID))
  if(verbose)
    cat(sprintf("and %sk anchors\n",
                round(nrow(cull4)/1e3, 1)))

  qstart <- qend <- qname <- qpos <- tname <- tstart <- tend <- tpos <- NULL
  cull4[,tpos := convert_pos2ord(
    start = tpos, end = tpos, group = blkID, to = 1)]
  cull4[,qpos := convert_pos2ord(
    start = qpos, end = qpos, group = blkID, to = 1)]

  tmp <- NULL
  cull4[,tmp := dbscan(frNN(
    x = cbind(tpos, qpos),
    eps = sqrt((finNgaps+2))),
    minPts = minBlkSize)$cluster,
    by = c("blkID")]
  cull5 <- subset(cull4, tmp != 0)
  cull5[,blkID := paste(tmp, regID)]

  nInGpd <- gpd <- blkID <- NULL
  cull5[,gpd := dbscan(frNN(
    x = cbind((tstart + tend)/2, (qstart + qend) / 2),
    eps = maxGapInBlk*2),
    minPts = 1)$cluster,
    by = c("blkID")]
  cull5[,nInGpd := .N, by = c("blkID", "gpd")]

  qstart <- qend <- qname <- tname <- tstart <- tend <- blkID <- NULL
  out <- subset(cull5, nInGpd >= minSplitSize)
  out[,blkID := paste(blkID, gpd)]
  setkey(out, tname, tstart, tend, qname, qstart, qend)
  out[,blkID := sprintf(
    "blk_%s", as.integer(factor(blkID, levels = unique(blkID))))]
  bc <- calc_pafBlkCoords(out)
  if(verbose)
    cat(sprintf("\tFinal synteny set: %s hits in %s blocks\n",
                nrow(out), nrow(bc)))
  return(list(hits = out[,names(bc), with = F], blkCoords = bc))
}
