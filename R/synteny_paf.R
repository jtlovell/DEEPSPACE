#' @title synteny inference from windowed alignments
#' @description
#' \code{synteny_paf} methods to constrain windowed hits stored in a paf to
#' syntenic regions
#' @name synteny_paf
#'
#' @param paf character string coercible to a file path specifying the
#' location of the query (to-be-windowed) fasta file
#' @param nGaps character string coercible to a file path specifying the
#' location of the target fasta file
#' @param blkSize character string to name the query (windowed) genome
#' @param MCScanX_hCall character string coercible to a file path pointing to
#' the installed executable for MCScanX_h
#' @param tmpDir character string coercible to a file path where the MCScanX
#' input and output files should be stored.
#' @param cleanTmpDir logical, should the tmpDir be cleaned out after running?
#' Always set to TRUE unless manually troubleshooting
#' @param verbose logical, should updates be printed to the console?
#' @param fafiles named vector of file.paths to assembly fasta files
#' @param outdir file.path to the output directory
#' @param minChrSize numeric, minimum chromosome size to be included
#' @param stripChrname regex passed to gsub to strip text of chromosome names
#' @param windowSize integer, sliding window width
#' @param stepSize integer, step between sliding window start positions\
#' @param mapq20nmatchQuantile whats the quantile of nmatches for map20 hits
#' @param mapq60nmatchQuantile whats the quantile of nmatches for map60 hits
#' @param round2 what value should positions be rounded to?
#' @param targetNhitsQuantile quantile to calculate the threshold
#' @param maxBestMaps integer, the maximum number of best-mapping hits to keep
#' @param topNhits2keep integer, the number of top hits to keep
#' @param anchorPaf paf containing the anchor hits
#' @param buffer distance in hit rank order to include in search
#' MCScanX_hCall character string coercible to a valid file path pointing
#' to the MCScanX_h executable
#' @param minimap2call character string coercible to a valid file path pointing
#' to the minimap2 program. If installed with conda or in the path, use
#' "minimap2".
#' @param makePlots integer, specifying how many dotplots should be made. O
#' means no plots. 3 means make one for each step of synteny curation
#' @param blkSize the minimum number of windowed hits to be called a block.
#' @param tmpDir character string coercible to a file.path where temporary
#' files should be written. Cleaned out after run finishes unless this
#' directory already existed
#' @param ppi integer, pixels per inch in the dotplot
#' @param minPlotDim numeric, narrowest dimension of the dotplot in inches
#' @param maxBestMaps integer, the maximum number of best-mapping hits to keep
#' @param topNhits2keep integer, the number of top hits to keep
#' @param buffer numeric, position-rank order buffer around anchor hits to
#' search for additional anchors.
#' @param maxDistBtwHitsInBlk numeric, the maximum distance between two
#' adjacent hits in the same block
#' @param keepEverything logical, should final synteny searches be conducted on
#' all hits (including masked repeats etc)?
#' @param targetNhitsQuantile quantile to choose repetitive regions in the
#' target
#' @param verbose logical, should updates be printed to the console?
#' \cr
#' If called, \code{synteny_paf} returns its own arguments.
#'
#'

#' @title choose_topNmatch
#' @description
#' \code{choose_topNmatch} choose_topNmatch
#' @rdname synteny_paf
#' @import data.table
#' @importFrom stats quantile
#' @export
choose_topNmatch <- function(paf,
                             mapq20nmatchQuantile = 0.1,
                             mapq60nmatchQuantile = 0){
  hits <- data.table(paf)
  mapq <- nmatch <- NULL
  q20t <- with(subset(hits, mapq >= 20), quantile(nmatch, mapq20nmatchQuantile))
  q60t <- with(subset(hits, mapq >= 60), quantile(nmatch, mapq60nmatchQuantile))
  hits <- subset(hits, nmatch >= min(c(q20t, q60t)))
  return(hits)
}

#' @title drop_repetitiveTargets
#' @description
#' \code{drop_repetitiveTargets} drop_repetitiveTargets
#' @rdname synteny_paf
#' @import data.table
#' @export
drop_repetitiveTargets <- function(paf,
                                   round2,
                                   targetNhitsQuantile){
  hits <- data.table(paf)
  hits[,twind := round_toInteger((tstart + tend) / 2, round2)]
  hits[,n := .N, by = c("tname", "twind")]
  nhitThresh <- quantile(hits$n[!duplicated(hits$n)], targetNhitsQuantile)
  hits <- subset(hits, n <= nhitThresh)
  return(hits)
}

#' @title choose_localUniqHits
#' @description
#' \code{choose_localUniqHits} choose_localUniqHits
#' @rdname synteny_paf
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
choose_localUniqHits <- function(paf, round2){
  hits <- data.table(paf)
  hits[,`:=`(
    qord = convert_pos2ord(
      start = qstart, end = qend, group = qname, to = round2),
    tord = convert_pos2ord(
      start = tstart, end = tend, group = tname, to = round2))]
  setorder(hits, -mapq, -nmatch)
  hits <- subset(hits, !duplicated(paste(qord, tord)))
  hits <- hits[,colnames(paf), with = F]
  return(hits)
}


#' @title choose_topWindowHits
#' @description
#' \code{choose_topWindowHits} choose_topWindowHits
#' @rdname synteny_paf
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
choose_topWindowHits <- function(paf,
                                     maxBestMaps = 5,
                                     topNhits2keep = 2){
  hits <- data.table(paf)
  if(!"windID" %in% names(hits)){
    warning("paf data.table doesn't have window IDs ... there may be a problem\n")
    paf[,windID := 1:.N]
  }

  # 1. get per-window rankings
  hits[,`:=`(tmp1 = -mapq, tmp2 = -nmatch)]
  hits[,maprank := frankv(
    hits, cols = c("windID", "tmp1", "tmp2"),
    ties.method = "dense")]
  hits[,minrank := min(maprank), by = "windID"]
  hits[,maprank := (maprank + 1) - minrank]
  hits[,`:=`(minrank = NULL, tmp1 = NULL, tmp2 = NULL)]
  setkey(hits, qname, windID, maprank)

  # 2. calculate the number of best hits per window
  hits[,ntop := sum(maprank == 1), by = "windID"]
  hits <- subset(hits, ntop <= maxBestMaps)
  hits <- subset(hits, maprank <= topNhits2keep)

  # 3. Calculate the nmatches for the top hit for each window
  hits[,topMatch := max(nmatch[maprank == 1]), by = "windID"]
  # hits[,`:=`(ntop = NULL, maprank = NULL)]
  return(hits)
}

#' @title project_anchorsOnPaf
#' @description
#' \code{project_anchorsOnPaf} project_anchorsOnPaf
#' @rdname synteny_paf
#' @import data.table
#' @export
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
  setkey(hits, qord, tord)
  chrs <- with(hits, unique(paste(qname, tname)))

  tmp <- subset(raw, paste(qname, tname) %in% chrs)
  tmp <- subset(tmp, !index %in% hits$index)
  tmp[,isAnchor := FALSE]

  hits <- rbind(hits, tmp, fill = T)

  setkey(hits, qname, qstart, qend)
  hits[,qiord := interp_xy(x = (qstart + qend)/2, y = qord), by = "qname"]
  setkey(hits, tname, tstart, tend)
  hits[,tiord := interp_xy(x = (tstart + tend)/2, y = tord), by = "tname"]

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

#' @title choose_localBestHits
#' @description
#' \code{choose_localBestHits} choose_localBestHits
#' @rdname synteny_paf
#' @import data.table
#' @export
choose_localBestHits <- function(paf, round2){
  hits <- data.table(paf)

  hits[,`:=`(
    tmp1 = round_toInteger((qstart + qend) / 2, round2),
    tmp2 = round_toInteger((tstart + tend) / 2, round2))]
  hits[,`:=`(
    qid = paste(qname, tmp1),
    tid = paste(tname, tmp2))]
  hits[,`:=`(
    qid = as.numeric(as.factor(qid)),
    tid = as.numeric(as.factor(tid)),
    tmp1 = NULL, tmp2 = NULL)]
  hits[,`:=`(tmp1 = -mapq, tmp2 = -nmatch)]
  hits[,maprank := frankv(
    hits, cols = c("qid", "tmp1", "tmp2"),
    ties.method = "dense")]
  hits[,minrank := min(maprank), by = "qid"]
  hits[,maprank := (maprank + 1) - minrank]
  hits <- subset(hits, maprank == 1)
  hits[,`:=`(minrank = NULL, tmp1 = NULL, tmp2 = NULL)]

  return(hits)
}

#' @title synteny_paf
#' @description
#' \code{synteny_paf} synteny_paf
#' @rdname synteny_paf
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
synteny_paf <- function(paf,
                            makePlots = 3,
                            pdfFile = file.path(getwd(), "dotplots.pdf"),
                            ppi = 96,
                            minPlotDim = 12,
                            maxBestMaps = 5,
                            topNhits2keep = 2,
                            stepSize,
                            windowSize,
                            tmpDir,
                            MCScanX_hCall,
                            blkSize,
                            buffer = 5,
                            keepEverything = TRUE,
                            targetNhitsQuantile = 0.95,
                            mapq20nmatchQuantile = .1,
                            mapq60nmatchQuantile = 0,
                            maxDistBtwHitsInBlk = (blkSize * stepSize) * 2,
                            verbose = TRUE){
  ##############################################################################
  # 1. Get parameters set up
  # -- plotting parameters
  paf[,index := 1:.N]
  qtot <- sum(paf$qlen[!duplicated(paf$qname)])
  ttot <- sum(paf$tlen[!duplicated(paf$tname)])
  qtrat <- qtot / ttot
  if(qtot > ttot){
    pltHt <- minPlotDim
    pltWt <- minPlotDim * qtrat
  }else{
    pltWt <- minPlotDim
    pltHt <- minPlotDim / qtrat
  }
  qbin <- qtot / (ppi * pltWt)
  tbin <- ttot / (ppi * pltWt)

  if(makePlots > 0){
    pdf(pdfFile, height = pltHt, width = pltWt)
    on.exit(expr = dev.off())
  }


  ##############################################################################
  # 2. Plot raw hits
  if(makePlots > 1){
    plt0 <- dotplot_paf(
      paf = paf,
      round2q = qbin,
      round2t = tbin,
      maxN2show = 100,
      minN2plot = 10) +
      ggtitle(sprintf("%sk raw hits", round(nrow(paf) / 1e3)))
    print(plt0)
  }

  if(verbose)
    cat(sprintf(" ... %sk total hits\n",
                round(nrow(paf) / 1e3)))

  ##############################################################################
  # 3. Initial parsing (paf1)
  # -- top hits by nmatch for each window
  paf1 <- choose_topNmatch(
    paf = paf,
    mapq20nmatchQuantile = mapq20nmatchQuantile,
    mapq60nmatchQuantile = mapq60nmatchQuantile)
  if(verbose)
    cat(sprintf("\t\t%sk (nmatch)", round(nrow(paf1) / 1e3)))
  if(makePlots > 2){
    plt1a <- dotplot_paf(
      paf = paf1,
      round2q = qbin,
      round2t = tbin,
      maxN2show = 50,
      minN2plot = 5) +
      ggtitle(sprintf(
        "%sk hits, dropping those with the very lowest n. matches",
        round(nrow(paf1)/1e3)))
    print(plt1a)
  }

  # -- choose window hits by whether they look repetitive
  paf1 <- choose_topWindowHits(
    paf = paf1,
    maxBestMaps = maxBestMaps,
    topNhits2keep = topNhits2keep)
  if(verbose)
    cat(sprintf(", %sk (bestWindows)", round(nrow(paf1) / 1e3)))
  if(makePlots > 2){
    plt1a <- dotplot_paf(
      paf = paf1,
      round2q = qbin,
      round2t = tbin,
      maxN2show = 30,
      minN2plot = 3) +
      ggtitle(sprintf(
        "%sk hits, keeping those with top nmatches per window",
        round(nrow(paf1)/1e3)))
    print(plt1a)
  }

  # -- Mask problematic regions on the target
  paf1 <- drop_repetitiveTargets(
    paf = paf1,
    round2 = stepSize * 2,
    targetNhitsQuantile = targetNhitsQuantile)
  paf1 <- paf1[,colnames(paf), with = F]

  if(verbose)
    cat(sprintf(", %sk (no target repeats)", round(nrow(paf1) / 1e3)))
  if(makePlots > 0){
    plt1a <- dotplot_paf(
      paf = paf1,
      round2q = qbin,
      round2t = tbin,
      maxN2show = 30,
      minN2plot = 3) +
      ggtitle(sprintf(
        "%sk hits, keeping those with top nmatches per window and dropping repetitive target regions",
        round(nrow(paf1)/1e3)))
    print(plt1a)
  }

  ##############################################################################
  # 4. Synteny parsing (paf2)
  # -- if the paf is getting big, speed things up by choosing local best hits
  paf2 <- choose_localBestHits(
    paf = paf1,
    round2 = stepSize + windowSize)
  paf2 <- choose_localUniqHits(
    paf = paf2,
    round2 = stepSize + windowSize)

  if(verbose)
    cat(sprintf("\n\t\t%sk (fed to MCScanX)", round(nrow(paf2) / 1e3)))
  if(makePlots > 2){
    plt1a <- dotplot_paf(
      paf = paf2,
      round2q = qbin,
      round2t = tbin,
      maxN2show = 10,
      minN2plot = 1) +
      ggtitle(sprintf(
        "%sk hits, initially passed to MCScanX",
        round(nrow(paf2)/1e3)))
    print(plt1a)
  }

  # -- initial MCScanX
  paf2 <- mcscanx_paf(
    paf = paf2,
    nGaps = 50,
    blkSize = 5,
    MCScanX_hCall = MCScanX_hCall,
    tmpDir = tmpDir)
  paf2 <- subset(paf2, !is.na(mcscanxBlkID))

  if(verbose)
    cat(sprintf(", %sk (initial anchors)", round(nrow(paf2) / 1e3)))
  if(makePlots > 1){
    plt1a <- dotplot_paf(
      paf = paf2,
      round2q = qbin,
      round2t = tbin,
      maxN2show = 10,
      minN2plot = 1) +
      ggtitle(sprintf(
        "%sk hits that are the initial MCScanX anchors",
        round(nrow(paf2)/1e3)))
    print(plt1a)
  }

  ##############################################################################
  # 5. Finalize the anchors/blocks
  # -- get orders of paf1 in coordinates of paf3
  if(keepEverything){
    raw <- data.table(paf)
  }else{
    raw <- data.table(paf1)
  }
  # paf2 <<- paf2
  # stepSize <<- stepSize
  # windowSize <<- windowSize
  # buffer <<- buffer
  paf3 <- project_anchorsOnPaf(
    anchorPaf = paf2,
    paf = raw,
    round2 = stepSize + windowSize,
    buffer = buffer)
  if(verbose)
    cat(sprintf("\n\t\t%sk (hits near anchors)", round(nrow(paf3) / 1e3)))
  if(makePlots > 2){
    plt1a <- dotplot_paf(
      paf = paf3,
      round2q = qbin,
      round2t = tbin,
      maxN2show = 100,
      minN2plot = 1) +
      ggtitle(sprintf(
        "%sk hits that are nearby the initial MCScanX anchors",
        round(nrow(paf3)/1e3)))
    print(plt1a)
  }

  # -- final mcscanx
  paf4 <- mcscanx_paf(
    paf = paf3,
    nGaps = 1,
    blkSize = blkSize,
    MCScanX_hCall = MCScanX_hCall,
    tmpDir = tmpDir,
    verbose = FALSE)
  paf4 <- subset(paf4, !is.na(mcscanxBlkID))

  if(verbose)
    cat(sprintf(", %sk (final anchors)", round(nrow(paf4) / 1e3)))
  if(makePlots > 0){
    plt1a <- dotplot_paf(
      paf = paf4,
      round2q = qbin,
      round2t = tbin,
      maxN2show = 100,
      minN2plot = 1) +
      ggtitle(sprintf(
        "%sk hits that are the final MCScanX anchors",
        round(nrow(paf4)/1e3)))
    print(plt1a)
  }

  # 12. final blocks
  paf4[,blkID := dbscan::dbscan(dbscan::frNN(
    x = cbind((qstart + qend) / 2, (tstart + tend) / 2),
    eps = maxDistBtwHitsInBlk),
    minPts = 1)$cluster,
    by = c("qname", "tname")]
  out <- paf4[,c(colnames(paf)[1:13], "blkID"), with = F]
  out[,blkID := as.numeric(as.factor(paste(qname, tname, blkID)))]
  if(verbose)
    cat(sprintf(", in %s blocks\n", uniqueN(out$blkID)))

  if(makePlots > 0){
    plt1a <- dotplot_paf(
      paf = out,
      round2q = qbin,
      round2t = tbin) +
      ggtitle(sprintf(
        "%sk hits that are the final MCScanX anchors, colored by %s blockIDs",
        round(nrow(paf4)/1e3), uniqueN(out$blkID)))
    print(plt1a)
  }

  return(out)
}


#' @title mcscanx_paf
#' @description
#' \code{mcscanx_paf} mcscanx_paf
#' @rdname synteny_paf
#' @import data.table
#' @export
mcscanx_paf <- function(paf,
                        nGaps,
                        blkSize,
                        MCScanX_hCall,
                        tmpDir,
                        verbose = FALSE){

  ##############################################################################
  # 1. Get paths in order
  tmpd <- file.path(
    tmpDir,
    paste0("tmp_", paste(
      sample(c(letters, LETTERS), 20, replace = T),
      collapse = "")))
  if(dir.exists(tmpd))
    unlink(tmpd, recursive = T)
  dir.create(tmpd)

  blFile <- file.path(tmpd, "mcs.homology")
  gfFile <- file.path(tmpd, "mcs.gff")
  colFile <- file.path(tmpd, "mcs.collinearity")

  on.exit(expr = unlink(tmpd, recursive = T))

  hits <- data.table(paf)

  ##############################################################################
  # 2. get a unique identifier (if it doesn't yet exist)
  if(!"qid" %in% names(hits))
    hits[,qid := as.numeric(as.factor(paste(qname, qstart)))]
  if(!"tid" %in% names(hits))
    hits[,tid := as.numeric(as.factor(paste(tname, tstart)))]
  hits[,tid := paste0(tid, "xxxxAddedxxxx")]

  ##############################################################################
  # 3. make gff-like mcscanx file
  u1 <- subset(hits, !duplicated(qid))
  u2 <- subset(hits, !duplicated(tid))

  setkey(u1, qname, qstart)
  setkey(u2, tname, tstart)

  u1[,qchr := paste0(
    "aa", as.numeric(factor(qname, levels = unique(qname))))]
  u2[,tchr := paste0(
    "bb", as.numeric(factor(tname, levels = unique(tname))))]

  mcs1 <- u1[,c("qchr", "qid", "qstart", "qend")]
  setnames(mcs1, c("chr", "id", "start", "end"))

  mcs2 <- u2[,c("tchr", "tid", "tstart", "tend")]
  setnames(mcs2, c("chr", "id", "start", "end"))
  mcsGffIn <- rbind(mcs1, mcs2)

  ##############################################################################
  # 4. make blast-like mcscanx mappings file
  mcsBlsIn <- hits[,c("qid", "tid", "nmatch")]
  setorder(mcsBlsIn, -nmatch)
  mcsBlsIn <- subset(mcsBlsIn, !duplicated(paste(tid, qid)))

  ##############################################################################
  # 5. write both files
  fwrite(
    mcsGffIn, file = gfFile, sep = "\t", quote = FALSE, col.names = FALSE,
    showProgress = FALSE, verbose = FALSE)
  fwrite(
    mcsBlsIn, file = blFile, sep = "\t", quote = FALSE, col.names = FALSE,
    showProgress = FALSE, verbose = FALSE)

  ##############################################################################
  # 6. run MCScanX_h
  mcsCom <- sprintf(
    "-a -b 2 -c 2 -m %s -s %s %s",
    nGaps, blkSize, file.path(tmpd, "mcs"))
  if(verbose)
    cat(sprintf("running MCScanX_h %s ...\n", mcsCom))
  comout <- system2(MCScanX_hCall, mcsCom, stdout = TRUE, stderr = TRUE)
  if(verbose)
    cat("####\n\t", comout, sep = "\n\t")

  ##############################################################################
  # 7. read output
  suppressWarnings(collin <- fread(
    cmd = sprintf("cat %s | grep xxxxAddedxxxx | grep :", colFile),
    col.names = c("blkID","gn1","gn2"),
    select = 1:3,
    showProgress = FALSE,
    header = FALSE))

  ##############################################################################
  # 8. combine with paf and return
  colv <- collin$blkID; names(colv) <- with(collin, paste(gn1, gn2))
  colv <- gsub("-.*", "", colv)
  hits[,mcscanxBlkID := colv[paste(qid, tid)]]
  if(verbose)
    cat(sprintf("%s of %s hits in paf are collinear\n",
                sum(!is.na(hits$mcscanxBlkID)), nrow(hits)))
  coln <- names(paf)
  return(hits[,c(coln, "mcscanxBlkID"), with = F])
}

#' @title calc_pafBlkCoords
#' @description
#' \code{calc_pafBlkCoords} calc_pafBlkCoords
#' @rdname synteny_paf
#' @import data.table
#' @export
calc_pafBlkCoords <- function(paf){
  hits <- data.table(paf)
  if(!"blkID" %in% colnames(hits))
    stop("can't calculate block coordinates without a blkID column\n")
  pafBlkCols <- c("qname", "qlen", "qstart", "qend", "strand",
                  "tname", "tlen", "tstart", "tend",
                  "nmatch", "alen", "mapq", "blkID")
  if(!all(pafBlkCols %in% colnames(hits)))
    stop("paf should have the following column names:",
         paste(pafBlkCols, collapse = ", "), "\n")

  # -- get the columns and complete observations for these
  bhits <- subset(hits, !is.na(blkID))
  if(nrow(bhits) > 0){
    # -- get the genome1 coordinates
    setkey(bhits, qstart)
    blks <- bhits[,list(
      qstart = as.integer(min(as.numeric(qstart))),
      qend = as.integer(max(as.numeric(qend))),
      tmin = as.integer(min(as.numeric(tstart))),
      tmax = as.integer(max(as.numeric(tend))),
      nmatch = as.integer(sum(as.numeric(nmatch))),
      alen = as.integer(sum(as.numeric(alen))),
      mapq = as.integer(median(as.numeric(mapq))),
      strand = ifelse(.N <= 1, "+",
                      ifelse(cor(jitter(qstart),
                                 jitter(tstart)) > 0,"+", "-"))),
      by = c("blkID", "qname","tname", "qlen", "tlen")]
    blks[,`:=`(tstart = ifelse(strand == "+", tmin, tmax),
               tend = ifelse(strand == "+", tmax, tmin),
               tmax = NULL, tmin = NULL)]
    setcolorder(
      blks,
      pafBlkCols)
  }else{
    blks <- NULL
  }

  return(blks)
}
