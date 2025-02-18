#' @title mcscanx_paf
#' @description
#' \code{mcscanx_paf} mcscanx_paf
#' @name mcscanx_paf
#'
#' @param paf DNAStringSet containing the fasta to be windowed
#' @param nGaps DNAStringSet containing the fasta to serve as the target
#' @param blkSize character string coercible to a file.path where the temporary
#' files should be written. Taken as the current working directory by default.
#' @param MCScanX_hCall character string coercible to a file path specifying the
#' location of the output (paf) file
#' @param tmpDir integer, the basepair size of the window in the query genome
#' @param verbose integer, the basepair size of the step between windows
#'
#' @details coming soon
#'
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
  qname <- qstart <- tname <- tstart <- qid <- tid <- qchr <- tchr <- NULL
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
  nmatch <- qid <- tid <- NULL
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
  mcscanxBlkID <- gn1 <- gn2 <- NULL
  if(nrow(collin) > 0){
    colv <- collin$blkID; names(colv) <- with(collin, paste(gn1, gn2))
    colv <- gsub("-.*", "", colv)
    hits[,mcscanxBlkID := colv[paste(qid, tid)]]
    if(verbose)
      cat(sprintf("%s of %s hits in paf are collinear\n",
                  sum(!is.na(hits$mcscanxBlkID)), nrow(hits)))
    coln <- names(paf)
    return(hits[,c(coln, "mcscanxBlkID"), with = F])
  }else{
    hits[,mcscanxBlkID := NA]
    coln <- names(paf)
    return(hits[0,c(coln, "mcscanxBlkID"), with = F])
  }
}
