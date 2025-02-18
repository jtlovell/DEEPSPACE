#' @title Calculate block coordinates
#' @description
#' \code{calc_pafBlkCoords} Convert hits data.table to blocks
#' @name synteny_windows
#'
#' @param paf data.table with the 12 standard paf columns and "blkID" specifying
#' the syntenic block for each hit.
#'
#' @details coming soon
#'
#' @import data.table
#' @importFrom stats cor median
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
  strand <- tmax <- tmin <- mapq <- nmatch <-
    tend <- tstart <- qend <- qstart <- blkID <- NULL
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
      alen = as.integer(diff(range(c(tstart, tend)))),
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
