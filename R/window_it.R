#' @title flexible sliding window constructor
#' @description
#' \code{window_it} Build a granges sliding window from a seqinfo, stringset or
#' granges object. Identical output to GenomicRanges::slidingWindows(), except
#' allowing non-overlapping windows.
#'
#' @param x GRanges, DNAStringSet or SeqInfo object. The start and end
#' coordinates for the windows are taken from this.
#' @param width integer of length 1, specifying the window size in basepairs
#' @param step integer of length 1, specifying the step size in basepairs
#' @param minWindowSize numeric of length 1, specifying the minimum window size
#' to be retained. Defaults to 50% of the width.
#'
#' @details Constructs the windows in data.table
#'
#' @return GRanges object with the window ranges without considering strand
#'
#' @examples
#' \dontrun{
#' # coming soon.
#' }
#'
#' @importFrom GenomicRanges GRanges start end makeGRangesFromDataFrame seqnames
#' @importFrom data.table data.table foverlaps CJ setnames setkey :=
#' @importFrom Biostrings seqinfo
#' @importFrom methods is
#' @export
window_it <- function(x,
                      width,
                      step,
                      minWindowSize = width / 2){

  ignoreMd <- TRUE
  # -- if stringset of seqinfo, convert to granges
  if(is(x, "DNAStringSet")){
    x <- GRanges(Biostrings::seqinfo(x))
  }else{
    if(is(x, "Seqinfo")){
      x <- GRanges(x)
    }
  }

  # -- check if granges and if so, split by region
  if(is(x, "GRanges")){
    # -- for a gr with no metadata, just add a filler column 'tmp'
    y <- as.data.frame(x)
    hasMetadata <- ncol(y) > 5
    if(!hasMetadata){
      x$tmp <- 1:length(x)
    }else{
      ignoreMd <- FALSE
      # -- for a gr with 1> md, join them and rename as 'tmp'
      if(ncol(y) > 6){
        warning(">1 metadata columns observed - these will be joined\n")
        x$tmp <- apply(elementMetadata(x), 1, paste, sep = "_XXX_")
      }else{
        x$tmp <- y[[6]]
      }
    }

    x <- makeGRangesFromDataFrame(data.table(as.data.frame(x))[,list(
      start = min(start), end = max(end)), by = c("seqnames", "tmp")],
      keep.extra.columns = TRUE)
  }else{
    # -- if not gr, error out
    stop("x can only be a DNAStringSet, Granges of Seqinfo object\n")
  }

  # -- check other parameters
  width <- as.integer(width)
  if(length(width) != 1 || is.na(width) || is.null(width) || width <= 2)
    stop("width must be a single integer with value > 2\n")

  step <- as.integer(step)
  if(length(step) != 1 || is.na(step) || is.null(step) || step < 1)
    stop("step must be a single integer with value > 0\n")

  minWindowSize <- as.numeric(minWindowSize)
  if(length(minWindowSize) != 1 || is.na(minWindowSize) ||
     is.null(minWindowSize) || minWindowSize < 1)
    stop("minWindowSize must be a single integer with value > 0\n")

  # -- build out a data.table of the chromosomes
  seqnames <- start <- end <- NULL
  dt <- data.table(
    seqnames = as.character(seqnames(x)),
    start = start(x),
    end = end(x),
    grp = x$tmp)

  # -- add the starts
  dt <- dt[,list(
    start = seq(from = min(start),
                to = max(end),
                by = step),
    maxp = max(end)),
    by = c("seqnames", "grp")]

  # -- add the ends
  dt[,end := start + (width - 1)]

  # -- subset ends so they are within the bounds of the ranges
  dt$end[dt$end > dt$maxp] <- dt$maxp[dt$end > dt$maxp]

  # -- remove any windows too small
  windWidth <- NULL
  dt[,windWidth := end - start]
  dt <- subset(dt, windWidth > minWindowSize)

  # -- return the granges object
  gr <- makeGRangesFromDataFrame(
    dt[,c("seqnames", "start", "end")])
  if(!ignoreMd)
    gr$region <- dt$grp
  return(gr)
}
