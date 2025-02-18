#' @title Convenience functions related to GRanges objects
#' @description
#' \code{gr_utils} Simplify common scripting goals on GRanges objects
#' @name utils
#'
#' @param x GRanges object
#' @param min.gapwidth passed to GenomicRanges::reduce specifying the minimum
#' gap across which ranges should be reduced
#' @param ... additional parameters, typically passed to reduce(), but depends.
#' \cr
#' If called, \code{gr_utils} returns its own arguments.
#'

#' @title join adjacent ranges
#' @description
#' \code{join_gappedRanges} Reduces ranges, returning a granges with proportion
#' of sequence and total sequence covered.
#' @rdname gr_utils
#' @importFrom GenomicRanges sort reduce width
#' @importFrom methods is
#' @export
join_gappedRanges <- function(x,
                              min.gapwidth = 0){
  # -- check input parameters
  if(!is(x, "GRanges"))
    stop("x must be a GRanges object\n")

  min.gapwidth <- as.integer(min.gapwidth)
  if(length(min.gapwidth) != 1 || is.null(min.gapwidth) || is.na(min.gapwidth))
    stop("min.gapwidth must be a single integer value\n")

  # -- sort and combine all overlapping ranges
  gr <- sort(x)
  gr <- reduce(
    gr,
    drop.empty.ranges = TRUE,
    ignore.strand = TRUE)

  # -- combine overlapping ranges with specified width, keeping mapping info
  out <- reduce(
    gr,
    min.gapwidth = min.gapwidth,
    drop.empty.ranges = TRUE,
    ignore.strand = TRUE,
    with.revmap = TRUE)

  # -- extract all ranges that were not joined with gaps
  out1 <- subset(out, sapply(out$revmap, length) == 1)
  out1$matchBases <- width(out1)
  out1$matchProp <- 1

  # -- extract ranges joined by gaps
  out2 <- subset(out, sapply(out$revmap, length) > 1)

  # -- calculate amount of sequence mapped in each range and combine with out1
  if(length(out2) > 0){
    rmap <- as.list(out2$revmap)
    wds <- sapply(rmap, function(x) sum(width(gr[x])))
    out2$matchBases <- wds
    out2$matchProp <- wds/width(out2)

    out1 <- c(out1, out2)
  }

  out1$revmap <- NULL
  return(out1)
}

#' @title Reduce granges by group
#' @description
#' \code{reduce_grByGroup} Join overlapping ranges within unique values of the
#' first elementMetadata column of the GRanges object
#' @rdname gr_utils
#' @importFrom GenomicRanges elementMetadata ranges reduce
#' @importFrom data.table uniqueN
#' @importFrom methods is
#' @export
reduce_grByGroup <- function(x, ...){

  # -- check input parameters
  if(!is(x, "GRanges"))
    stop("x must be a GRanges object\n")

  # -- if there is no element metadata just reduce
  if(length(names(attr(x, "elementMetadata"))) == 0){
    x <- reduce(x, ...)
  }else{

    # -- same if all same value
    grpDat <- elementMetadata(x)[,1]
    if(length(unique(grpDat)) == 1){
      x <- reduce(x, ...)
    }else{

      # -- for each unique element, reduce
      splb <- NULL
      tmp <- x
      spl <- split(x, splb)
      x <- unlist(reduce(spl, ...))

      # -- add elementMetadata back in
      ren <- names(elementMetadata(tmp))[1]
      x$tmp <- names(ranges(x))
      names(attr(x, "elementMetadata")) <- ren
      names(ranges(x)) <- NULL
    }
  }
  return(x)
}

#' @title Make granges hierarchical
#' @description
#' \code{make_hierarchical} Ensure ranges present in the previous group are
#' masked
#' @rdname gr_utils
#' @importFrom GenomicRanges elementMetadata makeGRangesFromDataFrame setdiff
#' @importFrom methods as is
#' @export
make_hierarchical <- function(x){

  # -- check input parameters
  if(!is(x, "GRanges"))
    stop("x must be a GRanges object\n")

  # -- if no elementMetadata or all same value, return
  if(ncol(elementMetadata(x)) == 0)
    return(x)
  grpDat <- elementMetadata(x)[,1]
  if(length(unique(grpDat)) == 1)
    return(x)

  # -- ensure we respect the order of data in grpDat
  if(is.factor(grpDat)){
    u <- levels(grpDat)
  }else{
    u <- unique(grpDat)
  }

  # -- clean out granges object
  x <- makeGRangesFromDataFrame(as.data.frame(x))

  # -- add group data back in
  x$grp <- as.character(grpDat)

  # -- split and order the values
  spl <- split(x, x$grp)[u]

  # -- toss anything with no observations
  spl <- spl[sapply(spl, length) > 0]

  # -- choose the first value as the start
  out <- spl[[1]]
  grps <- names(spl)

  # -- for each group, remove overlapping ranges already in the out file
  for(i in grps[-1]){
    xi <- spl[[i]]
    tmp <-  setdiff(xi, out)
    if(length(tmp) > 0){
      tmp$grp <- i
      out <- c(out, tmp)
    }
  }

  # -- ensure formatting is right and return
  out$grp <- factor(out$grp, levels = grps)
  out <- out[order(out$grp),]

  return(out)
}
