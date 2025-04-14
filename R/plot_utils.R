#' @title Generic internal functions used by DEEPSPACE
#' @description
#' \code{plot_utils} Convenience functions for DEEPSPACE, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name plot_utils
#'
#' @param paf data.table containing a windowed paf alignment file
#' @param blks a paf file read into R with the following additional fields:
#' query: name of the query genome, target: name of the target genome, qy: y
#' coordinate of the the query genome, ty: y coordinate of the target genome,
#' color: the name of the color for the alignment. This object should contain
#' all of the pafs for each pairwise contrast together.
#' @param chrordQuery vector of chromosome names in the order to be plotted. For
#' the query genome
#' @param chrordTarget vector of chromosome names in the order to be plotted. For
#' the target genome
#' @param gapProp numeric (0-1) giving the proportion of the x axis that should
#' be gaps between chromosomes
#' @param howRound numeric (0-Inf) specifying how rouded the chromosome
#' polygons should be.
#' @param yQuery numeric, y-value of the query genome
#' @param yTarget numeric, y-value of the target genome
#' @param braidOffset distance between the chromosome midpoint and the y-start
#' of each daisychain braid. If the same value as chrWidth, the braid starts
#' exactly at the edge of the chromosome polygon
#' @param xoffsetQuery numeric, how far to the left should the query genome be
#' offset?
#' @param xoffsetTarget numeric, how far to the left should the target genome be
#' offset?
#' @param sameXrange logical, should both the query and target genome have the
#' same x range?
#' @param colorColumn character, matching one of the paf columns that contains
#' the color for that block
#' @param chrWidth numeric, providing the thinkness of the chromosomes
#' @param xleft left x-position
#' @param ybottom bottom y-position
#' @param xright right x-position
#' @param ytop top y-position
#' @param plotWidth width of the plot
#' @param plotHeight height of the plot
#' @param xrange range of x-values
#' @param yrange range of y-values
#' @param npts number of points that make up the polygon
#' @param x value to transform
#' @param min numeric, minimum value
#' @param max numeric, maximum value
#' @param scale1toMean logical, should the mean be scaled to 1?
#' @param col color
#' @param alpha opacity

#' If called, \code{plot_utils} returns its own arguments.
#'

#' @title order_bySynteny
#' @description
#' \code{order_bySynteny} order_bySynteny
#' @rdname plot_utils
#' @importFrom GenomicRanges start end width slidingWindows
#' @import Biostrings
#' @importFrom Rsamtools scanFaIndex
#' @export
order_bySynteny <- function(paf,
                            chrordTarget){
  nh <- paf[,list(n = sum(alen)), by = c("tname", "qname")]
  nh[,tname := factor(tname, levels = chrordTarget)]
  setorder(nh, tname, -n)
  nh <- subset(nh, complete.cases(nh))
  nh[,o := as.numeric(tname)]
  wm <- nh[,list(wtmean = weighted.mean(o, n)), by = "qname"]
  setkey(wm, wtmean)
  return(wm$qname)
}

#' @title transform_paf2riparian
#' @description
#' \code{transform_paf2riparian} transform_paf2riparian
#' @rdname plot_utils
#' @importFrom GenomicRanges start end width slidingWindows
#' @import Biostrings
#' @importFrom Rsamtools scanFaIndex
#' @export
transform_paf2riparian <- function(paf,
                                   chrordQuery,
                                   chrordTarget,
                                   gapProp = 0.05,
                                   howRound = .5,
                                   yQuery = 0,
                                   yTarget = 1,
                                   braidOffset = .1,
                                   xoffsetQuery = 0,
                                   xoffsetTarget = 0,
                                   sameXrange = FALSE,
                                   colorColumn = NULL,
                                   chrWidth = braidOffset){
  blks <- data.table(paf)
  blks[, qname := as.character(qname)]
  blks[, tname := as.character(tname)]
  plWt <- dev.size()[1]
  plHt <- dev.size()[2]

  chrordQuery <- chrordQuery[!duplicated(chrordQuery)]
  chrordTarget <- chrordTarget[!duplicated(chrordTarget)]

  if(sameXrange)
    blks <- transform_paf2sameRange(blks, xrange = c(0,1))

  # -- get the linear positions of blocks and chromosomes
  qlens <- with(blks, tapply(qlen, qname, min))
  if(!all(names(chrordQuery) %in% names(qlens)))
    stop("supplied query chromosome order does not match qnames\n")
  chrordQuery <- chrordQuery[chrordQuery %in% names(qlens)]
  qlens <- qlens[chrordQuery]
  tlens <- with(blks, tapply(tlen, tname, min))
  if(!all(names(chrordTarget) %in% names(tlens)))
    stop("supplied target chromosome order does not match tnames\n")
  chrordTarget <- chrordTarget[chrordTarget %in% names(tlens)]
  tlens <- tlens[chrordTarget]

  qgap <- max(qlens) * gapProp
  tgap <- max(tlens) * gapProp

  qst <- calc_linStart(seqLens = qlens, gapSize = qgap) + xoffsetQuery
  tst <- calc_linStart(seqLens = tlens, gapSize = tgap) + xoffsetTarget

  # -- get chromosome coordinates
  qen <- as.numeric(qst + qlens)
  names(qen) <- names(qst)
  ten <- as.numeric(tst + tlens)
  names(ten) <- names(tst)

  xr <-  howRound * range(c(tst, ten, qst, qen))

  qchr <- rbindlist(lapply(names(qst), function(i)
    data.table(chrID = i, round_rect(
      xleft = qst[i], xright = qen[i],
      ybottom = yQuery - (chrWidth),
      ytop = yQuery + (chrWidth),
      yrange = c(0,1),
      xrange = xr,
      plotWidth = plWt,
      plotHeight = 1))))

  tchr <- rbindlist(lapply(names(tst), function(i)
    data.table(chrID = i, round_rect(
      xleft = tst[i], xright = ten[i],
      ybottom = yTarget - (chrWidth),
      ytop = yTarget + (chrWidth),
      yrange = c(0,1),
      xrange = xr,
      plotWidth = plWt,
      plotHeight = 1))))

  tp <- data.table(blks)
  tp[,`:=`(qstart = qstart + qst[qname],
           qend = qend + qst[qname],
           tstart = tstart + tst[tname],
           tend = tend + tst[tname])]
  if(is.null(colorColumn)){
    tp[,col := NA]
  }else{
    tp[,col := tp[[colorColumn]]]
  }
  braids <- rbindlist(lapply(1:nrow(tp), function(i)
    data.table(
      index = tp$index[i], qname = tp$qname[i], tname = tp$tname[i], col = tp$col[i],
      with(tp[i,], calc_curvePolygon(
        start1 = qstart, end1 = qend,
        start2 = tstart, end2 = tend,
        y1 = ifelse(yQuery < yTarget,
                    yQuery + braidOffset,
                    yQuery - braidOffset),
        y2 = ifelse(yQuery < yTarget,
                    yTarget - braidOffset,
                    yTarget + braidOffset))))))
  if(is.null(colorColumn)){
    braids[,col := NULL]
  }else{
    setnames(braids, "col", colorColumn)
  }
  return(list(queryChrs = qchr,
              targetChrs = tchr,
              braids = braids))
}

#' @title calculate coordinates for rounded rectange polygons
#' @description
#' \code{round_rect} from x-y coordinates, make a rounded rectangle
#' @rdname plot_utils
#' @importFrom graphics par
#' @importFrom grDevices dev.size
#' @export
round_rect <- function(xleft,
                       ybottom,
                       xright,
                       ytop,
                       plotWidth,
                       plotHeight,
                       xrange,
                       yrange,
                       npts = 50){

  if (ytop <= ybottom){
    tmp <- ytop
    ytop <- ybottom
    ybottom <- tmp
  }
  if (xright <= xleft){
    tmp <- xleft
    xleft <- xright
    xright <- tmp
  }

  # measure graphics device
  pars <- c(xrange, yrange)
  asp <- diff(pars[3:4]) / diff(pars[1:2])
  dev <- plotWidth / plotHeight

  # make a curve and split into left and right
  radius <- (ytop - ybottom) / 2
  centerY <- ytop - radius
  centerX <- mean(c(xleft, xright))
  theta <- seq(0, 2 * pi, length = npts)
  circX <- cos(theta)
  circY <- sin(theta)
  leftC <- which(circX <= 0)
  rightC <- which(circX >= 0)

  xR <- circX[rightC]
  yR <- circY[rightC]
  ordYR <- rev(order(yR))
  xR <- xR[ordYR]
  yR <- yR[ordYR]

  xL <- circX[leftC]
  yL <- circY[leftC]
  ordYL <- order(yL)
  xL <- xL[ordYL]
  yL <- yL[ordYL]

  # project onto graphics device and scale
  xRightS <- xright - (radius / asp / dev)
  xLeftS <- xleft + (radius / asp / dev)
  if (centerX < xLeftS)
    xLeftS <- centerX
  if (centerX > xRightS)
    xRightS <- centerX
  xLS <- scale_between(xL, xleft, xLeftS)
  xRS <- scale_between(xR, xRightS, xright)
  yLS <- scale_between(yL, ybottom, ytop)
  yRS <- scale_between(yR, ybottom, ytop)
  return(data.table(x = c(xRS,xLS), y = c(yRS, yLS)))
}

#' @title transform_paf2sameRange
#' @description
#' \code{transform_paf2sameRange} transform_paf2sameRange
#' @rdname plot_utils
#' @import data.table
#' @export
transform_paf2sameRange <- function(paf,
                                    xrange){
  pafOut <- data.table(paf)
  xrange <- xrange[!duplicated(xrange)]
  xrange <- xrange[order(xrange)]
  xdiff <- diff(xrange)

  qlens <- with(paf, tapply(qlen, qname, min))
  tlens <- with(paf, tapply(tlen, tname, min))
  qtots <- sum(qlens)
  ttots <- sum(tlens)

  qscl <- xdiff / qtots
  tscl <- xdiff / ttots

  xmin <- xrange[1]

  pafOut[,`:=`(
    qlen = (qlen * qscl) + xmin,
    qstart = (qstart * qscl) + xmin,
    qend = (qend * qscl) + xmin,
    tlen = (tlen * tscl) + xmin,
    tstart = (tstart * tscl) + xmin,
    tend = (tend * tscl) + xmin)]
  return(pafOut)
}

#' @title scale_between
#' @description
#' \code{scale_between} scale_between
#' @rdname plot_utils
#' @import data.table
#' @export
scale_between <- function(x,
                          min,
                          max,
                          scale1toMean = TRUE){
  if(length(unique(x)) > 1){
    return((x - min(x)) / (max(x) - min(x)) * (max - min) + min)
  }else{
    if(scale1toMean){
      return(mean(c(min, max)))
    }else{
      return(max)
    }
  }
}

#' @title add transparency
#' @description
#' \code{add_alpha} add transparency to a color
#' @rdname plot_utils
#' @importFrom grDevices col2rgb rgb
#' @export
add_alpha <- function(col, alpha = 1){

  if(missing(col) || !all(are_colors(col)))
    stop("Colors are misspecified\n")
  if(length(alpha) != 1 || alpha > 1 || alpha < 0)
    stop("alpha is misspecified\n")

  return(apply(sapply(col, col2rgb)/255, 2, function(x)
    rgb(x[1], x[2], x[3], alpha = alpha)))
}

#' @title check if a vector is coercible to R colors
#' @description
#' \code{are_colors} check if a vector is coercible to R colors
#' @rdname plot_utils
#' @importFrom grDevices col2rgb
#' @export
are_colors <- function(col) {
  sapply(col, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}

#' @title convert cosine points to polygon
#' @description
#' \code{calc_curvePolygon} from 2d coordinates, make a curve
#' @rdname plot_utils
#' @export
calc_curvePolygon <- function(start1,
                              end1 = NULL,
                              start2,
                              end2 = NULL,
                              y1,
                              y2,
                              npts = 250,
                              keepat = round(npts / 20)){
  cosine_points <- function(npts, keepat){
    # initial number of points
    # grid to keep always
    grid <- seq(from = 0, to = pi, length.out = npts) # grid
    x <- (1 - cos(grid)) / max((1 - cos(grid))) # scaled cosine
    y <- grid / max(grid) # scaled grid
    # calculate slope for each point
    x1 <- x[-1];  y1 <- y[-1]
    x2 <- x[-length(x)];  y2 <- y[-length(y)]
    s <-  (y1 - y2) / (x1 - x2)
    # choose points that capture changes in slope
    ds <- cumsum(abs(diff(s)))*5
    wh <- c(1,which(!duplicated(round(ds))), length(x))
    wh2 <- c(wh, seq(from = 0, to = length(x), by = round(keepat)))
    wh <- c(wh, wh2)[!duplicated(c(wh, wh2))]
    wh <- wh[order(wh)]
    return(cbind(x[wh], y[wh]))
  }

  scaledCurve <- cosine_points(npts = npts, keepat = keepat)
  # print(scaledCurve)
  if (!is.null(end1) | !is.null(end2)) {
    sc1 <- scaledCurve[,1]
    sc2 <- scaledCurve[,2]

    tp <- rbind(
      start1 = data.table(
        x = start1, y = y1),
      poly1 = data.table(
        x = scale_between(x = sc1, min = start1, max = start2),
        y = scale_between(x = sc2, min = y1, max = y2)),
      start2 = data.table(x = start2, y = y2),
      end2 = data.table(
        x = end2, y = y2),
      poly2 = data.table(
        x = scale_between(x = sc1, min = end2, max = end1),
        y = scale_between(x = sc2, min = y2, max = y1)),
      end1 = data.table(
        x = end1, y = y1))
  }else{
    tp <- data.table(
      x = scale_between(x = scaledCurve[,1], min = start1, max = start2),
      y = scale_between(x = scaledCurve[,2], min = y1, max = y2))
  }

  return(tp)
}

#' @title calculate coordinates for rounded rectange polygons
#' @description
#' \code{round_rect} from x-y coordinates, make a rounded rectangle
#' @rdname plot_utils
#' @importFrom graphics par
#' @importFrom grDevices dev.size
#' @export
round_rect <- function(xleft, ybottom, xright, ytop, plotWidth, plotHeight, xrange, yrange, npts = 50){

  if (ytop <= ybottom){
    tmp <- ytop
    ytop <- ybottom
    ybottom <- tmp
  }
  if (xright <= xleft){
    tmp <- xleft
    xleft <- xright
    xright <- tmp
  }

  # measure graphics device
  pars <- c(xrange, yrange)
  asp <- diff(pars[3:4]) / diff(pars[1:2])
  dev <- plotWidth / plotHeight

  # make a curve and split into left and right
  radius <- (ytop - ybottom) / 2
  centerY <- ytop - radius
  centerX <- mean(c(xleft, xright))
  theta <- seq(0, 2 * pi, length = npts)
  circX <- cos(theta)
  circY <- sin(theta)
  leftC <- which(circX <= 0)
  rightC <- which(circX >= 0)

  xR <- circX[rightC]
  yR <- circY[rightC]
  ordYR <- rev(order(yR))
  xR <- xR[ordYR]
  yR <- yR[ordYR]

  xL <- circX[leftC]
  yL <- circY[leftC]
  ordYL <- order(yL)
  xL <- xL[ordYL]
  yL <- yL[ordYL]

  # project onto graphics device and scale
  xRightS <- xright - (radius / asp / dev)
  xLeftS <- xleft + (radius / asp / dev)
  if (centerX < xLeftS)
    xLeftS <- centerX
  if (centerX > xRightS)
    xRightS <- centerX
  xLS <- scale_between(xL, xleft, xLeftS)
  xRS <- scale_between(xR, xRightS, xright)
  yLS <- scale_between(yL, ybottom, ytop)
  yRS <- scale_between(yR, ybottom, ytop)
  return(data.table(x = c(xRS,xLS), y = c(yRS, yLS)))
}
