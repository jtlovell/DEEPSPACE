#' @title Calculate points along a curved polygon
#' @description
#' \code{calc_curvePolygon} Use block coordinates to make a polygon (like
#' riparian plot braids).
#' @name calc_curvePolygon
#'
#' @param start1 numeric, the left value of the first coordinates
#' @param start2 numeric, the left value of the second coordinates
#' @param end1 numeric, the right value of the first coordinates
#' @param end2 numeric, the right value of the second coordinates
#' @param y1 numeric, the y value of the first coordinates
#' @param y2 numeric, the y value of the second coordinates
#' @param npts integer, number of points to use
#' @param keepat integer, the rounding of y to use
#'
#' @details coming soon
#'
#' @import data.table
#' @importFrom stats cor median
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
