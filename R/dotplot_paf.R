#' @title dotplot_paf genome-wide dotplotter
#' @description
#' \code{dotplot_paf} lightweight and fast genome-wide dotplotter.
#'
#' @param paf data.table containing the paf-formatted mappings
#' @param gapProp numeric [0-1] giving the size of the gaps between chromosomes
#' as the fraction of the size of the largest chromosome
#' @param round2q integer > 0 specifying the rounding factor to use for the
#' query bins
#' @param round2t integer > 0 specifying the rounding factor to use for the
#' target bins
#' @param maxN2show integer, specifying the top of the color scale. Any bins
#' with more than this number of hits are white.
#' @param minN2plot integer, specifying the minimum number of hits in a bin
#' for that bin to be plotted.
#' @param gridsEvery numeric, how many bases should be between gridlines?
#' @param forceHeatmap logical, should the heatmap be plotted no matter what?
#' @param pointAlpha numeric specifying the opacity of the points
#' @param fixedAxisCoords logical, should the coordinates be fixed?
#'
#' @details coming soon
#'
#' @import data.table
#' @import ggplot2
#' @export
dotplot_paf <- function(paf,
                        gapProp = .05,
                        round2q = 100,
                        round2t = round2q,
                        maxN2show = NULL,
                        pointAlpha = 1,
                        gridsEvery = 10e6,
                        minN2plot = 0,
                        forceHeatmap = FALSE,
                        fixedAxisCoords = TRUE){

  hits <- data.table(paf)
  if(!"blkID" %in% names(hits)){
    hasBlkID <- FALSE
    hits[,blkID := 1]
  }else{
    hasBlkID <- TRUE
  }

  # 1. Get the lengths of the sequences
  qlens <- with(hits, tapply(qlen, qname, min, na.rm = T))
  tlens <- with(hits, tapply(tlen, tname, min, na.rm = T))

  # 2. Calculate the gap sizes for each genome
  qgap <- max(qlens) * gapProp
  tgap <- max(tlens) * gapProp
  qgap <- tgap <- max(c(qgap, tgap))

  # 3. Calculate the linear start positions of chromosomes
  qst <- calc_linStart(seqLens = qlens, gapSize = qgap)
  tst <- calc_linStart(seqLens = tlens, gapSize = tgap)
  qen <- as.numeric(qst + qgap)
  ten <- as.numeric(tst + tgap)

  # 4. Convert paf to x/y data.table
  tp <- with(hits, data.table(
    x = ((qstart + qend) / 2) + qst[qname],
    y = ((tstart + tend) / 2) + tst[tname],
    blkID = blkID))

  # 5. Round the positions
  round2q <- as.integer(round2q)
  round2t <- as.integer(round2t)
  if(!is.integer(round2q))
    stop("round2q must be coercible to an integer\n")
  if(!is.integer(round2t))
    stop("round2t must be coercible to an integer\n")
  if(round2q > 0)
    tp[,x := round_toInteger(x, round2q)]
  if(round2t > 0)
    tp[,y := round_toInteger(y, round2t)]

  # 6. Calculate the number of observations within each bin
  tp <- tp[, list(n = .N), by = c("x", "y", "blkID")]

  # 7. Subset dataset to bins with > minN2plot hits
  tp <- subset(tp, n >= minN2plot)

  # 8. Transform n hits to a color palette
  if(is.null(maxN2show))
    maxN2show <- ceiling(quantile(tp$n[tp$n > 1], .99))
  tp[,col := ifelse(n > maxN2show, maxN2show, n)]
  setorder(tp, col)

  if(minN2plot <= 1){
    pltCols <- c("grey10", pal_deepspace(10)[-c(1:2)])
    pltCols[2] <- "darkslateblue"
  }else{
    pltCols <- pal_deepspace(10)[-c(1:2)]
    pltCols[1] <- "darkslateblue"
  }

  # 9. Get x and y axis label positions
  xlabShift <- qst + (qlens / 2)
  ylabShift <- tst + (tlens / 2)

  # 10. Get the gridline positions
  xgrd <- seq(from = gridsEvery, to = max(qst + qlens), by = gridsEvery)
  ygrd <- seq(from = gridsEvery, to = max(tst + tlens), by = gridsEvery)

  # 11. Get the chromosome bound breaks to make it look facetted
  brks <- data.table(
    xmin = c(qst[-1] - qgap,
             rep(0, length(tst) - 1)),
    xmax = c(qst[-1],
             rep(max(qst + qlens), length(tst) - 1)),
    ymin = c(rep(0, length(qst) - 1),
             tst[-1] - tgap),
    ymax = c(rep(max(tst + tlens), length(qst) - 1),
             tst[-1]))

  # 12. Make the plot
  plt <- ggplot()+

    # -- gridlines
    geom_hline(
      yintercept = ygrd,
      col = add_alpha("white", .5),
      linetype = 3,
      linewidth = .1)+
    geom_vline(
      xintercept = xgrd,
      col = add_alpha("white", .5),
      linetype = 3,
      linewidth = .1)+

    # -- gaps between chromosomes
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey20",
      color = NA,
      data = brks)+



    # -- axis labels and expansion
    scale_x_continuous(
      expand = c(0.001, 0.001),
      breaks = xlabShift,
      labels = names(xlabShift))+
    scale_y_continuous(
      expand = c(0.001, 0.001),
      breaks = ylabShift,
      labels = names(ylabShift))+

    # -- plot themes
    theme(
      panel.background = element_rect(fill = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text.x = element_text(
        angle= 90, family = "Helvetica", size = 5, vjust = .5, hjust = 1),
      axis.text.y = element_text(
        angle= 0, family = "Helvetica", size = 5, vjust = .5, hjust = 1),
      axis.title = element_text(family = "Helvetica", size = 6),
      axis.ticks = element_blank(),
      plot.title = element_text(family = "Helvetica", size = 7))+

    # -- axis labels
    labs(
      x = sprintf("Query genome (%sMbp), gridlines every %sMb",
                  round(sum(qlens)/1e6, 1),
                  format(round(gridsEvery/1e6, 2), scientific = F)),
      y = sprintf("Target genome (%sMbp), gridlines every %sMb",
                  round(sum(tlens)/1e6, 1),
                  format(round(gridsEvery/1e6, 2), scientific = F)),
      color = "n hits")

  if(forceHeatmap | !hasBlkID){
    # -- colored bins as points
    plt <- plt +
      geom_point(
        aes(x = x, y = y, col = col),
        pch = ".",
        alpha = pointAlpha,
        data = tp) +
      # -- color palette
      scale_color_gradientn(colors = pltCols)
  }else{
    # -- colored by blocks
    cols <- c("#C4645C", "#F5915C", "#FFC765",
              "#FCF8D5", "#BEF3F9", "#66B8FF", "#6666FF", "#9C63E1",
              "#F4BDFF")
    pal <- colorRampPalette(cols)(uniqueN(hits$blkID))
    plt <- plt +
      geom_point(
        aes(x = x, y = y, col = as.factor(blkID)),
        pch = ".",
        alpha = pointAlpha,
        data = tp) +
      # -- color palette
      scale_color_manual(values = sample(pal), guide = "none")
  }

  # 13. If necessary add fixed coordinate system to the plot
  if(fixedAxisCoords)
    plt <- plt + coord_fixed()

  return(plt)
}
