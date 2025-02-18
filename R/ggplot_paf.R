#' @title dot and riparian plots
#' @description
#' \code{ggplot_paf} lightweight and fast genome-wide dotplotter and riparian
#' mapper for two genomes
#' @name ggplot_paf

#' @param paf data.table containing the paf-formatted mappings
#' @param colorColumn character string matching the name of a column in the paf-
#' like data.frame. If paf contains a column 'blkID', then this is taken by
#' default and must be turned off by setting colorColumn = NULL.
#' @param chrColor value coercible to a color for the chromosome separater lines
#' @param gapProp numeric (typically close to 0) specifying the proportion of
#' the total plot that should be gaps between chromosomes
#' @param colorPalette R color palette function for variable coloring
#' @param pointAlpha numeric specifying the opacity of the points
#' @param randomizeColors logical, should the colors be randomized so adjancent
#' blocks (etc) are not similarly colored?
#' @param pointShape value coercible to pch (point shape)
#' @param pointAlpha numeric, passed to geom_point alpha
#' @param defaultColor value coercible to a color for the points
#' @param addThemes ggplot2 theme to be added to the plot. can also be added
#' with + theme(xx) to the output
#' @param chrWidth numeric, passed to linewidth for the chromosome bounds
#' @param minWidth4segment integer, specifying the block width value below which
#' hits are represented as points and above which as lines. By default must be
#' larger than 1/100th of the largest target chromosome.
#' @param segLineWidth for blocks with alen > minWidth4segment, the line
#' width of the segment.
#' @param chrBuffer numeric (typically positive and close to 0) specifying
#' how far from the chromosome midpoints the braids should start
#' @param chrLabBackground value coercible to color giving the background color
#' of the chromosome labels
#' @param chrLabSize numeric value giving the font size of the chromosome labels
#' @param chrLabColor value coercible to R color specying the text label color
#' @param chrLabOffset numeric (typically positive and close to 0) specifying
#' how far the chromosome labels should be from the midpoints
#' @param genomeLabOffset numeric (typically positive and close to 0) specifying
#' how far from the chromosome midpoints the genome labels should be
#' @param queryLab character string to append to the query genome label
#' @param targetLab character string to append to the target genome label
#' @param braidAlpha numeric [0-1] specifying the opacity of the braids
#' @param reorderByColumn character string pointing to a column name to order by
#' @param ... passed to the main geom
#' \cr
#' If called, \code{utils} returns its own arguments.
#'  Undocumented arguments in Rd file 'ggplot_paf.Rd'
#' @details coming soon
#'
#' @title paf dotplotter
#' @description
#' \code{dotplot_paf} make a dotplot from a paf-like data.frame or similar. If
#' regions of various size are given, small ones are represented by dots, while
#' larger ones are segments.
#' @rdname ggplot_paf
#' @import data.table
#' @importFrom ggplot2 ggplot element_rect scale_x_continuous geom_segment aes
#' scale_y_continuous scale_color_identity geom_point
#' @export
dotplot_paf <- function(paf,
                        colorColumn = "blkID",
                        chrColor = add_alpha("grey50", .25),
                        gapProp = 0,
                        chrWidth = .5,
                        colorPalette = pal_ds,
                        randomizeColors = TRUE,
                        pointShape = ".",
                        defaultColor = add_alpha("white", .5),
                        addThemes = theme(panel.background = element_rect(fill = "black")),
                        minWidth4segment = NULL,
                        segLineWidth = .5,
                        ...){

  # -- check parameters
  paf <- data.table(paf)
  if(nrow(paf) == 0)
    stop("paf input is empty")
  if(ncol(paf) < 12)
    stop("paf input appears to be malformed")

  colorColumn <- as.character(colorColumn[1])
  if(!colorColumn %in% names(paf) || length(colorColumn) != 1)
    colorColumn <- as.character(NULL)

  gapProp <- as.numeric(gapProp)[1]
  if(is.na(gapProp) || is.null(gapProp))
    stop("gapProp has to be a single numeric value\n")

  chrWidth <- as.numeric(chrWidth)[1]
  if(is.na(chrWidth) || is.null(chrWidth))
    stop("chrWidth has to be a single numeric value\n")

  if(!is.function(colorPalette))
    stop("colorPalette must be a function that produces a vector of colors\n")
  test <- colorPalette(2)
  if(!all(are_colors(test)))
    stop("colorPalette must be a function that produces a vector of colors\n")
  if(length(test) != 2)
    stop("colorPalette function must produce a color vector of length n\n")

  minWidth4segment <- as.integer(minWidth4segment)
  if(length(minWidth4segment) == 0 || is.na(minWidth4segment))
    minWidth4segment <- max(paf$tlen) / 1000

  chrColor <- chrColor[1]
  if(!are_colors(chrColor)){
    warning("chrColor must be coercible to a color. setting to lightgrey\n")
    chrColor <- "lightgrey"
  }
  # -- figure out how colors should be handled
  ptColorToUse <- NULL
  if(is.null(colorColumn) || is.na(colorColumn) || length(colorColumn) == 0){
    paf[,ptColorToUse := defaultColor]
  }else{
    if(!colorColumn %in% names(paf)){
      warning(sprintf("colorColumn %s is not one of the paf object names,
                      setting to NULL\n"))
      paf[,ptColorToUse := defaultColor]
    }else{
      cols <- paf[[colorColumn]]
      if(all(are_colors(cols))){
        paf[,ptColorToUse := cols]
      }else{
        colsu <- unique(cols[!is.na(cols)])
        colsn <- length(colsu)
        colsv <- colorPalette(colsn)
        if(randomizeColors)
          colsv <- sample(colsv, replace = FALSE)
        names(colsv) <- colsu
        paf[,ptColorToUse := colsv[cols]]
      }
    }
  }

  # -- calculate linear coordinates
  linCoords <- linearize_pafCoords(
    paf = paf, gapProp = gapProp)

  # -- pull plot limits
  xrng <- range(unlist(linCoords$chrCoords[["query"]][,c("xStart", "xEnd")]))
  yrng <- range(unlist(linCoords$chrCoords[["target"]][,c("xStart", "xEnd")]))

  # -- convert to vectors
  tp <- linCoords$paf
  qchr <- linCoords$chr$query
  tchr <- linCoords$chr$target

  qmin <- min(qchr$xStart)
  qmax <- max(qchr$xEnd)
  tmin <- min(tchr$xStart)
  tmax <- max(tchr$xEnd)

  # -- get mean positions per rows
  xqstart <- xtstart <- xtend <- xqend <- alen <- NULL
  tp[,`:=`(x = (xqstart + xqend) / 2,
           y = (xtstart + xtend) / 2)]

  tpSeg <- subset(tp, alen >= minWidth4segment)
  tpDot <- subset(tp, alen < minWidth4segment)

  # -- make the plot
  qmin <- qmax <- xStart <- tmax <- tmin <- x <- y <- ptColorToUse <- NULL
  plt <- ggplot()+
    geom_segment(
      data = tchr[-1,],
      aes(x = qmin, xend = qmax, y = xStart, yend = xStart),
      col = chrColor,
      linewidth = chrWidth)+
    geom_segment(
      data = qchr[-1,],
      aes(y = tmin, yend = tmax, x = xStart, xend = xStart),
      col = chrColor,
      linewidth = chrWidth)
  if(nrow(tpDot) > 0){
    plt <- plt +
      geom_point(
        data = tpDot,
        aes(x = x, y = y, color = ptColorToUse),
        shape = pointShape, ...)
  }
  xtstart <- xtend <- xqstart <- xqend <- ptColorToUse <- NULL
  if(nrow(tpSeg) > 0){
    plt <- plt +
      geom_segment(
        data = tpSeg,
        aes(y = xtstart, yend = xtend,
            x = xqstart, xend = xqend, color = ptColorToUse),
        linewidth = segLineWidth)
  }
  plt <- plt +
    scale_x_continuous(
      expand = c(0.01,0.01),
      limits = xrng,
      labels = qchr$chr,
      breaks = qchr$mid,
      name = "query chromosomes")+
    scale_y_continuous(
      expand = c(0.01,0.01),
      limits = yrng,
      labels = tchr$chr,
      breaks = tchr$mid,
      name = "target chromosomes")+
    scale_color_identity(guide = "none")+
    theme_dotplot()
  if(!is.null(addThemes))
    plt <- plt + addThemes
  return(plt)
}

#' @title Riparian plot from a paf
#' @description
#' \code{riparian_paf} two genome riparian plot
#' @rdname ggplot_paf
#' @importFrom ggplot2 theme element_rect element_blank element_text
#' geom_polygon geom_label annotate scale_fill_identity scale_alpha_identity
#' scale_linewidth_identity
#' @export
riparian_paf <- function(paf,
                         reorderByColumn = colorColumn,
                         colorColumn = "tname",
                         chrColor = "grey50",
                         gapProp = 0.01,
                         chrWidth = 2,
                         colorPalette = pal_ds,
                         randomizeColors = FALSE,
                         pointShape = ".",
                         defaultColor = add_alpha("white", .5),
                         addThemes = theme(panel.background = element_rect(fill = "black")),
                         chrBuffer = 0.01,
                         chrLabBackground = add_alpha("white", .8),
                         chrLabSize = 6,
                         chrLabColor = "black",
                         chrLabOffset = -chrBuffer,
                         genomeLabOffset = .075,
                         queryLab = "query genome",
                         targetLab = "target genome",
                         braidAlpha = 0.5,
                         ...){
  paf <- data.table(paf)
  # -- chromosome segment sizes
  chrWidth <- as.numeric(chrWidth)[1]
  if(is.na(chrWidth) || is.null(chrWidth))
    stop("chrWidth has to be a single numeric value\n")

  # -- chromosome gap sizes
  gapProp <- as.numeric(gapProp)[1]
  if(is.na(gapProp) || is.null(gapProp))
    stop("gapProp has to be a single numeric value\n")

  # -- color palette
  if(!is.function(colorPalette))
    stop("colorPalette must be a function that produces a vector of colors\n")
  test <- colorPalette(2)
  if(!all(are_colors(test)))
    stop("colorPalette must be a function that produces a vector of colors\n")
  if(length(test) != 2)
    stop("colorPalette function must produce a color vector of length n\n")

  # -- braid colors
  colorColumn <- as.character(colorColumn[1])
  if(!colorColumn %in% names(paf) || length(colorColumn) != 1)
    colorColumn <- as.character(NULL)

  # -- orderBy (so that colors get ordered correctly)
  reorderByColumn <- as.character(reorderByColumn[1])
  if(!reorderByColumn %in% names(paf) || length(reorderByColumn) != 1)
    reorderByColumn <- as.character(NULL)

  if(reorderByColumn %in% names(paf)){
    x <- tmpcol <- num <- char <- x <- NULL
    tmp <- data.table(x = unique(paf[[reorderByColumn]]))
    suppressWarnings(
      tmp[,`:=`(num = frank(as.numeric(x), ties.method = "average"),
                char = frank(as.character(x), ties.method = "average"))])
    setkey(tmp, num, char, x)
    paf[,tmpcol := factor(paf[[reorderByColumn]], levels = tmp$x)]
    setkey(paf, tmpcol)
    paf[,tmpcol := NULL]
    tmp <- NULL
  }

  # -- figure out how colors should be handled
  if(is.null(colorColumn) || is.na(colorColumn) || length(colorColumn) == 0){
    paf[,color2use := defaultColor]
  }else{
    if(!colorColumn %in% names(paf)){
      warning(sprintf("colorColumn %s is not one of the paf object names,
                      setting to NULL\n"))
      paf[,color2use := defaultColor]
    }else{
      cols <- paf[[colorColumn]]
      if(all(are_colors(cols))){
        paf[,color2use := cols]
      }else{
        colsu <- unique(cols[!is.na(cols)])
        colsn <- length(colsu)
        colsv <- colorPalette(colsn)
        if(randomizeColors)
          colsv <- sample(colsv, replace = FALSE)
        names(colsv) <- colsu
        paf[,color2use := colsv[cols]]
      }
    }
  }
  color2use <- NULL

  linCoords <- linearize_pafCoords(
    paf = paf, gapProp = gapProp)

  # -- pull plot limits
  xrng <- range(unlist(linCoords$chrCoords[["query"]][,c("xStart", "xEnd")]))
  yrng <- range(unlist(linCoords$chrCoords[["target"]][,c("xStart", "xEnd")]))

  # -- convert to vectors
  tp <- linCoords$paf
  tp[,index := 1:.N]
  qchr <- linCoords$chr$query
  tchr <- linCoords$chr$target

  qmin <- min(qchr$xStart)
  qmax <- max(qchr$xEnd)
  tmin <- min(tchr$xStart)
  tmax <- max(tchr$xEnd)


  # -- braids
  xqstart <- xqend <- xtstart <- xtend <- NULL
  braids <- rbindlist(lapply(1:nrow(tp), function(i)
    data.table(index = tp$index[i],
               with(tp[i,], calc_curvePolygon(
                 start1 = xqstart, end1 = xqend,
                 start2 = xtstart, end2 = xtend,
                 y2 = 1 - chrBuffer,
                 y1 = 0 + chrBuffer)))))
  braids <- merge(tp, braids, by = "index")

  # -- chromosome segments

  x <- y <- index <- color2use <- xStart <- xEnd <- mid <- chr <- NULL
  p <- ggplot()+
    geom_polygon(
      data = braids,
      aes(x = x, y = y, group = index, fill = color2use,
          alpha = braidAlpha, color = NA,
          linewidth = 0))+
    geom_segment(
      data = linCoords$chrCoords$query,
      aes(x = xStart, xend = xEnd, y = 0, yend = 0),
      linewidth = chrWidth, col = chrColor, lineend = "butt")+
    geom_segment(
      data = linCoords$chrCoords$target,
      aes(x = xStart, xend = xEnd, y = 1, yend = 1),
      linewidth = chrWidth, col = chrColor, lineend = "butt")+
    geom_label(
      data = linCoords$chrCoords$query,
      aes(x = mid, y = 0 - chrLabOffset, label = chr),
      fill = chrLabBackground, col = chrLabColor, size = chrLabSize/2.845276,
      vjust = 1)+
    geom_label(
      data = linCoords$chrCoords$target,
      aes(x = mid, y = 1 + chrLabOffset, label = chr),
      fill = chrLabBackground, col = chrLabColor, size = chrLabSize/2.845276,
      vjust = 0)+
    annotate(
      geom = "text",
      x = min(linCoords$chrCoords$target$xStart),
      y = 1+genomeLabOffset,
      label = sprintf("%s (%sMb)",
                      targetLab,
                      round(sum(linCoords$chrCoords$target$width) / 1e6, 2)),
      col = "white", size = 8/2.845276, hjust = 0, vjust = 1)+
    annotate(
      geom = "text",
      x = min(linCoords$chrCoords$query$xStart),
      y = 0-genomeLabOffset,
      label = sprintf("%s (%sMb)",
                      queryLab,
                      round(sum(linCoords$chrCoords$query$width) / 1e6, 2)),
      col = "white", size = 8/2.845276, hjust = 0, vjust = 0)+
    scale_color_identity(guide = "none") +
    scale_fill_identity(guide = "none") +
    scale_alpha_identity(guide = "none") +
    scale_linewidth_identity(guide = "none")+
    theme_pafriparian()+
    scale_x_continuous(
      expand = c(0.01, 0.01),
      name = "chromosomes, scaled by physical position")+
    scale_y_continuous(
      expand = c(0.01, 0.01),
      name = NULL)+
    # ggtitle("DEEPSPACE pairwise riparian plot\n")+
    addThemes

  return(p)
}

#' @title theme_dotplot
#' @description
#' \code{theme_dotplot} theme_dotplot
#' @rdname ggplot_paf
#' @importFrom ggplot2 theme element_rect element_blank element_text
#' @export
theme_dotplot <- function(...){
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(
          angle= 90, family = "Helvetica", size = 5, vjust = .5, hjust = 1),
        axis.text.y = element_text(
          angle= 0, family = "Helvetica", size = 5, vjust = .5, hjust = 1),
        axis.title = element_text(family = "Helvetica", size = 6),
        axis.ticks = element_blank(),
        plot.title = element_text(family = "Helvetica", size = 7))
}

#' @title theme_riparian
#' @description
#' \code{theme_dotplot} theme_riparian
#' @rdname ggplot_paf
#' @importFrom ggplot2 theme element_rect element_blank element_text
#' @export
theme_pafriparian <- function(...){
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(color = "black", size = 8),
        title = element_text(color = "black", size = 8))
}
