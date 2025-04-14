#' @title Generic internal functions used by DEEPSPACE
#' @description
#' \code{riparian_paf} Convenience functions for DEEPSPACE, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name riparian_paf
#'
#' @param pafFiles character vector coercible to a file.path pointing to the
#' paf files. These must be named as '$QUERYGENOMEID'__s__'$TARGETGENOMEID'.paf.
#' @param genomeIDs character vector specifying the top-to-bottom genome order.
#' Must specify a daisy chain where all edges in the vector have a corresponding
#' .paf.
#' @param refGenome one of the genomeIDs that is the reference to which all
#' alignments are phased. If absent, sameChrOnly should be TRUE
#' @param orderBySynteny logical, should chromosomes be ordered by synteny
#' (TRUE, default) or by the name of chromosomes (FALSE)
#' @param sameChrOnly logical, should the search only be conducted on alignments
#' where the query and target sequence names are identical? This speeds things
#' up, but should only be used in very similar genomes where inter-chromosomal
#' SVs are known to not exist.
#' @param round2join integer specifying the amount of fuzziness should be
#' present in the joining of alignments against the reference and within the
#' daisy chain. Larger numbers means less precision, but likely more complete
#' phasing.
#' @param labelTheseGenomes character string specifying which of the genomeIDs
#' should have chromosomes labeled in the plot
#' @param braidColors vector of colors used in the palette for the braid colors.
#' If a single value, all blocks are not phased against the reference genome.
#' @param chrWidth numeric, providing the thinkness of the chromosomes
#' @param sameScale logical, do all genomes have the same scale (TRUE) or all
#' have the same width (FALSE)
#' @param braidAlpha numeric (0-1) specifying the opacity of the braids
#' @param braidOffset distance between the chromosome midpoint and the y-start
#' of each daisychain braid. If the same value as chrWidth, the braid starts
#' exactly at the edge of the chromosome polygon
#' @param highlightInversions either NULL or a color specifying what color the
#' inverted sequence should be
#' @param braidOutline color (or NA) specifying the color of the braid outlines
#' @param chrOutlineColor color (or NA) specifying the color of the chr outlines
#' @param chrOutlineWidth numeric, giving the thickness of the chr outlines
#' @param braidOutlineWidth numeric, giving the thickness of the braid outlines
#' @param gapProp numeric (0-1) giving the proportion of the x axis that should
#' be gaps between chromosomes
#' @param howRound numeric (0-Inf) specifying how rouded the chromosome
#' polygons should be.
#'
#' \cr
#' If called, \code{riparian_paf} returns its own arguments.
#'
#'
#' @title riparian_paf
#' @description
#' \code{riparian_paf} riparian_paf
#' @rdname riparian_paf
#' @importFrom GenomicRanges start end width slidingWindows
#' @import ggplot2
#' @import data.table
#' @export
riparian_paf <- function(pafFiles,
                         labelTheseGenomes = NULL,
                         refGenome = NULL,
                         genomeIDs,
                         sameChrOnly = FALSE,
                         round2join = 10e3,
                         braidColors = c("#C4645C", "#F5915C", "#FFC765", "#FCF8D5", "#BEF3F9", "#66B8FF", "#6666FF", "#9C63E1", "#F4BDFF"),
                         chrWidth = .05,
                         sameScale = TRUE,
                         orderyBySynteny = TRUE,
                         braidAlpha = 0.5,
                         braidOffset = 0.06,
                         highlightInversions = FALSE,
                         braidOutline = NA,
                         chrOutlineColor = "grey",
                         chrOutlineWidth = 0.1,
                         braidOutlineWidth = 0,
                         gapProp = 0.05,
                         howRound = .2){
  ##############################################################################
  # 1. Get the paf files in order
  # -- get pafFiles in order including genomeIDs
  if(any(!file.exists(pafFiles)))
    stop("some pafFiles do not exist\n")
  md <- data.table(path = pafFiles[!duplicated(pafFiles)])
  md[,c("query", "target") := tstrsplit(gsub(
    ".paf$|.paf.gz$|.syn.paf|.syn.paf.gz", "", basename(path)),
    "_windows__vs__")]

  # -- check that the correct paf files exist
  if(!sameChrOnly){
    if(is.null(refGenome))
      stop("if !sameChrOnly, must supply a reference genome to phase blocks\n")
    rh <- genomeIDs[genomeIDs != refGenome]
    rhmd <- subset(md, target == refGenome)
    if(!all(rh %in% rhmd$query))
      stop("not all genomeIDs have paf files mapped to the reference\n")
  }

  # -- check that the daisy chain exists
  ng <- length(genomeIDs)
  dcmd <- data.table(
    query = c(genomeIDs[-1], genomeIDs[-ng]),
    target = c(genomeIDs[-ng], genomeIDs[-1]),
    index = rep(1:(ng-1), 2))
  dcmd <- merge(dcmd, md, by = c("query", "target"))
  setkey(dcmd, index)
  if(uniqueN(dcmd$index) != (ng - 1))
    stop("something is wrong: number of paf files does not match expected daisy chain\n")

  if(is.null(labelTheseGenomes)){
    labelTheseGenomes <- genomeIDs
  }else{
    labelTheseGenomes <- genomeIDs[genomeIDs %in% labelTheseGenomes]
    if(length(labelTheseGenomes) == 0){
      warning("none of labelTheseGenomes are in genomeIDs ... labeling all\n")
      labelTheseGenomes <- genomeIDs
    }
  }

  ##############################################################################
  # 2. Phase the hits by reference genome chromosomes
  # -- read in the ref hits
  pafNames <- c(
    "qname", "qlen", "qstart", "qend", "strand",
    "tname", "tlen", "tstart", "tend", "nmatch",
    "alen", "mapq", "windID", "blkID")

  if(!sameChrOnly){
    refh <- rbindlist(lapply(1:nrow(rhmd), function(i)
      data.table(
        fread(rhmd$path[i], col.names = pafNames),
        query = rhmd$query[i],
        target = rhmd$target[i])))
    refh[, qname := as.character(qname)]
    refh[, tname := as.character(tname)]

    # -- round query hit positions for a fuzzy join
    refh[,rnd := round_toInteger((qstart + qend) / 2, round2join)]

    # -- calculate the reference chromosomes for each rounded query position
    refd <- refh[,list(refChr = unique(tname)), by = c("query", "qname", "rnd")]
  }

  # -- for each member of the daisy chain ...
  dcblks <- rbindlist(lapply(1:nrow(dcmd), function(i){
    x <- dcmd[i,]
    # -- read in the daisy chain hits
    pafi <- fread(x$path, col.names = pafNames)
    pafi[, qname := as.character(qname)]
    pafi[, tname := as.character(tname)]

    if(sameChrOnly)
      pafi <- subset(pafi, tname == qname)

    pafi[,`:=`(query = x$query, target = x$target)]

    # -- round to the fuzzy join
    if(!sameChrOnly){
      pafi[,rnd := round_toInteger((qstart + qend) / 2, round2join)]
      pafi <- subset(pafi, !duplicated(paste(qname, rnd, windID, blkID)))
      # -- merge with the reference mappings
      pafm <- merge(pafi, refd, by = c("query", "qname", "rnd"))
      # -- rename the block ID to contain the reference chr names
      pafm[,blkID := sprintf("%s_refChrXX%s", blkID, refChr)]
      # -- calculate the block coordinates
    }else{
      pafm <- data.table(pafi)
    }

    bl <- calc_pafBlkCoords(pafm)

    if(!is.null(bl)){
      # -- add in the y positions and reference chr ID
      bl[,`:=`(qy = which(genomeIDs == x$query) * -1,
               ty = which(genomeIDs == x$target) * -1,
               query = x$query,
               target = x$target)]
      if(sameChrOnly){
        bl[,refChr := tname]
      }else{
        bl[,refChr := gsub(".*refChrXX", "", blkID)]
      }
      return(bl)
    }
  }))

  ##############################################################################
  # 3. integrate user-defined parameters
  suppressWarnings(dcblks[,tnamen := as.numeric(gsub('\\D+','', tname))])
  setorder(dcblks, target, query, tnamen, tname, qname, na.last = TRUE)

  # -- unique reference chromosomes
  if(sameChrOnly){
    refChrs <- unique(dcblks$tname)
    refChrs <- refChrs[order(refChrs)]
    pal <- colorRampPalette(braidColors)
    cls <- pal(length(refChrs))
    names(cls) <- refChrs
  }else{
    suppressWarnings(refh[,tnamen := as.numeric(gsub('\\D+','', tname))])
    refChrs <- subset(refh, !duplicated(tname))
    setorder(refChrs, tnamen, tname, na.last = TRUE)
    refChrs <- split(refChrs$tname, refChrs$target)
    pal <- colorRampPalette(braidColors)
    cls <- pal(length(refChrs[[1]]))
    names(cls) <- refChrs[[1]]
  }

  # -- colors
  dcblks[,cols := cls[refChr]]

  # -- chromosome order
  if(orderyBySynteny & !sameChrOnly){
    suppressWarnings(refh[,tnamen := as.numeric(gsub('\\D+','', tname))])
    setorder(refh, tnamen, tname, na.last = TRUE)
    spl <- split(refh, by = "query")
    chroList <- sapply(spl, USE.NAMES = T, simplify = F, function(x)
      order_bySynteny(paf = x, chrordTarget = refChrs[[1]]))
    chroList <- c(refChrs,  chroList)
  }else{

    suppressWarnings(dcblks[,qnamen := as.numeric(gsub('\\D+','', qname))])

    cl <- with(dcblks, data.table(chr = c(qname, tname), genome = c(query, target)))
    suppressWarnings(cl[,chrn := as.numeric(gsub('\\D+','', chr))])
    cl <- subset(cl, !duplicated(cl))
    setorder(cl, genome, chrn, chr, na.last = TRUE)
    chroList <- split(cl$chr, cl$genome)
  }

  # highlightInversions
  if(are_colors(highlightInversions)[1]){
    dcblks$cols[dcblks$strand == "-"] <- highlightInversions
  }

  ##############################################################################
  # 4. make the plot
  if(sameScale){
    pl <- riparian_blksSameScale(
      blks = data.table(dcblks),
      chrList = chroList,
      colorColumn = "cols",
      chrWidth = chrWidth,
      braidOffset = braidOffset,
      braidAlpha = braidAlpha,
      braidOutline = braidOutline,
      braidOutlineWidth = braidOutlineWidth,
      gapProp = gapProp,
      labelTheseGenomes = labelTheseGenomes,
      howRound = howRound)
  }else{
    pl <- riparian_blksSameWidth(
      blks = data.table(dcblks),
      chrList = chroList,
      colorColumn = "cols",
      chrWidth = chrWidth,
      braidAlpha = braidAlpha,
      braidOffset = braidOffset,
      braidOutline = braidOutline,
      braidOutlineWidth = braidOutlineWidth,
      gapProp = gapProp,
      labelTheseGenomes = labelTheseGenomes,
      howRound = howRound)
  }
  return(pl)
}


#' @title riparian_blksSameWidth
#' @description
#' \code{riparian_blksSameWidth} riparian_blksSameWidth
#' @rdname riparian_paf
#' @importFrom GenomicRanges start end width slidingWindows
#' @import ggplot2
#' @import data.table
#' @export
riparian_blksSameWidth <- function(blks,
                                   chrList,
                                   labelTheseGenomes,
                                   colorColumn = "color",
                                   braidAlpha = 0.5,
                                   braidOutline = NA,
                                   braidOutlineWidth = 0,
                                   braidPalette = gs_colors,
                                   chrOutlineColor = "grey",
                                   chrOutlineWidth = 0.1,
                                   gapProp = 0.05,
                                   howRound = .2,
                                   braidOffset = .12,
                                   chrWidth = .1){
  ##############################################################################
  # 1. get colors in order
  # -- if the color column isn't in the data table, use query name

  blks[,index := 1:.N]
  if(!colorColumn %in% names(blks))
    colorColumn <- "qname"

  # -- make a new column "color" to store the colors
  blks[,color := blks[[colorColumn]]]

  # -- if the color column is not a vector of colors, convert to colors
  if(any(!are_colors(blks$color))){
    pal <- braidPalette(n = uniqueN(blks$color))
    names(pal) <- unique(blks$color)
    blks[,color := pal[color]]
  }

  ##############################################################################
  # 2. Get blocks data in order
  # -- calculate query and target genome sizes
  tmp <- with(blks, data.table(
    genome = c(query, target), name = c(qname, tname), len = c(qlen, tlen)))
  tmp <- subset(tmp, !duplicated(paste(genome, name)))
  clens <- with(tmp, tapply(len, genome, sum))
  blks[,`:=`(qsize = clens[query], tsize = clens[target])]

  # -- get the length of the largest genome
  mxlen <- with(blks, max(c(qsize, tsize), na.rm = T))

  # -- determine scale as % of the difference between the largest genome
  blks[,`:=`(qscale = mxlen / qsize,
             tscale = mxlen / tsize)]
  blks[,`:=`(qstart = qscale * qstart,
             tstart = tscale * tstart,
             qlen = qscale * qlen,
             tlen = tscale * tlen,
             qend = qscale * qend,
             tend = tscale * tend)]
  ##############################################################################
  # 3. Calculate riparian coordinates
  # -- split the blocks files by the query and target
  spl <- split(blks, by = c("query", "target"))

  # -- for each pairwise combination, build the chromosome and braid polygons
  ripdat <- rbindlist(lapply(spl, function(x){
    tmp <- transform_paf2riparian(
      paf = x,
      chrordQuery = chrList[[x$query[1]]],
      chrordTarget = chrList[[x$target[1]]],
      yQuery = x$qy[1],
      yTarget = x$ty[1],
      sameXrange = FALSE,
      colorColumn = "color",
      xoffsetQuery = 0,
      xoffsetTarget = 0,
      gapProp = gapProp,
      howRound = howRound,
      braidOffset = braidOffset,
      chrWidth = chrWidth)
    q <- x$query[1]; t = x$target[1]
    tmp$queryChrs[,`:=`(type = "chrs", genome = q, index = 1)]
    tmp$targetChrs[,`:=`(type = "chrs", genome = t, index = 2)]
    tmp$braids[,`:=`(chrID = NA, type = "braids", genome = paste(q, t))]
    tmp <- rbindlist(tmp, use.names = T, fill = TRUE)
    return(tmp)
  }))

  # -- split out chromosome polygons, only keeping one set per genome
  chrdat <- subset(ripdat, type == "chrs")
  setkey(chrdat, genome, index)
  tk <- with(subset(chrdat, !duplicated(paste(genome, index))),
             paste(genome, index))
  chrdat <- subset(chrdat, paste(genome, index) %in% tk)

  # -- make chr polygons of unlabeled genomes half the height
  chrdat <<- chrdat
  spln <- split(subset(chrdat, !genome %in% labelTheseGenomes), by = "genome")
  spln <- rbindlist(lapply(spln, function(splg){
    yr <- range(splg$y)
    yadj <- diff(yr)/4
    splg[,y := scale_between(y, min = yr[1] + yadj, max = yr[2] - yadj)]
    return(splg)
  }))
  chrdat <- rbind(subset(chrdat, genome %in% labelTheseGenomes), spln)

  # -- split out the braid data
  braiddat <- subset(ripdat, type == "braids")

  # -- get the chromosome names and position of labels
  chrnames <- chrdat[,list(x = mean(range(x)), y = mean(range(y))),
                     by = c("genome", "chrID")]
  chrnames <- subset(chrnames, genome %in% labelTheseGenomes)

  # -- calcuate the genome names and sizes
  gnames <- with(blks, data.table(
    genome = c(query, target), bp = c(qsize, tsize), y = c(qy, ty)))
  gnames <- subset(gnames, !duplicated(genome))

  ##############################################################################
  # 4. make the plot object and return
  p <- ggplot()+
    # -- braid polygon
    geom_polygon(
      aes(x = x, y = y, fill = as.factor(index), group = index),
      alpha = braidAlpha,
      data = braiddat, color = braidOutline,
      linewidth = braidOutlineWidth)+
    scale_fill_manual(
      values = blks$color,
      guide = "none")+
    geom_polygon(aes(x = x, y = y, group = paste(chrID, genome)),
                 data = chrdat, fill = "white",
                 col = chrOutlineColor, linewidth = chrOutlineWidth)+
    geom_text(aes(x = x, y = y, label = chrID), data = chrnames, size = 5*0.36,
              hjust = "middle", vjust = "center")+
    scale_y_continuous(
      labels = sprintf("%s\n%sMbp", gnames$genome, round(gnames$bp/1e6, 1)),
      breaks = gnames$y,
      expand = c(0.01, 0.01))+
    scale_x_continuous(
      expand = c(0.01, 0.01))+
    theme(panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black", size = 7),
          axis.title.x = element_text(color = "black", size = 8),
          title = element_text(color = "black", size = 8))+
    labs(x = "chromosomes scaled to the same x-range across genomes",
         title = "DEEPSPACE synteny map")
  return(p)
}


#' @title riparian_blksSameScale
#' @description
#' \code{riparian_blksSameScale} riparian_blksSameScale
#' @rdname riparian_paf
#' @importFrom GenomicRanges start end width slidingWindows
#' @import ggplot2
#' @import data.table
#' @export
riparian_blksSameScale <- function(blks,
                                   chrList,
                                   labelTheseGenomes,
                                   chrOutlineColor = "grey",
                                   chrOutlineWidth = 0.1,
                                   colorColumn = "color",
                                   braidAlpha = 0.5,
                                   braidOutline = NA,
                                   braidOutlineWidth = 0,
                                   braidPalette = gs_colors,
                                   gapProp = 0.05,
                                   howRound = .2,
                                   braidOffset = .12,
                                   chrWidth = .1){
  ##############################################################################
  # 1. get colors in order
  # -- if the color column isn't in the data table, use query name
  blks[,index := 1:.N]

  if(!colorColumn %in% names(blks))
    colorColumn <- "qname"

  # -- make a new column "color" to store the colors
  blks[,color := blks[[colorColumn]]]

  # -- if the color column is not a vector of colors, convert to colors
  if(any(!are_colors(blks$color))){
    pal <- braidPalette(n = uniqueN(blks$color))
    names(pal) <- unique(blks$color)
    blks[,color := pal[color]]
  }

  ##############################################################################
  # 2. Get blocks data in order
  # -- calculate query and target genome sizes
  tmp <- with(blks, data.table(
    genome = c(query, target), name = c(qname, tname), len = c(qlen, tlen)))
  tmp <- subset(tmp, !duplicated(paste(genome, name)))
  clens <- with(tmp, tapply(len, genome, sum))
  blks[,`:=`(qsize = clens[query], tsize = clens[target])]

  # -- get the length of the largest genome
  mxlen <- with(blks, max(c(qsize, tsize)))

  # -- determine offset as 50% of the difference between the largest genome
  blks[,`:=`(qoffset = (mxlen - qsize) / 2,
             toffset = (mxlen - tsize) / 2)]

  ##############################################################################
  # 3. Calculate riparian coordinates
  # -- split the blocks files by the query and target
  spl <- split(blks, by = c("query", "target"))

  # -- for each pairwise combination, build the chromosome and braid polygons
  ripdat <- rbindlist(lapply(spl, function(x){
    tmp <- transform_paf2riparian(
      paf = x,
      chrordQuery = chrList[[x$query[1]]],
      chrordTarget = chrList[[x$target[1]]],
      yQuery = x$qy[1],
      yTarget = x$ty[1],
      sameXrange = FALSE,
      colorColumn = "color",
      xoffsetQuery = x$qoffset[1],
      xoffsetTarget = x$toffset[1],
      gapProp = gapProp,
      howRound = howRound,
      braidOffset = braidOffset,
      chrWidth = chrWidth)
    q <- x$query[1]; t = x$target[1]
    tmp$queryChrs[,`:=`(type = "chrs", genome = q, index = 1)]
    tmp$targetChrs[,`:=`(type = "chrs", genome = t, index = 2)]
    tmp$braids[,`:=`(chrID = NA, type = "braids", genome = paste(q, t))]
    tmp <- rbindlist(tmp, use.names = T, fill = TRUE)
    return(tmp)
  }))

  # -- split out chromosome polygons, only keeping one set per genome
  chrdat <- subset(ripdat, type == "chrs")
  setkey(chrdat, genome, index)
  tk <- with(subset(chrdat, !duplicated(paste(genome, index))),
             paste(genome, index))
  chrdat <- subset(chrdat, paste(genome, index) %in% tk)

  chrdat <<- chrdat
  spln <- split(subset(chrdat, !genome %in% labelTheseGenomes), by = "genome")
  spln <- rbindlist(lapply(spln, function(splg){
    yr <- range(splg$y)
    yadj <- diff(yr)/4
    splg[,y := scale_between(y, min = yr[1] + yadj, max = yr[2] - yadj)]
    return(splg)
  }))
  chrdat <- rbind(subset(chrdat, genome %in% labelTheseGenomes), spln)

  # -- split out the braid data
  braiddat <- subset(ripdat, type == "braids")

  # -- get the chromosome names and position of labels
  chrnames <- chrdat[,list(x = mean(range(x)), y = mean(range(y))),
                     by = c("genome", "chrID")]
  chrnames <- subset(chrnames, genome %in% labelTheseGenomes)

  # -- calcuate the genome names and sizes
  gnames <- with(blks, data.table(
    genome = c(query, target), bp = c(qsize, tsize), y = c(qy, ty)))
  gnames <- subset(gnames, !duplicated(genome))

  ##############################################################################
  # 4. make the plot object and return
  p <- ggplot()+
    # -- braid polygon
    geom_polygon(
      aes(x = x, y = y, fill = as.factor(index), group = index),
      alpha = braidAlpha,
      data = braiddat,
      color = braidOutline,
      linewidth = braidOutlineWidth)+
    scale_fill_manual(
      values = blks$color,
      guide = "none")+
    geom_polygon(aes(x = x, y = y, group = paste(chrID, genome)),
                 data = chrdat, fill = "white",
                 col = chrOutlineColor, linewidth = chrOutlineWidth)+
    geom_text(aes(x = x, y = y, label = chrID), data = chrnames, size = 5*0.36,
              hjust = "middle", vjust = "center")+
    scale_y_continuous(
      labels = sprintf("%s\n%sMbp", gnames$genome, round(gnames$bp/1e6, 1)),
      breaks = gnames$y,
      expand = c(0.01, 0.01))+
    scale_x_continuous(
      expand = c(0.01, 0.01))+
    theme(panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black", size = 7),
          axis.title.x = element_text(color = "black", size = 8),
          title = element_text(color = "black", size = 8))+
    labs(x = "chromosomes scaled by physical size", title = "DEEPSPACE synteny map")
  return(p)
}

