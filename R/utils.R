#' @title Generic internal functions used by DEEPSPACE
#' @description
#' \code{utils} Convenience functions for DEEPSPACE, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name utils
#'
#' @param fafiles named vector of file.paths to assembly fasta files
#' @param outdir file.path to the output directory
#' @param minChrSize numeric, minimum chromosome size to be included
#' @param stripChrname regex passed to gsub to strip text of chromosome names
#' @param windowSize integer, sliding window width
#' @param stepSize integer, step between sliding window start positions
#' @param verbose logical, should updates be printed to the console?
#' \cr
#' If called, \code{utils} returns its own arguments.
#'
#'
#' @title startup messages
#' @description
#' \code{.onAttach} startup messages
#' @rdname utils
#' @export
.onAttach <- function(...) {
  packageStartupMessage(paste(strwrap(
    "DEEPSPACE v0.1.0: genome-wide comparative genomics in R ... this is a
    development version and may be (definitely is) buggy\n",
    indent = 0, exdent = 8), collapse = "\n"))
}

#' @title prep_fafiles
#' @description
#' \code{prep_fafiles} prep_fafiles
#' @rdname utils
#' @importFrom GenomicRanges start end width slidingWindows
#' @import Biostrings
#' @importFrom Rsamtools scanFaIndex
#' @export
prep_fafiles <- function(fafiles,
                         outdir,
                         minChrSize = 0,
                         stripChrname = "\\s.*",
                         nWindows = 1e6,
                         windowSize = 1e3,
                         stepSize = 1e3,
                         verbose = TRUE){

  if(is.null(names(fafiles))){
    if(verbose)
      cat("faFiles vector is not named. Pulling genome IDs from the file name\n")
    names(fafiles) <- gsub(
      ".fa$|.fasta$|.fa.gz$|.fasta.gz", "", basename(fafiles))
  }
  gids <- names(fafiles)

  if(any(duplicated(gids)))
    stop("some fafile names are duplicated - not allowed\n")

  if(any(!file.exists(fafiles)))
    stop("some fafiles do not exist - provide a valid path to fastas\n")

  outfa <- file.path(outdir, sprintf("%s.fa", gids))
  names(outfa) <- gids
  outwi <- file.path(outdir, sprintf("%s.window.fa", gids))
  names(outwi) <- gids

  for(i in gids){
    if(verbose)
      cat(sprintf("%s ... ", i))
    ss <- readDNAStringSet(fafiles[i])

    if(verbose)
      cat(sprintf("%s Mbp / %s scaffolds",
                  round(sum(width(ss))/1e6, 2),
                  length(ss)))
    names(ss) <- gsub(stripChrname, "", names(ss))
    if(any(duplicated(names(ss))))
      stop("chromosome name stripping resulted in non-unique names; re-specify stripChrname\n")

    if(minChrSize > 0){
      ss <- ss[width(ss) > minChrSize]
      cat(sprintf(" (%s Mbp / %s scaffolds > %s Mb)",
                  round(sum(width(ss))/1e6, 2),
                  length(ss),
                  minChrSize/1e6))
      if(sum(width(ss)) == 0)
        stop(sprintf("no sequences < %s; use a smaller minChrSize", minChrSize))
    }

    writeXStringSet(ss, filepath = outfa[i])

    faif <- sprintf("%s.fai", outfa[i])
    if(file.exists(faif))
      unlink(faif)
    chrgr <- scanFaIndex(outfa[i])

    if(is.null(nWindows)){
      stepSize <- sum(width(ss)) / nWindows
    }


    wind <- slidingWindows(chrgr, width = windowSize, step = stepSize)


    wind <- BiocGenerics::do.call(c, wind)
    suppressMessages(windss <- ss[wind])

    if(verbose)
      cat(sprintf(" ... %s windows\n",
                  length(windss)))
    names(windss) <- sprintf(
      "%s_XstartX%s_XendX%s", seqnames(wind), start(wind), end(wind))
    writeXStringSet(windss, filepath = outwi[i])
  }

  # -- pull sizes
  szs <- sapply(outfa, USE.NAMES = T, simplify = T, function(x)
    sum(width(scanFaIndex(x))))
  return(list(genome = gids, ref = outfa, wind = outwi, totalbases = szs))
}

#' @title Read the first 12 columns of a .paf file
#' @description
#' \code{fread_paf} Use data.table::fread import the 12 columns of a .paf file
#' @rdname utils
#' @import data.table
#' @export
fread_paf <- function(x, verbose = FALSE){

  if(!file.exists(x))
    stop(sprintf("paf file: %s does not exist\n", x))
  rl <- length(strsplit(readLines(x, 1), "\t")[[1]])
  if(rl < 12)
    stop("data in paf file has < 12 columns\n")

  pafNames <- c(
    "qname", "qlen", "qstart", "qend", "strand",
    "tname", "tlen", "tstart", "tend", "nmatch",
    "alen", "mapq")
  pafClasses <- c(
    "character", "integer", "integer", "integer", "character",
    "character", "integer", "integer", "integer", "integer",
    "integer", "integer")

  if(rl > 12)
    pafClasses <- c(pafClasses, rep("character", rl - 12))

  if(verbose)
    cat(sprintf("\t%s", basename(x)))
  p <- fread(
    file = x,
    select = 1:12,
    col.names = pafNames,
    colClasses = pafClasses,
    header = F,
    showProgress = F,
    sep = "\t",
    fill = T)
  if(verbose)
    cat(sprintf("... read in %s paf entries\n", nrow(x)))
  if(any(!complete.cases(p)))
    warning("some data could not be coerced to the correct data class - there may be a problem with your paf file\n")

  return(p)
}

#' @title convert_pos2ord
#' @description
#' \code{convert_pos2ord} convert_pos2ord
#' @rdname utils
#' @importFrom data.table frankv
#' @export
convert_pos2ord <- function(start, end, group, to = 0){
  x <- (start + end)/2
  to <- as.integer(to)
  if(to > 0){
    x <- round_toInteger(x, to)
  }
  ord <- frankv(
    list(group, x),
    ties.method = "dense")
  return(ord)
}

#' @title round_toInteger
#' @description
#' \code{round_toInteger} round_toInteger
#' @rdname utils
#' @export
round_toInteger <- function(x, to){
  round(x/to, 0) * to
}

#' @title parse_windPaf
#' @description
#' \code{parse_windPaf} parse_windPaf
#' @rdname utils
#' @import data.table
#' @export
parse_windPaf <- function(paf,
                          queryFastaFile,
                          stripChrname){

  querySeqlens <- fasta.seqlengths(queryFastaFile)
  names(querySeqlens) <- gsub(stripChrname, "", names(querySeqlens))
  paf[,windID := as.integer(as.factor(qname))]

  paf[,windChr := gsub("_XstartX.*", "", qname)]
  paf[,winds := as.numeric(gsub(".*_XstartX(.*?)_XendX.*", "\\1", qname))]

  paf[,`:=`(
    qname = windChr,
    qstart = as.integer(qstart + winds),
    qend = as.integer(qend + winds),
    windChr = NULL, winds = NULL)]

  paf[,qlen := querySeqlens[qname]]
  setorder(paf, qname, qstart, -mapq, -nmatch)
  return(paf)
}

#' @title calc_linStart
#' @description
#' \code{calc_linStart} calc_linStart
#' @rdname utils
#' @import data.table
#' @import Biostrings
#' @export
calc_linStart <- function(fastaFile = NULL, seqLens = NULL, gapSize = 0){
  if(is.null(seqLens)){
    seqLens <- fasta.seqlengths(fastaFile)
  }

  if(length(seqLens) < 1)
    stop("there must be at least one sequence!\n")

  if(length(seqLens) == 1){
    x <- 0
    names(x) <- names(seqLens)
  }else{
    g <- data.table(
      name = names(seqLens),
      len = as.numeric(seqLens),
      gapLeft = c(0, rep(gapSize, length(seqLens) - 1)))

    g[,seqLeft := c(0, cumsum(len)[-.N])]
    g[,gapLeft := cumsum(gapLeft)]
    g[,st := gapLeft + seqLeft]

    x <- g$st; names(x) <- g$name
  }

  if(any(duplicated(names(x))))
    warning("some duplicated names ... something might be wrong\n")
  return(x)
}

#' @title parse_pafWindows
#' @description
#' \code{parse_pafWindows} parse_pafWindows
#' @rdname utils
#' @import data.table
#' @export
parse_pafWindows <- function(paf,
                             querySeqlens,
                             addMaprank = TRUE,
                             addOrd = TRUE,
                             flagRepeat = TRUE,
                             addnp = TRUE,
                             pidThresh = .9,
                             ntopThresh = 4,
                             npidThresh = ntopThresh * 2,
                             round2ord = max(paf$qend - paf$qstart)){
  p <- data.table(paf)
  p[,windChr := gsub("_XstartX.*", "", qname)]
  p[,winds := as.numeric(gsub(".*_XstartX(.*?)_XendX.*", "\\1", qname))]
  windsize <- with(p, max(qend - qstart))
  p[,`:=`(
    qstart = as.integer(qstart + winds),
    qend = as.integer(qend + winds))]

  if(addMaprank){
    p[,mapp := nmatch / windsize]
    p[,maprank := frank_multicolumn(
      dt = p,
      cols = c("qname", "mapq", "mapp", "nmatch"),
      order = c(1, -1, -1, -1),
      rankWithin = "qname")]

    if(flagRepeat){
      p[,isRepeat := sum(maprank == 1) >= ntopThresh |
            sum(mapp >= pidThresh) >= npidThresh, by = "qname"]
    }
  }

  if(addOrd){
    p[,tmp := round_toInteger(x = (qstart + qend) / 2, to = round2ord)]
    p[,qord := frank_multicolumn(
      dt = p,
      cols = c("windChr", "tmp"),
      order = c(1, 1))]
    p[,tmp := round_toInteger(x = (tstart + tend) / 2, to = round2ord)]
    p[,tord := frank_multicolumn(
      dt = p,
      cols = c("tname", "tmp"),
      order = c(1, 1))]
    p[,tmp := NULL]
  }

  p[,`:=`(qname = windChr, windChr = NULL)]
  if(addMaprank){
    setkey(p, qname, qstart, qend, maprank)
  }else{
    setkey(p, qname, qstart, qend)
  }

  p[,qlen := querySeqlens[qname]]
  return(p)
}

#' @title frank_multicolumn
#' @description
#' \code{frank_multicolumn} frank_multicolumn
#' @rdname utils
#' @import data.table
#' @export
frank_multicolumn <- function(dt, cols, order, rankWithin = NULL){
  x <- NULL
  x <- data.table(dt)
  if(!all(cols %in% names(x)))
    stop("all cols must be in the dt\n")
  x[,index := 1:.N]
  setorderv(x = x, cols = cols, order = order)
  x[,dupRank := duplicated(x[,cols, with = F])]
  if(!is.null(rankWithin)){
    if(!rankWithin %in% cols)
      stop("rankWithin must be a column in the dt matching one of cols\n")
    x[,mapRank := cumsum(!dupRank), by = rankWithin]
  }else{
    x[,mapRank := cumsum(!dupRank)]
  }
  setkey(x, index)
  return(x$mapRank)
}


#' @title pal_deepspace
#' @description
#' \code{pal_deepspace} pal_deepspace
#' @rdname utils
#' @importFrom grDevices colorRampPalette
#' @export
pal_deepspace <- function(n){
  cols <- c(
    "navy", "dodgerblue4", "darkcyan",
    "mediumspringgreen","greenyellow", "yellow",
    "white")
  pal <- colorRampPalette(cols)
  return(pal(n))
}


#' @title check_isDNAFasta
#' @description
#' \code{check_isDNAFasta} check_isDNAFasta
#' @rdname utils
#' @import data.table
#' @importFrom Biostrings readDNAStringSet width extractAt writeXStringSet
#' @export
check_isDNAFasta <- function(faFile){
  if(!file.exists(faFile))
    stop(sprintf("faFile %s does not exist\n", faFile))
  l <- readLines(faFile, n = 2)
  hasHead <- nchar(l[[1]]) > 1 & substr(l[[1]],1,1) == ">"
  hasDNA <- all(sapply(strsplit(substr(l[[2]], 1, 25), "")[[1]],
                       function(x) toupper(x) %in% c("N", "A", "T", "C", "G")))
  hasDNA <- hasDNA & nchar(l[[2]] > 1)
  if(!hasDNA)
    stop(sprintf("something is wrong with the fasta file `%s` ... doesn't have DNA sequences\n", faFile))
  if(!hasHead)
    stop(sprintf("something is wrong with the fasta file `%s` ... doesn't have a proper header\n", faFile))
  return(faFile)
}

#' @title window_ss
#' @description
#' \code{window_ss} window_ss
#' @rdname utils
#' @import data.table
#' @importFrom Biostrings readDNAStringSet width extractAt writeXStringSet
#' @export
window_ss <- function(ss, windowSize, windowStep, returnGranges = FALSE){
  dt <- data.table(
    seqnames = names(ss),
    start = 1,
    end = width(ss))
  dt <- dt[,list(
    start = seq(from = min(start),
                to = max(end),
                by = windowStep),
    maxp = max(end)),
    by = "seqnames"]
  dt[,end := start + (windowSize - 1)]
  dt$end[dt$end > dt$maxp] <- dt$maxp[dt$end > dt$maxp]
  dt[,windWidth := end - start]
  dt <- subset(dt, windWidth > (windowSize / 2))
  gr <- as(as.data.frame(dt[,c("seqnames", "start", "end")]), "GRanges")
  if(returnGranges){
    return(gr)
  }else{
    out <- ss[gr]

    names(out) <- sprintf(
      "%s_XstartX%s_XendX%s",
      seqnames(gr),
      gsub(" ", "", format(start(gr), scientific = FALSE)),
      gsub(" ", "", format(end(gr), scientific = FALSE)))
    return(out)
  }
}



#' @title set_params
#' @description
#' \code{set_params} set_params
#' @rdname utils
#' @import GenomicRanges
#' @import Biostrings
#' @importFrom Rsamtools scanFaIndex
#' @export
set_params <- function(faFiles,
                       wd,
                       minChrSize = 0,
                       stripChrname = "\\s.*",
                       nWindows = 1e6,
                       windowSize = 1e3,
                       stepSize = 1e3,
                       genomeIDs = NULL,
                       refGenome = NULL,
                       minimap2.call = "minimap2",
                       MCScanX_h.call = "MCScanX_h",
                       minimap2.args = "-B4 -O6,26 -U50,500 -g100 --no-long-join -r250,500 --rmq=no -A1 -E2,1 -s200 --frag=yes -F1500 -f0 -p0.5 -k19 -w19 -N10",
                       nCores = 1,
                       blkSize = 5,
                       verbose = T){

  #
  find_refWind <- function(genomeIDs,
                           genomeSizes){
    rmd <- data.table(ref = genomeIDs[1], wind = genomeIDs[-1])
    if(length(genomeIDs) > 2){
      dcmd <- data.table(
        g1 = genomeIDs[-c(1, length(genomeIDs))],
        g2 = genomeIDs[-c(1:2)])
      dcmd[,`:=`(s1 = genomeSizes[g1],
                 s2 = genomeSizes[g2])]
      dcmd[,`:=`(
        ref = ifelse(s1 > s2, g1, g2),
        wind = ifelse(s1 > s2, g2, g1))]
      rmd <- rbind(rmd, dcmd[,c("ref", "wind")])
    }
    return(rmd)
  }


  # -- check fasta files
  famd <- prep_fafiles(
    fafiles = fafiles,
    outdir = wd,
    minChrSize = minChrSize,
    stripChrname = stripChrname,
    nWindows = nWindows,
    windowSize = windowSize,
    stepSize = stepSize,
    verbose = verbose)

  # -- check general parameters
  genomeIDs <- check_genomeIDs(faFiles, genomeIDs, refGenome)
  refGenome <- check_refGenome(refGenome, genomeIDs)
  nonRefGenomes <- genomeIDs[genomeIDs != refGenome]
  genomeIDs <- c(refGenome, nonRefGenomes)

  # -- pull out the reference and window genomes
  mps <- with(famd, find_refWind(
    genomeIDs = genome,
    genomeSizes = totalbases))

  # -- add reference and window file.paths
  mps[,`:=`(
    refFile = famd$ref[ref],
    windFile = famd$wind[wind],
    outFile = file.path(wd, sprintf("%s_vs_%s.paf", ref, wind)))]

  # -- run minimap2
  for(i in 1:nrow(mps)){
    ers <- system2(
      command = minimap2.call,
      args = c(minimap2.args,
               sprintf("-o %s", mps$outFile[i]),
               mps$refFile[i],
               mps$windFile[i]),
               stdout = TRUE,
               stderr = TRUE)
    cat(ers, sep = "\n")
  }

}

#' @title interp_xy
#' @description
#' \code{interp_xy} interp_xy
#' @rdname utils
#' @importFrom stats approx
#' @export
interp_xy <- function(x, y){

  interp_approx <- function(x, y){
    y <- as.numeric(y)
    if(all(is.na(y)))
      stop("must have some non-NA values in y")
    x <- as.numeric(x)
    if(any(is.na(x)))
      stop("x cannot contain NAs")
    if(any(is.na(y))){
      newy <- approx(
        x = x,
        y = y,
        ties = mean,
        xout = x[is.na(y)])$y
      y[is.na(y)] <- newy
    }
    return(y)
  }

  if(is.na(y[1]))
    y[1] <- min(y, na.rm = T) - .49
  if(is.na(y[length(y)]))
    y[length(y)] <- max(y, na.rm = TRUE) + .49

  out <- interp_approx(x = x, y = y)
  return(out)
}

#' @title pull_proxHits
#' @description
#' \code{pull_proxHits} pull_proxHits
#' @rdname utils
#' @import data.table
#' @importFrom dbscan frNN
#' @export
pull_proxHits <- function(x, y, isAnchor, radius){
  ##############################################################################
  # -- parameter checking
  if(length(x) != length(y))
    stop("x and y must be vectors of the same length\n")
  if(sum(isAnchor) == 0)
    return(rep(FALSE, length(x)))
  if(sum(!isAnchor) == 0)
    return(rep(TRUE, length(x)))
  if(length(x) < 2)
    stop("x and y must have lengths > 1\n")
  x <- as.numeric(x)
  if(any(is.na(x)))
    stop("x values must all be numeric or integer\n")
  y <- as.numeric(y)
  if(any(is.na(y)))
    stop("y values must all be numeric or integer\n")

  isAnchor <- as.logical(isAnchor)
  if(any(is.na(isAnchor)))
    stop("values given to isAnchor (logical vector of isAnchor observations are anchors) must all be coercible to logical\n")

  radius <- as.numeric(radius[1])
  if(is.na(radius) || is.null(radius))
    stop("radius must be a single numeric value > 0\n")
  radius <- as.numeric(radius[1])
  if(radius <= 0)
    stop("radius must be a single numeric value > 0\n")

  ##############################################################################
  # -- function to do it
  # -- get fixed radius nearest neighbors
  if(all(isAnchor) || all(!isAnchor)){
    return(isAnchor)
  }else{
    xy <- data.frame(x = x, y = y, isAnchor = isAnchor, inBuffer = isAnchor)
    fn <- frNN(
      x = subset(xy, isAnchor)[c("x", "y")],
      query = subset(xy, !isAnchor)[c("x", "y")],
      eps = radius)

    # -- subset to positions with any nearest neighbors in radius
    hasAnch <- fn$id[sapply(fn$id, length) > 0]
    wh <- as.integer(names(hasAnch))
    # -- set these observations to true and return
    if(length(wh) > 0)
      xy$inBuffer[wh] <- TRUE
    return(xy$inBuffer)
  }
}

#' @title check_MCScanXhInstall
#' @description
#' \code{check_MCScanXhInstall} check_MCScanXhInstall
#' @rdname utils
#' @export
check_MCScanXhInstall <- function(filepath){
  filepath <- path.expand(filepath[1])
  if(!file.exists(filepath))
    stop(sprintf("could not find a valid MCScanX_h install at %s", filepath))
  chk <- suppressWarnings(grepl(
    "prefix_fn",
    system2(filepath, "-h", stdout = TRUE, stderr = FALSE)[1]))

  if(!chk)
    stop("could not find valid MCScanX_h install\n")
  return(filepath)
}

#' @title check_minimap2Install
#' @description
#' \code{check_minimap2Install} check_minimap2Install
#' @rdname utils
#' @export
check_minimap2Install <- function(filepath){
  filepath <- path.expand(filepath[1])
  chk <- suppressWarnings(grepl(
    "Usage",
    system2(filepath, "--help", stdout = TRUE, stderr = FALSE)[1]))

  if(!chk)
    stop("could not find valid minimap2 install\n")
  return(filepath)
}


