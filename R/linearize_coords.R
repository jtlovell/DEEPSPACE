#' @title Functions to convert genomic coordinates to linear
#' @description
#' \code{linearize_coords} Take within-chromosome coordinates and convert to
#' global positions, optionally including offsets and gaps
#' @name linearize_coords
#'
#' @param x granges, stringset or named vector of integers
#' @param paf paf-like data.table or data.frame
#' @param qnameOrder input order of qnames
#' @param tnameOrder input order of tnames
#' @param gapProp numeric, proportion of the genome to be gaps
#' @param qGapSize query gap size
#' @param tGapSize target gap size
#' @param qXoffset query offset
#' @param tXoffset target offset
#' @param gapWidth gap width
#' @param stripChrName regex to parse chromosome names
#' @param xoffset left offset of data
#' \cr
#' If called, \code{linearize_coords} returns its own arguments.
#'

#' @title join adjacent ranges
#' @description
#' \code{linearize_pafCoords} Reduces ranges, returning a granges with proportion
#' of sequence and total sequence covered.
#' @rdname linearize_coords
#' @importFrom GenomicRanges sort reduce width
#' @export
linearize_pafCoords <- function(paf,
                                qnameOrder = NULL,
                                tnameOrder = NULL,
                                gapProp = 0,
                                qGapSize = NULL,
                                tGapSize = NULL,
                                qXoffset = NULL,
                                tXoffset = NULL){

  # -- get chromosome orders
  num <- x <- char <- NULL
  if(is.null(qnameOrder)){
    tmp <- data.table(x = unique(paf$qname))
    suppressWarnings(
      tmp[,`:=`(num = frank(as.numeric(x), ties.method = "average"),
              char = frank(as.character(x), ties.method = "average"))])
    setkey(tmp, num, char, x)
    qnameOrder <- tmp$x

    tmp <- data.table(x = unique(paf$tname))
    suppressWarnings(
      tmp[,`:=`(num = frank(as.numeric(x), ties.method = "average"),
                char = frank(as.character(x), ties.method = "average"))])
    setkey(tmp, num, char, x)
    tnameOrder <- tmp$x
  }
  qnord <- c(paste("query", qnameOrder), paste("target", tnameOrder))
  qordv <- c(1:(length(qnameOrder) + length(tnameOrder)))
  names(qordv) <- qnord

  # -- build table of chromosome lengths
  qname <- qlen <- tname <- tlen <- NULL
  paf <- data.table(paf)
  chrLens <- with(paf, rbind(
    data.table(genome = "query", chr = qname, len = qlen),
    data.table(genome = "target", chr = tname, len = tlen)))

  # -- add in chromosome order
  genome <- chr <- ord <- NULL
  chrLens <- subset(chrLens, paste(genome, chr) %in% qnord)
  chrLens[,ord := qordv[paste(genome, chr)]]
  chrLens <- subset(chrLens, !duplicated(chrLens))
  setkey(chrLens, ord)

  # -- calculate gap size if not given
  len <- nGaps <- seqWidth <- NULL
  if(is.null(qGapSize) || is.null(tGapSize)){
    chrTot <- chrLens[,list(seqWidth = sum(len),
                            nGaps = .N - 1),
                      by = "genome"]
    mxWidth <- max(chrTot$seqWidth)
    gapSize <- mxWidth * gapProp
    gapSize <- c(gapSize, gapSize)
  }else{
    gapSize <- c(qGapSize, tGapSize)
  }
  names(gapSize) <- c("query", "target")

  # -- calculate left x offset if not given
  totWidth <- leftOffset <- NULL
  if(is.null(qXoffset) || is.null(tXoffset)){
    # -- get genome-specific gap sizes
    chrTot[,totWidth := (nGaps * gapSize) + seqWidth]
    chrTot[,leftOffset := (max(totWidth) - totWidth) / 2]
    xoffset <- chrTot$leftOffset; names(xoffset) <- chrTot$genome
  }else{
    xoffset <- c(qXoffset, tXoffset)
    names(xoffset) <- c("query", "target")
  }

  # -- make list of chrlen vectors
  chrlv <- as.integer(chrLens$len); names(chrlv) <- chrLens$chr
  chrlv <- split(chrlv, chrLens$genome)

  # -- calculate linear chromosome starts
  linChr <- sapply(names(chrlv), USE.NAMES = TRUE, simplify = F, function(i)
    linearize_chrCoords(
      chrlv[[i]],
      gapWidth = gapSize[i],
      xoffset = xoffset[i]))

  linSt <- lapply(linChr, function(x){
    y <- x$xStart; names(y) <- x$chr
    return(y)
  })

  # -- add to input data
  qname <- qstart <- qend <- tstart <- tname <- tend <- NULL
  paf[,`:=`(xqstart = qstart + linSt$query[qname],
           xqend = qend + linSt$query[qname],
           xtstart = tstart + linSt$target[tname],
           xtend = tend + linSt$target[tname])]
  # paf[,index := 1:.N]
  return(list(paf = paf, chrCoords = linChr))
}

#' @title Calculate the linear positions of chromosomes
#' @description
#' \code{linearize_chrCoords} Produce a data.table with the linear positions for
#' plotting
#' @rdname linearize_coords
#' @import GenomicRanges
#' @importFrom Biostrings fasta.seqlengths
#' @importFrom GenomeInfoDb seqlengths
#' @export
linearize_chrCoords <- function(x,
                                gapWidth = 0,
                                stripChrName = "",
                                xoffset = 0){

  # -- get the input data formats into a vector
  if(is(x, "DNAStringSet")){
    si <- width(x)
    names(si) <- names(x)
  }else{
    if(is(x, "Seqinfo")){
      si <- seqlengths(x)
    }else{
      if(is(x, "integer") && !is.null(names(x))){
        si <- x
      }else{
        if(is(x, "character") && length(x) == 1){
          if(file.exists(x)){
            si <- fasta.seqlengths(x)
          }else{
            stop("cannot find file: ", x, "\n")
          }
        }else{
          stop("x must be a DNAStringSet, Seqinfo, or path to a fasta file\n")
        }
      }
    }
  }

  # -- convert to data table
  dt <- data.table(chr = names(si), seqWidth = as.numeric(si))

  # -- get the starts and ends
  gp <- xStart <- seqWidth <- xEnd <- mid <- NULL
  dt[,gp := gapWidth]
  dt[,xStart := c(1, cumsum(seqWidth[-.N]) + cumsum(gp[-.N]) + 1)]
  dt[,`:=`(xEnd = xStart + seqWidth, gp = NULL,
           mid = (xStart + seqWidth + xStart)/2)]

  # -- rename chrs if necessary
  chr <- NULL
  if(stripChrName != ""){
    dt[,chr := gsub(stripChrName, "", chr)]
  }

  if(xoffset != 0){
    dt[,`:=`(xStart = xStart + xoffset,
             xEnd = xEnd + xoffset,
             mid = mid + xoffset)]
  }
  return(dt)
}

