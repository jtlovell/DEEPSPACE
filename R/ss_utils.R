#' @title Convenience functions related to StringSet objects
#' @description
#' \code{ss_utils} Simplify common scripting goals on StringSet objects
#' @name ss_utils
#'
#' @param x StringSet object
#' @param kmerSS StringSet object containing sequences to kmerize
#' @param min.gapwidth passed to GenomicRanges::reduce specifying the minimum
#' gap across which ranges should be reduced
#' @param dropDuplicates logical, should duplicate entries be dropped?
#' @param kmerSize integer, width of the interval
#' @param faFile file.path to the fasta file
#' @param minChrLen integer, minimum chromosome length to retain
#' @param regexs character vector of regular expressions to parse the seq names.
#' @param nCores integer, the number of parallel processes to run
#'

#' @title kmerize stringset
#' @description
#' \code{kmerize_ss} Split sequences into 1bp offset sliding windows
#' @rdname ss_utils
#' @importFrom GenomicRanges slidingWindows GRanges
#' @importFrom Biostrings seqinfo
#' @export
kmerize_ss <- function(x,
                       kmerSize = 10,
                       dropDuplicates = TRUE){
  # -- check parameters
  if(!is(x, "StringSet"))
    stop("x must be a DNA or AA StringSet\n")
  kmerSize <- as.integer(kmerSize)
  if(length(kmerSize) != 1 || is.null(kmerSize) || is.na(kmerSize))
    stop("kmerSize must be a single integer value\n")
  if(!is.logical(dropDuplicates) || length(dropDuplicates) != 1)
    stop("dropDuplicates must be a single logical value")

  # -- add names to the object
  if(is.null(names(x)))
    names(x) <- 1:length(x)

  # -- pull seqinfo for sliding windows
  gr <- GRanges(seqinfo(x))
  sw <- unlist(slidingWindows(gr, width = kmerSize, step = 1))

  # -- chop of ss into windows
  sso <- x[sw]

  # -- drop any partial windows
  sso <- sso[width(sso) == kmerSize]

  # -- drop anything identical (if specified)
  if(dropDuplicates)
    sso <- sso[!duplicated(sso)]

  # -- return
  names(sso) <- NULL
  return(sso)
}

#' @title pull out blocks of exact kmer hits
#' @description
#' \code{pull_kmerBlocks} exact kmer extraction and summarization
#' @rdname ss_utils
#' @importFrom GenomicRanges slidingWindows GRanges
#' @importFrom Biostrings seqinfo
#' @importFrom parallel mclapply
#' @importFrom BiocGenerics do.call
#' @export
pull_kmerBlocks <- function(kmerSS,
                            x,
                            min.gapwidth = 100,
                            nCores = 1){
  # -- check parameters
  if(!is(x, "StringSet"))
    stop("x must be a DNA or AA StringSet\n")
  width <- as.integer(width)
  if(length(width) != 1 || is.null(width) || is.na(width))
    stop("width must be a single integer value\n")

  # -- pdict only works on same length sequences, kmerize if necessary
  if(length(unique(width(kmerSS))) > 1){
    minSeqLen <- min(width(kmerSS))
    kmerSS <- kmerize_ss(x = kmerSS, kmerSize = minSeqLen)
    kmerSS <- kmerSS[!duplicated(kmerSS)]
  }

  # -- renamne and summarize input kmers
  names(kmerSS) <- 1:length(kmerSS)
  kmerLen <- width(kmerSS)[1]
  pdictKmer <- PDict(kmerSS)

  # -- loop through each target sequence and find exact matches
  kmerLocs <- parallel::mclapply(
    names(x),
    mc.preschedule = FALSE,
    mc.cores = nCores,
    function(i){
      mdict <- makeGRangesFromDataFrame(data.table(
        chr = i,
        as.data.frame(unlist(matchPDict(
          pdict = pdictKmer,
          subject = x[[i]])))[,1:2]))

      out <- join_gappedRanges(
        x = mdict,
        min.gapwidth = min.gapwidth)
      return(out)
    })
  names(kmerLocs) <- NULL
  suppressWarnings(out <- do.call(c, kmerLocs))
  return(out)
}

#' @title import fasta as DNAStringSet
#' @description
#' \code{import_dnass} read in and optionally parse a fasta file
#' @rdname ss_utils
#' @importFrom Biostrings readDNAStringSet width
#' @export
import_dnass <- function(faFile, minChrLen = 0, regexs = NULL){
  ss <- readDNAStringSet(faFile)

  ss <- ss[width(ss) >= minChrLen]

  regexs <- as.character(regexs)
  regexs <- regexs[!is.na(regexs)]
  if(length(regexs) > 0){
    for(i in regexs)
      names(ss) <- gsub(i,"", names(ss))
  }
  return(ss)
}
