#' @title Convenience functions related to StringSet objects
#' @description
#' \code{paf_utils} Simplify common scripting goals on StringSet objects
#' @name paf_utils
#'
#' @param x paf file or connection
#' @param paf data.table or data.frame coercible to a paf file.
#' @param convert1index logical, should the paf be converted into 1-indexed
#' formatting (to then be converted into GRanges etc)?
#' @param keepExtraCols logical, should columns beyond the first 12 be kept?
#' @param extraColNames character string vector with the length matching the
#' number of extra columns in the input file (only used if keepExtraCols = TRUE)
#' @param overwrite logical, if x exists, should it be overwritten?
#' @param keepColnames logical, should column names be written?
#' @param verbose logical, should updates be printed to the console?
#' \cr
#' If called, \code{paf_utils} returns its own arguments.
#'

#' @title Read the first 12 columns of a .paf file
#' @description
#' \code{fread_paf} Use data.table::fread import the 12 columns of a .paf file
#' @rdname paf_utils
#' @import data.table
#' @importFrom stats complete.cases
#' @export
fread_paf <- function(x,
                      verbose = FALSE,
                      convert1index = TRUE,
                      keepExtraCols = FALSE,
                      extraColNames = NULL){

  # -- check arguments
  verbose <- verbose[1]
  if(!is.logical(verbose))
    verbose <- TRUE

  if(!file.exists(x))
    stop(sprintf("paf file: %s does not exist\n", x))
  rl <- length(strsplit(readLines(x, 1), "\t")[[1]])
  if(rl < 12)
    stop("data in paf file has < 12 columns\n")
  keepExtraCols <- keepExtraCols[1]
  if(!is.logical(keepExtraCols))
    keepExtraCols <- FALSE
  convert1index <- convert1index[1]
  if(!is.logical(convert1index))
    convert1index <- TRUE

  # -- provide the names and classes
  pafNames <- c(
    "qname", "qlen", "qstart", "qend", "strand",
    "tname", "tlen", "tstart", "tend", "nmatch",
    "alen", "mapq")
  pafClasses <- c(
    "character", "numeric", "numeric", "numeric", "character",
    "character", "numeric", "numeric", "numeric", "numeric",
    "numeric", "numeric")

  # -- determine the number of additional columns
  if(rl > 12){
    pafClasses <- c(pafClasses, rep("character", rl - 12))
  }

  # -- read in data if we are keeping extra columns
  if(verbose)
    cat(sprintf("\t%s", basename(x)))
  if(keepExtraCols){
    if(rl > 12){
      if(!is.null(extraColNames)){
        if(length(extraColNames) != rl-12)
          stop("extraColNames length must match the number of extra columns\n")
        pafNames <- c(pafNames, extraColNames)
      }else{
        pafNames <- c(pafNames, sprintf("col%s", 12:rl))
      }
    }
    p <- fread(
      file = x,
      col.names = pafNames,
      colClasses = pafClasses,
      header = F,
      showProgress = F,
      sep = "\t",
      fill = T)
  }else{

    # -- read in only first 12 columns if not keeping extra columns
    p <- fread(
      file = x,
      select = 1:12,
      col.names = pafNames,
      colClasses = pafClasses,
      header = F,
      showProgress = F,
      sep = "\t",
      fill = T)
  }

  # -- convert to 1-index if we desired
  qstart <- tstart <- NULL
  if(convert1index)
    p[,`:=`(qstart = as.numeric(qstart) + 1,
            tstart = as.numeric(tstart) + 1)]
  if(verbose)
    cat(sprintf("... read in %s paf entries\n", nrow(x)))
  if(any(!complete.cases(p)))
    warning("some data could not be coerced to the correct data class - there may be a problem with your paf file\n")

  return(p)
}

#' @title Read the first 12 columns of a .paf file
#' @description
#' \code{fread_paf} Use data.table::fread import the 12 columns of a .paf file
#' @rdname paf_utils
#' @import data.table
#' @export
fwrite_paf <- function(paf,
                       x,
                       keepColnames = FALSE,
                       keepExtraCols = FALSE,
                       overwrite = FALSE){
  # -- check the paf object
  if(!is.data.frame(paf))
    stop("paf object must be a data.frame or data.table\n")
  if(ncol(paf) < 12)
    stop("paf object appears to be malformed\n")
  if(nrow(paf) == 0)
    stop("paf object has no data\n")

  # -- check output file path
  x <- as.character(x)
  if(is.null(x) || is.na(x))
    stop("x is not a valid path to write paf file\n")
  x <- path.expand(x)
  if(!file.exists(dirname(x)))
    stop("parent directory of x does not exist\n")

  # -- check other parameters
  keepColnames <- match.arg(keepColnames, choices = c(TRUE, FALSE))
  keepExtraCols <- match.arg(keepExtraCols, choices = c(TRUE, FALSE))
  overwrite <- match.arg(overwrite, choices = c(TRUE, FALSE))

  # -- check overwrite
  if(file.exists(x) && !overwrite){
    warning("x file exists and !overwrite, so not writing\n")
  }else{

    # -- choose the columns
    if(!keepExtraCols){
      out <- data.table(paf[,1:12])
    }else{
      out <- data.table(paf)
    }

    # -- write it
    fwrite(
      out,
      file = x,
      quote = FALSE,
      sep = "\t",
      col.names = keepColnames,
      showProgress = FALSE,
      verbose = FALSE)
  }
}
