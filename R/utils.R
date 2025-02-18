#' @title Generic internal functions used by DEEPSPACE
#' @description
#' \code{utils} Convenience functions for DEEPSPACE, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name utils
#'
#' @param col color or column
#' @param x depends on function
#' @param to end of result
#' @param dt input data.table
#' @param cols columns to consider
#' @param order what to order by
#' @param rankWithin logical
#' @param alpha transparency
#' @param n number to keep
#' @param faFile file.path to a fasta file
#' @param y depends on function
#' @param filepath path to file
#' @param min min value
#' @param max max value
#' @param scale1toMean logical
#' @param ... passed on
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
    "DEEPSPACE v0.2.1: genome-wide comparative genomics in R ... this is a
    development version and may be (definitely is) buggy\n",
    indent = 0, exdent = 8), collapse = "\n"))
}

#' @title are these colors?
#' @description
#' \code{are_colors} check if a vector of numbers or character strings can be
#' coerced into a vector of colors
#' @rdname utils
#' @export
are_colors <- function(col) {
  sapply(col, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}

#' @title DEEPSPACE colors
#' @description
#' \code{cols_ds} the nine colors that make up the deepspace palette
#' @rdname utils
#' @export
cols_ds <- function(){
  return(c("#C4645C", "#F5915C", "#FFC765",
           "#FCF8D5", "#BEF3F9", "#66B8FF",
           "#6666FF", "#9C63E1", "#F4BDFF"))
}

#' @title DEEPSPACE palette
#' @description
#' \code{pal_ds} colorRampPalette on DEEPSPACE colors
#' @rdname utils
#' @export
pal_ds <- function(x){
  colorRampPalette(colors = cols_ds())(x)
}

#' @title round_toInteger
#' @description
#' \code{round_toInteger} round_toInteger
#' @rdname utils
#' @export
round_toInteger <- function(x, to){
  as.numeric(round(x/to, 0) * to)
}


#' @title frank_multicolumn
#' @description
#' \code{frank_multicolumn} frank_multicolumn
#' @rdname utils
#' @import data.table
#' @export
frank_multicolumn <- function(dt, cols, order, rankWithin = NULL){
  x <- dupRank <- mapRank <- NULL
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
  index <- NULL
  setkey(x, index)
  return(x$mapRank)
}

#' @title add transparency
#' @description
#' \code{add_alpha} add transparency to a color
#' @rdname utils
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

#' @title pal_deepspace
#' @description
#' \code{pal_deepspace} pal_deepspace
#' @rdname utils
#' @importFrom grDevices colorRampPalette
#' @export
pal_heatmap <- function(n){
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

#' @title scale_between
#' @description
#' \code{scale_between} scale_between
#' @rdname utils
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

