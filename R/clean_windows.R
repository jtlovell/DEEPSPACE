#' @title Fast window-based synteny analysis
#' @description
#' \code{clean_windows} The main functionality for DEEPSPACE. Checks, parses,
#' and aligns genomic fasta files. Then subsets the resulting .paf files to
#' hits in synteny and produces dotplot and riparian plot visualizations
#' @name clean_windows
#'
#' @param faFiles character vector coercible to file.paths pointing to the
#' genome assembly fasta files to use in the run. If faFiles is a named vector,
#' genomeIDs are taken as the names.
#' @param wd character string pointing to the directory where deepspace output
#' should be stored
#' @param genomeIDs if faFiles are not named, supply genome names here
#' @param refGenome the genomeID to be used as the reference. Ignored if
#' sameChrOnly = TRUE.
#' @param preset more simple specification of parameters to run fast or search
#' for repeats and if the genomes are really close or distantly related.
#' Options are: "close", "dist", "fast", "reps", "close fast", "close reps",
#' "dist fast", "dist reps".
#' @param sameChrOnly logical, should the search only be conducted on alignments
#' where the query and target sequence names are identical? This speeds things
#' up, but should only be used in very similar genomes where inter-chromosomal
#' SVs are known to not exist.
#' @param stripChrname character string coercible to a regular expression,
#' passed to gsub to strip off text from chromosome names
#' @param minChrLen integer, the minimum sequence length. All sequences in the
#' fasta files shorter than this are ignored.
#' @param stepSize integer, the base-pair distance between adjacent window
#' mid points. Alternative specification to nWindows
#' @param windowSize integer, the base-pair size of each window
#' @param nWindows integer, the number of windows in the query genome. If
#' specified, supercedes stepSize
#' @param filterTopReps integer, passed to minimap2 -f. Setting to 0 can
#' significantly slow operations but will allow visualization of all repeats
#' @param maxTopHits integer, passed to minimap2 -N
#' @param minSecRat numeric, passed to minimap2 -p
#' @param mismatchPenalty numeric, passed to minimap2 -B
#' @param kmerSize integer, passed to minimap2 -k
#' @param kmerStep integer, passed to minimap2 -w
#' @param nCores integer, passed to minimap2 -t
#' @param basicParams character string of additional parameters passed to
#' minimap2.
#' @param mm2verbose logical, should minimap2 output be printed to the console?
#' @param MCScanX_hCall character string coercible to a valid file path pointing
#' to the MCScanX_h executable
#' @param minimap2call character string coercible to a valid file path pointing
#' to the minimap2 program. If installed with conda or in the path, use
#' "minimap2".
#' @param makePlots integer, specifying how many dotplots should be made. O
#' means no plots. 3 means make one for each step of synteny curation
#' @param blkSize the minimum number of windowed hits to be called a block.
#' @param tmpDir character string coercible to a file.path where temporary
#' files should be written. Cleaned out after run finishes unless this
#' directory already existed
#' @param ppi integer, pixels per inch in the dotplot
#' @param minPlotDim numeric, narrowest dimension of the dotplot in inches
#' @param maxBestMaps integer, the maximum number of best-mapping hits to keep
#' @param topNhits2keep integer, the number of top hits to keep
#' @param buffer numeric, position-rank order buffer around anchor hits to
#' search for additional anchors.
#' @param maxDistBtwHitsInBlk numeric, the maximum distance between two
#' adjacent hits in the same block
#' @param keepEverything logical, should final synteny searches be conducted on
#' all hits (including masked repeats etc)?
#' @param targetNhitsQuantile quantile to choose repetitive regions in the
#' target
#' @param overwrite logical, should all output be overwritten?
#' @param overwriteInput logical, should windowed fastas be overwritten?
#' @param overwriteMinimap2 logical, should minimap2 output be overwritten?
#' @param overwriteSynpaf logical, should synteny output be overwritten?

#' @param orderyBySynteny logical, should chromosomes be ordered by synteny
#' (true) or by chromosome name order (false)?
#' @param sameScale logical, should all genomes have the same scale (TRUE) or
#' the same plotted size on the x-axis (FALSE)
#' @param braidOffset numeric, the amount (in inches) that the braids should
#' start along the vertical axis from the midpoint of a genome
#' @param highlightInversions logical, NULL or a single color, specifying
#' what to do with inversions. If a color, all inversions are colored as that.
#' @param howRound numeric, (0-inf) where lower values are more square polygons
#' for the chromosomes
#' @param gapProp numeric (0-1) the proportion of the x-axis size that should
#' be made up of gaps between chromosomes
#' @param braidAlpha numeric (0-1) specifying the opacity of the braids
#' @param verbose logical, should updates be printed to the console?
#' @param ... additional arguments passed to riparian_paf
#'
#' @details This is an annotation-free version of GENESPACE. Feed it the genome
#' assemblies and it produces the standard GENESPACE graphical outputs
#'  (dotplots, riparian). DEEPSPACE does not produce pan-gene sets. There are
#'  four main steps to the pipeline:
#'
#'  1. Parse the fasta files. The following alignments are run: (a) all
#'  non-reference genomes windowed against the reference genome ... by default,
#'  the first genome provided is the "reference", but this can be changed by the
#'   user; (b) a 'daisy chain' of all genomes in the order given in genomeIDs so
#'    that for genomes A, B, C, alignments are done A --> B and B --> C where
#'    the larger genomes are windowed and aligned to the smaller genome
#'
#' 2. Run minimap2. The windowed genome is always the 'query' and the other
#' genome is the target. The user can specify all possible parameters to
#' minimap2 with a series of parameters in `clean_windows`. See the help file
#' for more details.
#'
#' 3. Parse the minimap2 results (.paf) into syntenic blocks. This generally
#' follows the GENESPACE synteny pipeline, treating minimap2 windowed alignments
#'  as it would the peptide-peptide alignments.

#'
#' @importFrom GenomicRanges start end width slidingWindows
#' @import ggplot2
#' @import data.table
#' @export
clean_windows <- function(faFiles,
                          wd,
                          preset = NULL,
                          genomeIDs = names(faFiles),
                          refGenome = genomeIDs[1],
                          sameChrOnly = FALSE,
                          stripChrname = "\\s.*",
                          minChrLen = 10e6,
                          stepSize = NULL,
                          windowSize = 1e3,
                          nWindows = 500e3,

                          filterTopReps = 1e-04,
                          maxTopHits = 50,
                          minSecRat = 0.9,
                          mismatchPenalty = 4,
                          kmerSize = 25,
                          kmerStep = 20,
                          nCores = 1,
                          basicParams = "-A1 -U50,500 --no-long-join -r0,0 --rmq=no -E2,1 -n1 -m0 --frag=yes",
                          mm2verbose = FALSE,

                          MCScanX_hCall = "MCScanX_h",
                          minimap2call = "minimap2",

                          makePlots = 3,
                          blkSize = 10,
                          tmpDir = file.path(wd, "DEEPSPACEtmpwd"),
                          ppi = 96,
                          minPlotDim = 8,
                          maxBestMaps = 5,
                          topNhits2keep = 4,
                          buffer = 5,
                          maxDistBtwHitsInBlk = 500e3,
                          keepEverything = FALSE,
                          targetNhitsQuantile = 0.95,
                          minMapq = 0,

                          overwrite = FALSE,
                          overwriteInput = FALSE,
                          overwriteMinimap2 = FALSE,
                          overwriteSynpaf = FALSE,

                          orderyBySynteny = TRUE,
                          sameScale = TRUE,
                          braidOffset = .05,
                          highlightInversions = FALSE,
                          howRound = .5,
                          gapProp = .1,
                          braidAlpha = .75,
                          ripHeightPerGenome = .9,
                          ripWidth = 8,

                          verbose = TRUE,
                          ...){

  if(overwrite){
    overwriteInput = TRUE
    overwriteMinimap2 = TRUE
    overwriteSynpaf = TRUE
  }

  if(!is.null(preset)){
    prs <- c("close", "dist", "fast", "reps",
             "close fast", "close reps", "dist fast", "dist reps")
    preset <- match.arg(preset, choices = prs)

    if(preset == "close fast"){
      topNhits2keep <- 1
      maxTopHits <- 4
    }
    if(grepl("close", preset)){
      windowSize <- 500
    }

    if(grepl("dist", preset)){
      blkSize <- 5
      windowSize <- 5000
    }

    if(grepl("fast", preset)){
      filterTopReps <- 0.01
      kmerSize <- 19
      kmerStep <- 19
      nWindows <- 200e3
      minMapq <- 12
    }

    if(grepl("reps", preset)){
      keepEverything <- TRUE
      filterTopReps <- 0
      maxTopHits <- Inf
      minSecRat <- 0.5
    }
  }

  wd <- as.character(wd)
  if(is.null(wd) || is.na(wd))
    stop("outDir must be a character string coercible to a file.path\n")
  wd <- path.expand(wd)

  if(!dir.exists(wd)){
    if(!dir.exists(dirname(wd)))
      stop("parent directory of outDir does not exist. specify a valid output directory\n")
    cat(sprintf("output directory %s does not exists, creating it...\n", wd))
    dir.create(wd)
  }

  # 0. Check paramters
  # -- dependency installs
  MCScanX_hCall <- check_MCScanXhInstall(MCScanX_hCall)
  minimap2call <- check_minimap2Install(minimap2call)

  # -- basic run parameters
  minChrLen <- as.integer(minChrLen[1])
  if(!is.integer(minChrLen))
    stop("minChrLen needs to be coercible to an integer")

  windowSize <- as.integer(windowSize[1])
  if(!is.integer(windowSize))
    stop("windowSize needs to be coercible to an integer")

  # -- n windows
  nWindows <- as.integer(nWindows[1])
  if(is.null(minChrLen) | is.na(minChrLen))
    stop("minChrLen needs to be coercible to an integer")

  # -- minimap2 parameters are check in minimap2_windows
  # -- riparian plotting parameters are checked in riparian_paf
  windowSize <- as.integer(windowSize[1])
  if(is.null(windowSize) | is.na(windowSize))
    stop("windowSize needs to be coercible to an integer")

  if(is.null(names(faFiles)) && length(genomeIDs) != length(faFiles))
    stop("faFiles must be a named vector\n")
  if(any(!file.exists(faFiles)))
    stop("one or more of the faFiles does not exist\n")

  # input faFiles
  faFiles <- sapply(faFiles, check_isDNAFasta)

  ##############################################################################
  # 1. Setup input file parameters

  if(verbose)
    cat("1. Initializing run ... ")

  # -- 1.1 check genome IDs
  genomeIDs <- names(faFiles)

  # -- 1.2 calculate genome sizes
  genomeSizes <- sapply(faFiles, function(x)
    sum(fasta.seqlengths(x)))

  # -- 1.3 raw input files against the reference
  if(!sameChrOnly){
    nonRefGenomes <- genomeIDs[genomeIDs != refGenome]
    againstRef <- data.table(
      target = refGenome, query = nonRefGenomes)
    againstRef[,`:=`(rawTargetFa = faFiles[target],
                     rawQueryFa = faFiles[query],
                     type = "againstRef")]
  }else{
    againstRef <- NULL
  }

  # -- 1.4 raw input files daisy chained
  dchain <- data.table(
    genome1 = genomeIDs[-1],
    genome2 = genomeIDs[-length(genomeIDs)])
  dchain[,`:=`(size1 = genomeSizes[genome1],
               size2 = genomeSizes[genome2])]
  dchain[,`:=`(query = ifelse(size1 > size2, genome1, genome2),
               target = ifelse(size1 > size2, genome2, genome1))]
  dchain[,`:=`(rawTargetFa = faFiles[target],
               rawQueryFa = faFiles[query],
               type = "daisyChain")]

  if(!sameChrOnly){
    dchain <- dchain[,colnames(againstRef), with = F]
    rawFiles <- rbind(againstRef, dchain)
  }else{
    rawFiles <- dchain
  }

  rawFiles <- subset(rawFiles, !duplicated(paste(query, target)))
  rawFiles[,`:=`(
    prepTargetFile = file.path(wd, sprintf("%s.fa", target)),
    prepQueryFile = file.path(wd, sprintf("%s.window.fa", query)),
    pafFile = file.path(wd, sprintf("%s_windows__vs__%s.paf.gz", query, target)),
    synFile = file.path(wd, sprintf("%s_windows__vs__%s.syn.paf.gz", query, target)),
    dotplotFile = file.path(wd, sprintf("%s_windows__vs__%s.dotplot.pdf", query, target)))]

  if(verbose)
    cat(sprintf("will conduct %s windowed alignments\n", nrow(rawFiles)))

  ##############################################################################
  # 2. Window and copy fasta files
  if(verbose)
    cat("2. Preparing minimap2 input files\n")

  # -- for each row in the input metadata
  paraml <- lapply(1:nrow(rawFiles), function(i){
    x <- rawFiles[i,]

    # -- 2.1 read in the query
    if(!file.exists(x$prepQueryFile) || overwriteInput){
      if(verbose)
        cat(sprintf("\tquery: %s", x$query))
      qss <- readDNAStringSet(x$rawQueryFa)

      # -- 2.2 subset to large enough chrs
      qss <- qss[width(qss) > minChrLen]

      # -- 2.3 parse the names
      names(qss) <- gsub(stripChrname, "", names(qss))

      # -- 2.4 calculate window size
      nbp <- sum(width(qss))
      stepSize <- round(nbp / nWindows)
      if(stepSize < 1)
        stepSize <- 1

      # -- 2. 5 chop into windows
      qssw <- window_ss(
        ss = qss,
        windowSize = windowSize,
        windowStep = stepSize,
        returnGranges = FALSE)

      # -- 2.6 write query
      writeXStringSet(qssw, filepath = x$prepQueryFile)
      if(verbose)
        cat(sprintf(" (%sMb in %s windows)\n",
                    round(sum(width(qssw)) / 1e6, 1), length(qssw)))
    }

    # -- 2.7 read in the target
    if(!file.exists(x$prepTargetFile) || overwriteInput){
      if(verbose)
        cat(sprintf("\ttarget: %s", x$target))
      tss <- readDNAStringSet(x$rawTargetFa)

      # -- 2.8  parse names and keep long enuf chrs
      names(tss) <- gsub(stripChrname, "", names(tss))
      tss <- tss[width(tss) > minChrLen]

      # -- 2.9 write uncompressed
      writeXStringSet(tss, filepath = x$prepTargetFile)
      cat(sprintf(" (%sMb total in %s seqeunces longer than %sMb)\n",
                  round(sum(width(tss)) / 1e6, 1), length(tss),
                  round(minChrLen/1e6, 1)))
    }
  })

  if(verbose)
    cat("\tDone!\n")

  ##############################################################################
  # 3. Run minimap2 window --> target
  if(verbose)
    cat("3. Running minimap2 on windows --> targets\n")
  mml <- lapply(1:nrow(rawFiles), function(i){
    x <- rawFiles[i,]
    if(verbose)
      with(x, cat(sprintf("\t%s (window) vs. %s (target) ... ", query, target)))
    if(!file.exists(x$pafFile) || overwriteMinimap2){
      outmm2 <- minimap2_windows(
        queryFasta = x$prepQueryFile,
        targetFasta = x$prepTargetFile,
        pafFile = x$pafFile,
        maxTopHits = maxTopHits,
        kmerSize = kmerSize,
        kmerStep = kmerStep,
        filterTopReps = filterTopReps,
        nCores = nCores,
        verbose = mm2verbose,
        overwrite = TRUE,
        windowSize = windowSize)
      if(verbose)
        cat("Done\n")
      return(outmm2)
    }else{
      cat("file exists and not overwriting\n")
    }
  })

  ##############################################################################
  # 4. Synteny constraint
  if(verbose)
    cat("4. Parsing synteny ... \n")

  synl <- lapply(1:nrow(rawFiles), function(i){
    x <- rawFiles[i,]
    if(verbose)
      cat(sprintf("\t%s (window) vs. %s (target)", x$query, x$target))
    if(!file.exists(x$synFile) || overwriteSynpaf){

      # -- 4.1 read in the paf
      h1 <- fread_paf(x$pafFile)
      setorder(h1, qname, -nmatch, tlen)

      # -- 4.2 parse input
      h1 <- parse_windPaf(
        h1,
        queryFastaFile = x$rawQueryFa,
        stripChrname = stripChrname)

      if(sameChrOnly)
        h1 <- subset(h1, tname == qname)

      if(minMapq > 0)
        h1 <- subset(h1, mapq >= minMapq)

      if(!dir.exists(tmpDir)){
        dir.create(tmpDir)
        on.exit(expr = unlink(tmpDir, recursive = T))
      }

      # -- 4.3 synteny and dotplotting
      synHits <- synteny_paf(
        paf = h1,
        makePlots = makePlots,
        blkSize = blkSize,
        tmpDir = tmpDir,
        pdfFile = x$dotplotFile,
        ppi = ppi,
        minPlotDim = minPlotDim,
        maxBestMaps = maxBestMaps,
        MCScanX_hCall = MCScanX_hCall,
        topNhits2keep = topNhits2keep,
        stepSize = 1000,
        windowSize = windowSize,
        buffer = buffer,
        maxDistBtwHitsInBlk = maxDistBtwHitsInBlk,
        keepEverything = keepEverything,
        targetNhitsQuantile = targetNhitsQuantile)

      # -- write it
      fwrite(
        synHits,
        file = x$synFile,
        sep = "\t", quote = FALSE, col.names = FALSE)
    }else{
      cat("... file exists and not overwriting\n")
    }
    return(x)
  })

  if(verbose)
    cat("5. Generating riparian plot ... ")
  pl1 <- riparian_paf(
    pafFiles = rawFiles$synFile,
    refGenome = refGenome,
    genomeIDs = genomeIDs,
    sameChrOnly = sameChrOnly,
    orderyBySynteny = TRUE,
    sameScale = FALSE,
    braidOffset = .075,
    chrWidth = .05,
    highlightInversions = highlightInversions,
    braidAlpha = braidAlpha,
    ...)
  pl2 <- riparian_paf(
    pafFiles = rawFiles$synFile,
    refGenome = refGenome,
    genomeIDs = genomeIDs,
    sameChrOnly = sameChrOnly,
    orderyBySynteny = TRUE,
    sameScale = TRUE,
    braidOffset = .075,
    chrWidth = .05,
    highlightInversions = highlightInversions,
    braidAlpha = braidAlpha,
    ...)
  if(verbose)
    cat(sprintf("\n\t wrote to: riparian_phasedBy%s.pdf", refGenome))

  ripFile <- file.path(wd, sprintf("riparian_phasedBy%s.pdf", refGenome))
  pdf(ripFile, height = length(genomeIDs) * ripHeightPerGenome, width = 8)
  print(pl1)
  print(pl2)
  dev.off()

  return(list(mapFilePaths = rawFiles,
              genomeIDs = genomeIDs,
              sameChrOnly = sameChrOnly,
              refGenome = refGenome))
}
