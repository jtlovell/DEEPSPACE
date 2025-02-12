

binContigs_byChr <- function(refGenomeFasta,
                             contigsFasta,
                             tmpDir,
                             maxMarkersPerContig = 1000,
                             minMarkersPerContig = 50,
                             markerSize = 100,
                             verbose = TRUE,
                             overwrite = FALSE,
                             minChrSize = 0,
                             nCores = 1){

  ##############################################################################
  if(!dir.exists(tmpDir))
    dir.create(tmpDir)

  qfa <- file.path(tmpDir, "contigs_windows.fa")
  tfa <- file.path(tmpDir, "ref.fa")
  paff <- file.path(tmpDir, "contigs_chr.paf")

  if(overwrite && file.exists(paff))
    unlink(paff)

  ##############################################################################
  # 1. Modified window fasta, with fewer windows on large contigs
  if(!file.exists(qfa) || !file.exists(tfa) || overwrite){
    # -- 1.1 Read in and window the contigs fasta
    qf <- "/Users/jlovell/Desktop/agaveCSP/phaseMX/MX_100window_step200.fasta"
    if(verbose)
      cat("Reading in contigs fasta file ... ")
    ss <- readDNAStringSet(contigsFasta)

    # -- 1.2 Window the contigs fasta
    if(verbose)
      cat("Done!\nWindowing contigs ... ")
    dt <- data.table(
      seqnames = names(ss),
      start = 1,
      end = width(ss))

    # -- 1.3 adaptively choose the number of windows based on above rules
    dt[,n := floor(end / markerSize), by = "seqnames"]
    dt[,n := ifelse(n > maxMarkersPerContig, maxMarkersPerContig,
                    ifelse(n < minMarkersPerContig, minMarkersPerContig, n))]

    # -- 1.4 make the windows
    dt <- dt[,list(
      start = seq(from = min(start),
                  to = max(end),
                  length.out = n),
      maxp = max(end)),
      by = "seqnames"]
    dt[,start := floor(start)]
    dt[,end := start + (markerSize - 1)]
    dt <- subset(dt, end <= maxp)
    dt[,maxp := NULL]

    # -- 1.5 convert to granges
    gr <- as(as.data.frame(dt[,c("seqnames", "start", "end")]), "GRanges")

    # -- 1.6 subset stringset
    out <- ss[gr]

    # -- 1.7 rename the stringset to be parsed later
    names(out) <- sprintf(
      "%s_XstartX%s_XendX%s",
      seqnames(gr),
      gsub(" ", "", format(start(gr), scientific = FALSE)),
      gsub(" ", "", format(end(gr), scientific = FALSE)))

    # -- 1.8 write the output
    writeXStringSet(out, filepath = qfa)

    # -- 1.9 copy the reference over
    if(verbose)
      cat("Done!\nCopying reference to tmp directory ... ")
    ss <- readDNAStringSet(refGenomeFasta)
    ss <- ss[width(ss) >= minChrSize]
    writeXStringSet(ss, filepath = tfa)

    if(verbose)
      cat("Done!\n")
  }

  ##############################################################################
  # 2. Run minimap
  if(!file.exists(paff) || overwrite){
    if(verbose)
      cat(sprintf(
        "Running minimap2 on %s windows mapped to %s chromosomes ... \n",
        length(out), length(ss)))
    out <- minimap2_windows(
      queryFasta = qfa,
      targetFasta = tfa,
      pafFile = paff,
      windowSize = markerSize,
      filterTopReps = .01,
      maxTopHits = 2,
      minSecRat = .99,
      kmerSize = 19,
      kmerStep = 19,
      nCores = nCores)
  }

  ##############################################################################
  # 3. Parse the results
  if(verbose)
    cat("\nParsing paf results ... ")
  p <- fread_paf(paste0(paff, ".gz"))
  p <- subset(p, mapq == 60)
  p <- subset(p, alen >= (markerSize * .8))
  setorder(p, -nmatch, -alen)
  p <- subset(p, !duplicated(qname))
  out <- parse_windPaf(
    paf = p, queryFastaFile = contigsFasta, stripChrname = "")
  out[,nHits := .N, by = "qname"]
  outl <- out[,list(n = .N, min = min(qstart), max = max(qend)),
              by = c("qname", "tname", "nHits", "qlen", "tlen")]
  outl[,prop := n / nHits]
  outl[,flag := prop < .9 | n < (minMarkersPerContig / 2)]
  t1 <- subset(outl, !flag)
  t2 <- subset(subset(outl, flag), !duplicated(qname))
  if(verbose)
    cat(sprintf(
      "Done!\n\t%s contigs (%s Mb) assigned to chrs\n\t%s unassigned contigs (%s Mb)\n",
      nrow(t1), round(sum(t1$qlen)/1e6, 1), nrow(t2), round(sum(t2$qlen)/1e6, 1)))
  return(outl)
}



parse_scafPaf <- function(pafFile, minMapq = 55){
  y <- fread_paf(pafFile)
  y[,hap := ifelse(grepl("pri", qname), "Primary", "Alt")]
  y <- subset(y, mapq > minMapq)
  setorder(y, tstart, qlen)
  y[,maxLen := max(alen), by = "qname"]
  y <- subset(y, alen == maxLen)
  y[,qord := factor(qname, levels = unique(qname))]
  print(ggplot(y, aes(x = tstart/1e6, xend = tend/1e6, y = qord, yend = qord, col = hap))+
          geom_segment()+
          theme(axis.text.y = element_text(size = 5))+
          labs(title = y$tname[1],
               x = "target physical position (Mb)", col = "haplotype"))
  return(y)
}

flag_internalContigs <- function(paf){
  pafgr <- with(paf, as(
    data.frame(seqnames = tname, start = tstart, end = tend),
    "GRanges"))
  ovlp <- findOverlaps(pafgr, pafgr, type = "within")
  win <- data.table(subset(as.data.frame(ovlp), queryHits != subjectHits))
  win[,`:=`(qname1 = paf$qname[queryHits],
            qname2 = paf$qname[subjectHits],
            qsize1 = paf$alen[queryHits],
            qsize2 = paf$alen[subjectHits])]
  smallWithin <- unique(with(win, ifelse(qsize1 < qsize2, qname1, qname2)))
  paf[,smallWithin := qname %in% smallWithin]
  p1 <- ggplot(paf, aes(x = tstart/1e6, xend = tend/1e6, y = qord, yend = qord, col = smallWithin))+
    geom_segment()+
    theme(axis.text.y = element_text(size = 5))+
    labs(title = paf$tname[1],
         x = "target physical position (Mb)", col = "is smallWithin")
  p2 <- ggplot(subset(paf, !smallWithin), aes(x = tstart/1e6, xend = tend/1e6, y = qord, yend = qord, col = hap))+
    geom_segment()+
    theme(axis.text.y = element_text(size = 5))+
    labs(title = paf$tname[1],
         x = "target physical position (Mb)", col = "haplotype")
  print(gridExtra::grid.arrange(p1, p2, nrow = 1))
  return(paf)
}

cut_ovlpContigs <- function(paf,
                            minGap2fill = 1e3,
                            plotIt = FALSE,
                            verbose = TRUE){
  cut_engine <- function(pafChr, plotIt, minGap2fill){
    x <- subset(pafChr, !smallWithin)
    xsv <- data.table(x)
    x[,`:=`(dist2r = qlen - qend)]
    setorder(x, tstart, -qlen)

    x[,tendb := ifelse(strand == "+", tend + dist2r, tend + qstart)]
    x[,tstartb := ifelse(strand == "-", tstart - dist2r, tstart - qstart)]

    for(i in 1:(nrow(x) - 1)){
      x1 <- x[i,]
      x2 <- x[i+1,]
      if(x1$tend > x2$tstart){
        if(x1$qlen > x2$qlen){
          x[i+1,"tstart"] <- x$tend[i] + 1
        }else{
          x[i,"tend"] <- x$tstart[i+1] - 1
        }
      }
    }

    x[,bigEnough := abs(tend - tstart) > minGap2fill]
    qns <- unique(x$qname[x$bigEnough])
    xsv1 <- subset(xsv, qname %in% qns)
    x1 <- subset(x, bigEnough)
    if(plotIt)
      print(ggplot(xsv1, aes(x = tstart/1e6, xend = tend/1e6, y = qord, yend = qord, col = hap))+
              geom_segment(data = xsv1, color = "grey", linewidth = 2)+
              geom_segment(data = xsv1, color = "black")+
              geom_segment(data = x1)+
              theme(axis.text.y = element_text(size = 5))+
              labs(title = xsv1$tname[1],
                   x = "target physical position (Mb)", col = "haplotype"))
    return(x)
  }

  out <- rbindlist(lapply(split(paf, by = "tname"),
                          cut_engine, plotIt = plotIt, minGap2fill = minGap2fill))

  if(verbose)
    with(subset(out, bigEnough), cat(
      sprintf("%s contigs (%sMb N50) covering %sMb (%s%%) of the genome\n",
              uniqueN(qname),
              round(N50(abs(qstart - qend))/1e6, 1),
              round(sum(tend-tstart)/1e6, 1),
              round(sum(tend-tstart)/sum(tlen[!duplicated(tname)]), 3)*100)))
  return(out)
}
