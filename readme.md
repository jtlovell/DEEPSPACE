
# This is DEEPSPACE v0.2

## Release notes

Currently only the alpha v0.2 release is available. Testing has only been completed on a limited set of genomes. **DO NOT USE THIS WITHOUT CAREFULLY VETTING THE RESULTS**. 

To help with development, if you do use < v0.2, please contact JT Lovell and let him know what you have used it for and how it went. If you find issues, please report them. 

## Overview 

The DEEPSPACE R package is an extension of GENESPACE to the whole genome (not just genes). 

In addition to synteny/riparian mapping functionalities, many of the utility functions and genome-wide visualization tools that were initially incorporated into GENESPACE have migrated here. 
DEEPSPACE synteny mapping functionality is designed to quickly summarize and plot genome-wide complexity-reduced assembly-based alignments across fairly closely related organisms. Depending on settings and genomes, DEEPSPACE may function 2-1000X faster than minimap2-SyRI. For example, with the 'fast' setting, we can produce a synteny map across 10 2Gb+ mammalian genomes in XXX minutes with 8 cores (YMMV). 

DEEPSPACE synteny mapping is NOT a substitute for genome-wide variant detection. The smallest SVs that can be detected are a function of the alignment window size and the gap join size parameters (see below). The defaults for these in the human genome result in a minimum INDEL size of XXX and a minimum inversion/translocation size of XXX; any variant smaller than these will be missed. 

## Functionality

DEEPSPACE offers fast sequence summaries and comparisons in R with the following:

1. `mm2_windows`: Short sequence ('windowed') alignments between paird of genomes
2. `dotplot_paf`: Dot or segment plots of paf-formatted alignments
3. `syntenic_windows`: Gene-free syntenic hits and blocks from windowed alignments. 
4. `riparian_paf`: One-line riparian plotting from a set of paf-formatted alignments.
5. `annotated_riparian`: Highly customizable riparian plotting allowing for tracks with sliding windows and genome categorization blocks.
6. `categorize_genome`: Neural network genome binning and categorization into 'telomere', 'arm', 'centromere', and 'pericentromere'. 
7. `summarize_genome`: Calculate contiguity statistics and the positions/sizes of telomeres
8. `plot_contigs`: Contig maps annotated with the positions of gaps, contigs, and telomeres
9. Many utilities to easily script out sliding windows, kmer positions, gene/repeat densities, etc including: `count_overlapsByGroup`, `window_gr`, `pull_kmerBlocks`, `find_telomeres`, `find_gaps`, `find_contigs`, ....
10. Utility functions that operate on or convert to GenomicRanges objects: `read_gffAsGr`, `make_hierarchical`, `join_gappedRanges`, ....
11. Utility functions to convert from GENESPACE

Like GENESPACE, there are high-level pipelines that simply takes a vector of file paths and produces graphics: `clean_windows` (riparian and dot plotting of synteny), and `qc_genomes` (telomere/contig/gap statistics and plotting). Details regarding these pipelines are below. 

## Installation

To install from scratch, you need R, two 3rd party programs (optionally to run synteny mapping), and a few R packages. Here is how to install these:

#### Install R

DEEPSPACE is meant to be run interactively in the R environment for statistical computing. So, you need to have R installed. See [CRAN](https://www.r-project.org/) for the most recent release. 

#### Install minimap2

minimap2 is most simply installed via conda (in the shell, not R).

```
conda create -n minimap2
conda activate minimap2
conda install -c bioconda minimap2 
```

If conda is not available on your machine, you can install minimap2 from a number of other sources. See minimap2 documentation for details.
Regardless of how minimap2 is installed, ensure that you have minimap2 version >= 2.28.

#### Install MCScanX

`MCScanX` should be installed from [github](https://github.com/wyp1125/MCScanX). 

#### Install DEEPSPACE

Once the above 3rd party dependencies are installed, get into R. If you made a conda environment, its useful to open R directly from that environment so that minimap2 stays in the path.

```
conda activate minimap2
open -na rstudio # if using rstudio, otherwise, simply `R`
```

Once in R, install R dependencies: if they are not yet installed, install_github will install a few dependencies directly (ggplot2, dbscan, R.utils, data.table). However, you will need to install the bioconductor packages separately:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer", "GenomicRanges", "Rsamtools"))
```

The easiest way to install DEEPSPACE uses the package devtools (which may need to be installed separately):

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/DEEPSPACE")

library(DEEPSPACE)
```

## one-line multi-genome synteny maps

**See below for pipeline details, or running the steps independently.**

For any set of genome assemblies, `clean_windows` builds a daisy chain (G1-->G2, G2-->G3, ... Gn-1-->Gn) and does pairwise windowed alignments. If the genomes are from groups with intrachromosomal translocations (which is assumed by default), then a reference genome is also aligned to all genomes and any windows in the other genomes are phased by the chromosome they hit in the reference. The windowed alignments are then parsed so that only syntenic hits are retained and the syntenic hits are converted into blocks (and optionally phased by the reference). Finally, the blocks are plotted in a multiple-genome synteny map called a riparian plot. 

To run this, both `minimap2` must be in the path and `MCScanX` must be installed. 

```
genomeFastas <- c("path/genome1.fa", "path/genome2.fa")
names(genomeFastas) <- c("genomeID1", "genomeID2")
workHere <- "path/to/working/directory"
blks <- clean_windows(
    faFiles = genomeFastas, 
    wd = workHere,
    MCScanX_hCall = "path/to/MCScanX/MCScanX_h")
```

This will produce a riparian plot in the working directory and output an annotated paf-like data.table that can then be passed into `riparian_paf` for further customization. 

#### considerations when running `clean_windows`

The three most important parameters are:

- `speed` (options are: "slow", "fast" (default), "veryFast") ... this alters the mm2 word and kmer gap size parameters and the number of windows
- `divergence` (options are: 5, 10, 20 (default))` ... this alters the mm2 XXX and DEEPSPACE window size 
- `repeats` (options are: "ignore", "some" (default), "all")`... this alters the mm2 masking parameter and the number of top hits retained. 

The pipeline is built for single-copy moderately closely related genomes. It is very difficult to track paralogs between adjacent genomes except in the case of very recent whole genome duplications. If seeing all paralogs is a goal, we strongly recommend GENESPACE instead. Adjacent genomes with variable copy number (so that both paralogs in duplicated genomes are orthologs to the single copy genome) often work just fine in DEEPSPACE. These contrasts can be improved by adjusting the 'ploidy' parameter which adjusts the number of syntenic hits that can be seeds for any given window. 

## one-line genome QC and contig maps

COMING SOON



## Step-by-step synteny maps for a single pair of genomes

The options below should be used for complex situations where genomes are very differently related

#### 1. `mm2_windows`: align genomic windows (query) to a full genome (target)

This applies a sliding window to the query genome, write the uncompressed windowed query and uncompressed full assembly target to file, then align the two sequences using minimap2 with the user-specified arguments. This outputs a paf-like data.table object with just the 12 standard paf columns. We can then write this to file if desired. 

```
rawPaf <- mm2_windows(query, target, width, step, mm2args, ...)
dotplot_paf(rawPaf, ...)
fwrite_paf(paf = rawPaf, file = "/path/to/raw.paf", keepExtraColumns = FALSE)
```

#### 2. `synteny_windows`: convert windowed alignments into syntenic blocks

This take a paf-formatted "hits" data.table and parses it into syntenic regions. 

```
synPaf <- synteny_windows(paf, ...)
dotplot_paf(synPaf, ...)
fwrite_paf(paf = synPaf, file = "/path/to/syn.paf", keepExtraColumns = TRUE)
```

#### 3. `phase_windows`: determine groupings for windows from a paf

For a query (windowed) genome, determine positions against a reference. Windows may or may not be allowed to hit multiple places

```
windowPhases <- phase_windows(paf, ...)
phasePaf <- add_phase2paf(synPaf, windowPhases)
riparian_paf(phasePaf, colorColumn = "phase")
```
