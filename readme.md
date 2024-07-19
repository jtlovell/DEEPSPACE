
# This is DEEPSPACE v0.1

## Release notes

Currently only the alpha v0.1 release is available. Testing has only been completed on a limited set of genomes. **DO NOT USE THIS WITHOUT CAREFULLY VETTING THE RESULTS**. 

To help with development, if you do use <v0.2, please contact JT Lovell and let him know what you have used it for and how it went. If you find issues, please report them. 

DEEPSPACE v0.1 is only to be used in single-copy genomes. It will generally ignore homeologs (in polyploids) and other forms of large-scale syntenic paralogs (e.g. segmental duplicates). This is run with `ploidy = 1`.

**NOTE: known issue** If a genome is a recent polyploid and one subgenome is missing a large segment of a chromosome, the alternative subgenome may align in the '1x' implementation. No work-around is planned for this ... so, either carefully vet polyploid genome alignments when running them with `ploidy = 1`, or wait until ~v0.5 for ploidy > 1 to be allowed. This may be partially resolved using a larger window and a higher `minPid` value.


## Functionality 

DEEPSPACE has three basic functionalities:

1. Genome visualization and quality control. Numerous functions and tools are presented in [coming soon].
2. Genome binning and categorization: `categorize_genome()`
3. Gene-free synteny maps: `clean_windows()`

### `categorize_genome()` pipeline overview

[coming soon]


## `clean_windows` 

#### Pipeline overview

DEEPSPACE is an annotation-free version of GENESPACE. Feed it the genome assemblies and it produces the standard GENESPACE graphical outputs (dotplots, riparian). DEEPSPACE does not produce pan-gene sets. There are four main steps to the pipeline:

1. Parse the fasta files. The following alignments are run: (a) all non-reference genomes windowed against the reference genome ... by default, the first genome provided is the "reference", but this can be changed by the user; (b) a 'daisy chain' of all genomes in the order given in genomeIDs so that for genomes A, B, C, alignments are done A --> B and B --> C where the larger genomes are windowed and aligned to the smaller genome
2. Run minimap2. The windowed genome is always the 'query' and the other genome is the target. The user can specify all possible parameters to minimap2 with a series of parameters in `clean_windows`. See the help file for more details. 
3. Parse the minimap2 results (.paf) into syntenic blocks. This generally follows the GENESPACE synteny pipeline, treating minimap2 windowed alignments as it would the peptide-peptide alignments. 

#### How to run it

`clean_windows` takes two required inputs, the directory in which to conduct the run and the genome assembly fasta files to use. The fasta files should be a named vector. The names of the files are used for labels in the output files. 

As an example, download three reference genomes from NCBI (or your favorite repo). 

```{r}
# -- increase timeout since we are downloading three big files
options(timeout = 1000) 

# -- make a fresh directory to put the files in
workHere <- "~/Desktop/tryDeepspace" 
dir.create(workHere)

# -- provide urls to the genome assembly files
urls <- c(
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/281/585/GCF_029281585.2_NHGRI_mGorGor1-v2.0_pri/GCF_029281585.2_NHGRI_mGorGor1-v2.0_pri_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz")

# -- provide names of the files to download to
fastaFiles <- file.path(workHere, basename(urls))

# -- give both the urls and the fasta files unique IDs, these will be the labels
genomeIDs <- c("human", "gorilla", "mouse")
names(urls) <- genomeIDs
names(fastaFiles) <- genomeIDs

# -- download the files
for(i in genomeIDs)
  download.file(urls[i], fastaFiles[i])

```

Once the files are downloaded, we can start the run. We have to add a few user-defined parameters to this run because (a) mouse-[human/gorrila] divergence is deep enough that synteny is not super clear at the basepair level, and (b) NCBI genome fasta headers need some parsing. 

Here, in addition to the faFiles (named paths to the assembly fastas) and the wd (path to the working directory), we give the following:

- preset = "dist fast", which indicates that the genomes are diverged and we want to run the pipeline as fast is we can without losing much sensitivity. The distant setting decreases the minimum blkSize: from 10 to 5, increase the windowSize from from 1kb to 5kb and increases the maxDistBtwHitsInBlk from 500kb to 1Mb. The fast setting switches the minimap2 kmer size and kmer step from 25/20 to 19/19, decreases the number of windows from 500k to 200k, increases the minimum mapq from 0 to 12 and increase the repetitive kmer screening from 0.0001 to 0.001. 
- stripChrname: NCBI genomes store lots of info in the chromosome names (fasta headers) to get just the chr number, strip off everything before 'chromosome ' and after the first ',' with this regex ".*chromosome |,.*"
- minChrLen: NCBI keeps lots of ALT and unanchored sequences in some assemblies. You can ignore these by setting this number to a bit below the  size of the smallest chromosome you are interested in. 

You also need to give some parameters that are specific to your system: the path to the MCScanX_h executable and minimap2. If these are in the path, just the executable name is fine. 

Also, be sure that there is >15Gb of space free on your machine. 
Depending on the number of cores you can give minimap2, each alignment should take 1-5min and the synteny search should take about 1min per pair. In total, the run should take 5-20mins. 

```{r}
test <- clean_windows(
  faFiles = fastaFiles,
  genomeIDs = c("gorilla", "human", "mouse"),
  wd = workHere,
  
  preset = "dist fast",
  
  stripChrname = ".*chromosome |,.*",
  minChrLen = 20e6, 
  
  nCores = 8,
  MCScanX_hCall = "/path/to/MCScanX/MCScanX_h",
  minimap2call = "/path/to/minimap2")
```

#### `clean_windows` output

The main outputs the users will be interested in are:

- synteny-constrained windowed hits stored in [workingDirectory]/[QUERY]_window__vs__[TARGET].syn.paf 
- dotplots: [QUERY]_window__vs__[TARGET]_dotplots.pdf
- syntenic block break points: [QUERY]_window__vs__[TARGET]_blocks.paf
- riparian plots: riparian_phasedBy[referenceGenome].pdf

#### Customizing the riparian plot

Once the synteny-constrained windowed hits (.paf) files have been produced, you can customize riparian plots in a similar manner to those in GENESPACE. `riparian_paf` is called internally in `clean_windows`, but it can be called directly to customize the output. For example, if you wanted to make a plot with only dark grey for collinear regions and orange to highlight all inverted sequence, and order by the chromosome names (and not synteny to the reference genome), you could run ...

```{r}
riparian_paf(
    pafFiles = test$mapFilePaths$synFile,
    refGenome = test$refGenome,
    genomeIDs = test$genomeIDs[c(2,1,3)],
    orderyBySynteny = FALSE,
    braidOffset = .075,
    braidColors = "grey20", 
    highlightInversions = "darkorange",
    braidAlpha = .9)
```

There are many other ways to customize the riparian plot. Not quite as flexible as the options in GENESPACE, but it will get there eventually. 



