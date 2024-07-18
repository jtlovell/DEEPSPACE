
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


### `clean_windows` pipeline overview

DEEPSPACE is an annotation-free version of GENESPACE. Feed it the genome assemblies and it produces the standard GENESPACE graphical outputs (dotplots, riparian). DEEPSPACE does not produce pan-gene sets. There are four main steps to the pipeline:

1. Parse the fasta files. The following alignments are run: (a) all non-reference genomes windowed against the reference genome ... by default, the first genome provided is the "reference", but this can be changed by the user; (b) a 'daisy chain' of all genomes in the order given in genomeIDs so that for genomes A, B, C, alignments are done A --> B and B --> C where the larger genomes are windowed and aligned to the smaller genome
2. Run minimap2. The windowed genome is always the 'query' and the other genome is the target. The user can specify all possible parameters to minimap2 with a series of parameters in `clean_windows`. See the help file for more details. 
3. Parse the minimap2 results (.paf) into syntenic blocks. This generally follows the GENESPACE synteny pipeline, treating minimap2 windowed alignments as it would the peptide-peptide alignments. 

### `clean_windows` output

The main outputs the users will be interested in are:

- synteny-constrained windowed hits stored in [workingDirectory]/[QUERY]_window__vs__[TARGET].syn.paf 
- dotplots: [QUERY]_window__vs__[TARGET]_dotplots.pdf
- syntenic block break points: [QUERY]_window__vs__[TARGET]_blocks.paf
- riparian plots: riparian_phasedBy[referenceGenome].pdf

Once the synteny-constrained windowed hits (.paf) files have been produced, 
you can customize riparian plots in a similar manner to those in GENESPACE. 

### `clean_windows` methods and options

There are many parameters for fine-tuning your run (see `?clean_windows`) and three basic settings that help choose the best settings for your run, which is specified in `setting=c("exploreRepeats", "fast", "default")`

1. The `default` clean_windows run is meant to optimize speed and sensitivity. Given the lightweight nature of the pipeline, its easy to re-run many times for a subset of your genomes and figure out what parameter combinations work best for your goals (look at the dotplots). The main parameters for the users to explore is .. 

    - `nWindows` the number of windows to chop the query genome into
    - `windowSize` the basepair size of each window: larger windows are better for genomes with lots of paralogs or significant divergence
    - `filterTopReps` the `-p` parameter for minimap2 which controlls alignment in repetitive regions
    - `blkSize` like the parameter of the same name in GENESPACE: how many hits need to be syntenic for a region to be called?
    - `maxDistBtwHitsInBlk` the maximum distance between two syntenic anchors within the same block. Smaller values give more apparent INDELs (often at repeat regions). 
    
2. `fast` runs minimap keeping only the primary alignment for each query window and only mapq=60 alignments. This will produce riparian plots with apparent gaps anywhere the target genome, and potentially the query windowed genome, is repetitive. Minimap2 kmer sensitivity is also reduced. The subsequent smaller input .paf file for synteny searches massively speeds up all downstream analysis. 
3. `exploreRepeats` runs minimap with -p0 and keeps all hits regardless of whether they map multiple times. This is much slower than the default and should generally be used only for a handfull of genomes. 


