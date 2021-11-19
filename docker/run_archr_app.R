#!/usr/bin/env Rscript
library(optparse)
# https://github.com/GreenleafLab/ArchR/discussions/1044#discussioncomment-1405648
# Needed when using R 4.1.1
library(parallel)
library(magick)

option_list = list(
  make_option(
    c("-b", "--project_directory"),
    type="character",
    default=NULL,
    help="path to directory which contains ArrowFiles folder."
  )
)
#  make_option(
#    c("-b", "--bam_file"),
#    type="character",
#    default=NULL,
#    help="BAM file to use."
#  ),
#  make_option(
#    c("-t", "--threads"),
#    type="integer",
#    default=2,
#    help="Number of subprocesses/threads to use."
#  ),
#  make_option(
#    c("-e", "--minTSS"),
#    type="double",
#    default=1.5,
#    help="The minimum numeric transcription start site (TSS) enrichment score required to pass filtering. E.g. 1.5"
#  ),
#  make_option(
#    c("-g", "--minFrags"),
#    type="integer",
#    default=2000,
#    help="The minimum number of mapped ATAC-seq fragments required per cell to pass filtering. E.g. 2000"
#  )
#)
#
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#
if (is.null(opt$project_directory)){
  print_help(opt_parser)
  stop("--project_directory must be specified.", call.=FALSE)
}
#
## First, we load the ArchR library. If this fails, you have not properly installed
## ArchR and should revisit the installation instructions. We also recommend setting
## and remembering a known seed to facilitate replication of operations requiring randomization.
library(ArchR)
#
## Next, we set the default number of threads for parallelized operations in ArchR
## functions. You should change the value passed to threads to match the
## specifications of your local machine.
#addArchRThreads(threads = opt$threads) 
#
#script.dir <- getwd()
#message(paste("directory used is:",script.dir))
#
#inputFiles <- c(opt$bam_file)
#names(inputFiles) <- c("BAM_data")
#inputFiles
#
## Before we begin, we need add a reference genome annotation for ArchR to have
## access to chromosome and gene information. ArchR natively supports hg19, hg38, mm9, and mm10.
#addArchRGenome("hg38")
#
## Creating Arrow Files
## Now we will create our Arrow files. For each sample, this step will:
#
## Read accessible fragments from the provided input files.
## Calculate quality control information for each cell (i.e. TSS enrichment scores and nucleosome info).
## Filter cells based on quality control parameters.
## Create a genome-wide TileMatrix using 500-bp bins.
## Create a GeneScoreMatrix using the custom geneAnnotation that was defined when we called addArchRGenome().
#ArrowFiles <- createArrowFiles(
#  inputFiles = inputFiles,
#  sampleNames = names(inputFiles),
#  minTSS = opt$minTSS, # Dont set this too high because you can always increase later
#  minFrags = opt$minFrags,
#  addTileMat = TRUE,
#  addGeneScoreMat = TRUE,
#  bamFlag = list(isMinusStrand = FALSE, isProperPair = TRUE, isDuplicate = FALSE),
#  bcTag = "CB" # We added this tag to the SAM file and then converted it to a BAM
#)
#
## Inferring Doublets
## After Arrow file creation, we can infer potential doublets (a single droplet
## containing multiple cells) that can confound downstream results. This is
## done using the addDoubletScores() function.
#doubletScores <- addDoubletScores(
#  input = ArrowFiles,
#  k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
#  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
#)
#
#projSci <- ArchRProject(
#  ArrowFiles = ArrowFiles, 
#  outputDirectory = "ArchRProjFiles",
#  copyArrows = TRUE # This is recommened so that if you modify the Arrow files you have an original copy for later usage.
#)
#projSci
#
## We can check how much memory is used to store the ArchRProject in memory within R:
#message(paste0("Memory Size = ", round(object.size(projSci) / 10^6, 3), " MB"))
### [1] “Memory Size = 37.135 MB”
#
## We can also ask which data matrices are available within the ArchRProject which
## will be useful downstream once we start adding to this project:
#getAvailableMatrices(projSci)
### [1] “GeneScoreMatrix” “TileMatrix”
#
#
## We can access the cell names associated with each cell:
#head(projSci$cellNames)
### [1] “scATAC_BMMC_R1#TTATGTCAGTGATTAG-1” “scATAC_BMMC_R1#AAGATAGTCACCGCGA-1”
### [3] “scATAC_BMMC_R1#GCATTGAAGATTCCGT-1” “scATAC_BMMC_R1#TATGTTCAGGGTTCCC-1”
### [5] “scATAC_BMMC_R1#TCCATCGGTCCCGTGA-1” “scATAC_BMMC_R1#AGTTACGAGAACGTCG-1”
#
## We can access the sample names associated with each cell:
#head(projSci$Sample)
### [1] “scATAC_BMMC_R1” “scATAC_BMMC_R1” “scATAC_BMMC_R1” “scATAC_BMMC_R1”
### [5] “scATAC_BMMC_R1” “scATAC_BMMC_R1”
#
##We can access the TSS Enrichment Scores for each cell:
#quantile(projSci$TSSEnrichment)
### 0% 25% 50% 75% 100%
### 4.027 13.922 16.832 19.937 41.782
#
#
#saveArchRProject(ArchRProj = projSci, outputDirectory = "ArchRProjFiles", load = FALSE)
#
### Now we can filter putative doublets based on the previously determined
### doublet scores using the filterDoublets() function. This doesn’t physically
### remove data from the Arrow files but rather tells the ArchRProject to ignore
## these cells for downstream analysis.
##
#projSci <- filterDoublets(ArchRProj = projSci)
#
### Dimensionality Reduction and Clustering
### ArchR implements an iterative LSI dimensionality reduction via the addIterativeLSI() function.
#projSci <- addIterativeLSI(
#    ArchRProj = projSci,
#    useMatrix = "TileMatrix",
#    name = "IterativeLSI",
#    iterations = 5,
#    clusterParams = list(
#        resolution = c(2),
#        sampleCells = 10000,
#        maxClusters = 6, 
#        n.start = 10),
#    varFeatures = 25000
#    )
##
## To call clusters in this reduced dimension sub-space, we use the addClusters()
## function which uses Seurat’s graph clustering as the default clustering method.
#projSci <- addClusters(input = projSci, reducedDims = "IterativeLSI")
#
#
#message(paste("Adding UMAP"))
#projSci <- addUMAP(ArchRProj = projSci, reducedDims = "IterativeLSI")
#
#message(paste("Getting embedding"))
#projSciEmbeddingWClustersDF = getEmbedding(ArchRProj = projSci, embedding = "UMAP", returnDF = TRUE)
##write.csv(projSciEmbeddingWClustersDF, file='umap_embedding.csv')
#
#
#
### First, we load the ArchR library. If this fails, you have not properly installed ArchR and should revisit the installation instructions. We also recommend setting and remembering a known seed to facilitate replication of operations requiring randomization.
##
##library(ArchR)
##set.seed(1)
### Next, we set the default number of threads for parallelized operations in ArchR functions. You should change the value passed to threads to match the specifications of your local machine.
##
##addArchRThreads(threads = 4)
#### Setting default number of Parallel threads to 16.
##
###The Hematopoeisis tutorial data can be downloaded using the getTutorialData() function. The tutorial data is approximately 0.5 GB in size. If you have already downloaded the tutorial in the current working directory, ArchR will bypass downloading.
##
##inputFiles <- c("/mnt/ArchR/app_test/HemeFragments/scATAC_BMMC_R1.fragments.tsv.gz", "/mnt/ArchR/app_test/HemeFragments/scATAC_CD34_BMMC_R1.fragments.tsv.gz", "/mnt/ArchR/app_test/HemeFragments/scATAC_PBMC_R1.fragments.tsv.gz")
##names(inputFiles) <- c("scATAC_BMMC_R1", "scATAC_CD34_BMMC_R1", "scATAC_PBMC_R1")
###inputFiles <- getTutorialData("Hematopoiesis")
##inputFiles
##
### Before we begin, we need add a reference genome annotation for ArchR to have access to chromosome and gene information. ArchR natively supports hg19, hg38, mm9, and mm10.
##
##addArchRGenome("hg19")
#### Setting default genome to Hg19.
##
### Creating Arrow Files
### Now we will create our Arrow files which will take 10-15 minutes. For each sample, this step will:
##
### Read accessible fragments from the provided input files.
### Calculate quality control information for each cell (i.e. TSS enrichment scores and nucleosome info).
### Filter cells based on quality control parameters.
### Create a genome-wide TileMatrix using 500-bp bins.
### Create a GeneScoreMatrix using the custom geneAnnotation that was defined when we called addArchRGenome().
##ArrowFiles <- createArrowFiles(
##  inputFiles = inputFiles,
##  sampleNames = names(inputFiles),
##  minTSS = 4, #Dont set this too high because you can always increase later
##  minFrags = 1000,
##  addTileMat = TRUE,
##  addGeneScoreMat = TRUE,
###  bamFlag = list(isMinusStrand = FALSE, isProperPair = TRUE, isDuplicate = FALSE),
##)
#### Using GeneAnnotation set by addArchRGenome(Hg19)!
#### Using GeneAnnotation set by addArchRGenome(Hg19)!
#### ArchR logging to : ArchRLogs/ArchR-createArrows-69ef6ba2e1c7-Date-2020-04-21_Time-16-18-35.log
#### If there is an issue, please report to github with logFile!
#### Cleaning Temporary Files
#### 2020-04-21 16:18:35 : Batch Execution w/ safelapply!, 0 mins elapsed.
#### ArchR logging successful to : ArchRLogs/ArchR-createArrows-69ef6ba2e1c7-Date-2020-04-21_Time-16-18-35.log
##
### We can inspect the ArrowFiles object to see that it is actually just a character vector of Arrow file paths.
##
#
## In addition to plotting gene scores per cell as a UMAP overlay, we can browse the local chromatin accessibility at these marker genes on a per cluster basis with genome browser tracks. To do this, we use the plotBrowserTrack() function which will create a list of plots, one for each of the genes specified by markerGenes.
#pBrowserTrack <- plotBrowserTrack(
#    ArchRProj = projSci,
#    groupBy = "Clusters",
#    geneSymbol = markerGenes,
#    upstream = 50000,
#    downstream = 50000
#)
## Last but certainly not least, ArchR natively supports an interactive and dynamic genome browser that can be launched locally via a shiny app. To do this, we use the ArchRBrowser() function.
#
#ArchRBrowser(ArchRProj = proj)
#
## Saving and Loading an ArchRProject
## To easily save an ArchRProject for later use or for sharing with collaborators, we use the saveArchRProject() function. This copies the current ArchRProject object and all of the Arrow files to a specified directory. If we don’t specify an output directory (as below), saveArchRProject() uses the output directory that we specified upon creation of our ArchRProject. In this case that is the folder “HemeTutorial”.
#
#projSci <- saveArchRProject(ArchRProj = projSci)
### Saving ArchRProject…
### Loading ArchRProject…
### Successfully loaded ArchRProject!
#
## When we are ready to load this saved ArchRProject we use the loadArchRProject() object and provide the path to the folder containing the saved ArchRProject object.
#
message(paste("Loading ArchR project in project directory:", opt$project_directory))
proj <- loadArchRProject(path = opt$project_directory)
## Successfully loaded ArchRProject!

message(paste("Calling ArchR browser"))
ArchRBrowser(ArchRProj = proj, logFile = createLogFile(name = "ArchRBrowser", logDir = paste(opt$project_directory,"ArchRLogs", sep = "/")))
# This launches a dynamic genome browser session with a whole host of features including export of vectorized tracks for publication.



# Session Information
# This tutorial was run on the date specified below.

Sys.Date()
## [1] “2020-04-21”

# The sessionInfo() at run time was:

sessionInfo()
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
##
## Matrix products: default
## BLAS/LAPACK: /share/software/user/open/openblas/0.2.19/lib/libopenblasp-r0.2.19.so
##
## locale:
## [1] LC_CTYPE=en_US.UTF-8 LC_NUMERIC=C
## [3] LC_TIME=en_US.UTF-8 LC_COLLATE=en_US.UTF-8
## [5] LC_MONETARY=en_US.UTF-8 LC_MESSAGES=en_US.UTF-8
## [7] LC_PAPER=en_US.UTF-8 LC_NAME=C
## [9] LC_ADDRESS=C LC_TELEPHONE=C
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
##
## attached base packages:
## [1] parallel stats4 stats graphics grDevices utils datasets
## [8] methods base
##
## other attached packages:
## [1] ArchR_0.9.1 magrittr_1.5
## [3] rhdf5_2.30.1 Matrix_1.2-17
## [5] data.table_1.12.8 SummarizedExperiment_1.16.1
## [7] DelayedArray_0.12.2 BiocParallel_1.20.1
## [9] matrixStats_0.56.0 Biobase_2.46.0
## [11] GenomicRanges_1.38.0 GenomeInfoDb_1.22.1
## [13] IRanges_2.20.2 S4Vectors_0.24.3
## [15] BiocGenerics_0.32.0 ggplot2_3.2.1
##
## loaded via a namespace (and not attached):
## [1] Rcpp_1.0.4 pillar_1.4.3 compiler_3.6.1
## [4] XVector_0.26.0 tools_3.6.1 bitops_1.0-6
## [7] zlibbioc_1.32.0 lifecycle_0.1.0 tibble_2.1.3
## [10] gtable_0.3.0 lattice_0.20-38 pkgconfig_2.0.3
## [13] rlang_0.4.5 GenomeInfoDbData_1.2.2 withr_2.1.2
## [16] dplyr_0.8.4 grid_3.6.1 tidyselect_1.0.0
## [19] glue_1.4.0 R6_2.4.1 Rhdf5lib_1.8.0
## [22] purrr_0.3.3 scales_1.1.0 assertthat_0.2.1
## [25] colorspace_1.4-1 RCurl_1.98-1.1 lazyeval_0.2.2
## [28] munsell_0.5.0 crayon_1.3.4
