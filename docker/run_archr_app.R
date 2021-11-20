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
  ),
  make_option(
    c("-t", "--threads"),
    type="integer",
    default=2,
    help="Number of subprocesses/threads to use."
  )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#
if (is.null(opt$project_directory)){
  print_help(opt_parser)
  stop("--project_directory must be specified.", call.=FALSE)
}
message(paste("Project directory used is:",opt$project_directory))

## First, we load the ArchR library. If this fails, you have not properly installed
## ArchR and should revisit the installation instructions. We also recommend setting
## and remembering a known seed to facilitate replication of operations requiring randomization.
library(ArchR)

## Next, we set the default number of threads for parallelized operations in ArchR
## functions. You should change the value passed to threads to match the
## specifications of your local machine.
addArchRThreads(threads = opt$threads) 

script.dir <- getwd()
message(paste("Working directory used is:",script.dir))

# When we are ready to load this saved ArchRProject we use the loadArchRProject()
# object and provide the path to the folder containing the saved ArchRProject object.
message(paste("Loading ArchR project in project directory:", opt$project_directory))
proj <- loadArchRProject(path = opt$project_directory)
## Successfully loaded ArchRProject!

message(paste("Calling ArchR browser"))
ArchRBrowser(ArchRProj = proj, logFile = createLogFile(name = "ArchRBrowser", logDir = paste(opt$project_directory,"ArchRLogs", sep = "/")))
# This launches a dynamic genome browser session with a whole host of features including export of vectorized tracks for publication.

# Session Information
# This tutorial was run on the date specified below.
Sys.Date()

sessionInfo()


