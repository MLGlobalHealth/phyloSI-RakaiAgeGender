cat("\n\n===== check_trees.R =====\n\n")
# The phyloscanner run for TSI estimates - make trees if not exist

# Preamble
# This script aims to make trees using each alignment 
# if the trees were not built in make_trees.R due to the time constrain.

# Load the required packages
library(data.table)
library(seqinr)
library(tidyverse)
library(dplyr)

#
# Define input arguments that can be changed by users
#
option_list <- list(
  optparse::make_option(
    c("-v", "--verbose"),
    action = "store_true",
    default = FALSE,
    help = "Print extra output [default]",
    dest = "verbose"
  ),
  optparse::make_option(
    "--seed",
    type = "integer",
    default = 42L,
    help = "Random number seed [default %default]",
    dest = "seed"
  ),
  optparse::make_option(
    "--save_data",
    action = "store_true",
    default = TRUE,
    help = "Save data [default]",
    dest = 'if_save_data'
  ),
  optparse::make_option(
    "--env_name",
    type = "character",
    default = 'phylor4',
    help = "Conda environment name [default]",
    dest = 'env_name'
  ),
  optparse::make_option(
    "--iqtree_method",
    type = "character",
    default = "GTR+F+R6",
    help = "Methods in IQTREE [default]",
    dest = 'iqtree_method'
  ),
  optparse::make_option(
    "--iqtree_root",
    type = "character",
    default = NULL,
    help = "Root in IQTREE [default]",
    dest = 'iqtree_root'
  ),
  optparse::make_option(
    "--pkg_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to package directory, used as long we don t build an R package [default]",
    dest = 'prj.dir'
  ),
  optparse::make_option(
    "--out_dir_base",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all output is stored [default]",
    dest = 'out.dir'
  ),
  optparse::make_option(
    "--prog_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where phyloscanner is stored [default]",
    dest = 'prog.dir'
  ),
  optparse::make_option(
    "--date",
    type = 'character',
    default = '2022-02-04',
    metavar = '"YYYY-MM-DD"',
    help = 'As of date to extract data from.  Defaults to today.',
    dest = 'date'
  ),
  optparse::make_option(
    "--walltime_idx",
    type = "integer",
    default = 2,
    help = "Indicator for amount of resources required by job. Values ranging from 1 (lala) to 3 (lala)",
    dest = "walltime_idx"
  )
  
)

args <-  optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# test
#

if(0){
  args <- list(
    out.dir="/rds/general/project/ratmann_deepseq_analyses/live/seroconverters3_alignXX",
    pkg.dir="/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software",
    iqtree_method="GTR+F+R6",
    env_name = 'phylostan',
    date = '2022-07-20',
    seed = 42,
    walltime_idx=1
  )
}

#
# Add constants that should not be changed by the user
#
max.per.run <- 950

# if prj.dir and out.dir are not manually set, default to here()
if (is.na(args$prj.dir))
{
  args$prj.dir <- here::here()
  args$out.dir <- here::here()
  args$prog.dir <- here::here()
}

# Set default output directories relative to out.dir
args$date <- gsub('-','_',args$date)
.f <- function(x)
{
        out <- file.path(args$out.dir, paste0(args$date, x))
        stopifnot(file.exists(out))
        out
}
args$out.dir.data <- .f('_phsc_input')
args$out.dir.work <- .f('_phsc_work')
args$out.dir.output <- .f('_phsc_output')

# Source functions
source(file.path(args$prj.dir, "utility.R"))

#
#	produce trees
# 


# iqtree option
iqtree.pr <- 'iqtree'
iqtree.args <- paste0('-m ',args$iqtree_method)

if(!is.null(args$iqtree_root)){
  iqtree.args	<- paste0(iqtree.args, ' -o ', args$iqtree_root)
}
if(!is.na(args$seed)){
  iqtree.args	<- paste0(iqtree.args, ' -seed ', args$seed)
}

# Load alignments
infiles <- list.files(args$out.dir.output, 
                      pattern='InWindow_(.*)?_v2.fasta$',
                      full.names=TRUE, recursive=TRUE)

# Extract info from name
infiles	<- data.table(FI=infiles)
infiles[, FO:= gsub('.fasta$','',FI)]
infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(FI)))]
infiles[, W_FROM:= as.integer(gsub('.*InWindow_([0-9]+)_.*','\\1',basename(FI)))]		
infiles[is.na(W_FROM),W_FROM:= as.integer(gsub('.*PositionsExcised_([0-9]+)_.*','\\1',basename(FI)))]
infiles[,PositionsExcised:=grepl('PositionsExcised',FI)]

# get files with excision if exist
setkey(infiles, PTY_RUN, W_FROM)
tmp <- infiles[,list(NUM=length(PositionsExcised)),by=c('PTY_RUN', 'W_FROM')]
infiles <- merge(infiles, tmp, by=c('PTY_RUN', 'W_FROM'))
infiles <- infiles[!(PositionsExcised==F & NUM==2),]

# check tree completed: treefiles.
infiles[,FO_NAME:=paste0(FO,'.treefile')]
infiles[,FO_EXIST:=file.exists(FO_NAME)]
infiles <- infiles[FO_EXIST==FALSE]
if(args$verbose){
  cat('These trees were not built sucessfully...\n')
  print(head(infiles))
  cat('Build again...\n')
}

# Set up jobs (/1 is an artefact of XX's code, who divided by 4)
df <- infiles[, list(CMD=cmd.iqtree(FI, outfile=FO, pr=iqtree.pr, pr.args=iqtree.args)), by=c('PTY_RUN','W_FROM')]
df[, ID:=ceiling(seq_len(nrow(df))/1)]
df <- df[, list(CMD=paste(CMD, collapse='\n',sep='')), by='ID']

# Create PBS job array
if(nrow(df) > max.per.run){
  df[, CASE_ID:= rep(1:max.per.run,times=ceiling(nrow(df)/max.per.run))[1:nrow(df)]]
  df[, JOB_ID:= rep(1:ceiling(nrow(df)/max.per.run),each=max.per.run)[1:nrow(df)]]
  df[, ID:=NULL]
}else{
  setnames(df, 'ID', 'CASE_ID')
  df[, JOB_ID:=1]
}

# Chose PBS specifications according to PBS index
# Idea is that this script can be run multiple times with increasing res reqs

list(
  `1`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 4 , hpc.mem = "2gb" ,hpc.q = NA),
  `2`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 23, hpc.mem = "2gb" ,hpc.q = NA),
  `3`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 71, hpc.mem = "63gb",hpc.q = NA)
) -> pbs_headers

tmp <-pbs_headers[[args$walltime_idx]]
list2env(tmp,globalenv())

#
# Write cmd and submit
#

# header
pbshead	<- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=hpc.select, hpc.walltime=hpc.walltime, hpc.q=hpc.q, hpc.mem=hpc.mem,  hpc.nproc=hpc.nproc, hpc.load=NULL)

for (i in 1:df[, max(JOB_ID)]) {

  tmp<-df[JOB_ID==i,]

  hpc.array <- nrow(tmp)

  # Here no pbsJ_bool?
  cmd<-tmp[, list(CASE=paste0(CASE_ID,')\n',CMD,';;\n')), by='CASE_ID']
  cmd<-cmd[, paste0('case $PBS_ARRAY_INDEX in\n',paste0(CASE, collapse=''),'esac')]		
  tmp <- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')	
  tmp <- paste0(tmp,'\n module load anaconda3/personal \n source activate ',args$env_name)
  cmd<-paste(tmp,cmd,sep='\n')	
  
  #	Prepare sh file  
  date <- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
  date <- gsub(':','_',date)
  outfile <- paste("srx",paste0('job',i), date,'sh',sep='.')
  outfile <- file.path(args$out.dir.work, outfile)
  cat(cmd, file=outfile)

  # Change directory...
  cmd <- paste('cd', dirname(outfile))
  cat(system(cmd, intern=TRUE))

  # ... and submit job!
  cmd<-paste("qsub", outfile)
  cat(cmd)
  cat(system(cmd, intern= TRUE))
}

