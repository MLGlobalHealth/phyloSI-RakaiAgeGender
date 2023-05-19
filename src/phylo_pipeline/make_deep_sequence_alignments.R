# The phyloscanner run - make alignments

# Preamble
# The set of scripts aims to run phyloscanner.
# This script aims to make alignments in each potential transmission network.
# TODO: check whether problematic_windows is fixed...
# Also qsub.next...

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
    default = TRUE,
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
    "--windows_start",
    type = "integer",
    default = 800,
    help = "Left extremity of leftmost window[default]",
    dest = "windows_start"
  ),
  optparse::make_option(
    "--windows_end",
    type = "integer",
    default = 9175,
    help = "Left extremity of rightmost windows[default]",
    dest = "windows_end"
  ),
  optparse::make_option(
    "--sliding_width",
    type = "integer",
    default = NA_integer_,
    help = "Sliding width [default %default]",
    dest = "sliding_width"
  ),
  optparse::make_option(
    "--window_size",
    type = "integer",
    default = 250L,
    help = "Window size [default %default]",
    dest = "window_size"
  ),
  optparse::make_option(
    "--window_cutoff",
    type = "numeric",
    default = NA,
    help = "Cutoff of proportion of windows indicating close relationship [default %default]",
    dest = "window_cutoff"
  ),
  optparse::make_option(
    "--n_control",
    type = "integer",
    default = 0,
    help = "Number of controls added [default %default]",
    dest = "n_control"
  ),
  optparse::make_option(
    "--cluster_size",
    type = "integer",
    default = 100L,
    help = "Maximum cluster size [default %default]",
    dest = "cluster_size"
  ),
  optparse::make_option(
    "--save_data",
    action = "store_true",
    default = TRUE,
    help = "Save data [default]",
    dest = 'if_save_data'
  ),
  optparse::make_option(
    "--pkg_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to package directory, used as long we don t build an R package [default]",
    dest = 'pkg.dir'
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
    "--reference",
    type = "character",
    default = NA_character_,
    help = "path to reference consensus sequences .fasta file required for the analysis. Basename is sufficient if file is 'standard'.",
    dest = 'reference'
  ),
  optparse::make_option(
    "--date",
    type = 'character',
    default = as.character(Sys.Date()),
    metavar = '"YYYY-MM-DD"',
    help = 'As of date to extract data from.  Defaults to today.',
    dest = 'date'
  ),
  optparse::make_option(
    "--tsi_analysis",
    type = 'logical',
    default = FALSE,  
    help = 'Indicator on whether we want to perform a Time Since Infection analysis[default]. Close to being deprecated',
    dest = 'tsi_analysis'
  ),
  optparse::make_option(
    "--mafft",
    type = 'character',
    default = '--globalpair --maxiterate 1000',  
    help = 'options for alignment program mafft', 
    dest = 'mafft.opt'
  ),
  optparse::make_option(
    "--rm_vloops",
    action = "store_true",
    default = TRUE,
    help = "Indicator on whether to avoid alignments on vloop region.[default]",
    dest = 'rm_vloops'
  ),
  optparse::make_option(
    "--controller",
    type = "character",
    default = NA_character_, 
    help = "Path to sh script irecting the full anysis",
    dest = 'controller'
  ),
  optparse::make_option(
    "--walltime_idx",
    type = "integer",
    default = 2,
    help = "Indicator for amount of resources required by job. Values ranging from 1 (lala) to 3 (lala)",
    dest = "walltime_idx"
  )
)
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# stop if arguments are not provided:
if( is.na(args$sliding_width) ) stop('No sliding_width provided')

#
# Helpers
#

.percent <- function(x){paste0(format(round(x*100, 2), nsmall=2), '%') }

.check.existing.outputs <- function(regex, outdir, pattern='')
{
  # returns a vector of TRUE or FALSE based on whether regex is found in outdir
  if( ! file.exists(outdir) ){return(rep(NA_character_, length(regex)))}
  
  files <- list.files(outdir, pattern=pattern)
  # .f <- function(rgx) any(grepl(x=files,rgx))
  .f <- function(rgx) grep(x=files,rgx, value=T)[1]
  sapply(regex, .f) -> tmp
  unlist(tmp)
} 

.remove_problematic_windows <- function(DT)
{
  # check if file problematic_windows.csv exists
  tmp <- list.files(args$out.dir.output, pattern='problematic_windows.csv', 
                    recursive = T, full.names = T)
  if( length(tmp) == 0 )
  {
    cat('No problematic windows files found\n')
    return(DT)
  }
  names(tmp) <- basename(dirname(tmp))
  dprob <- lapply(tmp, fread)
  dprob <- rbindlist(dprob, idcol = 'PTY_RUN')
  names(dprob) <- c('PTY_RUN', 'FROM', 'TO', 'PBS_ARRAY_INDEX', 'SH', 'QUEUE')
  dprob[,  PTY_RUN:=as.integer(gsub('[A-z]|_', '', PTY_RUN))]
  
  # dprob <- dprob[!SH %in% c('readali.job1.Sat_Jul_23_205803_2022.sh', 'readali.job1.Sun_Jul_24_133928_2022.sh')]
  
  # This below is more for debugging purposes
  # names of logfiles
  .f <- function(x)
  {
    # finds prefix for LOG, needs only PBS_ARRAY_INDEX
    dir <- file.path(args$out.dir.work, gsub('.sh$','',x))
    fil <- list.files(dir,  full.names = T)
    pre <- unique(gsub('\\.[0-9]+$','.',fil))
    if(length(pre) == 1) return(pre)
    NA_character_
  }
  dprob[, LOG:=paste0(.f(SH), PBS_ARRAY_INDEX), by='SH']
  
  # compare with outputs 
  merge(
    DT[, .(W_FROM, PTY_RUN, OUT1, OUT2)],
    dprob[, .(W_FROM = FROM, PTY_RUN, LOG)],
    all.x=TRUE,
    by=c('W_FROM', 'PTY_RUN')
  ) -> dlogs
  
  tmp <- dlogs[ is.na(LOG), .(W_FROM, PTY_RUN, OUT2) ]
  # stopifnot(tmp[, all(is.na(OUT2))])
  tmp$OUT2 <- NULL
  
  merge(DT, tmp, by=c('W_FROM', 'PTY_RUN'))
}

.modify_cmd_for_existing_v1 <- function(cmd, outdir, out1)
{
  if(length(out1) == 0){ stop('error: empty input to .modify_cmd_for_existing_v1')}
  # cmd <- pty.c[, CMD[1]]
  out1 <- basename(out1)
  out2 <- gsub('\\.fasta', '\\_v2.fasta', out1)
  lines <- strsplit(cmd, '\n')[[1]]
  
  # remove useless parts
  rm1_start <- grep('phyloscanner_make_trees', lines)
  rm1_end <- grep('Performing realignment', lines) - 1 
  rm1 <- rm1_start:rm1_end
  rm2 <- grep("(?=.*mv.*)(?=.*AlignedReads.*)", lines, perl = TRUE)
  rm2 <- (rm2-1):(rm2+1)
  rm3 <- grep('problematic_windows.csv', lines)
  lines <- lines[ -c(rm1, rm2, rm3)]
  
  # copy out1 to work dir
  cp_pos <- grep('^cp', lines)[1]
  cp <- lines[cp_pos]
  cp <- gsub('^cp "(.*?)" (.*?)$', paste0('cp "', file.path(outdir, out1)  ,'" \\2'), cp)
  lines[cp_pos] <- cp
  
  # substitute AlignedReads* with the name of our file
  mafft_pos <- grep('mafft', lines)[1]
  mafft <- lines[mafft_pos]
  mafft2 <- gsub('\\\t', '', mafft)
  mafft2 <- gsub('\\$file', out1, mafft2)
  mafft2 <- gsub('\\$\\{file//.fasta/_v2.fasta\\}', out2, mafft2)
  lines[mafft_pos] <- mafft2
  rm4 <- c(mafft_pos + 1, mafft_pos -1)
  lines <- lines[ -rm4]
  
  cmd <- paste0(lines, collapse='\n')
  return(cmd)
}

.write.job <- function(DT)
{
  # Define PBS header for job scheduler
  pbshead <- cmd.hpcwrapper.cx1.ic.ac.uk(
    hpc.select = hpc.select,
    hpc.nproc = hpc.nproc,
    hpc.walltime = hpc.walltime,
    hpc.q = hpc.q,
    hpc.mem = hpc.mem,
    hpc.array = DT[, max(CASE_ID)],
    hpc.load = "module load intel-suite/2015.1 mpi raxml/8.2.9 mafft/7 anaconda/2.3.0 samtools"
  )
  
  cmd <- DT[, list(CASE = paste0(CASE_ID, ')\n', CMD, ';;\n')), by = 'CASE_ID']
  cmd <-    cmd[, paste0('case $PBS_ARRAY_INDEX in\n',
                         paste0(CASE, collapse = ''),
                         'esac')]
  cmd <- paste(pbshead, cmd, sep = '\n')
  cmd
}

#
# test
#
if(0){
  args <- list(
    verbose = T,
    seed = 42,
    windows_start=550L,
    windows_end=9500L,
    sliding_width = 25L,
    window_size = 250L,
    window_cutoff = NA,
    n_control = 0,
    cluster_size = 50,
    if_save_data = T,
    date = '2022-08-22',
    out.dir = "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI/",
    pkg.dir = "/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software",
    prog.dir = "/rds/general/user/ab1820/home/git/phyloscanner",
    reference = 'ConsensusGenomes.fasta',
    tsi_analysis=FALSE,
    rm_vloops=FALSE,
    mafft.opt='--globalpair --maxiterate 1000',
    controller='/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/runall_TSI_pairs2.sh',
    walltime_idx = 2
  )
}

# Source functions
source(file.path(args$pkg.dir, "utility.R"))

# I can run at most 1000 simultaneous jobs on the short q.
list(
  `1`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 7 , hpc.mem = "4gb",hpc.q = NA, max.per.run=1000),
  # `2`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 7, hpc.mem = "4gb",hpc.q = NA, max.per.run=1000),
  `2`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 24, hpc.mem = "64gb",hpc.q = NA, max.per.run=500),
  `3`= list(hpc.select = 1, hpc.nproc = 1, hpc.walltime = 71, hpc.mem = "64gb",hpc.q = NA, max.per.run=500)
) -> pbs_headers
tmp <- pbs_headers[[args$walltime_idx]]
cat('selected the following PBS specifications:\n')
print(tmp)
invisible(list2env(tmp,globalenv()))

#
# Add constants that should not be changed by the user
#
dir.data <-  '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
dir.analyses <- '/rds/general/project/ratmann_deepseq_analyses/live'
dir.net <- file.path(args$out.dir, "potential_network")

ifelse(
  !is.na(args$window_cutoff),
  paste0('_n_control_', args$n_control), ''
) -> tmp

infile.runs <- file.path(
  args$out.dir,
  paste0(
    'phscinput_runs_clusize_', args$cluster_size,
    '_ncontrol_', args$n_control,
    tmp,'.rds'
  )
)
stopifnot(file.exists(infile.runs))

# Set default output directories relative to out.dir
args$date <- gsub('-','_',args$date)

.f <- function(x)  
{
  dir <- file.path(args$out.dir, paste0(args$date, x))
  
  # If date has been passed as par, (ie. args$date != today's date) 
  # check that everything exists...
  cnd <-  as.Date(args$date, format='%Y_%m_%d') != as.character(Sys.Date()) 
  
  if(!dir.exists(dir))
  {
    if ( cnd )
      stop('No existing directories prefixed by the input --date')
    dir.create(dir)
  }
  dir
}
args$out.dir.data <- .f('_phsc_input')
args$out.dir.work <- .f('_phsc_work')
args$out.dir.output <- .f('_phsc_output')


# Look for RDS file containing subjobs CMDS, if can't find, KEEP GOING
cmds.path <- file.path(args$out.dir.work, 'align_commands.rds')

if(file.exists(cmds.path))
{
  # Load commands
  cat('Cleaning work directory...\n')
  move.logs(args$out.dir.work)
  pty.c <- readRDS(cmds.path)
  
}else{
  
  cat('Writing .rds of commands...\n')
  # Write commands
  
  # Copy files into input folder
  tmp  <- file.path(dir.data, 'PANGEA2_RCCS/phyloscanner_input_data')
  tmp1 <- file.path(dir.data, 'PANGEA2_RCCS/220419_reference_set_for_PARTNERS_mafft.fasta')
  tmp <- c(list.files( tmp, full.name=T), tmp1)
  file.copy(
    tmp,  args$out.dir.data,
    overwrite = T,
    recursive = FALSE,
    copy.mode = TRUE
  )
  
  # Set consensus sequences.
  # (default consensus/reference for tsi and pair analyses if no arg is passed)
  # not really sure whether oneeach is needed anywhere
  
  infile.consensus <- args$reference 
  
  if(is.na(infile.consensus))
  {
    # if NA, look in args$out.dir.data
    tmp <- ifelse(args$tsi_analysis,
                  yes='220419_reference_set_for_PARTNERS_mafft.fasta',
                  no='ConsensusGenomes.fasta')
    infile.consensus <- file.path(args$out.dir.data, tmp)
  }else{
    # If does not exists, look within args$out.dir.data
    if(! file.exists(infile.consensus) )
    {
      infile.consensus <- file.path(args$out.dir.data, basename(infile.consensus))
    }
    
  }
  stopifnot(file.exists(infile.consensus))
  
  # Load sequences and remove duplicates if existing
  pty.runs <- data.table(readRDS(infile.runs))
  if ('ID_TYPE' %in% colnames(pty.runs)) {
    setorder(pty.runs, PTY_RUN, -ID_TYPE, UNIT_ID)
  } else{
    setorder(pty.runs, PTY_RUN, UNIT_ID)
  }
  
  tmp <- pty.runs[, list(idx=duplicated(SAMPLE_ID)), by = 'PTY_RUN']
  tmp <- tmp[, which(idx)]
  if(length(tmp)!=0){
    pty.runs <- pty.runs[-tmp, ]
  }
  
  tmp <- pty.runs[, uniqueN(SAMPLE_ID) == .N, by='PTY_RUN' ]
  stopifnot( all(tmp$V1) )
  
  # Load backgrounds and extract HXB2 for pairwise MSA. 
  consensus_seq <- seqinr::read.fasta(infile.consensus)
  consensus_seq_names <- names(consensus_seq)
  hxb2 <- grep('HXB2', names(consensus_seq), value = T)
  hxb2_seq <- consensus_seq[[hxb2]]
  # removed root_seq as it seemed useless
  
  # Change the format
  if ('ID_TYPE' %in% colnames(pty.runs)) {
    pty.runs[ID_TYPE == 'control', UNIT_ID := paste0('CNTRL-', UNIT_ID)]
    pty.runs[ID_TYPE == 'control', RENAME_ID := paste0('CNTRL-', RENAME_ID)]
  }
  pty.runs[, BAM := paste0(dir.data, SAMPLE_ID, '.bam')]
  pty.runs[, REF := paste0(dir.data, SAMPLE_ID, '_ref.fasta')]
  setkey(pty.runs, PTY_RUN, RENAME_ID)
  
  # 
  # Set the alignment options
  #
  
  # MAFFT: reformat options so they are readily pasted in sh command.
  # args$mafft.opt <- '"mafft --globalpair --maxiterate 1000"'
  args$mafft.opt <- gsub('mafft|"', '', args$mafft.opt)
  args$mafft.opt <- paste0('"mafft ', args$mafft.opt, '"')
  # cat(args$mafft.opt)
  
  # GENOMIC WINDOWS
  # excision.default will excise more positions, atm I group together with remove vloops
  
  stopifnot(args$windows_start <= args$window_end)
  ptyi <- seq(args$windows_start, args$windows_end, by=args$sliding_width)
  
  # by remove loops, I also perform different 
  if(args$rm_vloops)
  {
    ptyi <- c(ptyi[ptyi <= 6615 - args$window_size], 6825, 6850, ptyi[ptyi >= 7636])
    excision.default.bool <- TRUE
  }else{
    excision.default.bool <- FALSE
  }
  
  # Now write command
  write.pty.command <- function(w_from)
  {
    pty.args <- list(
      prog.pty = file.path(args$prog.dir, "phyloscanner_make_trees.py"),
      prog.mafft = args$mafft.opt,
      data.dir = args$out.dir.data,
      work.dir = args$out.dir.work,
      out.dir = args$out.dir.output,
      alignments.file = infile.consensus,
      # alignments.root = root_seq, # is this even doing anything?
      alignments.pairwise.to = hxb2,
      window.automatic = '',
      merge.threshold = 0,
      min.read.count = 1,
      quality.trim.ends = 23,
      min.internal.quality = 23,
      merge.paired.reads = TRUE,
      discard.improper.pairs = TRUE,
      no.trees = TRUE,
      dont.check.duplicates = FALSE,
      dont.check.recombination = TRUE,
      num.bootstraps = 1,
      all.bootstrap.trees = TRUE,
      strip.max.len = 350,
      min.ureads.individual = NA,
      win = c(w_from, w_from + args$window_size, args$sliding_width, args$window_size),
      keep.overhangs = FALSE,
      mem.save = 0,
      verbose = TRUE,
      select = NA,
      default.coord = excision.default.bool,
      realignment = TRUE
    )
    pty.c <- phsc.cmd.phyloscanner.multi(pty.runs, pty.args)
    pty.c[, W_FROM := w_from ]
    pty.c
  }
  
  pty.c	<- lapply(ptyi, write.pty.command)
  pty.c	<- do.call('rbind', pty.c)
  setkey(pty.c, PTY_RUN, W_FROM)
  stopifnot(pty.c[, .N, by=c('PTY_RUN', 'W_FROM')][, all(N == 1)])
  
  pty.c[ , OUTDIR := file.path(args$out.dir.output, paste0('ptyr', PTY_RUN, '_trees')) ]
  pty.c[, OUT_REGEX_V1 := paste0( W_FROM, '_to_', W_FROM + args$window_size - 1,  '.fasta') ]
  pty.c[, OUT_REGEX_V2 := paste0( W_FROM, '_to_', W_FROM + args$window_size - 1,  '_v2.fasta') ]
  
  saveRDS(pty.c, cmds.path)
}

# Now check which outputs are not ready yet and need re-running
cols <- c('OUT1', 'OUT2')
pty.c[, (cols) := lapply(.SD, .check.existing.outputs, outdir=OUTDIR, pattern='fasta' ) , by=OUTDIR, .SDcols= names(pty.c) %like% 'REGEX']

pty.c[, cat(sum(!is.na(OUT2)),
            '(', .percent(mean(!is.na(OUT2))),')', 
            'of desired outputs were generated in the previous runs\n')]
pty.c[, cat(sum(!is.na(OUT1)),
            '(', .percent(mean(!is.na(OUT1))),')', 
            'of intermediary outputs were generated in the previous runs\n')]

# Select jobs with missing final outcome.
# remove problematic windows as informed by the shell scripts output
pty.c <- pty.c[ is.na(OUT2)]
pty.c <- .remove_problematic_windows(pty.c)

# if intermediary output is there, change CMD.
if(pty.c[ !is.na(OUT1), .N > 0])
{
  pty.c[ ! is.na(OUT1) ,  .modify_cmd_for_existing_v1(cmd=CMD, outdir=OUTDIR, out1=OUT1)  , by='OUT1' ]
}

#
# IF ALL OUTPUTS WERE CREATED, GO TO NEXT STEP IN THE ANALYSIS
#

if( nrow(pty.c) == 0 )
{
  # ISN'T THIS BEAUTIFUL?
  qsub.next.step(file=args$controller,
                 next_step='btr', 
                 res=1,
                 redo=0
  )
  stop('Alignment step completed, submitted the following task')
}

#
# NOW THAT WE HAVE THE IND CMDs to RUN, WE NEED TO AGGREGATE THEM
#

n_jobs <- ceiling( nrow(pty.c) / max.per.run)
idx <- 1:nrow(pty.c)

pty.c[, CASE_ID := rep(1:max.per.run, times = n_jobs)[idx] ]
pty.c[, JOB_ID := rep(1:n_jobs, each = max.per.run)[idx] ]

# Write and submit:
djob <- pty.c[, .(CMD=.write.job(.SD)), by=JOB_ID]
if(args$walltime_idx == 3){
  # Assign one every 5 jobs to pqeelab 
  djob[, Q := {
    idx <- seq(5, .N, 5)
    z <- rep('', .N)
    z[idx] <- 'pqeelab'
    z
    }]
}


ids <- djob[, list(ID=.store.and.submit(.SD, prefix='readali')), by=JOB_ID, .SDcols=names(djob)]
ids <- as.character(ids$ID)
cat('Submitted job ids are:', ids, '...\n')

# qsub alignment step again, to check whether everything has run...
qsub.next.step(file=args$controller,
               ids=ids, 
               next_step='ali', 
               res=args$walltime_idx + 1, 
               redo=0
)
