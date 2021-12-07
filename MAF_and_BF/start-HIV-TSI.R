# command line parsing if any
args_line <-  as.list(commandArgs(trailingOnly=TRUE))
if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-l')	
  stopifnot(args_line[[3]]=='-m')
  stopifnot(args_line[[5]]=='-f')
  stopifnot(args_line[[7]]=='-outdir')
  file.path.locs <- args_line[[2]]
  # file.path.maf <- args_line[[4]]
  file.path.patstats <- as.integer(args_line[[6]])
  outdir <- as.integer(args_line[[8]])
}else{
 
  git.dir <- "~/git/phyloflows"
  outdir <- '~/Documents/2021/phyloTSI'
  file.path.locs <- file.path(git.dir, 'MAF_and_BF/filelocs_with_selsamples.csv')
  file.path.patstats <- file.path(outdir, 'PatStats', 'patstats_02_05_30_min_read_100_max_read_noPTY.csv')
}

if(!dir.exists('/home/andrea') )
{
  usr <- '/rds/general/user/ab1820/home'
  git.dir <- file.path(usr, "/git/phyloflows")
  outdir <- file.path(usr, '/projects/2021/phyloTSI')
  file.path.locs <- file.path(outdir, 'MAF_and_BF/filelocs_with_selsamples.csv')
  file.path.patstats <- file.path(outdir, 'PatStats', 'patstats_02_05_30_min_read_100_max_read_noPTY.csv')
}


date <- format(Sys.time(), "%y%m%d")


########################################
# Header Function
########################################

make.PBS.header <- function(hpc.walltime=47, hpc.select=1, hpc.nproc=1, hpc.mem= "6gb", hpc.load= "module load anaconda3/personal\nsource activate hivphylotsi", hpc.q=NA, hpc.array=1 )
{	
  pbshead <- "#!/bin/sh"
  
  if(!dir.exists('/home/andrea'))
  {
    tmp <- paste("#PBS -l walltime=", hpc.walltime, ":59:00", sep = "")
    pbshead <- paste(pbshead, tmp, sep = "\n")
    tmp <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":ompthreads=", hpc.nproc,":mem=", hpc.mem, sep = "")	
    pbshead <- paste(pbshead, tmp, sep = "\n")
    pbshead <- paste(pbshead, "#PBS -j oe", sep = "\n")
  }else{
    hpc.load="conda activate hivphylotsi"
  }
  
  if(hpc.array>1)
  {
    pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep='')
  }				
  if(!is.na(hpc.q))
  {
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  }		
  pbshead	<- paste(pbshead, hpc.load, '', sep = "\n")
  pbshead
}

########################################
# Body Function
########################################

make_single_cmd <- function(file.path.patstats){
  
  # Run Lele's script
  file.path.maf <- gsub('^(.*)read_(.*).csv$','\\2', basename(file.path.patstats))
  file.path.maf <- paste0('MAF_',file.path.maf, '_', format(Sys.time(), "%y%m%d"), '.csv')
  file.path.maf <- file.path(outdir, 'phyloTSIinput', file.path.maf)
  
  
  cmd <- paste0("\necho ----- STEP 1: run Lele s script -----\n" )
  cmd <- paste0(cmd, 'Rscript generate_sample_MAF.R')
  cmd <- paste0(cmd, " '", file.path.locs, "' '", file.path.maf, "' '", file.path.patstats, "'\n")
  
  
  # 3  - Run Tanya's program
  file.path.output <- gsub('^(.*)read_(.*).csv$','\\2', basename(file.path.patstats))
  file.path.output <- paste0('TSI_',file.path.output, '_', format(Sys.time(), "%y%m%d"), '.csv')
  file.path.output <- file.path(outdir, 'phyloTSIoutput', file.path.output)
  
  cmd <- paste0(cmd, "\necho  ----- STEP 2: run Tanya s script  -----\n")
  cmd <- paste0(cmd, "python /rds/general/user/ab1820/home/git/HIV-phyloTSI/HIVPhyloTSI.py ",
                "-d /rds/general/user/ab1820/home/git/HIV-phyloTSI/Model ",
                "-p '", file.path.patstats, "' ",
                "-m '", file.path.maf, "' ",
                "-o '", file.path.output, "'\n")
  
  return(cmd)
}

########################################
#	put everything together and submit job
########################################

cmd <- ''
if(dir.exists(file.path.patstats))
{
  tmp <- list.files(file.path.patstats)
  N <- length(tmp)
  cmd <- paste0(cmd, 'case $PBS_ARRAY_INDEX in\n')
  
  for(i in 1:N){ # CHANGE HERE
    cmd <- paste0(cmd, i, ')\n')
    cmd <- paste0(cmd, make_single_cmd(tmp[i]))
    cmd <- paste0(cmd, ';;\n')
  }
  
  cmd <- paste0(cmd, 'esac')
  
}else{
  N <- 1
  cmd <- paste0(cmd, make_single_cmd(file.path.patstats)) 
}

cmd <- paste(make.PBS.header(hpc.array = N), cmd)

jobfile <- gsub(':','',paste("phyTSI",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),'sh',sep='.'))
jobfile <- file.path(outdir, jobfile)
cat(cmd, file=jobfile)

cmd <- ifelse(!dir.exists('/home/andrea'),
              paste("qsub", jobfile),
              paste("source", jobfile))

cat(cmd)
cat(system(cmd, intern= TRUE))


