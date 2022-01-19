# WANT: bash script that runs Tanya s script for each individual PTY run
# This means:
#     - Get the PatStats file. 
#     - Use R to subset it according to PTY_identifier 
#     - Use subsetted PS + MAF to run Tanya's
#     - Delete subsetted PS

args_line <-  as.list(commandArgs(trailingOnly=TRUE))

git.phy <- "~/git/phyloflows/MAF_and_BF"
git.tsi <- "~/git/HIV-phyloTSI-main"

if(dir.exists('/home/andrea')){LOCAL <- 1}else{LOCAL <- 0}

if(length(args_line) > 0) 
{
  stopifnot(args_line[[1]]=='-p')	
  stopifnot(args_line[[3]]=='-m')
  stopifnot(args_line[[5]]=='-o')
  file.path.patstats <- args_line[[2]]
  file.path.maf <- args_line[[4]]
  outdir <- args_line[[6]]

}else{
  
  if(LOCAL)
  {
    outdir <- '~/Documents/2021/phyloTSI/outputs'
    file.path.patstats <- "~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/patstats_02_05_30_min_read_100_max_read.csv"
    file.path.maf <- "~/Documents/Box/2021/phyloTSI/MAF_20220118.csv"
  }else{
    outdir <- '~/Documents/2021/phyloTSI/outputs'
    file.path.patstats <- "~/Documents/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/patstats_02_05_30_min_read_100_max_read.csv"
    file.path.maf <- "~/Documents/2021/phyloTSI/MAF_20220118.csv"
  }
  
}



########################################
# Body function
########################################
# Loop assuminig PTY_idx does not have 'holes' from 1 to max element

PTY_idx = 1:409 #max(patstats$PTY_RUN) 
output_pty <- file.path(outdir,paste0('tsi_', PTY_idx, '.csv'))

# Subset Patstats according to PTY_RUN
tmp_pat <- file.path(outdir, paste0('tmp_pat',PTY_idx,'.csv'))
tmp_maf <- file.path(outdir, paste0('tmp_maf',PTY_idx,'.csv'))

cmd <- 'echo "----- Create patstats ------" \n'
cmd <- paste0(cmd, 'head -q -n 1 ',file.path.patstats,"| cut -d , -f 1,2,3 --complement | tr -d \'\"\' > ", tmp_pat, '\n')
cmd <- paste0(cmd, 'cut -d , -f 1,2 --complement ', file.path.patstats, "| tr -d \'\"\' ") 
# cmd <- paste0(cmd, "| awk -F , '$1 == ", PTY, "' | cut -d , -f 1 --complement >> ",tmp_pat," \n")
cmd <- paste0(cmd, "| awk -F , '$1 == ", PTY_idx, "' | cut -d , -f 1 --complement | sed 's/\\.5//' >> ",tmp_pat," \n")
# there is a strange mistake in the maf, where x.positions are NOT integers, rather they end by.5 ...

cmd <- paste0(cmd, 'echo "------ Create MAF -------" \n')
cmd <- paste0(cmd, "Rscript ", file.path(git.phy, 'run_HIV_TSI.R'),
              ' -idx ', PTY_idx,
              ' -p ', tmp_pat,
              ' -m ', file.path.maf, '\n')

# Run Tanya
cmd <- paste0(cmd, 'echo "------ Run Tanya -------" \n')
cmd <- paste0(cmd, paste('python',file.path(git.tsi,'HIVPhyloTSI.py'),
                         '-d',file.path(git.tsi,'Model'),
                         '-p', tmp_pat,
                         '-m', tmp_maf,
                         '-o', output_pty, '\n'))

cmd <- paste0(cmd, 'rm ', tmp_pat, '\n')
cmd <- paste0(cmd, 'rm ', tmp_maf, '\n')

if(length(PTY_idx>1))
{
  if(!LOCAL){cmd <- paste0(PTY_idx, ')',cmd)}
  
  cmd <- paste0('\necho ',PTY_idx, '\n', cmd, ';;\n')
  cmd <- paste0(cmd, collapse = '\n')
  
  if(!LOCAL){cmd <- paste0('case $PBS_ARRAY_INDEX in\n', cmd, '\nesac\n')}else{cmd <- gsub(';;','',cmd)}
}


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
}

########################################
# Combine
########################################

cmd0 <- make.PBS.header()
cmd <- paste0(cmd0, cmd)
date <- format(Sys.Date(),'%Y%m%d')
cat(cmd, file = paste0('run_HIV_TSI_',date,'.sh'))
