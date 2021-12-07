# Script originally written by Lele Zhao at BDI
# Produces MAF.csv files necessary to run Tanya's HIV-phyloTSI-main program
# INPUT: 
# - csv containing path/to/sample/basefreq1, sampleID
# - output filename.
# USAGE:  Rscript generate_sample_MAF.R path_to_basefreq.csv MAF.csv (patstats)

  
# args <- c('~/git/phyloflows/MAF_and_BF/filelocs_with_selsamples.csv',
#          '~/Documents/2021/phyloTSI/phyloTSIinput/MAF_noPTY_211207.csv',
#          '~/Documents/2021/phyloTSI/PatStats/patstats_02_05_30_min_read_100_max_read_noPTY.csv')

args = commandArgs(trailingOnly=TRUE)

# 
path_file<-read.csv(args[1],header = F)
outfilename<-args[2]
if(!is.na(args[3]))
{
  # Do I want to process patstats file?
  
  # Take patstats file, check which individuals are in it, and 'grab' these from path_file
  file.path.patstats <- args[3]
  patstats.AIDs <- read.csv(file.path.patstats)
  patstats.AIDs <- unique(patstats.AIDs$host.id)
  # path_file <- path_file[tmp,]
}

if(ncol(path_file) > 2) # Can add the option of using FILE when running locally
{
  cat("More than 2 rows provided\nProcessing assuming the file 'filelocs_with_selsamples.csv'\n")
  if(!grepl('filelocs_with_selsamples',args[1])){warning('Base Frequency file is not filelocs_with_selsamples.csv')}
  
  if(is.na(args[3]))
  {
    stop('Error: Processing needs patstats file')
  }
  
  rownames(path_file) <- NULL
  colnames(path_file) <- path_file[1,]
  path_file <- path_file[-1,]
  path_file <- path_file[, c('FILE_HPC', 'AID')]
  path_file <- path_file[which(!is.na(path_file$AID)),]
  tmp <- which(path_file$AID %in% patstats.AIDs)
  tmp1 <- which(paste0('CNTRL-', path_file$AID) %in% patstats.AIDs)
  tmp1 <- path_file[tmp1,]
  tmp1$AID <- paste0('CNTRL-', tmp1$AID)
  path_file <- rbind(tmp1, path_file[tmp,])
  
  colnames(path_file) <- c('V1', 'V2')
  # What about entries AID that are not in path file?
  # sum(!patstats.AIDs %in% path_file$AID)
  # sum(! gsub('CNTRL-', '', patstats.AIDs) %in% gsub('CNTRL-', '', path_file$AID))
  # There are 360. Why are there more than the 351 we saw previously?
}


# Lele's part

MAF_matrix<-matrix(0,ncol=10000,nrow=(nrow(path_file)+1))
MAF_matrix[1,]<-seq(1,10000)
rownames(MAF_matrix)<-c('pos',path_file$V2)

N <- nrow(path_file)
for(i in 1:N){

  cat(paste0(i, ' out of ',N,': ', path_file$V2[i],'\n'))
  basefile <- path_file$V1[i]
  
  if (file.exists(basefile))
  {
    basefile<-read.csv(basefile,header = T)
    indexes_HXB2_pos<-which(basefile$Position.in.B.FR.83.HXB2_LAI_IIIB_BRU.K03455!='-')
    sample_HXB2_pos<-as.numeric(basefile[indexes_HXB2_pos,1])
    sample_MAFs<-apply(basefile[indexes_HXB2_pos,4:7], 1, function(x) 1-(max(x,na.rm = T)/sum(x,na.rm = T)))
    sample_MAFs<-as.numeric(gsub(NaN, 0,sample_MAFs))
    MAF_matrix[(1+i),sample_HXB2_pos]<-sample_MAFs 
  }else{
    MAF_matrix[(1+i),]<-NA 
  }
}

write.table(MAF_matrix,file = outfilename,quote = F,sep = ',',col.names = F) # consider col.names=T
outfilename <- gsub('.csv$','_noNA.csv',outfilename) # Consider removing
tmp <- which(is.na(MAF_matrix[, 2]))
write.table(MAF_matrix[-tmp, ],file = outfilename,quote = F,sep = ',',col.names = F)
