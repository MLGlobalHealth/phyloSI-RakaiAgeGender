# Script originally written by Lele Zhao at BDI
# Produces MAF.csv files necessary to run Tanya's HIV-phyloTSI-main program
# INPUT: 
# - csv containing path/to/sample/basefreq1, sampleID
# - output filename.
# USAGE:  Rscript generate_sample_MAF.R path_to_basefreq.csv MAF.csv

args = commandArgs(trailingOnly=TRUE)

path_file<-read.csv(args[1],header = F)
outfilename<-args[2]

MAF_matrix<-matrix(0,ncol=10000,nrow=(nrow(path_file)+1))
MAF_matrix[1,]<-seq(1,10000)
rownames(MAF_matrix)<-c('pos',path_file$V2)

N <- nrow(path_file)
for(i in 1:N){
  cat(paste0(i, ' out of ',N,'\n'))
  basefile<-read.csv(path_file$V1[i],header = T)
  sample_name<-path_file$V2[i]
  indexes_HXB2_pos<-which(basefile$Position.in.B.FR.83.HXB2_LAI_IIIB_BRU.K03455!='-')
  sample_HXB2_pos<-as.numeric(basefile[indexes_HXB2_pos,1])
  sample_MAFs<-apply(basefile[indexes_HXB2_pos,4:7], 1, function(x) 1-(max(x,na.rm = T)/sum(x,na.rm = T)))
  sample_MAFs<-as.numeric(gsub(NaN, 0,sample_MAFs))
  MAF_matrix[(1+i),sample_HXB2_pos]<-sample_MAFs
 
}

write.table(MAF_matrix,file = outfilename,quote = F,sep = ',',col.names = F)
