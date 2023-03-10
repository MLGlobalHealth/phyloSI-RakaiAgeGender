cat('\n\n=====  TSI_run_predictions.R =====\n\n')

library(lubridate)
library(data.table)

option_list <- list(
  optparse::make_option(
    "--out_dir_base",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all output is stored [default]",
    dest = 'out.dir'
  ),
  optparse::make_option(
    "--pkg_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to package directory, used as long we don t build an R package [default]",
    dest = 'pkg.dir'
  ),
  optparse::make_option(
    "--relationship_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to directory containing analyse_trees.R results", 
    dest = 'rel.dir'
  ),
  optparse::make_option(
    "--TSI_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to HIV-phylo-TSI-main repository", 
    dest = 'TSI.dir'
  ),
  optparse::make_option(
    "--env_name",
    type = "character",
    default = 'hivphylotsi',
    help = "Conda environment name to run HIV-phylo-TSI analyses [default] ",
    dest = 'env_name'
  ),
  optparse::make_option(
    "--input_samples",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to input samples rds containing PANGEA_IDs and RENAME_IDs", 
    dest = 'phsc.samples'
  ),
  optparse::make_option(
    "--walltime",
    type = "integer",
    default = 4L,
    help = "Job walltime (hours) [default]",
    dest = 'walltime'
  ),
  optparse::make_option(
    "--memory",
    type = "character",
    default = "2gb",
    help = "Job memory (GB) [default]",
    dest = 'memory'
  ),
  optparse::make_option(
    "--date",
    type = 'character',
    default = NA,
    metavar = '"YYYY-MM-DD"',
    help = 'As of date to extract data from.  Defaults to today.',
    dest = 'date'
  ),
  optparse::make_option(
    "--controller",
    type = "character",
    default = NA_character_, # Think about adding the controller in the software directory
    help = "Path to sh script directing the full analysis",
    dest = 'controller'
  )

)

args <-  optparse::parse_args(optparse::OptionParser(option_list = option_list))


###############
# helper f's
###############

.unzip.patstats <- function(x)
{
        csv.name <- unzip(x, list = TRUE)$Name
        csv.name <- grep('_patStats.csv$',csv.name,value = T)
        patstat <- as.data.table(read.csv(unz(x, csv.name), header = TRUE, sep = ","))
        csv.name <- file.path(dirname(x), csv.name)
        write.csv(patstat, file=csv.name)
        return(csv.name)
}

generate.sample <- function(maf) 
{
        cat('\n', 'Running generate.sample()...', '\n')
        # Lele's function to compute Minor Allele Frequencies. 
        MAF_matrix<-matrix(0,ncol=10000,nrow=(nrow(maf)+1))
        MAF_matrix[1,]<-seq(1,10000)
        rownames(MAF_matrix)<-c('pos', sort(maf$SAMPLE_ID))
        names_nobf <- maf[HXB2_EXISTS==FALSE, SAMPLE_ID] 
        idx_nobf <- which(rownames(MAF_matrix) %in% names_nobf)
        MAF_matrix[idx_nobf, ] <- NA
        spls <- names(which(!is.na(MAF_matrix[,1])))
        spls <- spls[spls != 'pos']

        for(sp in spls){
                sp_file <- maf[SAMPLE_ID == sp, HXB2_PATH]
                # cat(paste0('Computing MAF for: ', sp_file,'\n'))
                basefile<-read.csv(sp_file, header=T, stringsAsFactors=FALSE)
                indexes_HXB2_pos<-which(basefile$Position.in.B.FR.83.HXB2_LAI_IIIB_BRU.K03455!='-')
                sample_HXB2_pos<-as.numeric(basefile[indexes_HXB2_pos,1])
                sample_MAFs<-apply(basefile[indexes_HXB2_pos,4:7], 1, function(x) 1-(max(x,na.rm = T)/sum(x,na.rm = T)))
                sample_MAFs<-as.numeric(gsub(NaN, 0,sample_MAFs))
                MAF_matrix[sp,sample_HXB2_pos]<-sample_MAFs 
        }

        return(MAF_matrix)
}

get.sampling.dates <- function(phsc.samples = args$phsc.samples)
{
        # Find files containing all sample collection dates 
        db.sharing.path.rccs <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv'
        db.sharing.path.mrc <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv'

        tmp <- c(db.sharing.path.rccs,db.sharing.path.mrc, phsc.samples)
        stopifnot(all(file.exists(tmp)))

        # Convert PANGEA_IDs to RENAME_ID (sample IDs)
        dsamples <- setDT(readRDS(phsc.samples))
        dsamples <- unique(dsamples[, .(PANGEA_ID, RENAME_ID)])
        dsamples[, PANGEA_ID:=gsub('^.*?_','',PANGEA_ID)]

        # Get sampling dates 
        ddates <- setDT(read.csv(db.sharing.path.mrc))
        ddates <- unique(ddates[, .(pangea_id, visit_dt)])
        tmp <- fread(db.sharing.path.rccs)
        tmp <- unique(tmp[, .(pangea_id, visit_dt=as.character(visit_dt))])
        ddates <- rbind(tmp, ddates)
        ddates[, visit_dt:=as.Date(visit_dt, format="%Y-%m-%d")]
        stopifnot(ddates[, anyDuplicated(pangea_id) == 0,])

        # Subset to pop of interest
        ddates <- merge(dsamples, ddates, all.x=TRUE,
                        by.x='PANGEA_ID', by.y='pangea_id')
        ddates[, PANGEA_ID := NULL]
        
        # Order based on sampling dates 
        setnames(ddates, 'RENAME_ID', 'SAMPLE_ID')
        ddates[, AID := gsub('-fq.*?$','', SAMPLE_ID)]
        setorder(ddates, AID, -visit_dt)

        return(ddates)
}


# TODO:
write.allmafs.from.samples.with.bf <- function(path.samples.with.bf,
                                               bf.prefix='/rds/general/project/ratmann_pangea_deepsequencedata/live')
{
    # Load files 
    dsamples <- readRDS(path.samples.with.bf)

    # select samples for which we can compute base frequencies.
    idx <- dsamples[!is.na(BF), .(RENAME_ID, BF), ]

    # function to read bf files and compute Minore allele freqs.
    .get.minor.allele.freqs <- function(path)
    {

        # prepare output 
        out <- rep(NA, 10000)

        # read input base frequency files and rename columns
        # path <- idx[file.exists(file.path(prefix, BF)),BF]
        dbfs <- fread(path, header=TRUE)
        names(dbfs)[names(dbfs) %like% 'HXB2'] <- 'Position_HXB2'
        dbfs <- dbfs[Position_HXB2 != '-', ]

        # Compute MAFs for given file
        names(dbfs) <- gsub(' ', '_', names(dbfs))
        dbfs[,{ 
             z <- c(A_count, C_count, G_count, T_count)
             list(maf=1 - (max(z, na.rm=TRUE)/sum(z, na.rm=TRUE)))
        }, by='Position_HXB2'] -> dmafs

        # put into output
        out[as.integer(dmafs$Position_HXB2)] <- dmafs$maf
        return(out)
    }

    # Apply to every file and combine
    idx <- idx[file.exists(file.path(bf.prefix, BF)), ]
    idx[, BF := file.path(prefix, BF)]
    maf_all <- lapply(idx$BF, .get.minor.allele.freqs)
    maf_all <- Reduce('rbind', maf_all)
    maf_all <- as.data.table(maf_all)
    names(maf_all) <- as.character(1:10000)
    maf_all[, pos := idx$RENAME_ID]
    setcolorder(maf_all, 'pos')

    # save 
    stopifnot(nrow(maf_all) == nrow(idx))
    path.output <- gsub('samples_with_bf.rds$', 'samples_maf.csv', path.samples.with.bf)
    fwrite(maf_all, file = path.output)
}


# Should break this into two!
write.mafs.and.cmds <-function(pty_idx)
{
        # For each PTY index, computes the MAF matrix and stores it in the phscrel directory
        # Also, writes bash command to set up HIV-phylo-TSI

        cat('Processing PTY index: ', pty_idx, '...\n')

        # Load patstats
        files_pty <- as.vector(dfiles[pty==pty_idx]) 
        patstats <- fread(files_pty$pat.path, header=TRUE, sep=",", stringsAsFactors=F)
        if(any((patstats$xcoord %% 1) != 0)){
                patstats[, xcoord:=ceiling(xcoord)]
                write.csv(patstats, files_pty$pat.path)
        }
        ph.input <- fread(files_pty$phi.path, header=FALSE, sep=",",stringsAsFactors=F)
        colnames(ph.input) <- c('BAM_PATH','FASTA_PATH', 'SAMPLE_ID')
        ph.input[, AID:=gsub('-fq.*?$','', SAMPLE_ID)]
        stopifnot( all(patstats[, unique(host.id)] %in% ph.input[, unique(AID)]) )

        # (Checked that BAM and FASTA files are consistent in terms of namings and locs)

        # Find BAM_PATH then get the MAF
        ph.input[, HXB2_PATH := gsub('.bam$','_BaseFreqs_WithHXB2.csv', basename(BAM_PATH))]
        ph.input[, HXB2_PATH := file.path(dirname(BAM_PATH), HXB2_PATH)]
        ph.input[, HXB2_EXISTS := file.exists(HXB2_PATH)]
        maf <- ph.input[, .(SAMPLE_ID, HXB2_PATH, HXB2_EXISTS)]
        maf_mat <- generate.sample(maf)
        # cat(maf_mat[2, 1:10], '\n')

        # If there are multiple sequences associated to one AID:
        # take sequence with associated BaseFreq file ("HXB2")
        # with latest collection date  
        tmp1 <- gsub('-fq.*?$','',rownames(maf_mat))
        tmp1 <- unique(tmp1[duplicated(tmp1)])
        tmp1 <- data.table(AID = tmp1)
        if(tmp1[, .N>1])
        {
                tmp1 <- tmp1[, list(FQ=grep(AID, rownames(maf_mat), value=T)),by=AID]
                # which do not have HXB2?
                tmp1 <- tmp1[, list(HXB2 = !is.na(maf_mat[FQ, 1])), by=c("AID","FQ") ]
                tmp1 <- merge(tmp1, ddates, by.x=c('AID', 'FQ'), by.y=c("AID", "SAMPLE_ID"))
                setorder(ddates, AID, -visit_dt)
                tmp1 <- tmp1[, {
                        z <- which(HXB2 == TRUE)[1];
                        z <- ifelse(is.na(z), 1, z)
                        list(FQ=FQ[z], visit_dt=visit_dt[z])
                }, by='AID']


                rows_to_del <- rownames(maf_mat)[which(grepl(paste0(tmp1$AID, collapse='|'), rownames(maf_mat) ))]
                rows_to_del <- rows_to_del[! rows_to_del %in% tmp1$FQ]
                
                # Store -fq used so we can check exact sampling date
                maf_mat <- maf_mat[! rownames(maf_mat) %in% rows_to_del, ]
        }else{
                rows_to_del <- c()
        }
        filename=paste0('ptyr', pty_idx, '_basefreqs_used.csv')

        write.csv(rownames(maf_mat), 
                  file=file.path(dirname(files_pty$pat.path), filename))

        rownames(maf_mat) <- gsub('-fq[0-9]$', '', rownames(maf_mat))
        rm(tmp1, rows_to_del)
        maf_file=paste0('ptyr', pty_idx, '_maf.csv')
        maf_file=file.path(dirname(files_pty$pat.path), maf_file)
        cat(maf_mat[2, 1:10], '\n')
        write.table(maf_mat, file=maf_file, quote=F, col.names=FALSE, sep=',')


        # writing command
        #________________
        cmd <- paste0('python ',args$TSI.dir,'/HIVPhyloTSI.py \\\n',
                      ' -d ', args$TSI.dir, '/Model \\\n',
                      ' -p ', files_pty$pat.path,' \\\n',
                      ' -m ', files_pty$maf.path,' \\\n',
                      ' -o ', files_pty$tsi.path, '\n'
                     )
        dfiles[pty==pty_idx, CMD:=cmd]

        # guess I may need to free up some memory here
        rm(maf_mat, ph.input, patstats)
        gc()
}

write.cmds <- function(pty_idx, maf.path)
{
    # writes bash command to set up HIV-phylo-TSI
    cat('Processing PTY index: ', pty_idx, '...\n')
    files_pty <- dfiles[]

    # writing command
    #________________
    cmd <- paste0('python ',args$TSI.dir,'/HIVPhyloTSI.py \\\n',
                  ' -d ', args$TSI.dir, '/Model \\\n',
                  ' -p ', files_pty$pat.path,' \\\n',
                  ' -m ', maf.path,' \\\n',
                  ' -o ', files_pty$tsi.path, '\n'
    )
    dfiles[pty==pty_idx, CMD:=cmd]

    # guess I may need to free up some memory here
    gc()
}

###############
# testing
###############

if(0){
  args <- list(
    out.dir = "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI",
    pkg.dir = "/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software",
    rel.dir="/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI/2022_08_22_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_001_rla_T_zla_T",
    phsc.samples="/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI/220331_RCCSUVRI_phscinput_samples_with_bf.rds",
    TSI.dir="/rds/general/user/ab1820/home/git/HIV-phyloTSI",
    date = '2022-08-22',
    controller='/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/runall_TSI_pairs2.sh',
    env_name = 'hivphylotsi'
  )
}


################
# main
################

source(file.path(args$pkg.dir, "utility.R"))

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

stopifnot(dir.exists(args$rel.dir))
stopifnot(file.exists(args$phsc.samples))

# get mafs file if it doesn't exist already
path.all.maf <- gsub('samples_with_bf.rds$', 'samples_maf.csv', args$phsc.samples)
if(! file.exists(path.all.maf) )
    write.allmafs.from.samples.with.bf(path.samples.with.bf=args$phsc.samples)

# dates of collection for samples.
ddates <- get.sampling.dates(phsc.samples=args$phsc.samples)

# Collect pty files allowing to run Tanya's algorithm
patstats_zipped <- list.files(args$rel.dir, pattern='zip$', full.name=TRUE)
tmp <- gsub('^.*?ptyr|_otherstuff.zip','',patstats_zipped)
phsc_inputs <- file.path(args$out.dir.output, 
                         paste0('ptyr', tmp, '_trees'), 
                         paste0('ptyr', tmp, '_input.csv' ))


if(! file.exists(path.all.maf) )
{
    dfiles <- data.table(pty=as.integer(tmp), 
                         zip.path=patstats_zipped, 
                         phi.path=phsc_inputs)
    setkey(dfiles, pty)
    dfiles <- dfiles[ file.exists(zip.path) & file.exists(phi.path)]
    dfiles[, pat.path:=.unzip.patstats(zip.path), by='pty']
    dfiles[, maf.path:=gsub('patStats.csv$','maf.csv',pat.path)]
    dfiles[, tsi.path:=gsub('patStats.csv$','tsi.csv',pat.path)]
    dfiles[, CMD:=NA_character_]

    lapply(dfiles$pty, write.mafs.and.cmds)

}else{

    dfiles <- data.table(pty=as.integer(tmp), 
                         zip.path=patstats_zipped)
    setkey(dfiles, pty)
    dfiles <- dfiles[ file.exists(zip.path)]
    dfiles[, pat.path:=.unzip.patstats(zip.path), by='pty']
    dfiles[, tsi.path:=gsub('patStats.csv$','tsi.csv',pat.path)]
    dfiles <- dfiles[!file.exists(tsi.path)]
                         
    lapply(dfiles$pty, write.cmds, maf.path=path.all.maf)
}

##################################
# submit jobs 
##################################

dfiles[, IDX := 1:.N, ]
dfiles[, CMD := paste0(IDX, ')\n',CMD, ';;\n')]
cmd <- paste0(dfiles$CMD, collapse='\n')
cmd <- paste0('case $PBS_ARRAY_INDEX in\n', cmd, 'esac\n')

header <- cmd.hpcwrapper.cx1.ic.ac.uk(hpc.select=1, 
                                      hpc.walltime=args$walltime,
                                      hpc.mem=args$memory,
                                      hpc.nproc=1, 
                                      hpc.q=NA, 
                                      hpc.load=paste0('module load anaconda3/personal\nsource activate ', args$env_name, '\n'),
                                      hpc.array=nrow(dfiles))

# Patch together and write 
cmd <- paste0(header, cmd, sep='\n')

djob <- data.table(JOB_ID=1, CMD=cmd)
ids <- djob[, list(ID=.store.and.submit(.SD, prefix='tsi')), by=JOB_ID, .SDcols=names(djob)]

# queue the next step
qsub.next.step(file=args$controller,
               ids=ids, 
               next_step='dti')
