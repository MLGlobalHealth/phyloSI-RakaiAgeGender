cat("\n\n===== TSI_initialise.R =====\n\n")

require(data.table)

option_list <- list(
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
    "--cluster_size",
    type = "integer",
    default = 50L,
    help = "Minimum cluster size for host ids groupings[default %default]",
    dest = "cluster_size"
  ),
  optparse::make_option(
    "--controller",
    type = "character",
    default = NA_character_, 
    help = "Path to sh script directing the full analysis",
    dest = 'controller'
  ),
  optparse::make_option(
    "--transmission_chains",
    type = "character",
    default = NA_character_,
    help = "Optional: absolute file path to `phscnetwork.rda` containing individuals in potential transmission pairs",
    dest = 'file.path.chains'
  ),
  optparse::make_option(
    "--include_input",
    type = "character",
    default = NA_character_,
    help = "Optional: path to phscinput*rds file of individuals to include in the analysis (eg seroconverters)",
    dest = 'include.input'
  ),
  optparse::make_option(
    "--include_least_recent_only",
    type = "logical",
    default = FALSE,
    help = "Optional: Logical, indicating whether we want to only include the most recent blood sample for each individual, or not",
    dest = 'include.least.recent.only'
  ),
  optparse::make_option(
    "--include_twosample_individuals_only",
    type = "logical",
    default = FALSE,
    help = "Optional: Logical, indicating whether we want to only include participants with exactly two collected samples(test)", 
    dest = 'include.twosample.individuals.only'
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

###
# Helpers
###

get.sampling.dates <- function(phsc.samples = args$phsc.samples)
{
        # Find files containing all sample collection dates 
        .f <- function(x) file.path(indir.deepsequencedata, x)
        db.sharing.path.rccs <- .f('PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv')
        db.sharing.path.mrc  <- .f('PANGEA2_MRC/200319_pangea_db_sharing_extract_mrc.csv')

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

###
# Paths
###

usr <- Sys.info()[['user']]
if (usr == 'andrea')
{
        args$pkg.dir <- '~/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/'
        indir.deepsequence.analyses <- '~/Documents/Box/ratmann_deepseq_analyses/live'
        indir.deepsequence.xiaoyue <- '~/Documents/Box/ratmann_xiaoyue_jrssc2022_analyses/live'
        indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata'
}else{
        args$pkg.dir <- '~/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/'
        indir.deepsequence.analyses <- '/rds/general/project/ratmann_deepseq_analyses/live'
        indir.deepsequence.xiaoyue <- "/rds/general/project/ratmann_xiaoyue_jrssc2022_analyses/live"
        indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
}

indir.deepsequence.analyses.old <- file.path(indir.deepsequence.xiaoyue, 'PANGEA2_RCCS1519_UVRI')
file.phsc.input.samples.bf <- file.path(indir.deepsequence.analyses.old, '220331_RCCSUVRI_phscinput_samples_with_bf.rds' )

tmp <- c(indir.deepsequence.analyses.old,
         file.phsc.input.samples.bf)
if( ! is.na(args$file.path.chains) ) tmp <- c(tmp, args$file.path.chains)
if( ! is.na(args$include.input) ) tmp <- c(tmp, args$include.input)
stopifnot(all(file.exists(tmp)))


###
# Main 
###

# Create ouput directory
if(! dir.exists(args$out.dir))
{
        cat('Creating output directory\n')
        dir.create(args$out.dir)
        dir.create(file.path(args$out.dir, 'potential_network'))
}else{
        cat('Warning: output directory already exists\n')
}

# load individuals in potential transmission pairs if given
#__________________________________________________________

if( ! is.na(args$file.path.chains) )
{
        tmp <- new.env()
        load(args$file.path.chains, envir=tmp )
        dchain <- as.data.table(tmp$dchain)
        rm(tmp)
        include_pairs_aid <- dchain[, unique(c(H1, H2))]
}

if( ! is.na(args$include.input))
{
        tmp <- readRDS(args$include.input)
        include_rename_id <- unique(tmp$RENAME_ID)
        rm(tmp)
}



# Load old phsc input samples, and subset as required
#____________________________________________________

phsc_samples <- readRDS(file.phsc.input.samples.bf)
phsc_samples[, AID := gsub('-fq[0-9]+$', '', RENAME_ID)]
phsc_samples[, INCLUDE := TRUE]

if( args$include.least.recent.only)
{
        cat('Only use least recently/first collected sample per individual...\n')
        ddates <- get.sampling.dates(phsc.samples=file.phsc.input.samples.bf)
        setorder(ddates, AID, visit_dt)
        .f <- function(x) x[[1]]
        ddates <- ddates[, lapply(.SD, .f), by='AID']
        phsc_samples[ ! RENAME_ID %in% ddates$SAMPLE_ID, INCLUDE:=FALSE ]
}

if( ! is.na(args$file.path.chains))
{
        cat('Subsetting to individuals in chains...\n')
        stopifnot( all(include_pairs_aid %in% phsc_samples$AID ))
        phsc_samples[! AID %in% include_pairs_aid, INCLUDE := FALSE]
}

if( ! is.na(args$include.input))
{
        cat('Including seroconverters samples... \n')
        stopifnot( all(include_rename_id %in% phsc_samples$RENAME_ID ))
        phsc_samples[ RENAME_ID %in% include_rename_id, INCLUDE := TRUE]
}

# if( ! is.na(args$file.path.chains) | ! is.na(args$include.input) )
# {
#         tmp <- basename(filename)        
#         date <- format(Sys.Date(), '%y%m%d')
#         tmp <- gsub('^[0-9]+', date, tmp)
#         tmp <- gsub('\\.rds$', '_subset.rds', tmp)
#         filename <- file.path(dirname(filename), tmp)
# }

if( args$include.twosample.individuals.only)
{
        # Extract individuals with only two samples, and differing visit dates
        idx <- phsc_samples[, uniqueN(RENAME_ID) == 2, by='AID']
        idx <- idx[V1==TRUE, AID]
        phsc_samples[! AID %in% idx, INCLUDE := FALSE]

        # Get sampling dates, so we can make sure the samples were taken at different dates. 
        ddates <- get.sampling.dates(phsc.samples=file.phsc.input.samples.bf)
        tmp <- phsc_samples[INCLUDE == TRUE]
        tmp <- merge(tmp, ddates[,.(RENAME_ID=SAMPLE_ID, visit_dt)], by='RENAME_ID')
        setkey(tmp, AID, visit_dt)

        idx <- tmp[, visit_dt[2] - visit_dt[1] == 0, by='AID'][V1 == TRUE, AID]
        tmp <- tmp[! AID %in% idx, ]

        tmp <- rbind(
                tmp[, lapply(.SD, function(x) x[1]) ,  by='AID'][1:100],
                tmp[, lapply(.SD, function(x) x[2]) ,  by='AID'][1:100]
        ) 

        tmp[101:200,UNIT_ID := paste0(UNIT_ID, '+')]

        phsc_samples <- copy(tmp)
}

phsc_samples <- phsc_samples[INCLUDE == TRUE]
phsc_samples[, `:=` (AID=NULL, INCLUDE=NULL)]
filename=file.path(args$out.dir, basename(file.phsc.input.samples.bf))
saveRDS(phsc_samples, filename)


# Make clusters.rds
# ______________________________
set.seed(42)

make.clusters <- function(DT)
{
        idx <- unique(DT$UNIT_ID)
        NCLU <- DT[, floor(length(idx) / args$cluster_size) ]
        dclus <- DT[, list(ID=sample(idx), IDCLU=1:NCLU)]
        dclus[, CLU_SIZE:=.N, by='IDCLU']
        setkey(dclus, IDCLU)

        if(args$include.twosample.individuals.only)
        {
                cat('Samples from same individual run in different clusters\n',
                    '\tCheck whether this done on purpose.\n')
                dclus1 <- dclus[, list(
                                       ID = paste0(ID, '+'),
                                       IDCLU = IDCLU + max(IDCLU),
                                       CLU_SIZE = CLU_SIZE
                                       )]
                dclus <- rbind(dclus, dclus1)
        }
                dclus
}

dclus <- make.clusters(phsc_samples[! UNIT_ID %like% '\\+', ])


filename=file.path(args$out.dir, 'potential_network', 'clusters.rds')
saveRDS(dclus, filename)


# make phscinput_runs_clusize...
# ______________________________
suffix <- phsc_samples[, .(UNIT_ID,PANGEA_ID, RENAME_ID, SAMPLE_ID)]
dclus[, `:=` (PTY_RUN=IDCLU,  PTY_SIZE=CLU_SIZE) ]
dclus <- merge(dclus, suffix, by.x='ID', by.y='UNIT_ID')
setnames(dclus, 'ID', 'UNIT_ID')
setkey(dclus, IDCLU)
filename=file.path(args$out.dir,
                   paste0('phscinput_runs_clusize_', args$cluster_size,'_ncontrol_0.rds'))
saveRDS(dclus, filename)


# Source functions
source(file.path(args$pkg.dir, "utility.R"))
qsub.next.step(file=args$controller,
               next_step='ali', 
               res=1, 
               redo=0
)
