# AIMS:
# compute generation intervals for pairs in a network.
# needs phylo_dir as an input

################
# DEPENDENCIES #
################ 
library(data.table)
library(ggplot2)
library(lubridate)
library(xtable)
library(igraph)

################
#   OPTIONS    #
################

option_list <- list(
    optparse::make_option(
        "--path-chains",
        type = "character",
        default = NULL,
        help = "path to main output of source recipients pairs analyses",
        dest = 'path.chains.data'
    ),
    optparse::make_option(
        "--path-tsi",
        type = "character",
        default = NULL,
        help = "path to main output of time since infection analyses",
        dest = 'path.tsiestimates'
    ),
    optparse::make_option( 
        "--confidential",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to use the confidential data (if access is granted) [Defaults to TRUE]", 
        dest = 'confidential'
    ),
    optparse::make_option(
        "--rerun",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to re-run rejection sampling [Defaults to FALSE]", 
        dest = 'rerun'
    ),
    optparse::make_option(
        "--get-round-probabilities",
        type = "logical",
        default = TRUE,
        help = "Computes probability of infection occurring in different rounds", 
        dest = 'get.round.probabilities'
    ),
    optparse::make_option( "--RH-infectiousness",
        type = "numeric",
        default = 1,
        help = "Relative Hazard of transmission during the acute versus chronic infection phase [Defaults to 5]", 
        dest = 'RH.infectiousness'
    ),
    optparse::make_option(
        "--RH-duration",
        type = "numeric",
        default = 1,
        help = "Duration of acute phase, in years [Defaults to 2 months]", 
        dest = 'RH.duration'
    ),
    optparse::make_option( "--sensitivity-no-refinement",
        type = "logical",
        default = FALSE,
        help = "Relative Hazard of transmission during the acute versus chronic infection phase [Defaults to 5]", 
        dest = 'sensitivity.no.refinement'
    ),
    optparse::make_option(
        "--intermediate",
        type = "logical",
        default = FALSE,
        help = "Save intermediate outputs to outdir.", 
        dest = 'save.intermediate'
    ),
    optparse::make_option(
        "--outdir",
        type = "character",
        default = NULL,
        help = "Absolute file path to output directory to save the scripts outputs in", 
        dest = 'outdir'
    )
)

args <-  optparse::parse_args(optparse::OptionParser(option_list = option_list))
print(args)

catn <- function(x) cat("\n----", x, "----\n")

################
#    PATHS     #
################

usr <- Sys.info()[['user']]

gitdir <- here::here()
gitdir.data <- file.path(gitdir, 'data')

# if output directory is null, set it to gitdir.data
if( is.null(args$outdir) )
    args$outdir <- gitdir.data

# randomized version of the paths used for the analysis
path.meta <- file.path(gitdir.data,"Rakai_Pangea2_RCCS_Metadata_randomized.RData" )
path.sequence.dates <- file.path(gitdir.data, "sequences_collection_dates_randomized.rds")

if( args$confidential)
{
    if(usr == 'andrea')
    {
        args$phylo_dir <- '/home/andrea/HPC/project'
        intermed.out.dir<- '/home/andrea/HPC/ab1820/home/projects/2022/genintervals'
        indir.deepsequencedata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
    }

    # paths that should not be available to everyone, but were used for the analysis
    path.meta <- file.path(indir.deepsequencedata, 'RCCS_R15_R18/Rakai_Pangea2_RCCS_Metadata_20221128.RData')
    path.sequence.dates <- file.path(indir.deepsequencedata, "RCCS_R15_R18/sequences_collection_dates.rds")

}

# intermediate output directory:
if( args$save.intermediate )
{
    intermed.out.dir <- args$outdir
}else{
    intermed.out.dir <- tempdir()
    warning("output directory for intermediary results not set.\n Setting a temporary directory through `tempdir()`")
}

path.round.timeline <- file.path(gitdir.data, 'RCCS_round_timeline_220905.RData' )
path.chains.data <- fifelse(
    is.null(args$path.chains.data),
    yes=file.path(gitdir.data, 'Rakai_phscnetworks_ruleo_sero.rda'),
    no=args$path.chains.data
)
path.tsiestimates <- fifelse(
    is.null(args$path.tsiestimates),
    file.path(gitdir.data, 'TSI_estimates.csv'),
    no=args$path.tsiestimates
)

file.exists(
    path.meta,
    path.sequence.dates,
    path.chains.data,
    path.round.timeline,
    # file.anonymisation.keys,
    path.tsiestimates) |> all() |> stopifnot()


################
#    HELPERS   #
################

source(file.path(gitdir, 'functions', 'utils.R'))
source(file.path(gitdir, 'scripts_for_confidential_data/utils', 'functions_tsi_attribution.R'))
source(file.path(gitdir, 'functions', 'plotting_functions.R'))
source(file.path(gitdir, 'functions', 'summary_functions.R'))
find_palette_round()

################
#     MAIN     #
################

# settings 
threshold.likely.connected.pairs <- 0.5
threshold.direction <- 0.5
get.sero.extra.pairs <- FALSE
build.network.from.pairs <- TRUE
postpone.samesex.removal <- TRUE
is.metadata.randomized <- path.meta %like% 'randomized'
stopifnot(is.metadata.randomized==!args$confidential)

# initialise overleaf substitute-expressions.
overleaf_expr <- list()

# Load anon. keys
# aik <- fread(file.anonymisation.keys, header = TRUE, select=c('PT_ID', 'AID'))

meta <- load.meta.data(path.meta)
chain <- build.phylo.network.from.pairs(path.chains.data)
aids_of_interest <- chain[, unique(c(SOURCE,RECIPIENT))]
meta <- meta[aid %in% aids_of_interest]

drange <- get.infection.range.from.testing() |> 
    check.range.consistency()

double.merge(chain, meta[, .(AID=aid, SEX=sex)])[, table(SEX.RECIPIENT, SEX.SOURCE)] |>
    knitr::kable(caption="DOI algorithm is run for a total number of pairs equal to:") |> print()

dancestors <- get.ancestors.from.chain(chain)

# get plausible infection ranges.
chain <- check.inconsistent.testing(drange, switch_if_no_other_src = TRUE)
drange <- shrink.intervals(drange)

# update using TSI estimates
drange_tsi <- get.infection.range.from.tsi(path.tsiestimates, path.sequence.dates )
drange_tsi <- check.inconsistent.testing(drange_tsi, switch_if_no_other_src = FALSE)
drange_tsi <- shrink.intervals(drange_tsi)
setnames(drange_tsi, c('MIN', 'MAX'), c('TSI.MIN', 'TSI.MAX'))

if(args$confidential) # Count number of removed individuals/pairs
{
    tmp <- double.merge(chain, dcomms[, .(AID, COMM, SEX)]) |>
        subset( COMM.SOURCE %like% 'inland' & COMM.RECIPIENT %like% 'inland', select = c('SOURCE','RECIPIENT')) 

    noprobs <- drange_tsi$AID
    overleaf_expr[['N_noTSI']] <- tmp[ ! RECIPIENT %in% noprobs, .N ]
}

# check no contradictions after shrinkning
idx <- double.merge(chain, drange_tsi, by_col = "AID")[TSI.MAX.RECIPIENT - TSI.MIN.SOURCE < 0, unique(IDCLU)]
stopifnot(nrow(idx) == 0)

catn('Finish establishing "prior" infection ranges, and merge to pairs')
if(! args$sensitivity.no.refinement)
{
    drange <- intersect.serohistory.and.predictions.withrefinement(drange, drange_tsi)
}else{
    cat('\nSENSITIVITY ANALYSIS: Using phyloTSI only if not-incoherent\n')
    drange <- intersect.serohistory.and.predictions.withoutrefinement(drange, drange_tsi)
}
dpairs <- double.merge(chain, drange, by_col = "AID")
dpairs[MAX.RECIPIENT - MIN.SOURCE < 0, stopifnot(.N == 0)]

# get transmission cluster ids. 
dclus <- get.transmission.cluster.ids(dpairs, check=FALSE)

# initialise data.table summarising geometric properties
dclus[, { 
    g <- GROUP;
    out <- dclus[GROUP==g, .(SOURCE=SOURCE, RECIPIENT=RECIPIENT)]
    out <- unique(out)
    N_out <- out[, uniqueN(c(SOURCE,RECIPIENT))]
    out <- list(out)
    out <- list(IDS=out, N_IDS=N_out)
}, by=GROUP] -> dcohords

catn('specifying relative infectiousness as a function of time since infection')
dinfectiousness <- specify.relative.infectiousness(args)


catn(" ** RUN MCMC ** ")

# gi: generation intervals
range_gi <- c(.5, 1:10)
tmp <- dinfectiousness[!is.infinite(END), paste(round(END, 2), sep='')]
tmp <- gsub('\\.', '', tmp)

filename_drange <- 'networks_individualDOIrange'
filename_net <- 'networks_GICentroids'
filename_overleaf <- 'pairsinfo_overleaf'

suffix <- set.mcmc.outputs.suffix()

filename_drange <- file.path(intermed.out.dir, paste0(filename_drange, suffix,  '.rds'))
filename_net <- file.path(intermed.out.dir, paste0(filename_net, suffix, '.rds'))
filename_overleaf <- file.path(intermed.out.dir, paste0(filename_overleaf, suffix, '.rds'))

if(args$get.round.probabilities)
    df_round_gi <- get.round.dates(path.round.timeline)

if( file.exists(filename_net) & ! args$rerun )
{
    dcohords <- readRDS(filename_net)
}else{

    # get A LOT of Uniform Samples for importance sampling MC
    provide_samples <- lapply( dcohords[, seq(max(N_IDS))], function(x){runif(1e6, 0, 1)})
    n_iter <- nrow(dcohords)

    dcohords[ , {
        cat(GROUP, '/', n_iter, '\n')
        DT <- IDS[[1]]
        out <- get.volume.genints.multid(DT,
                                         N_IDS=N_IDS,
                                         range_gi=range_gi,
                                         verbose=T,
                                         importance_weights=dinfectiousness,
                                         get_volume=T,
                                         provide_samples=provide_samples,
                                         df_round=df_round_gi)
        list(out, NAME=c('CENTROID', 'VOLUME', 'GI', 'RPROBS'))
    } , by=GROUP] -> tmp

    Reduce('merge',
           list(tmp[NAME == 'CENTROID', .(GROUP, CENTROIDS=out)],
                tmp[NAME == 'VOLUME', .(GROUP, VOLUME=unlist(out))],
                tmp[NAME == 'GI', .(GROUP, GENINTS=out)]),
                tmp[NAME == 'RPROBS', .(GROUP, RPROBS=out)]
           ) -> tmp
    dcohords <- merge( dcohords[, .(GROUP, IDS, N_IDS)], tmp, by='GROUP')

    rm(provide_samples, n_iter)
    saveRDS(dcohords, filename_net)
}

# checks only really make sense with the real data
if(args$confidential)   
{   # compute 95% ranges
    dpred_ranges <- dcohords[, CENTROIDS[[1]] ,by=GROUP]
    stopifnot( dpred_ranges[, .N == uniqueN(ID)] )
    stopifnot( dpred_ranges[, all(IL <= M & M <= IU) ] )
    stopifnot( dpred_ranges[, all(CL <= IL & IU <= CU) ] )

    tmp <- merge(dpred_ranges, drange, by.x='ID', by.y='AID')
    stopifnot( tmp[, all(IL >= MIN)] )
    stopifnot( tmp[, all(IU <= MAX)] )
}

if( nrow(df_round_gi) )
{   # compute probability of assigning source to different period
    dprobs_roundallocation <- dcohords[ , RPROBS[[1]], by='GROUP']
    dprobs_roundallocation <- merge(
                                    df_round_gi[, .(ROUND, INDEX_TIME)],
                                    dprobs_roundallocation,
                                    by='ROUND', all.y=TRUE)
    dprobs_roundallocation <- dprobs_roundallocation[, 
            .(P=sum(P)),
            by=c('SOURCE','RECIPIENT','INDEX_TIME')]
    
    setkey(dprobs_roundallocation, RECIPIENT, SOURCE, INDEX_TIME)
    dprobs_roundallocation[is.na(INDEX_TIME), INDEX_TIME := 0 ]

    sr <- c('SOURCE','RECIPIENT')
    check_probs_sum_to_1 <- dprobs_roundallocation[, abs(sum(P) - 1) <= .00001, by=sr][,all(V1)]
    stopifnot(check_probs_sum_to_1)

    dprobs_roundallocation <- dcast(dprobs_roundallocation, SOURCE+RECIPIENT ~ INDEX_TIME)
    setnames(dprobs_roundallocation, c('0', '1', '2'), c('BeforeR10', 'R10_R15', 'R16_18'))
    cols <- c('BeforeR10', 'R10_R15', 'R16_18')
    dprobs_roundallocation[, (cols) := lapply(.SD, function(x){x[is.na(x)] <- 0; x} ) , .SDcols=cols]
    
    filename <- file.path(intermed.out.dir, paste0('probs_roundallocations', suffix, '.rds'))
    saveRDS(dprobs_roundallocation, filename)
}

# Final results
filename <- file.path(args$outdir, paste0('pairsdata_toshare', suffix, '.rds')) 

if(  file.exists(filename) & ! args$rerun == TRUE )
{
    dresults <- readRDS(filename)
}else{

    # prepare input for Bayesian model
    centroids <- dcohords[, CENTROIDS[[1]], by='GROUP']
    dresults <- prepare.pairs.input.for.bayesian.model(centroids)
    dresults[, table(DIRECTION, useNA='ifany')]

    # add homosexual pairs to results
    dhomosexualpairs[, DIRECTION := 'phyloscanner']
    dhomosexualpairs[, `:=` (SCORE=NULL, SCORE_DIR=NULL, SEX.SOURCE=NULL, SEX.RECIPIENT=NULL)]
    dhomosexualpairs <-  double.merge(dhomosexualpairs,meta[, .(AID=aid, SEX=sex)])
    rbind(
        dresults,
        dhomosexualpairs,
        fill=TRUE
    ) -> dresults
    dresults <- get.community.type.at.infection.date(dresults, comm_number = FALSE)

    # denote which directions were fixed through serohistory
    dresults <- rbind(
        serohistory_impact[, .(SOURCE=H1, RECIPIENT=H2, CHANGED_DIR)],
        serohistory_impact[, .(RECIPIENT=H1, SOURCE=H2, CHANGED_DIR)]  
    ) |> merge(x = dresults, all.x=TRUE)
    dresults[ CHANGED_DIR == TRUE, DIRECTION := 'sero-adjusted']
    dresults[, CHANGED_DIR := NULL]

    if(get.sero.extra.pairs & FALSE) # deprecated build.network.from.pairs
    {
        if(nrow(additional_pairs_from_serohistory))
        {
            setkey(dresults, SOURCE, RECIPIENT)
            idx <- additional_pairs_from_serohistory[, .(SOURCE, RECIPIENT)]
            dresults[idx, DIRECTION := 'serohistory']

            # 
            idx2 <- idx[! RECIPIENT %in% dresults$RECIPIENT]
            dadditional <- meta[, .(SOURCE=aid, SEX.SOURCE=sex, DIRECTION='serohistory')]
            dadditional <- merge(idx2, dadditional)
            tmp1 <- meta[, .(RECIPIENT=aid, SEX.RECIPIENT=sex)]
            dadditional <- merge(dadditional, tmp1, by='RECIPIENT')
            setcolorder(dadditional, c('SOURCE','RECIPIENT', 'SEX.SOURCE', 'SEX.RECIPIENT', 'DIRECTION'))
            rbind(
                dresults,
                dadditional,
                fill=TRUE
            ) -> dresults
        }
    }

    # get sample collection dates
    if(args$confidential)
    {
        tmp <- get.sample.collection.dates(get_first_visit=TRUE)
        names(tmp) <- toupper(gsub('_','.',names(tmp)))
        dresults <- double.merge(dresults, tmp)
        if(! args$sensitivity.no.refinement)
            stopifnot(dresults[, all(CU < DATE.COLLECTION.RECIPIENT, na.rm=TRUE)])
    }

    if(exists('overleaf_expr'))
    {
        overleaf_expr
        idx_inland <- dresults[ COMM.SOURCE == 'inland' & COMM.RECIPIENT == 'inland']
        idx_inland_hetero <- idx_inland[ SEX.SOURCE != SEX.RECIPIENT ]
        overleaf_expr[['N_final_inland']] <- idx_inland_hetero[, .N]

        overleaf_expr[['N_final_inland_instudy_diff']] <- idx_inland_hetero[ is.na(ROUND.M), .N ]
        overleaf_expr[['N_final_inland_instudy']] <- idx_inland_hetero[ ! is.na(ROUND.M), .N ]
        overleaf_expr[['N_final']] <- idx_inland_hetero[ !is.na(ROUND.M), .N ]
    }

    # summarise 
    dresults[!is.na(M) &
             SEX.SOURCE != SEX.RECIPIENT &
             COMM.SOURCE == 'inland' & COMM.RECIPIENT == 'inland', {
                 cat('There are', .N, 'pairs with median DOI estimate(', sum(is.na(ROUND.M)), 'outside of range).\n'); .SD
             } ][is.na(ROUND.M), .(SOURCE, RECIPIENT, M)] -> idx
    idx <- double.merge(idx, meta[, .(AID=aid, DFP=date_first_positive)])
    setkey(dprobs_roundallocation, SOURCE,RECIPIENT)
    setkey(idx, SOURCE,RECIPIENT)

    # add column statying whether source and recipient both were participants
    dcomms[, PARTICIPATED := COMM %like% 'inland']
    dresults <- double.merge(dresults, dcomms[, .(AID, PARTICIPATED)])
    dresults[, `:=` (
        BOTH_PARTICIPATED=PARTICIPATED.SOURCE & PARTICIPATED.RECIPIENT, 
        PARTICIPATED.SOURCE=NULL, 
        PARTICIPATED.RECIPIENT=NULL
    )]

    saveRDS(dresults, filename)    
    saveRDS(overleaf_expr, filename_overleaf)
}
