# AIMS:
# compute generation intervals for pairs in a network.

# TODO: atm, the specification of the min-max infection ranges is done somewhere else...
# TODO: MIN shouldn't be more than 15 year before first positive

################
# DEPENDENCIES #
################ 
library(data.table)
library(ggplot2)
library(lubridate)
library(xtable)

################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{

    out.dir <- '/home/andrea/Documents/Box/2022/genintervals/'
    indir.deepsequence_xiaoyue   <- '/home/andrea/Documents/Box/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI'
    indir.deepsequencedata <- '/home/andrea/Documents/Box/ratmann_pangea_deepsequencedata'
    indir.deepsequence_analyses   <- '~/Documents/Box/ratmann_deepseq_analyses/live'
    
    indir <- '/home/andrea/git/phyloflows'

}else{
    indir.deepsequence_xiaoyue  #  <- TODO
    indir.deepsequencedata #  <- TODO
    indir # <- TODO
}

.fp <- function(C, x)
{
    if(C=='X')
        indir <- indir.deepsequence_xiaoyue
    if(C=='A')
        indir <- indir.deepsequence_analyses
    if(C=='D')
        indir <- indir.deepsequencedata
    file.path(indir, x)
}

file.path.meta <- .fp('D', 'RCCS_R15_R18/Rakai_Pangea2_RCCS_Metadata_20220329.RData')
file.path.chains.data <- .fp('X','211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')
file.path.round.timeline <- .fp('D', 'RCCS_data_estimate_incidence_inland_R6_R18/220903/RCCS_round_timeline_220905.RData')
file.path.chains.data <- .fp('X','211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')
file.anonymisation.keys <- .fp('X','important_anonymisation_keys_210119.csv')
file.path.tsiestimates <- .fp('A', '/PANGEA2_RCCS_MRC_UVRI_TSI/2022_08_22_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_00/aggregated_TSI_with_estimated_dates.csv') 

################
#   OPTIONS    #
################
# eg: Rscript get_generation_intervals.R --rerun TRUE --RH-infectiousness 1 --RH-duration 1
# eg: Rscript get_generation_intervals.R --rerun TRUE --RH-infectiousness 5 --RH-duration 1/6
option_list <- list(
    optparse::make_option(
        "--rerun",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to re-run rejection sampling [Defaults to FALSE]", 
        dest = 'rerun'
    ),
    optparse::make_option(
        "--only-hetero",
        type = "logical",
        default = TRUE,
        help = "Excludes heterosexual couples from transmission network", 
        dest = 'include.only.heterosexual.pairs'
    ),
    optparse::make_option(
        "--get-round-probabilities",
        type = "logical",
        default = TRUE,
        help = "Computes probability of infection occurring in different rounds", 
        dest = 'get.round.probabilities'
    ),
    optparse::make_option(
        "--only-rccs",
        type = "logical",
        default = TRUE,
        help = "Excludes non-participants of the RCCS", 
        dest = 'include.only.rccs'
    ), optparse::make_option( "--RH-infectiousness",
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
    )
)
args <-  optparse::parse_args(optparse::OptionParser(option_list = option_list))
print(args)

################
#    HELPERS   #
################

source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'functions', 'gi_analysis_functions.R'))
source(file.path(indir, 'functions', 'plotting_functions.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))
# source(file.path(indir, 'functions', 'statistics_functions.R'))
# source(file.path(indir, 'functions', 'stan_utils.R'))
# source(file.path(indir, 'functions', 'check_potential_TNet.R'))


plot.rectangles <- function(DT, values=c('', 'TSI', 'INTERSECT'), idx=data.table())
{
    .p <- function(...) 
    {
        out <- paste(..., sep='.')
        out <- gsub('^\\.', '', out)
    }
    .c <- function(x)
        fcase(x=='', 'blue', x=='TSI', 'red', x=='INTERSECT','orange')
    .s <- 'SOURCE'; .r <- 'RECIPIENT'

    plot_dt <- double.merge(chain, DT)
    setkey(plot_dt, SOURCE,RECIPIENT)
    setkey(idx, SOURCE,RECIPIENT)
    values[.p(values, 'MIN', .s) %in% names(DT)]

    if(nrow(idx))
        plot_dt <- plot_dt[idx]


    add.rect <- function(pre)
    {
        geom_rect( aes_string(xmin=.p(pre,'MIN',.s),
                              xmax=.p(pre,'MAX',.s),
                              ymin=.p(pre,'MIN',.r),
                              ymax=.p(pre,'MAX',.r)), fill=NA, color=.c(pre) )
    }
    
    g <- ggplot(plot_dt)

    for(v in values)
        g <- g + add.rect(v)

    g <- g + 
        geom_abline(aes(intercept=0, slope=1), linetype='dashed', color='black')  +
        theme_bw() + 
        labs(x='DOI source', y='DOI recipient')
    g
}
################
#     MAIN     #
################

# Load anon. keys
aik <- fread(file.anonymisation.keys, header = TRUE, select=c('PT_ID', 'AID'))

if(0)
{
    filename <- paste0('network_attributed_doi',
                       '_chains', .assign.code.meta(file.path.meta),
                       '_meta',  .assign.code.chain(file.path.chains.data),
                       '_onlyhetero', as.integer(include.only.heterosexual.pairs),
                       '.rds')
    filename <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', filename)
}

# get sample collection dates.

# load chains 
load(file.path.chains.data)
threshold.likely.connected.pairs <- 0.5
dchain <- as.data.table(dchain)
chain <- keep.likely.transmission.pairs(as.data.table(dchain), threshold.likely.connected.pairs)
chain <- subset(chain, select=c('SOURCE', 'RECIPIENT', 'IDCLU'))
setkey(chain, IDCLU, SOURCE)

# load meta data
meta_env <- new.env()
load(file.path.meta, envir=meta_env)
meta <- subset(meta_env$meta_data,
               select=c('aid', 'sex', 'date_birth', 'date_first_positive', 'date_last_negative'))
meta <- unique(meta[!is.na(aid)])
stopifnot(meta[, uniqueN(aid) == .N])

# exclude non-RCCS
if( args$include.only.rccs )
{
    tmp <- nrow(chain)
    rccs_ids <- meta[, unique(aid)]
    chain <- chain[ SOURCE %in% rccs_ids & RECIPIENT %in% rccs_ids ]
    cat('Excluding', tmp - nrow(chain), 'of', tmp, 'pairs outside of RCCS\n')
    cat(nrow(chain), 'pairs remaining\n')

}

# exclude homosexuals
if(args$include.only.heterosexual.pairs)
{
    dsex <- meta[, .(aid, sex)] |> unique()
    idx <- merge(chain, dsex[, .(SOURCE=aid, SEX.SOURCE=sex)], all.x=T, by='SOURCE')
    idx <- merge(idx, dsex[, .(RECIPIENT=aid, SEX.RECIPIENT=sex)], all.x=T, by='RECIPIENT')
    # idx[, table(SEX.SOURCE,SEX.RECIPIENT, useNA = 'ifany'),]
    dhomosexualpairs <- idx[ (SEX.SOURCE==SEX.RECIPIENT) , ]
    dhomosexualpairs[, IDCLU:=NULL]
    idx <- idx[!(SEX.SOURCE==SEX.RECIPIENT), .(SOURCE, RECIPIENT)]

    tmp <- nrow(chain) - nrow(idx)
    cat('Excluding', tmp, 'of', nrow(chain), 'pairs of homosexual or unknown sex\n')
    setkey(chain, SOURCE,RECIPIENT)
    setkey(idx,SOURCE,RECIPIENT)
    chain <- chain[idx]
    cat(nrow(chain), 'pairs remaining\n')
}

# get plausible infection ranges.
drange <- get.infection.range.from.testing()
chain <- check.inconsistent.testing(drange, switch = TRUE)
tmp0 <- drange[, as.numeric(mean(MAX - MIN)/365.25)]
shrink.intervals(drange); tmp0

# update using TSI estimates
drange_tsi <- get.infection.range.from.tsi()
drange_tsi <- check.inconsistent.testing(drange_tsi, switch = FALSE)
tmp0 <- drange_tsi[, as.numeric(mean(MAX - MIN)/365.25)]
shrink.intervals(drange_tsi); tmp0
setnames(drange_tsi, c('MIN', 'MAX'), c('TSI.MIN', 'TSI.MAX'))

# intersect the 2
tmp <- merge(drange, drange_tsi)
tmp[,{
    rng_intersect <- (MAX>TSI.MIN & TSI.MAX>MIN)

    fifelse(rng_intersect,
            yes=pmax(MIN, TSI.MIN),
            no=MIN
    ) -> r_min
    fifelse(rng_intersect,
            yes=pmin(MAX, TSI.MAX),
            no=MAX
    ) -> r_max
    list(AID, INTERSECT.MIN=r_min, INTERSECT.MAX=r_max)
    }] -> tmp1


if(0)
{   # VISUALISE
    for(i in 1:100)
    {
        idx <- chain[i, .(SOURCE,RECIPIENT)]
        p <- plot.rectangles(tmp, idx=idx)
        print(p)
        Sys.sleep(1.5)
    }
}

cols <- c('MIN', 'MAX')
setnames(tmp1, paste0('INTERSECT.', cols), cols)
drange <- rbind( drange[! (AID %in% tmp1$AID)],tmp1)


# Specify relative likelihood of infectiousness
# _____________________________________________

cat('specifying relative infectiousness...\n')

# Bellan 2015: 5.3 EHM_\text{acute} d_\text{acute}=1.7months
with(args,
     data.table( START = c(0, RH.duration), VALUE = c(RH.infectiousness, 1))
) -> dinfectiousness
dinfectiousness[, END := c(START[-1], Inf)]
setcolorder(dinfectiousness, c('START', 'END', 'VALUE'))

# Load tranmsission pairs and run MC
# __________________________________

dpairs <- double.merge(chain, drange, by_col = "AID")
# dpairs[, sort(MAX.RECIPIENT - MIN.SOURCE)]
# dpairs[, sort(MAX.SOURCE-MIN.SOURCE)]
# dpairs[, sort(MAX.RECIPIENT-MIN.RECIPIENT)]

# get transmission cluster ids. 
dclus <- get.transmission.cluster.ids(dpairs, check=FALSE)

# initialise data.table summarising geometric properties
dclus[, { g <- GROUP;
         out <- dclus[GROUP==g, unique(.(SOURCE=SOURCE, RECIPIENT=RECIPIENT))]
         N_out <- out[, uniqueN(c(SOURCE,RECIPIENT))]
         out <- list(out)
         out <- list(IDS=out, N_IDS=N_out)
         }, by=GROUP] -> dcohords

# specify pdf through infectiousness

# RUN MCMC
# ________

range_gi <- c(.5, 1:10)
tmp <- dinfectiousness[!is.infinite(END), paste(round(END, 2), sep='')]
tmp <- gsub('\\.', '', tmp)
filename_net <- 'networks_GICentroids'
suffix <- ''
if( ! is.null(dinfectiousness))
    suffix <- paste0('_d', tmp, '_w', paste0( dinfectiousness$VALUE, collapse=''))
filename_net <- file.path(out.dir, paste0(filename_net, suffix, '.rds'))

if(args$get.round.probabilities)
{
    tmp_env <- new.env()
    load(file.path.round.timeline, envir = tmp_env)

    df_round_inland <- tmp_env$df_round_inland
    df_round_inland[, `:=` (min_sample_date = as.Date(min_sample_date), max_sample_date = as.Date(max_sample_date))]
    start_observational_period_inland <- df_round_inland[round == 'R010', min_sample_date] 
    stop_observational_period_inland <- df_round_inland[round == 'R018', max_sample_date] 
    cutoff_date <- df_round_inland[round == 'R016', min_sample_date] 

    stopifnot(start_observational_period_inland <= cutoff_date & stop_observational_period_inland >= cutoff_date)

    df_period <- make.df.period(start_observational_period_inland, stop_observational_period_inland, 
                                cutoff_date)

    df_round <- make.df.round(df_round_inland, df_period)
    df_round <- subset(df_round, select= ! names(df_round) %like% 'LABEL|ORIGINAL')
    rm(tmp_env)
}

if(file.exists(filename_net) & ! args$rerun )
{
    dcohords <- readRDS(filename_net)
}else{

    # get A LOT of Uniform Samples for MCMC
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
                                         df_round=df_round)
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

if(1)   
{   # compute 95% ranges
    dpred_ranges <- dcohords[, CENTROIDS[[1]] ,by=GROUP]
    stopifnot( dpred_ranges[, .N == uniqueN(ID)] )
    stopifnot( dpred_ranges[, all(IL <= M & M <= IU) ] )
    stopifnot( dpred_ranges[, all(CL <= IL & IU <= CU) ] )

    tmp <- merge(dpred_ranges, drange, by.x='ID', by.y='AID')
    stopifnot( tmp[, all(IL >= MIN)] )
    stopifnot( tmp[, all(IU <= MAX)] )
}

if(nrow(df_round))
{   # compute probability of assigning source to different period
    dprobs_roundallocation <- dcohords[ , RPROBS[[1]], by='GROUP']
    dprobs_roundallocation <- merge(
                                    df_round[, .(ROUND, INDEX_TIME)],
                                    dprobs_roundallocation,
                                    by='ROUND', all.y=TRUE)
    dprobs_roundallocation <- dprobs_roundallocation[, 
            .(P=sum(P)),
            by=c('SOURCE','RECIPIENT','INDEX_TIME')]
    
    setkey(dprobs_roundallocation, RECIPIENT, SOURCE, INDEX_TIME)

    dprobs_roundallocation[is.na(INDEX_TIME), INDEX_TIME := 0 ]

    stopifnot(dprobs_roundallocation[, abs(sum(P) - 1) <= .00001, by=c('SOURCE', 'RECIPIENT')][,all(V1)])

    dprobs_roundallocation <- dcast(dprobs_roundallocation, SOURCE+RECIPIENT ~ INDEX_TIME)
    setnames(dprobs_roundallocation, c('0', '1', '2'), c('BeforeR10', 'R10_R15', 'R16_18'))
    cols <- c('BeforeR10', 'R10_R15', 'R16_18')
    dprobs_roundallocation[, (cols) := lapply(.SD, function(x){x[is.na(x)] <- 0; x} ) , .SDcols=cols]
    
    filename <- file.path(out.dir, paste0('probs_roundallocations', suffix, '.rds'))
    saveRDS(dprobs_roundallocation, filename)
}

# prepare input for Bayesian model
centroids <- dcohords[, CENTROIDS[[1]], by='GROUP']
dresults <- prepare.pairs.input.for.bayesian.model()
dhomosexualpairs[, DIRECTION := 'phyloscanner']
rbind(
    dresults,
    dhomosexualpairs,
    fill=TRUE
) -> dresults

filename <-file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('pairsdata_toshare', suffix, '.rds')) 
saveRDS(dresults, filename)    


# make table
cols <- c('SOURCE', 'RECIPIENT', 'SEX.SOURCE', 'SEX.RECIPIENT', 'M', 'AGE_INFECTION.SOURCE', 'AGE_INFECTION.RECIPIENT', 'DIRECTION')
dtable <- subset(dresults,select=cols)
dtable[, table(SEX.SOURCE,SEX.RECIPIENT)]
setnames(dtable, 'M', 'INFECTION_DATE')
filename <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('pairsdata_table_toshare', suffix, '.csv'))
fwrite(dtable, filename)


if(0)
{
    # Compare against previous ones
    setkey(centroids, ID)
    setkey(dpairs, RECIPIENT)
    dpairs2 <- merge(dpairs,  centroids[, .(DOI2.RECIPIENT=M, RECIPIENT=ID)], by='RECIPIENT')
    dpairs2 <- merge(dpairs2, centroids[, .(DOI2.SOURCE=M, SOURCE=ID)], by='SOURCE')

    drounds <- load.rounds()[COMM=='aggregated']
    drounds[, max_sample_date := shift(min_sample_date, -1)]
    drounds[round=='R018', max_sample_date := as.Date('2018-06-01') ]

    dclass <- data.table(
                         uniform=dpairs2[, classify.transmission.round.given.centroid(DOI.RECIPIENT)],
                         weighted=dpairs2[, classify.transmission.round.given.centroid(DOI2.RECIPIENT)] 
    )
# knitr::kable(table(dclass))

    plot.scatter.comparison.vs.uniform(dpairs2)

    p1 <- plot.phylopair.dates.scores(dpairs2, doi.center.var = "DOI" )
    p2 <- plot.phylopair.dates.scores(dpairs2, doi.center.var = "DOI2" )
    p <- ggpubr::ggarrange(p1, p2)
    filename <- file.path(out.dir,'compare_uniform_weighted_dois',suffix,'.png')
    ggsave(filename, p, width=12, height=8)
}

if(1)   # study generation intervals
{       
    dgen <- dcohords[, GENINTS[[1]], by='GROUP']
    dgen <- dgen[, mean(V1), by=GI]
    dgen[GI %in% c(0.5, 1), {
            cat('We estimate', round(100*V1, 2), '% generation intervals to be less than', GI, 'years\n')
            V1
    }, by='GI']
}


if(0)   # Check whether MCMC performs well in 2d case: YES
{
        dcohords[ N_IDS == 2, {
                cat(GROUP, '\n')
                DT <- IDS[[1]]
                out <- get.volume.genints.multid(DT, N_IDS=N_IDS, range_gi=range_gi, verbose=T)
                list(out, NAME=c('CENTROID', 'VOLUME', 'GI'))
        } , by=GROUP] -> tmp
        Reduce('merge',
                list(tmp[NAME == 'CENTROID', .(GROUP, CENTROIDS=out)],
                     tmp[NAME == 'VOLUME', .(GROUP, VOLUME=unlist(out))],
                     tmp[NAME == 'GI', .(GROUP, GENINTS=out)])
        ) -> tmp

        cols <- c('GROUP', 'CENTROIDS', 'VOLUME', 'GENINTS')
        tmp0 <- subset(dcohords[ N_IDS == 2], select=cols)

        dcomp <- rbind(tmp[, TYPE:='MCMC'], tmp0[, TYPE:='EXACT'])
        # compuate volume errors
        dcomp[, (VOLUME[TYPE=='MCMC'] - VOLUME[TYPE=='EXACT'])/VOLUME[TYPE=='EXACT'] , by=GROUP]
        # compute GENINTS errors
        dcomp[GROUP == 1, (GENINTS[[1]]$PR - GENINTS[[2]]$PR),  by='GROUP']
        dcomp[, {
                z0 <- GENINTS[[1]]$PR
                z1 <- GENINTS[[2]]$PR
                max((z0 - z1)/z1)
        },  by='GROUP']
        # compute CENTROIDS errors
        dcomp[, {
                z0 <- unlist(CENTROIDS[[1]])
                z1 <- unlist(CENTROIDS[[1]])
                norm(z0-z1, type='2')
        }, by='GROUP']
}

if(0)   
{   # interesting case with source RK-B110064
    range_gi
    tmp0 <- dcohords[ N_IDS == 6 & GROUP==129, GENINTS[[1]] ]
    tmp0[, list(X=GI-.5, Y=diff(c(0, PR))), by=c('SOURCE','RECIPIENT')] |> 
            ggplot(aes(X, Y, color=RECIPIENT)) +
                    geom_line() + 
                    theme_bw()

    setkey(dpairs, SOURCE, RECIPIENT)
    idx <- dcohords[GROUP == 129, IDS[[1]]]
    dpairs[idx, ..cols0]

    cols <- names(dpairs)[names(dpairs) %like% 'REC|date_first_positive.SOURCE']
    cols0 <-names(dpairs)[names(dpairs) %like% '^SOURCE|^RECIPIENT|MAX|MIN']
    dpairs[RECIPIENT %like% 'G109883' , ]
    dpairs[SOURCE %like% 'B110064' , ..cols0]
    dpairs[RECIPIENT %like% 'G109883', ..cols]
    dpairs[RECIPIENT %like% 'G109883', ..cols0]
}

if(0)   
{   # Look at what is the minimum upper bound for generation time
    dgen[PR == 0, range(GI)]
    dgen[PR == 0 & GI >= 8 ] 

    idx <- dgen[PR == 0 & GI >= 8, unique(.(SOURCE=SOURCE, RECIPIENT=RECIPIENT))]

    setkey(dpairs, SOURCE, RECIPIENT)
    cols0 <-names(dpairs)[names(dpairs) %like% '^SOURCE|^RECIPIENT|MAX|MIN']
    dpairs[idx, ..cols0]

    dpairs[ , range(MIN.RECIPIENT - MAX.SOURCE)/365.25]
}
